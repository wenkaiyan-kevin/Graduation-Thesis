import numpy as np
import pandas as pd
import scanpy as sc

data_dir = '/gss1/home/yanwk/ywk_SingleCell_work/10_other_Analysis/01_Tian-Qi-Zhang_Rice_Root/01_Pre-Processing/scRNA_lib2/outs/filtered_feature_bc_matrix'

adata = sc.read_10x_mtx(data_dir, var_names='gene_ids', cache=True)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)

adata.var['mt'] = adata.var_names.str.startswith('gene-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save = 'violin.pdf')

adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.n_genes_by_counts > 500, :]
adata = adata[adata.obs.total_counts < 40000, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# 保存数据
adata.write('scRNA-RAM.h5ad')

####
import anndata
import networkx as nx
import scanpy as sc
import pandas as pd

# 读取表观组
# 需要对peak文件修改：awk '{print $1":"$2"-"$3}' peaks.bed > peaks.tsv
adata = sc.read_mtx("matrix.mtx").T
bc = pd.read_table("barcodes.tsv", header = None)
peak = pd.read_table("peaks.tsv", header = None)
adata.obs_names = bc[0]
adata.var_names = peak[0]

# 保存数据
adata.write('scATAC-RAM.h5ad')


###开始处理数据
# import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams

scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

rna = ad.read_h5ad("../01-Preprocessing/01-RNA/scRNA-RAM.h5ad")
rna

atac = ad.read_h5ad("../01-Preprocessing/02-ATAC/scATAC-RAM.h5ad")
atac

## Preprocess scRNA-seq data -------------------------------------------------
rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")
## Optionally
sc.pp.neighbors(rna, metric="cosine")
sc.tl.leiden(rna, resolution=0.4)
sc.tl.umap(rna)
sc.pl.umap(rna, color="leiden", save = 'scRNA-RAM-UMAP.pdf')

# Preprocess scATAC-seq data -------------------------------------------------
scglue.data.lsi(atac, n_components=100, n_iter=15)
## Optionally
sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
sc.tl.leiden(atac, resolution=0.6)
sc.tl.umap(atac)
sc.pl.umap(atac, color="leiden", save = 'scATAC-RAM-UMAP.pdf')

# Construct prior regulatory graph -------------------------------------------------
# Obtain genomic coordinate
scglue.data.get_gene_annotation(
    rna, gtf="/home/wkyan/ywk_lab/04-scATAC-analysis/00-seqlib/Oryza_sativa.IRGSP-1.0.53.gtf",
    gtf_by="gene_id")
rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()

# scATAC 数据集
atac.var_names[:5]
split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
atac.var.head()

# Graph construction -------------------------------------------------
guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
guidance
scglue.graph.check_graph(guidance, [rna, atac])
atac.var.head()

# Save preprocessed data files -------------------------------------------------
rna.write("rna-pp.h5ad", compression="gzip")
atac.write("atac-pp.h5ad", compression="gzip")
nx.write_graphml(guidance, "guidance.graphml.gz")


### 构建模型
from itertools import chain

import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams

scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

# Read preprocessed data -------------------------------------------------
rna = ad.read_h5ad("../02-Preprocessing/rna-pp.h5ad")
atac = ad.read_h5ad("../02-Preprocessing/atac-pp.h5ad")
guidance = nx.read_graphml("../02-Preprocessing/guidance.graphml.gz")

rna.obs['domain'] = 'rna'
atac.obs['domain'] = 'atac'

# Configure data -------------------------------------------------
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca")

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi")

guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()

# Train GLUE model -------------------------------------------------
glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory": "glue"})

glue.save("glue.dill")
# glue = scglue.models.load_model("glue.dill")

# Check integration diagnostics -------------------------------------------------
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, guidance_hvf
)
dx

# Apply model for cell and feature embedding -------------------------------------------------
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

combined = ad.concat([rna, atac])

# 聚类分群 -------------------------------------------------
sc.pp.neighbors(combined, use_rep="X_glue", n_pcs = 50, metric="cosine")
sc.tl.leiden(combined, resolution=0.7)
sc.tl.umap(combined, min_dist = 0.3)
sc.pl.umap(combined, color=["leiden", "domain"], wspace=0.65, save = 'Combined-RAM-UMAP.pdf')

sc.tl.embedding_density(combined, groupby='domain')
sc.pl.embedding_density(combined, groupby='domain', fg_dotsize=10, save = 'domain-distribution.pdf')

sc.pl.umap(combined, color='domain', groups=['rna'], save = 'rna.pdf')
sc.pl.umap(combined, color='domain', groups=['atac'], save = 'atac.pdf')

feature_embeddings = glue.encode_graph(guidance_hvf)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
feature_embeddings.iloc[:5, :5]

rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()

rna.write("rna-emb.h5ad", compression="gzip")
atac.write("atac-emb.h5ad", compression="gzip")
nx.write_graphml(guidance_hvf, "guidance-hvf.graphml.gz")


