import scvelo as scv

adata = scv.read("scRNA-SAM.h5ad")
adata

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)

# Consider installing `tqdm` as `pip install tqdm` and `ipywidgets` as `pip install ipywidgets`
scv.tl.velocity_graph(adata, n_jobs=18)
#scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters", 
#                                 legend_loc='none',
#                                 save='velocity_embedding_stream.pdf' )
#颜色：https://matplotlib.org/stable/tutorials/colors/colormaps.html

scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters",
                                 palette='Set2', legend_loc='right margin',
                                 save='velocity_embedding_stream.svg',figsize=(5,5))



#scv.tl.velocity_pseudotime(adata)
#scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')

scv.tl.recover_dynamics(adata, n_jobs=18)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color="latent_time", palette='BuGn', size=40,
               save='velocity_recover_dynamics.pdf',figsize=(5,5))



# plot velocity of a selected gene
#scv.pl.velocity(adata, var_names=['Os01g0914300'],  save='test.pdf', palette='Blues')

# ------------------------------作图---------------------------------

import scanpy as sc
import pandas as pd
from matplotlib.pyplot import rc_context

sc.set_figure_params(dpi=100, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()

adata = sc.read_h5ad('scRNA-SAM.h5ad')

adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')

with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color='seurat_clusters', add_outline=True, legend_loc='right margin',
               legend_fontsize=12, legend_fontoutline=2, frameon=False,size=50,
               title='clustering of cells', palette='Set2', save = 'test.pdf')

with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.tsne(adata, color='seurat_clusters', add_outline=True, legend_loc='on data',
               legend_fontsize=12, legend_fontoutline=2,frameon=False,size=80,
               title='clustering of cells', palette='Set1', save = 'test1.pdf')
