library(Seurat)
library(Signac)
library(rtracklayer)
library(tidyverse)
library(patchwork)
library(dplyr)
library(BSgenome.Osativa.IRGSP.IRGSP1)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 4)
set.seed(1234)

# 读取数据
scATAC <- readRDS("../01-CallPeak-QC/scATAC-SAM-QC.rds")

scATAC <- RunTFIDF(scATAC)
scATAC <- FindTopFeatures(scATAC, min.cutoff = 'q0')
scATAC <- RunSVD(object = scATAC, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

# 查看测序深度与数据特征之间的相关性
p1 <- DepthCor(scATAC)
ggsave("DepthCor.pdf", plot = p1, width = 8, height = 4)

scATAC <- FindNeighbors(object = scATAC, reduction = 'lsi', dims = 2:20)
scATAC <- FindClusters(object = scATAC, resolution = 0.8, algorithm = 3)
scATAC <- RunUMAP(object = scATAC, reduction = 'lsi', dims = 2:20) 

colors <- c("#E95C59", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3")
p2 <- DimPlot(scATAC, reduction = "umap", pt.size = 0.8, cols = colors)
ggsave("scATAC-SAM-Clusters-Seurat.pdf", plot = p2, width = 6, height = 6)

#> table(scATAC$seurat_clusters)
#   0    1    2    3    4
#2242 1841 1610 1508  799

# 保存数据
saveRDS(scATAC, file = "scATAC-SAM-Cluster.rds")

# 转换数据格式，用于作图
library(SeuratDisk)

SaveH5Seurat(scATAC, filename = "scATAC.h5Seurat")
Convert("scATAC.h5Seurat", dest = "h5ad")


########## 聚类图形在Python中重做 ###########
import scanpy as sc
from matplotlib.pyplot import rc_context

adata = sc.read_h5ad("scATAC.h5ad")
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')


#colors = ["#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3"]
#colors = ["#476D87", "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A"]
colors = ["#E95C59", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3"]
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color='seurat_clusters', add_outline=True, legend_loc='right margin',
           legend_fontsize=12, legend_fontoutline=2, frameon=False, size=50,
           title='clustering of cells', palette=colors, save = 'scATAC-SAM-Clusters-Scanpy.pdf')

##################################