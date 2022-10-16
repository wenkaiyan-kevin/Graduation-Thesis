library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

# ldat是一个列表：ldat$spliced、ldat$unspliced、ldat$ambiguous
ldat <- ReadVelocity(file = "/gss1/home/yanwk/ywk_Graduation_Project/01-Preprocessing/02-Rice/02-Rice-Shoot/01-Mydata/scRNA_Rice_SAM/velocyto/scRNA_Rice_SAM.loom")

# 导入之前整合后的Seurat对象
scRNA <- readRDS("../../02-Clustering/02-Rice/03-SAM/02-Cluster/scRNA-SAM-Cluster.rds")
#DefaultAssay(scRNA) <- 'RNA'
features <- rownames(scRNA)
cells <- colnames(scRNA)

# 更换spliced矩阵中细胞的名字
colnames(ldat$spliced) <- gsub("x","-1",colnames(ldat$spliced))
colnames(ldat$spliced) <- gsub("scRNA_Rice_SAM:","scRNA-SAM_",colnames(ldat$spliced))
# 更换unspliced矩阵中细胞的名字
colnames(ldat$unspliced) <- colnames(ldat$spliced)
# 更换ambiguous矩阵中细胞的名字
colnames(ldat$ambiguous) <- colnames(ldat$spliced)

bm <- as.Seurat(x = ldat)
bm <- subset(bm, features = features)
bm <- subset(bm, cells = cells)
bm[["RNA"]] <- bm[["spliced"]]

DefaultAssay(bm) <- "RNA"

bm@reductions <- scRNA@reductions
bm@reductions$pca@assay.used <- 'RNA'
bm@reductions$tsne@assay.used <- 'RNA'
bm@reductions$umap@assay.used <- 'RNA'
bm$seurat_clusters <- scRNA$seurat_clusters

SaveH5Seurat(bm, filename = "scRNA-SAM.h5Seurat")
Convert("scRNA-SAM.h5Seurat", dest = "h5ad")