# ------------------------------Create a gene activity matrix---------------------------------
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

scATAC <- readRDS("../02-Clustering/scATAC-SAM-Cluster.rds")

# quantify gene activity
gene.activities <- GeneActivity(scATAC, extend.upstream = 1000, extend.downstream = 200)

# add gene activities as a new assay
scATAC[["RNA"]] <- CreateAssayObject(counts = gene.activities)

scATAC <- NormalizeData(object = scATAC, assay = 'RNA',
                        normalization.method = 'LogNormalize',
                        scale.factor = median(scATAC$nCount_RNA))
# 保存数据
saveRDS(scATAC, file = "scATAC-SAM-GeneActivity.rds")

# ------------------------------Integrating with scRNA-seq data---------------------------------
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

#加载数据
scATAC <- readRDS("scATAC-SAM-GeneActivity.rds")
# 非常重要
DefaultAssay(scATAC) <- 'RNA'

# Load the pre-processed scRNA-seq data
scRNA <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/04-Shoot-Project/01-RNA/02-Cluster-Annotation/02-Cluster-rename/scRNA-SAM-annotation.rds")
scRNA <- FindVariableFeatures(object = scRNA, nfeatures = 3000)

transfer.anchors <- FindTransferAnchors(reference = scRNA, query = scATAC, reduction = 'cca', dims = 1:30)

predicted.labels <- TransferData(anchorset = transfer.anchors, refdata = scRNA$celltype,
                                 weight.reduction = scATAC[['lsi']], dims = 2:30)

scATAC <- AddMetaData(object = scATAC, metadata = predicted.labels)

table(scATAC$seurat_clusters, scATAC$predicted.id)


# 保存数据
saveRDS(scATAC, file = "scATAC-SAM-Integrated.rds")

#plot1 <- DimPlot(scRNA, group.by = 'celltype', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
#plot2 <- DimPlot(scATAC, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
#p <- plot1 + plot2
#ggsave("test.pdf", plot = p, width = 12, height = 6)

# ------------------------------Co-embedding---------------------------------

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
scATAC <- readRDS("scATAC-SAM-GeneActivity.rds")
scRNA <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/04-Shoot-Project/01-RNA/02-Cluster-Annotation/02-Cluster-rename/scRNA-SAM-annotation.rds")
scRNA$lib <- "scRNA"

scRNA <- FindVariableFeatures(object = scRNA, nfeatures = 3000)
genes.use <- VariableFeatures(scRNA)
refdata <- GetAssayData(scRNA, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
DefaultAssay(scATAC) <- 'RNA'
transfer.anchors <- FindTransferAnchors(reference = scRNA, query = scATAC, reduction = 'cca', dims = 1:30)

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata,
                           weight.reduction = scATAC[["lsi"]], dims = 2:30)

# this line adds the imputed data matrix to the pbmc.atac object
scATAC[["RNA"]] <- imputation
coembed <- merge(x = scRNA, y = scATAC)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)


p1 <- DimPlot(coembed, group.by = "lib")
p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
p <- p1 + p2
ggsave("test.pdf", plot = p, width = 18, height = 6)



