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

scATAC <- readRDS("../../03-scATAC-Analysis/01-Clustering/scATAC-RAM-Cluster.rds")

# quantify gene activity
gene.activities <- GeneActivity(scATAC, extend.upstream = 1000)

# add gene activities as a new assay
scATAC[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(scATAC) <- "ACTIVITY"
scATAC <- NormalizeData(scATAC)
scATAC <- ScaleData(scATAC, features = rownames(scATAC))

# 保存数据
saveRDS(scATAC, file = "scATAC-RAM-GeneActivity.rds")

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
scATAC <- readRDS("scATAC-RAM-GeneActivity.rds")

# Load the pre-processed scRNA-seq data
scRNA <- readRDS("../../04-scRNA-integrated/03-annotation/scRNA-RAM-annotation.rds")

scRNA_ <- subset(x = scRNA, subset = orig.ident == "scRNA-lib2")
cellinfo <- subset(scRNA_@meta.data, select = c("orig.ident", "nCount_RNA", "nFeature_RNA", "cellType"))
scRNA <- CreateSeuratObject(scRNA_[["RNA"]]@counts, meta.data = cellinfo)
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, rownames(scRNA))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = scRNA, query = scATAC,
                                        features = VariableFeatures(object = scRNA),
                                        reference.assay = "RNA", query.assay = "ACTIVITY",
                                        normalization.method = "LogNormalize",
                                        reduction = "cca", dims = 1:20)

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = scRNA$cellType,
                                     weight.reduction = scATAC[["lsi"]], dims = 2:20)

scATAC <- AddMetaData(scATAC, metadata = celltype.predictions)

table(scATAC$seurat_clusters, scATAC$predicted.id)

# 保存数据
saveRDS(scATAC, file = "RAM-rna-atac-Integrated.rds")