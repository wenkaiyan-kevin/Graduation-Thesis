library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 10)
rm(list = ls())
set.seed(1234)
#1.1. 读取数据
## 重新创建没有处理的经过降维等处理的数据

scRNA_Zhang <- readRDS("~/ywk_Graduation_Project/05-Root-project/02-scRNA-Seurat2AnnData/01-Zhang/scRNA_annotation.rds")
scRNA_Zhang$cellType <- Idents(scRNA_Zhang)
scRNA_Zhang$lib <- "scRNA_Zhang"
cellinfo <- subset(scRNA_Zhang@meta.data, select = c("orig.ident", "nCount_RNA", "nFeature_RNA", "lib", "cellType"))
scRNA_Zhang <- CreateSeuratObject(scRNA_Zhang[["RNA"]]@counts, meta.data = cellinfo)
scRNAlist_Zhang <- SplitObject(scRNA_Zhang, split.by = "orig.ident")

scRNA_Liu <- readRDS("~/ywk_Graduation_Project/05-Root-project/02-scRNA-Seurat2AnnData/02-Liu/scRNA-Nip-Cluster-Annotation.rds")
scRNA_Liu$cellType <- Idents(scRNA_Liu)
scRNA_Liu$lib <- "scRNA_Liu"
cellinfo <- subset(scRNA_Liu@meta.data, select = c("orig.ident", "nCount_RNA", "nFeature_RNA", "lib", "cellType"))
scRNA_Liu <- CreateSeuratObject(scRNA_Liu[["RNA"]]@counts, meta.data = cellinfo)

scRNAlist <- c(scRNAlist_Zhang, scRNA_Liu)

# 数据标准化与整合
scRNAlist <- lapply(X = scRNAlist, FUN = SCTransform, vst.flavor = "v2")

scRNA.features <- SelectIntegrationFeatures(scRNAlist, nfeatures = 3000)

scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = scRNA.features)

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,
                                        normalization.method = "SCT",
                                        anchor.features = scRNA.features)
scRNA.combined.sct <- IntegrateData(scRNA.anchors, normalization.method = "SCT")

# 降维聚类分群
scRNA.combined.sct <- RunPCA(scRNA.combined.sct, npcs = 100, verbose = FALSE)
scRNA.combined.sct <- FindNeighbors(scRNA.combined.sct, dims = 1:60)
scRNA.combined.sct <- FindClusters(scRNA.combined.sct, resolution = 0.5)
scRNA.combined.sct <- RunUMAP(scRNA.combined.sct, dims = 1:30)

# 作图
# 定义颜色
colors <- c('#D5EEBB','#7FC8A9','#5F7A61','#444941',
            '#F67280','#C06C84','#6C5B7B','#355C7D',
            '#F5E8C7','#DEBA9D','#9E7777','#6F4C5B',
            '#FBB448','#E3670C','#04009A','#77ACF1',
            '#3EDBF0','#867AE9','#FFF5AB','#C449C2',
            '#822659','#B34180','#E36BAE','#FF0000')

p1 <- DimPlot(scRNA.combined.sct, reduction = "umap", label = F, cols = c('#9ecae1','#3182bd'), group.by = 'lib', pt.size=0.05)
p2 <- DimPlot(scRNA.combined.sct, reduction = "umap", label = T, cols = colors, pt.size=0.05)
p <- p1 + p2
ggsave("scRNA_SCT_UMAP_Cluster.pdf", p, width = 14, height = 6)

p3 <- DimPlot(scRNA.combined.sct, reduction = "umap", cols = c('#9ecae1','#3182bd'), group.by = "lib", 
              split.by = "lib", ncol = 2, pt.size=0.01)
ggsave("scRNA_SCT_UMAP_Cluster_Split.pdf", p3, width = 10, height = 5)


p3 <- DimPlot(scRNA, reduction = "umap", cols = colors, group.by = "cellType", 
              split.by = "lib", ncol = 2, pt.size=0.01)
ggsave("test.pdf", p3, width = 10, height = 5)



####
p1 <- DimPlot(scRNA.combined.sct, reduction = "umap", label = F, cols = c('#9ecae1','#3182bd'), group.by = 'lib', pt.size=0.05)
p2 <- DimPlot(scRNA.combined.sct, reduction = "umap", label = F, cols = colors, pt.size=0.05)
p <- p1 + p2
ggsave("scRNA_SCT_UMAP_Cluste_.pdf", p, width = 14, height = 6)

# 保存数据
saveRDS(scRNA.combined.sct, file = "scRNA-RAM-Integrate.rds")

