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
# Sample-SAM------------------------------------------------------------
scRNA <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/02-Clustering/02-Rice/03-SAM/01-Quality-Control/scRNA-SAM-QC.rds")
Features <- grep("gene", rownames(scRNA), value=T, invert = T)
scRNA <- subset(scRNA, features = Features)

###
scRNA <- SCTransform(scRNA, vst.flavor = "v2", vars.to.regress = c("nCount_RNA","nFeature_RNA"))

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase 
scRNA <- RunPCA(scRNA, features = VariableFeatures(object=scRNA), verbose = FALSE)

# 聚类、分群、降维 ----------------------------------------------------------------------------
scRNA <- FindNeighbors(scRNA, dims = 1:10, annoy.metric='cosine', l2.norm = T, n.trees = 100)
scRNA <- FindClusters(scRNA, resolution = 0.2)
scRNA <- RunUMAP(scRNA, dims = 1:5, min.dist = 0.1, n.neighbors = 120, n.epochs=1000)
scRNA <- RunTSNE(scRNA, dims = 1:5)

# 作图 ----------------------------------------------------------------------------
cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
mycolors <- colorRampPalette(cluster_cols,space="rgb")(length(levels(scRNA)))

p1 <- DimPlot(scRNA, reduction = "tsne", cols=alpha(mycolors,0.7), pt.size=1)
p2 <- DimPlot(scRNA, reduction = "umap", cols=alpha(mycolors,0.7), pt.size=1)
p <- p1 + p2
ggsave("scRNA-SAM-Clusters.pdf", plot = p, width = 12, height = 6)

# 保存数据
saveRDS(scRNA, file = "scRNA-SAM-Cluster.rds")


## 细胞类群重新命名 ----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 10)
rm(list = ls())
set.seed(1234)

scRNA <- readRDS("../../../02-Clustering/02-Rice/03-SAM/02-Cluster/scRNA-SAM-Cluster.rds")
DefaultAssay(scRNA) <- 'RNA'
scRNA <- NormalizeData(scRNA, verbose = FALSE)

new.cluster.ids <- c("mesophyll cell", "meristerm epidermis cell", "vascular cell", "shoot meristerm cell",
                     "proliferating cell", "leaf epidermal cell", "undefined cell")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA$celltype <- Idents(scRNA)

## save seurat object
saveRDS(scRNA, "scRNA-SAM-annotation.rds")