library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 4)

# 读取数据-------------------------------------------------------------------------------------
dir <- "/gss1/home/yanwk/ywk_Graduation_Project/01-Preprocessing/02-Rice/02-Rice-Shoot/01-Mydata/scRNA_Rice_SAM/outs/filtered_feature_bc_matrix"
samples_name <- "scRNA-SAM"

# 创建Seurat对象-------------------------------------------------------------------------------
counts <- Read10X(data.dir = dir, gene.column = 1)
scRNA <- CreateSeuratObject(counts, project = samples_name, min.cells = 5, min.features = 200)
scRNA <- RenameCells(scRNA, add.cell.id = samples_name)
scRNA[["percent.Mt"]] <- PercentageFeatureSet(scRNA, pattern = '^gene-')

# 数据质控 ------------------------------------------------------------------------------------
# 设置绘图元素
plot.features <- c("nFeature_RNA", "nCount_RNA", "percent.Mt")
theme.set2 <- theme(axis.title = element_blank())
p <- VlnPlot(scRNA, pt.size = 0, features = plot.features) + theme.set2 + NoLegend()
dir.create("QC")
ggsave("QC/vlnplot_before_qc.pdf", plot = p, width = 6, height = 4)

### 设置质控标准-----------------------------------------------------------------------------------
maxGene = 7177.306
minGene = 696.694
maxUMI = 39751.56
minUMI = 4534
pctMt = 5

### 数据质控并绘制小提琴图
scRNA <- subset(scRNA, subset = nCount_RNA < maxUMI & 
                                nCount_RNA > minUMI &
                                nFeature_RNA < maxGene &
                                nFeature_RNA > minGene &
                                percent.Mt < pctMt)

# scRNA
An object of class Seurat 
26892 features across 9264 samples within 1 assay 
Active assay: RNA (26892 features, 0 variable features)

# 质控后小提琴图
p <- VlnPlot(scRNA, pt.size = 0, features = plot.features) + theme.set2 + NoLegend()
ggsave("QC/vlnplot_after_qc.pdf", plot = p, width = 6, height = 4)

# 保存数据
saveRDS(scRNA, file = "scRNA-SAM-QC.rds")