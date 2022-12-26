### 数据预处理
sort -k1,1 -k2,2n -k3,3n fragments.tsv > atac_fragments.sort.bed

sort -k 1,1 -k2,2n -k3,3n gene.bed > genes.sort.bed

awk '{if($6 == "+"){print $1"\t"$2-2000"\t"$2"\t"$4"\t"$4"\t"$6}else{print $1"\t"$3"\t"$3+2000"\t"$4"\t"$4"\t"$6}}' \
genes.sort.bed | sort -k 1,1 -k2,2n -k3,3n > promoters.sort.bed

bedmap --ec --delim "\t" --echo --echo-map-id promoters.sort.bed atac_fragments.sort.bed > atac_promoters_bc.bed
bedmap --ec --delim "\t" --echo --echo-map-id genes.sort.bed atac_fragments.sort.bed > atac_genes_bc.bed

###
library(Seurat)
library(rliger)
library(Signac)
library(tidyverse)
library(patchwork)
library(SeuratWrappers)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 4)
set.seed(1234)

##===scRNA与scATAC的整合===##

## 处理scATAC数据
genes.bc <- read.table(file = "../01-Preprocessing/atac_genes_bc.bed", sep = "\t", header = FALSE)
promoters.bc <- read.table(file = "../01-Preprocessing/atac_promoters_bc.bed", sep = "\t", header = FALSE)

bc <- genes.bc[,7]
bc_split <- strsplit(bc,";")
bc_split_vec <- unlist(bc_split)
bc_unique <- unique(bc_split_vec)
bc_counts <- table(bc_split_vec)
bc_filt <- names(bc_counts)[bc_counts > 800] #剩余的细胞数为 10047
barcodes <- bc_filt

gene.counts <- makeFeatureMatrix(genes.bc, barcodes)
promoter.counts <- makeFeatureMatrix(promoters.bc, barcodes)

gene.counts <- gene.counts[order(rownames(gene.counts)),]
promoter.counts <- promoter.counts[order(rownames(promoter.counts)),]
atac_count <- gene.counts + promoter.counts
colnames(atac_count) <- paste0("atac_",colnames(atac_count))

scATAC <- CreateSeuratObject(atac_count, min.cells=10, min.features=200)

#> scATAC
#An object of class Seurat
#31011 features across 10047 samples within 1 assay
#Active assay: RNA (31011 features, 0 variable features)

## 处理scRNA数据
scRNA_orig <- readRDS("../../../04-scRNA-integrated/03-annotation/scRNA-RAM-annotation.rds")
scRNA_ <- subset(x = scRNA_orig, subset = orig.ident == "scRNA-lib2")
cellinfo <- subset(scRNA_@meta.data, select = c("orig.ident", "nCount_RNA", "nFeature_RNA", "cellType", "seurat_clusters"))
colnames(cellinfo) <- c("orig.ident","nCount_RNA","nFeature_RNA","cellType","RNA_clusters")
scRNA <- CreateSeuratObject(scRNA_[["RNA"]]@counts, meta.data = cellinfo)
scRNA$orig.ident <- 'rna'

#> scRNA
#An object of class Seurat
#28444 features across 19658 samples within 1 assay
#Active assay: RNA (28444 features, 0 variable features)

int.data <- merge(scATAC, scRNA)

##liger整合流程
#数据标准化
int.data <- NormalizeData(int.data)
int.data <- FindVariableFeatures(int.data)
int.data <- ScaleData(int.data, split.by="orig.ident", do.center=FALSE) # do.center=FALSE 必须的
#设置矩阵分解的因子数，一般取值20-40 
#因式分解
int.data <- RunOptimizeALS(int.data, k=20, split.by="orig.ident")
#多样本整合
int.data <- RunQuantileNorm(int.data, split.by="orig.ident")
#整理因子顺序
int.data$clusters <- factor(int.data$clusters, levels=1:length(levels(int.data$clusters)))
#聚类
int.data <- FindNeighbors(int.data, reduction="iNMF", dims=1:20) %>% FindClusters(resolution = 0.8)
int.data <- RunUMAP(int.data, dims=1:20, reduction="iNMF")

# meta.data中的"clusters"为liger的聚类结果，seurat_clusters为seurat的结果
Idents(int.data) <- "clusters"

##可视化
p <- DimPlot(int.data, group.by = "orig.ident", cols = c("#2EB872", "#A3DE83"), pt.size = 0.01)
ggsave("atac-rna_integrated_liger.pdf", p, width = 5, height = 4)

p <- DimPlot(int.data, group.by = "orig.ident", split.by = "orig.ident",
             cols = c("#2EB872", "#A3DE83"), pt.size = 0.01)
ggsave("atac-rna_integrated_liger-split.pdf", p, width = 8, height = 4)


colors <- c('#D5EEBB','#7FC8A9','#5F7A61','#444941',
            '#F67280','#C06C84','#6C5B7B','#355C7D',
            '#F5E8C7','#DEBA9D','#9E7777','#6F4C5B',
            '#FBB448','#E3670C','#04009A','#77ACF1',
            '#3EDBF0','#867AE9','#FFF5AB','#C449C2',
            '#822659','#B34180','#E36BAE','#FF0000')

p <- DimPlot(int.data, group.by="clusters", label=T, label.size=3, cols = colors) + ggtitle("Clustered by liger") 
ggsave("atac-rna_clustered_liger.pdf", p, width=7.5, height=6)

p <- DimPlot(int.data, group.by="clusters", label=F, label.size=3, cols = colors) + ggtitle("Clustered by liger") 
ggsave("atac-rna_clustered_liger-onlabel.pdf", p, width=7.5, height=6)


#table(int.data$orig.ident, int.data$clusters)
        1    2    3    4    5    6    7    8    9   10   11   12   13   14
atac  559  744  629  545  484  334  241  458  637  444  725  670  736  679
rna   396  678 2523  919  310  859 2902  363  784  792   38  680 1091  411

       15   16   17   18   19   20
atac  254  194  469  312  692  241
rna   637  406 2633 2287  684  265

# 保存数据
saveRDS(int.data, file = "atac-rna_integrated_liger.rds")

### 数据分布作图
library(rliger)

data.int <- readRDS("../02-Intergrated/atac-rna_integrated_liger.rds")

data <- prop.table(table(data.int$orig.ident, data.int$clusters), 2)

write.csv(data, "source_liger_clusters.csv",row.names=T, quote=F)

## 注意 增加列名 Sample

####------------------------------ 作图 ------------------------------####
library(reshape2)
library(dplyr)
library(ggpubr)
library(ggalluvial)

data <- read.csv('source_liger_clusters.csv',header = T, check.names = F)

melted <- melt(data) %>% rename(Count = value, Cluster = variable) %>% collect()

color_clusters <- c("#FFD124", "#00AFC1")

p <- ggplot(melted, aes(x = Cluster, fill = Sample, stratum = Sample,
            alluvium = Sample, y = Count, label = Sample)) +
     ylab('Proportion') + 
     geom_flow() + geom_stratum() +
     scale_fill_manual(values = color_clusters) +
     theme_pubr(base_size = 16, legend = "right", border = T) + rotate_x_text(0)

ggsave("source_liger_clusters.pdf", p, width = 8, height = 4)

