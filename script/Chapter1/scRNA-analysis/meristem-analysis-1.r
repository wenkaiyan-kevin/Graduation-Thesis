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

scRNA <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/04-Root-Shoot-Project/02-Cluster-Annotation/02-Cluster-rename/scRNA-SAM-annotation.rds")
# 作图
mycolors <- c("grey50","grey50","grey50","#A6D854","grey50","grey50","grey50")
p <- DimPlot(scRNA, reduction = "umap", cols=mycolors, pt.size=1)
ggsave('scRNA-SMC-highlight.pdf', p, width = 8, height = 5)


#### 提取分生组织的细胞 #######
scRNA_SMC <- subset(scRNA, idents = 'shoot meristerm cell')

scRNA_SMC <- SCTransform(scRNA_SMC, vst.flavor = "v2", vars.to.regress = c("nCount_RNA","nFeature_RNA"))

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase 
scRNA_SMC <- RunPCA(scRNA_SMC, features = VariableFeatures(scRNA_SMC), verbose = FALSE)

# 聚类、分群、降维 ----------------------------------------------------------------------------
scRNA_SMC <- FindNeighbors(scRNA_SMC, dims = 1:5, annoy.metric='cosine', l2.norm = T, n.trees = 100)
scRNA_SMC <- FindClusters(scRNA_SMC, resolution = 0.1)
#scRNA_SMC <- RunUMAP(scRNA_SMC, dims = 1:5, min.dist = 0.1, n.neighbors = 120, n.epochs=1000)
scRNA_SMC <- RunUMAP(scRNA_SMC, dims = 1:5)
scRNA_SMC <- RunTSNE(scRNA_SMC, dims = 1:5)

# 作图 ----------------------------------------------------------------------------
mycolors <- c("#E64B35B2","#4DBBD5B2","#00A087B2")

p1 <- DimPlot(scRNA_SMC, reduction = "tsne", cols=alpha(mycolors,1), pt.size=1)
p2 <- DimPlot(scRNA_SMC, reduction = "umap", cols=alpha(mycolors,1), pt.size=1)
p <- p1 + p2
ggsave("scRNA-SMC-Cluster.pdf", plot = p, width = 12, height = 6)

# 保存数据
saveRDS(scRNA_SMC, file = "scRNA-SMC-Cluster.rds")


#-------------------------------- 基因集富集表达分析 ----------------------------------------
library(ggplot2)
library(Seurat)

mytheme <-theme_bw() +  
          theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
                plot.title = element_text(size = rel(2),hjust = 0.5),
                panel.grid.minor = element_blank(),
		            panel.grid.major = element_blank(),
		            panel.border = element_blank(),
		            panel.background = element_blank(),
                axis.line = element_line(color = 'black',size = 1),
								axis.ticks.length=unit(0.2, "cm"),
								axis.text = element_text(size=rel(1.5)),
								axis.title = element_text(size = rel(1.5)))


scRNA <- readRDS("scRNA-SMC-Cluster.rds")
DefaultAssay(scRNA) <- 'RNA'

###### PZ
PZ_Genes <- c('Os12g0274700','Os01g0760900','Os01g0914300','Os10g0554800','Os11g0150400',
              'Os10g0159700','Os01g0638000','Os03g0845000','Os02g0729700')
scRNA <- AddModuleScore(scRNA, features = list(PZ_Genes), ctrl = 100, name = "PZ")
df <- data.frame(scRNA@meta.data, scRNA@reductions$umap@cell.embeddings)
df$Enrichment <- df$PZ1

p1 <- ggplot(df, aes(UMAP_1, UMAP_2))  +
     geom_point(aes(colour  = Enrichment), size=2) + 
     ggtitle("Peripheral Zone meristerm cell (PZ)") +
     viridis::scale_color_viridis(option="A") + mytheme

###### CZ
genes <- c('Os08g0532800','Os03g0181500','Os07g0194300','Os02g0814200')
scRNA <- AddModuleScore(scRNA, features = list(genes), ctrl = 100,name = "CZ")
df<- data.frame(scRNA@meta.data, scRNA@reductions$umap@cell.embeddings)
df$Enrichment <- df$CZ1
p2 <- ggplot(df, aes(UMAP_1, UMAP_2))  +
     geom_point(aes(colour  = Enrichment), size=2) + 
     ggtitle("Central Zone meristerm cell (CZ)") +
     viridis::scale_color_viridis(option="A") + mytheme

## 合并作图
library(patchwork)
p <- p1 / p2
ggsave('scRNA-SMC-CZ-PZ-Enrichment.pdf', p, width = 8, height = 15, limitsize = FALSE)