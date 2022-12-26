library(Seurat)
library(tidyverse)
library(harmony)
library(patchwork)
library(dplyr)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 10)
rm(list = ls())
set.seed(1234)

#1.1. 读取数据
# Sample-Zhang------------------------------------------------------------
scRNA_Zhang <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/02-Clustering/01-Ath/01-Ath-Root/01-ZhangTianQi/scRNA-Zhang-QC.rds")
scRNA_Zhang$source <- "scRNA-Zhang"

# Sample-Ryu------------------------------------------------------------
scRNA_Ryu <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/02-Clustering/01-Ath/01-Ath-Root/02-Ryu/scRNA-Ryu-QC.rds")
scRNA_Ryu$source <- "scRNA-Ryu"
scRNAlist_Ryu <- SplitObject(scRNA_Ryu, split.by = "orig.ident")

# Sample-Denyer------------------------------------------------------------
scRNA_Denyer <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/02-Clustering/01-Ath/01-Ath-Root/03-Denyer/scRNA-Denyer-QC.rds")
scRNA_Denyer$source <- "scRNA_Denyer"

# Sample-Farmer------------------------------------------------------------
scRNA_Farmer <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/02-Clustering/01-Ath/01-Ath-Root/04-Andrew-Farmer/scRNA-Farmer-QC.rds")
scRNA_Farmer$source <- "scRNA-Farmer"
scRNAlist_Farmer <- SplitObject(scRNA_Farmer, split.by = "orig.ident")

# Sample-Wendrich------------------------------------------------------------
scRNA_Wendrich <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/02-Clustering/01-Ath/01-Ath-Root/05-Wendrich/scRNA-Wendrich-QC.rds")
scRNA_Wendrich$source <- "scRNA-Wendrich"
scRNAlist_Wendrich <- SplitObject(scRNA_Wendrich, split.by = "orig.ident")

scRNAlist <- c(scRNA_Zhang, scRNAlist_Ryu, scRNA_Denyer, scRNAlist_Farmer, scRNAlist_Wendrich)

names(scRNAlist)[1] <- 'scRNA-Zhang'
names(scRNAlist)[5] <- 'scRNA-Denyer'

saveRDS(scRNAlist, "scRNAlist.rds")

##==harmony整合多样本==##
library(Seurat)
library(tidyverse)
library(harmony)
library(patchwork)
library(dplyr)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 10)
rm(list = ls())
set.seed(1234)

scRNAlist <- readRDS("scRNAlist.rds")

##数据标准化与PCA降维
scRNA_harmony <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], scRNAlist[[4]], scRNAlist[[5]], 
                                           scRNAlist[[6]], scRNAlist[[7]], scRNAlist[[8]], scRNAlist[[9]], 
                                           scRNAlist[[10]], scRNAlist[[11]], scRNAlist[[12]], scRNAlist[[13]]))

scRNA_harmony <- SCTransform(scRNA_harmony, vst.flavor = "v2", verbose = FALSE)
scRNA_harmony <- RunPCA(scRNA_harmony, npcs = 100, verbose = FALSE)

##整合
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident", assay.use="SCT", max.iter.harmony = 20)

#降维聚类
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.5)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30)

# 保存结果
saveRDS(scRNA_harmony, file = "scRNA-RAM-Ath-Harmony-Integrate.rds")

##作图
my_color <- c("#EF3500","#FF5E47","#F7D098","#B7BCF3","#4E9FAA"
             ,"#AA91AA","#957AB9","#844AA7","#405DAC","#00909B"
             ,"#37C7CE","#477E63","#95B1D9","#71C89C","#50AE7A"
             ,"#00ED94","#BDE29B","#D6D18D","#B49685","#6F7E52"
             ,"#F0EFD0","#FF98A8","#D473BB","#DA4293")

p1 <- DimPlot(scRNA_harmony, reduction = "umap", label = F, cols = my_color, group.by='source', pt.size=0.05)
p2 <- DimPlot(scRNA_harmony, reduction = "umap", label = T, cols = my_color, pt.size=0.05)
p <- p1 + p2
ggsave("scRNA-UMAP-Ath-Harmony.pdf", p, width = 16, height = 6)

p1 <- DimPlot(scRNA_harmony, reduction = "tsne", label = F, cols = my_color, group.by='source', pt.size=0.05)
p2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = T, cols = my_color, pt.size=0.05)
p <- p1 + p2
ggsave("scRNA-TSNE-Ath-Harmony.pdf", p, width = 16, height = 6)

## 按照数据源进行作图
cells_ <- Cells(subset(x = scRNA_harmony, subset = source == "scRNA-Wendrich"))
p <- DimPlot(scRNA_harmony, reduction = "tsne", cells = cells_,label = F, cols = my_color, pt.size=0.05)
ggsave("scRNA-TSNE-Ath-Harmony-Wendrich.pdf", p, width = 8, height = 6)

library(Seurat)

scRNA <- readRDS("scRNA-RAM-Ath-Harmony-Integrate.rds")

data <- prop.table(table(scRNA$source, scRNA$seurat_clusters), 2)

write.csv(data, "source_seurat_clusters.csv",row.names=T, quote=F)

## 注意 增加列名 Sample

## ---------------------- 作图 ------------------------
library(reshape2)
library(dplyr)
library(ggpubr)
library(ggalluvial)

data <- read.csv('source_seurat_clusters.csv',header = T, check.names = F)

melted <- melt(data) %>% rename(Count = value, Cluster = variable) %>% collect()

color_clusters <- c("#844AA7","#405DAC","#00909B","#37C7CE","#477E63")

p <- ggplot(melted, aes(x = Cluster, fill = Sample, stratum = Sample,
            alluvium = Sample, y = Count, label = Sample)) +
     ylab('Proportion') + 
     geom_flow() + geom_stratum() +
     scale_fill_manual(values = color_clusters) +
     theme_pubr(base_size = 16, legend = "right", border = T) + rotate_x_text(0)

ggsave("source_seurat_clusters.pdf", p, width = 12, height = 4)

library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(RColorBrewer)
rm(list = ls())
set.seed(1234)

scRNA <- readRDS("../01-Clustering/scRNA-RAM-Ath-Harmony-Integrate.rds")

genes <- read.table("MarkerGene.txt")[,1]
genes <- as.vector(CaseMatch(search=genes, match=rownames(scRNA)))

for(i in seq(1,length(genes),9)){
   j <- i + 8
   if(j > length(genes)){
     j <- length(genes)
   }
   
   p <- FeaturePlot(scRNA, features = genes[i:j], reduction = 'umap', 
                    pt.size = 0.01, ncol = 3,
                    cols = c('#ccccca', '#e61a2e'), 
                    order = T)
   file_name <- paste(i,i+8,'heat','tsne.tiff',sep='-')
   ggsave(file_name, p, width = 18, heigh = 16)
}


for(i in seq(1,length(genes),6)){
   j <- i + 5
   if(j > length(genes)){
     j <- length(genes)
   }
   
   p <- FeaturePlot(scRNA, features = genes[i:j], reduction = 'umap', 
                    pt.size = 0.01, ncol = 3,
                    cols = c('#ccccca', '#e61a2e'), 
                    order = T)
   file_name <- paste(i,i+5,'heat','tsne.pdf',sep='-')
   ggsave(file_name, p, width = 14, heigh = 8)
}



p <- DotPlot(scRNA, features = genes) + RotatedAxis() +
        scale_x_discrete("") + scale_y_discrete("") +
        scale_color_gradientn(colours = brewer.pal(9,"BuGn")) +
        #coord_flip() +
        theme_bw() +
        theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "inches"))+
        theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45),
              axis.text.y = element_text(size=15))
ggsave("test.pdf", plot = p, width = 20, height = 10)

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

scRNA <- readRDS("../01-Clustering/scRNA-RAM-Ath-Harmony-Integrate.rds")

new.cluster.ids <- c("meristem (Atha)","vascular cylinder (Atha)","vascular cylinder (Atha)",
                     "root cap (Atha)","endodermis (Atha)",
                     "root cap (Atha)","epidermis/trichoblast (Atha)",
                     "cortex (Atha)","cortex (Atha)","vascular cylinder (Atha)",
                     "undefined cell (Atha)","epidermis/atrichoblast (Atha)",
                     "phloem (Atha)","undefined cell (Atha)",
                     "root cap (Atha)","xylem (Atha)","endodermis (Atha)",
                     "epidermis/trichoblast (Atha)","meristem (Atha)",
                     "undefined cell (Atha)","meristem (Atha)","phloem (Atha)",
                     "endodermis (Atha)","vascular cylinder (Atha)")

names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA$cellType <- Idents(scRNA)

## save seurat object
saveRDS(scRNA, "scRNA-RAM-Ath-Harmony-Integrate-annotation.rds")




