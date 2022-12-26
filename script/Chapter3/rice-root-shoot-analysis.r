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

# 读取水稻茎尖
scRNA_Rice_SAM <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/04-Shoot-Project/01-RNA/02-Cluster-Annotation/02-Cluster-rename/scRNA-SAM-annotation.rds")

new.cluster.ids <- c("mesophyll cell (SAM)", "meristem epidermis cell (SAM)","vascular cell (SAM)",
"shoot meristematic cell (SAM)","proliferating cell (SAM)","leaf epidermis cell (SAM)", "undefined cell (SAM)")
names(new.cluster.ids) <- levels(scRNA_Rice_SAM)
scRNA_Rice_SAM <- RenameIdents(scRNA_Rice_SAM, new.cluster.ids)
# 保存数据
saveRDS(scRNA_Rice_SAM, file = "scRNA-Rice-SAM-annotation.rds")

# 读取水稻根尖
scRNA_Rice_RAM <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/05-Root-project/04-scRNA-integrated/03-annotation/scRNA-RAM-annotation.rds")

new.cluster.ids <- c("endodermis (RAM)","cortex (RAM)","meristem (RAM)","exodermis (RAM)","dividing cell (RAM)",
"vascular cylinder (RAM)","epidermis/atrichoblast (RAM)","epidermis/meristem (RAM)","epidermis/trichoblast (RAM)",
"root cap (RAM)","pericycle (RAM)","phloem (RAM)","xylem (RAM)")
names(new.cluster.ids) <- levels(scRNA_Rice_RAM)
scRNA_Rice_RAM <- RenameIdents(scRNA_Rice_RAM, new.cluster.ids)

# 保存数据
saveRDS(scRNA_Rice_RAM, file = "scRNA-Rice-RAM-annotation.rds")

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

# 读取数据
# -------------------## Rice SAM##----------------------
scRNA_Rice_SAM <- readRDS("../01-data/scRNA-Rice-SAM-annotation.rds")
scRNA_Rice_SAM$cellType <- Idents(scRNA_Rice_SAM)
scRNA_Rice_SAM$source <- "SAM"
scRNA_Rice_SAM <- subset(scRNA_Rice_SAM, idents = "undefined cell (SAM)", invert = TRUE)

# -------------------## Rice leaf##----------------------
scRNA_Rice_Leaf <- readRDS("../01-data/scRNA-Rice-Leaf-annotation.rds")
scRNA_Rice_Leaf$cellType <- Idents(scRNA_Rice_Leaf)
scRNA_Rice_Leaf$source <- "SAM"
scRNA_Rice_Leaf_list <- SplitObject(scRNA_Rice_Leaf, split.by = "orig.ident")

# -------------------## Rice RAM##----------------------
scRNA_Rice_RAM <- readRDS("../01-data/scRNA-Rice-RAM-annotation.rds")
scRNA_Rice_RAM$cellType <- Idents(scRNA_Rice_RAM)
scRNA_Rice_RAM$source <- "RAM"
scRNA_Rice_RAM_list <- SplitObject(scRNA_Rice_RAM, split.by = "orig.ident")

scRNAlist <- c(scRNA_Rice_SAM, scRNA_Rice_Leaf_list, scRNA_Rice_RAM_list)
# 数据标准化与整合
scRNAlist <- lapply(X = scRNAlist, FUN = SCTransform, vst.flavor = "v2")

scRNA.features <- SelectIntegrationFeatures(scRNAlist, nfeatures = 2000)

scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
    x <- RunPCA(x, features = scRNA.features, verbose = FALSE)})

scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = scRNA.features)

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,
                                        normalization.method = "SCT",
                                        anchor.features = scRNA.features,
                                        reduction = "rpca")

scRNA.combined.sct <- IntegrateData(scRNA.anchors, normalization.method = "SCT")

scRNA.combined.sct <- RunPCA(scRNA.combined.sct, npcs = 100, verbose = FALSE)

scRNA.combined.sct <- FindNeighbors(scRNA.combined.sct, dims = 1:60)
scRNA.combined.sct <- FindClusters(scRNA.combined.sct, resolution = 0.5)

#table(scRNA.combined.sct$source, scRNA.combined.sct$seurat_clusters)
#table(scRNA.combined.sct$cellType, scRNA.combined.sct$seurat_clusters)
scRNA.combined.sct <- RunUMAP(scRNA.combined.sct, dims = 1:50)
scRNA.combined.sct <- RunTSNE(scRNA.combined.sct, dims = 1:50)

# 保存数据
saveRDS(scRNA.combined.sct, file = "scRNA-SAM-Ath-Rice-combined.rds")

# 定义颜色
cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
mycolors <- colorRampPalette(cluster_cols,space="rgb")(length(levels(scRNA.combined.sct)))

p1 <- DimPlot(scRNA.combined.sct, reduction = "umap", label = F, cols = c('#65C4CB','#1276B8'), group.by = 'source', pt.size=0.05)
p2 <- DimPlot(scRNA.combined.sct, reduction = "umap", label = T, cols = mycolors, pt.size=0.05)
p <- p1 + p2
ggsave("scRNA_SCT_UMAP_Cluster.pdf", p, width = 14, height = 6)

###
p1 <- DimPlot(scRNA.combined.sct, reduction = "tsne", label = F, cols = c('#65C4CB','#1276B8'), group.by = 'source', pt.size=0.05)
p2 <- DimPlot(scRNA.combined.sct, reduction = "tsne", label = T, cols = mycolors, pt.size=0.05)
p <- p1 + p2
ggsave("scRNA_SCT_TSNE_Cluster.pdf", p, width = 14, height = 6)

library(Seurat)

scRNA <- readRDS("../02-Clustering/scRNA-SAM-Ath-Rice-combined.rds")

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

color_clusters <- c('#65C4CB','#1276B8')

p <- ggplot(melted, aes(x = Cluster, fill = Sample, stratum = Sample,
            alluvium = Sample, y = Count, label = Sample)) +
     ylab('Cell Proportion') + 
     geom_flow() + geom_stratum() +
     scale_fill_manual(values = color_clusters) +
     theme_pubr(base_size = 16, legend = "right", border = T) + rotate_x_text(0)

ggsave("source_seurat_clusters.pdf", p, width = 10, height = 4)

library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(cowplot)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 5)
rm(list = ls())
set.seed(1234)

combined <- readRDS("../02-Clustering/scRNA-SAM-Ath-Rice-combined.rds")

data <- as.matrix(table(combined$cellType, combined$seurat_clusters))
# 输出数据对行名和行序进行重新的排列
write.csv(data, 'scRNA-five-tissue-celltype-cluster.csv',quote=F)

# 作图
data_sorted <- read.csv('scRNA-five-tissue-celltype-cluster.csv', header=T, row.names=1, check.names=F)

library(pheatmap)
my_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(101)
#my_palette <- colorRampPalette(c("white", "firebrick3"))(101)
pheatmap(data_sorted, 
         cluster_row = F, cluster_col = F,
         show_rownames = T, show_colnames = T,
         color = my_palette, 
         scale = "row",
         cellwidth=11, cellheight=10,
         border_color='black',
         annotation_legend = T,
				 fontsize = 8, angle_col = '0',
         filename = 'scRNA-SAM-Cluster-Celltype.pdf',width = 10, height = 5)

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

scRNA <- readRDS("../02-Clustering/scRNA-SAM-Ath-Rice-combined.rds")

#data <- prop.table(table(scRNA$lib, scRNA$seurat_clusters), 2)
table(scRNA$cellType, scRNA$seurat_clusters)

0 <- cortex (RAM)
1 <- mesophyll cell (SAM)
2 <- mesophyll cell (SAM)
3 <- dividing cell (RAM)
4 <- exodermis cell (RAM)
5 <- meristem cell (RAM)
6 <- endodermis cell (RAM)
7 <- endodermis cell (RAM)
8 <- epidermis cell (SAM)
9 <- exodermis cell (RAM)
10 <- vascular cell (RAM)
11 <- endodermis cell (RAM)
12 <- epidermis/atrichoblast (RAM)
13 <- dividing cell (SAM)
14 <- epidermis/meristem (RAM)
15 <- epidermis/trichoblast (RAM)
16 <- undefined cell
17 <- meristem cell (RAM)
18 <- xylem cell (SAM/RAM)
19 <- vascular cell (SAM)
20 <- meristem cell (SAM/RAM)
21 <- undefined cell
22 <- phloem cell (SAM/RAM)
23 <- root cap (RAM)
24 <- procambium cell (SAM)
25 <- epidermis cell (SAM)



new.cluster.ids <- c("cortex (RAM)","mesophyll cell (SAM)","mesophyll cell (SAM)","dividing cell (RAM)",
"exodermis cell (RAM)","meristem cell (RAM)","endodermis cell (RAM)","endodermis cell (RAM)","epidermis cell (SAM)",
"exodermis cell (RAM)","vascular cell (RAM)","endodermis cell (RAM)","epidermis/atrichoblast (RAM)","dividing cell (SAM)",
"epidermis/meristem (RAM)","epidermis/trichoblast (RAM)","undefined cell","meristem cell (RAM)","xylem cell (SAM/RAM)",
"vascular cell (SAM)","meristem cell (SAM/RAM)","undefined cell","phloem cell (SAM/RAM)","root cap (RAM)",
"procambium cell (SAM)","epidermis cell (SAM)")

names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)

mycolors <- c("#EF3500","#FF5E47","#F7D098","#B7BCF3","#4E9FAA"
                 ,"#AA91AA","#957AB9","#844AA7","#405DAC","#00909B"
                 ,"#37C7CE","#477E63","#95B1D9","#71C89C","#50AE7A"
                 ,"#00ED94","#BDE29B","#D6D18D","#B49685","#6F7E52"
                 ,"#F0EFD0","#FF98A8","#D473BB","#DA4293","#db7093")
p <- DimPlot(scRNA, reduction = "tsne", label = F, cols = mycolors, pt.size=0.01)
ggsave("scRNA-Rice-SAM-RAM-annotation.pdf", p, width = 7, height = 4)

saveRDS(scRNA, file = "scRNA-Rice-SAM-RAM-annotation.rds")


