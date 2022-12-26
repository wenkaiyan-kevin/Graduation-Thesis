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

new.cluster.ids <- c("mesophyll cell (Osat)", "meristem epidermis cell (Osat)","vascular cell (Osat)",
"shoot meristematic cell (Osat)","proliferating cell (Osat)","leaf epidermis cell (Osat)", "undefined cell (Osat)")
names(new.cluster.ids) <- levels(scRNA_Rice_SAM)
scRNA_Rice_SAM <- RenameIdents(scRNA_Rice_SAM, new.cluster.ids)
# 保存数据
saveRDS(scRNA_Rice_SAM, file = "scRNA-Rice-SAM-annotation.rds")

# 读取水稻叶片
scRNA_Rice_Leaf <- readRDS("/gss1/home/yanwk/ywk_SingleCell_work/06_leaf_scRNA_Analysis/04_Clusters_Annotation/scRNA-Leaf-annotation.rds")

new.cluster.ids <- c("mesophyll cell (Osat)","leaf epidermis cell (Osat)","procambium cell (Osat)","epidermal initial cell (Osat)",
"Mesophyll precursor (Osat)","large parenchyma/MO (Osat)","vascular initial cell (Osat)","preprocambium cell (Osat)",
"phloem cell (Osat)","mesophyll initial cell (Osat)","large parenchyma/PO (Osat)","mestome sheath cell (Osat)","xylem cell (Osat)")
names(new.cluster.ids) <- levels(scRNA_Rice_Leaf)
scRNA_Rice_Leaf <- RenameIdents(scRNA_Rice_Leaf, new.cluster.ids)
# 保存数据
saveRDS(scRNA_Rice_Leaf, file = "scRNA-Rice-Leaf-annotation.rds")

# 读取拟南芥茎尖+叶片
scRNA_Ath_SAM <- readRDS("/gss1/home/yanwk/ywk_Graduation_Project/02-Clustering/01-Ath/02-Ath-Shoot/05-Annotation/scRNA-SAM-Ath-Harmony-Integrate-Annotation.rds")

new.cluster.ids <- c("mesophyll cell (Atha)","shoot meristematic cell (Atha)","epidermis cell (Atha)",
"vascular cell (Atha)","proliferating cell (Atha)","guard cell (Atha)","companion cell (Atha)","undefined cell (Atha)")
names(new.cluster.ids) <- levels(scRNA_Ath_SAM)
scRNA_Ath_SAM <- RenameIdents(scRNA_Ath_SAM, new.cluster.ids)
# 保存数据
saveRDS(scRNA_Ath_SAM, file = "scRNA-Ath-SAM-annotation.rds")

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
scRNA_Rice_SAM <- readRDS("../01-data/scRNA-Rice-SAM-annotation.rds")
DefaultAssay(scRNA_Rice_SAM) <- "RNA"
scRNA_Rice_SAM$source <- "rice"
scRNA_Rice_SAM$cellType <- Idents(scRNA_Rice_SAM)

# 读取水稻叶片
scRNA_Rice_Leaf <- readRDS("../01-data/scRNA-Rice-Leaf-annotation.rds")
DefaultAssay(scRNA_Rice_Leaf) <- "RNA"
scRNA_Rice_Leaf$source <- "rice"
scRNA_Rice_Leaf$cellType <- Idents(scRNA_Rice_Leaf)

# 读取拟南芥数据
scRNA_Ath_SAM <- readRDS("../01-data/scRNA-Ath-SAM-annotation.rds")
DefaultAssay(scRNA_Ath_SAM) <- "RNA"
scRNA_Ath_SAM$source <- "ath"
scRNA_Ath_SAM$cellType <- Idents(scRNA_Ath_SAM)

# 获取原始数据中存在的基因name
df <- read.csv('One-to-One_orthologs.csv',header=T)
## 获取三个数据集的共同的基因集
scRNA_Rice_SAM_genes <- as.vector(CaseMatch(search=df$RAP, match=rownames(scRNA_Rice_SAM)))
scRNA_Rice_Leaf_genes <- as.vector(CaseMatch(search=df$RAP, match=rownames(scRNA_Rice_Leaf)))
scRNA_Ath_SAM_genes <- as.vector(CaseMatch(search=df$Atha, match=rownames(scRNA_Ath_SAM)))
scRNA_Ath_SAM_genes_ <- subset(df, df$Atha %in% scRNA_Ath_SAM_genes)$RAP

genes <- intersect(intersect(scRNA_Rice_SAM_genes, scRNA_Rice_Leaf_genes), scRNA_Ath_SAM_genes_) #9060
df <- subset(df, df$RAP %in% genes)
gene_rice <- df$RAP
gene_ath  <- df$Ath

# 创建新的SeuratObject
scRNA_Ath_SAM <- subset(scRNA_Ath_SAM, features = gene_ath)

cellinfo <- subset(scRNA_Ath_SAM@meta.data, select = c("orig.ident","nCount_RNA","nFeature_RNA",'cellType','source'))
scRNA_Ath_SAM_ <- CreateSeuratObject(scRNA_Ath_SAM[["RNA"]]@counts, meta.data = cellinfo)

# 修改名字
RenameGenesSeurat <- function(obj, newnames){
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

scRNA_Ath_SAM_chname <- RenameGenesSeurat(scRNA_Ath_SAM_, gene_rice)

# 提取水稻茎尖
scRNA_Rice_SAM <- subset(scRNA_Rice_SAM, features = gene_rice)
cellinfo <- subset(scRNA_Rice_SAM@meta.data, select = c("orig.ident","nCount_RNA","nFeature_RNA",'cellType','source'))
scRNA_Rice_SAM_ <- CreateSeuratObject(scRNA_Rice_SAM[["RNA"]]@counts, meta.data = cellinfo)

# 提取水稻叶片
scRNA_Rice_Leaf <- subset(scRNA_Rice_Leaf, features = gene_rice)
cellinfo <- subset(scRNA_Rice_Leaf@meta.data, select = c("orig.ident","nCount_RNA","nFeature_RNA",'cellType','source'))
scRNA_Rice_Leaf_ <- CreateSeuratObject(scRNA_Rice_Leaf[["RNA"]]@counts, meta.data = cellinfo)


# 保存数据
saveRDS(scRNA_Rice_SAM_, file = "scRNA-Rice-SAM.rds")
saveRDS(scRNA_Rice_Leaf_, file = "scRNA-Rice-Leaf.rds")
saveRDS(scRNA_Ath_SAM_chname, file = "scRNA-Ath-SAM.rds")

## 整合
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

#1.1. 读取数据
# -------------------## Rice SAM##----------------------
scRNA_Rice_SAM <- readRDS("../02-Preprocessing/scRNA-Rice-SAM.rds")
Idents(scRNA_Rice_SAM) <- 'cellType'
scRNA_Rice_SAM <- subset(scRNA_Rice_SAM, idents = "undefined cell (Osat)", invert = TRUE)

# -------------------## Rice SAM##----------------------
scRNA_Rice_Leaf <- readRDS("../02-Preprocessing//scRNA-Rice-Leaf.rds")
scRNAlist_Rice_Leaf <- SplitObject(scRNA_Rice_Leaf, split.by = "orig.ident")

# -------------------## Ath ##----------------------
scRNA_Ath <- readRDS("../02-Preprocessing/scRNA-Ath-SAM.rds")
Idents(scRNA_Ath) <- 'cellType'
scRNA_Ath <- subset(scRNA_Ath, idents = "undefined cell (Atha)", invert = TRUE)
scRNA_Ath <- subset(x = scRNA_Ath, downsample = 5000)
scRNAlist_Ath <- SplitObject(scRNA_Ath, split.by = "orig.ident")

scRNAlist <- c(scRNA_Rice_SAM, scRNAlist_Rice_Leaf, scRNAlist_Ath)

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

saveRDS(scRNA.combined.sct, file = "scRNA-SAM-Ath-Rice-combined-temp.rds")


scRNA.combined.sct <- FindNeighbors(scRNA.combined.sct, dims = 1:30)
scRNA.combined.sct <- FindClusters(scRNA.combined.sct, resolution = 0.5)

table(Idents(scRNA.combined.sct))
table(scRNA.combined.sct$source, scRNA.combined.sct$seurat_clusters)
table(scRNA.combined.sct$cellType, scRNA.combined.sct$seurat_clusters)



scRNA.combined.sct <- RunUMAP(scRNA.combined.sct, dims = 1:30)
scRNA.combined.sct <- RunTSNE(scRNA.combined.sct, dims = 1:30)

# 保存数据
saveRDS(scRNA.combined.sct, file = "scRNA-SAM-Ath-Rice-combined.rds")

# 定义颜色
cluster_cols <- c("#EF3500","#FF5E47","#F7D098","#B7BCF3","#4E9FAA"
                 ,"#AA91AA","#957AB9","#844AA7","#405DAC","#00909B"
                 ,"#37C7CE","#477E63","#95B1D9","#71C89C","#50AE7A"
                 ,"#00ED94","#BDE29B","#D6D18D","#B49685","#6F7E52"
                 ,"#F0EFD0","#FF98A8","#D473BB","#DA4293","#db7093")

mycolors <- colorRampPalette(cluster_cols,space="rgb")(length(levels(scRNA.combined.sct)))
p1 <- DimPlot(scRNA.combined.sct, reduction = "umap", label = F, cols = c('#65C4CB','#1276B8'), group.by = 'source', pt.size=0.05)
p2 <- DimPlot(scRNA.combined.sct, reduction = "umap", label = T, cols = mycolors, pt.size=0.05)
p <- p1 + p2
ggsave("scRNA_SCT_UMAP_Cluster.pdf", p, width = 14, height = 6)

mycolors <- colorRampPalette(cluster_cols,space="rgb")(26)
p3 <- DimPlot(scRNA.combined.sct, reduction = "umap", group.by = "cellType", 
              split.by = "source", ncol = 2, cols = mycolors, label = T, repel = T)
ggsave("scRNA-Integrate-UMAP-CCA-SCT-Cluster-Split.pdf", p3, width = 14, height = 6)


mycolors <- colorRampPalette(cluster_cols,space="rgb")(length(levels(scRNA.combined.sct)))
p1 <- DimPlot(scRNA.combined.sct, reduction = "tsne", label = F, cols = c('#65C4CB','#1276B8'), group.by = 'source', pt.size=0.05)
p2 <- DimPlot(scRNA.combined.sct, reduction = "tsne", label = T, cols = mycolors, pt.size=0.05)
p <- p1 + p2
ggsave("scRNA_SCT_TSNE_Cluster.pdf", p, width = 14, height = 6)

mycolors <- colorRampPalette(cluster_cols,space="rgb")(26)
p3 <- DimPlot(scRNA.combined.sct, reduction = "tsne", group.by = "cellType", 
              split.by = "source", ncol = 2, cols = mycolors, label = T)
ggsave("scRNA-Integrate-TSNE-CCA-SCT-Cluster-Split.pdf", p3, width = 16, height = 6)

library(Seurat)

scRNA <- readRDS("../03-Clustering/scRNA-SAM-Ath-Rice-combined.rds")

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

ggsave("source_seurat_clusters.pdf", p, width = 8, height = 4)

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

scRNA <- readRDS("../03-Clustering/scRNA-SAM-Ath-Rice-combined.rds")

#data <- prop.table(table(scRNA$lib, scRNA$seurat_clusters), 2)
table(scRNA$cellType, scRNA$seurat_clusters)

0 <- epidermis cell (Atha/Osat)
1 <- shoot meristematic cell (Atha)
2 <- mesophyll cell (Osat)
3 <- mesophyll cell (Atha)
4 <- mesophyll cell (Atha/Osat)
5 <- proliferating cell (Atha/Osat)
6 <- epidermis cell (Osat)
7 <- vascular cell (Atha/Osat)
8 <- shoot meristematic cell (Atha/Osat)
9 <- mesophyll cell (Atha)
10 <- vascular cell (Atha)
11 <- epidermis cell (Atha)
12 <- mesophyll cell (Osat)
13 <- procambium cell (Osat)
14 <- vascular cell (Atha)
15 <- proliferating cell (Atha)
16 <- epidermis cell (Atha)
17 <- mesophyll cell (Osat)
18 <- companion cell (Atha)
19 <- guard cell (Atha)
20 <- vascular cell (Osat)
21 <- epidermis cell (Osat)

new.cluster.ids <- c("epidermis cell","shoot meristematic cell","mesophyll cell","mesophyll cell","mesophyll cell",
"proliferating cell","epidermis cell","vascular cell","shoot meristematic cell","mesophyll cell","vascular cell",
"epidermis cell","mesophyll cell","procambium cell","vascular cell","proliferating cell","epidermis cell","mesophyll cell",
"companion cell","guard cell","vascular cell","epidermis cell")

names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)

mycolors <- c("#7372B8","#E66470","#83C875","#C4C26F","#8666A9",
              "#7AD3B5","#7AA3B7","#FF9822","#CD63AB","#C4A0C4",
              "#b15928")
p <- DimPlot(scRNA, reduction = "tsne", label = F, cols = mycolors, pt.size=0.01)
ggsave("scRNA-SAM-Ath-Rice-annotation.pdf", p, width = 6.5, height = 4)

saveRDS(scRNA, file = "scRNA-SAM-Ath-Rice-combined-annotation.rds")


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