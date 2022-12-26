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

# 读取拟南芥原始数据
scRNA_Ath <- readRDS("../01-Ath-Integrated/03-annotation/scRNA-RAM-Ath-Harmony-Integrate-annotation.rds")
DefaultAssay(scRNA_Ath) <- "RNA"

# 读取水稻原始数据
scRNA_Rice <- readRDS("../../../05-Root-project/04-scRNA-integrated/03-annotation/scRNA-RAM-annotation.rds")
DefaultAssay(scRNA_Rice) <- "RNA"
scRNA_Rice$source <- scRNA_Rice$lib

# 获取原始数据中存在的基因name
df <- read.csv('Ont-to-One_orthologs.csv',header=T)
## 获取拟南芥与水稻同源基因后，存在于数据集的基因
gene_ath <- as.vector(CaseMatch(search=df$Atha, match=rownames(scRNA_Ath)))
df <- subset(df, df$Atha %in% gene_ath)
# 获取水稻基因后，存在于数据集的基因
gene_rice <- as.vector(CaseMatch(search=df$RAP, match=rownames(scRNA_Rice)))
# 确定最后的水稻和拟南芥基因  9174
df <- subset(df, df$RAP %in% gene_rice)

# 创建新的SeuratObject
scRNA_Ath <- subset(scRNA_Ath, features = df$Atha)

cellinfo <- subset(scRNA_Ath@meta.data, select = c("orig.ident","nCount_RNA","nFeature_RNA",'cellType','source'))
scRNA_Ath_ <- CreateSeuratObject(scRNA_Ath[["RNA"]]@counts, meta.data = cellinfo)
#scRNA_Ath_ <- CreateSeuratObject(counts = GetAssayData(scRNA_Ath, slot = "counts"))

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

scRNA_Ath_chname <- RenameGenesSeurat(scRNA_Ath_, df$RAP)

# 提取水稻
scRNA_Rice <- subset(scRNA_Rice, features = df$RAP)
cellinfo <- subset(scRNA_Rice@meta.data, select = c("orig.ident","nCount_RNA","nFeature_RNA",'cellType','source'))
scRNA_Rice_ <- CreateSeuratObject(scRNA_Rice[["RNA"]]@counts, meta.data = cellinfo)

# 保存数据
saveRDS(scRNA_Ath_chname, file = "scRNA-RAM-Ath.rds")

saveRDS(scRNA_Rice_, file = "scRNA-RAM-Rice.rds")


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
# -------------------## Rice ##----------------------
scRNA_Rice <- readRDS("../02-Preprocessing/scRNA-RAM-Rice.rds")
scRNA_Rice$Species <- "Rice"
scRNAlist_Rice <- SplitObject(scRNA_Rice, split.by = "orig.ident")

# -------------------## Ath ##----------------------
scRNA_Ath <- readRDS("../02-Preprocessing/scRNA-RAM-Ath.rds")
scRNA_Ath$Species <- "Ath"
scRNAlist_Ath <- SplitObject(scRNA_Ath, split.by = "orig.ident")
# 去掉拟南芥Denyer和Farmer数据，为了保持数据等量
scRNAlist <- c(scRNAlist_Rice, scRNAlist_Ath[c(1,2,3,4,11,12,13)])

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

scRNA.combined.sct <- FindNeighbors(scRNA.combined.sct, dims = 1:50)
scRNA.combined.sct <- FindClusters(scRNA.combined.sct, resolution = 0.3)

scRNA.combined.sct <- RunUMAP(scRNA.combined.sct, dims = 1:50, umap.method = 'umap-learn', metric = 'correlation')
scRNA.combined.sct <- RunTSNE(scRNA.combined.sct, dims = 1:50)

# 保存数据
saveRDS(scRNA.combined.sct, file = "scRNA-RAM-Ath-Rice-combined.rds")

# 定义颜色
cluster_cols <- c("#EF3500","#FF5E47","#F7D098","#B7BCF3","#4E9FAA"
                 ,"#AA91AA","#957AB9","#844AA7","#405DAC","#00909B"
                 ,"#37C7CE","#477E63","#95B1D9","#71C89C","#50AE7A"
                 ,"#00ED94","#BDE29B","#D6D18D","#B49685","#6F7E52"
                 ,"#F0EFD0","#FF98A8","#D473BB","#DA4293")

mycolors <- colorRampPalette(cluster_cols,space="rgb")(length(levels(scRNA.combined.sct)))
p1 <- DimPlot(scRNA.combined.sct, reduction = "umap", label = F, cols = c('#65C4CB','#1276B8'), group.by = 'Species', pt.size=0.05)
p2 <- DimPlot(scRNA.combined.sct, reduction = "umap", label = T, cols = mycolors, pt.size=0.05)
p <- p1 + p2
ggsave("scRNA_SCT_UMAP_Cluster.pdf", p, width = 14, height = 6)

p3 <- DimPlot(scRNA.combined.sct, reduction = "umap", group.by = "cellType", 
              split.by = "Species", ncol = 2, cols = cluster_cols, label = T)
ggsave("scRNA-Integrate-UMAP-CCA-SCT-Cluster-Split.pdf", p3, width = 14, height = 6)


mycolors <- colorRampPalette(cluster_cols,space="rgb")(length(levels(scRNA.combined.sct)))
p1 <- DimPlot(scRNA.combined.sct, reduction = "tsne", label = F, cols = c('#65C4CB','#1276B8'), group.by = 'Species', pt.size=0.05)
p2 <- DimPlot(scRNA.combined.sct, reduction = "tsne", label = T, cols = mycolors, pt.size=0.05)
p <- p1 + p2
ggsave("scRNA_SCT_TSNE_Cluster.pdf", p, width = 14, height = 6)

p3 <- DimPlot(scRNA.combined.sct, reduction = "tsne", group.by = "cellType", 
              split.by = "Species", ncol = 2, cols = cluster_cols, label = T)
ggsave("scRNA-Integrate-TSNE-CCA-SCT-Cluster-Split.pdf", p3, width = 16, height = 6)

library(Seurat)

scRNA <- readRDS("../03-Clustering/scRNA-RAM-Ath-Rice-combined.rds")

data <- prop.table(table(scRNA$Species, scRNA$seurat_clusters), 2)

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

scRNA <- readRDS("../03-Clustering/scRNA-RAM-Ath-Rice-combined.rds")
scRNA_Rice <- subset(x = scRNA, subset = Species == "Rice")
scRNA_Ath <- subset(x = scRNA, subset = Species == "Ath")

# 开始水稻作图
Idents(scRNA_Rice) <- "cellType"
mycolors <- c("#FC6E52","#E2C9B3","#86ADD2","#8596AA","#854DA8","#3B5EAB","#0B9BA4",
              "#56B37F","#06E993","#CDC28B","#988B6F","#F9B1B1","#D94292")

p <- DimPlot(scRNA_Rice, reduction = "tsne", label = F, cols = mycolors, pt.size=0.01)
ggsave("scRNA-Rice-annotation.pdf", p, width = 7, height = 4)

# 开始拟南芥作图
Idents(scRNA_Ath) <- "cellType"
mycolors <- c("#854DA8","#DA79B8","#FC6E52","#D94292","#CDC28B","#56B37F","#0B9BA4","#F9B1B1","#86ADD2","#E2C9B3")

p <- DimPlot(scRNA_Ath, reduction = "tsne", label = F, cols = mycolors, pt.size=0.01)
ggsave("scRNA-Ath-annotation.pdf", p, width = 7, height = 4)

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

scRNA <- readRDS("../03-Clustering/scRNA-RAM-Ath-Rice-combined.rds")

#data <- prop.table(table(scRNA$lib, scRNA$seurat_clusters), 2)
table(scRNA$cellType, scRNA$seurat_clusters)

0 <- vascular cylinder (Atha)
1 <- root cap (Atha)
2 <- exodermis (Osat)/ endodermis (Atha)
3 <- endodermis (Osat)
4 <- cortex (Osat)
5 <- endodermis (Osat)
6 <- meristem (Atha)
7 <- meristem (Osat)
8 <- phloem (Atha/Osat)
9 <- meristem (Atha/Osat)
10 <- vascular cylinder (Atha)
11 <- epidermis/atrichoblast (Atha)
12 <- epidermis/trichoblast (Atha/Osat)
13 <- epidermis/atrichoblast (Osat)
14 <- vascular cylinder (Osat)
15 <- cortex (Atha/Osat)
16 <- undefined cell (Atha)
17 <- xylem (Atha/Osat)
18 <- endodermis (Atha)
19 <- root cap (Atha)
20 <- endodermis (Atha/Osat)

new.cluster.ids <- c("vascular cylinder","root cap","exodermis","endodermis","cortex",
                     "endodermis","meristem","meristem","phloem","meristem","vascular cylinder",
                     "epidermis/atrichoblast","epidermis/trichoblast","epidermis/atrichoblast",
                     "vascular cylinder","cortex","undefined cell","xylem","endodermis","root cap","endodermis")

names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)

mycolors <- c("#7372B8","#E66470","#83C875","#C4C26F","#8666A9",
              "#7AD3B5","#7AA3B7","#FF9822","#CD63AB","#C4A0C4",
              "#b15928")
p <- DimPlot(scRNA, reduction = "tsne", label = F, cols = mycolors, pt.size=0.01)
ggsave("scRNA-RAM-Ath-Rice-annotation.pdf", p, width = 6, height = 4)

## 差异基因分析
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

scRNA <- readRDS("../03-Clustering/scRNA-RAM-Ath-Rice-combined.rds")

scRNA$celltype.species <- paste(scRNA$seurat_clusters, scRNA$Species, sep = "_")
Idents(scRNA) <- "celltype.species"

scRNA <- PrepSCTFindMarkers(scRNA)

## 木质部（Cluster 17）差异基因分析
# 水稻
riec.xylem.response <- FindMarkers(scRNA, assay = "SCT", ident.1 = "17_Rice", ident.2 = "17_Ath", only.pos = TRUE)
riec.xylem.response <- subset(riec.xylem.response, subset = p_val_adj < 0.05)
write.table(riec.xylem.response, 'DEGs-xylem-rice-DEGs.txt', sep = '\t', quote=F)
# 拟南芥
Ath.xylem.response <- FindMarkers(scRNA, assay = "SCT", ident.1 = "17_Ath", ident.2 = "17_Rice", only.pos = TRUE)
Ath.xylem.response <- subset(Ath.xylem.response, subset = p_val_adj < 0.05)
write.table(Ath.xylem.response, 'DEGs-xylem-ath-DEGs.txt', sep = '\t', quote=F)
# Identify conserved cell type markers  57
Idents(scRNA) <- "seurat_clusters"
xylem.markers <- FindConservedMarkers(scRNA, assay = "SCT", ident.1 = "17", grouping.var = "Species")
write.table(xylem.markers, 'DEGs-xylem-conserved.txt', sep='\t', quote=F)
## 求平均表达量
scRNA_ <- subset(scRNA, idents = '17')
features <- c(rownames(xylem.markers), rownames(riec.xylem.response), rownames(Ath.xylem.response))
avg <- AverageExpression(scRNA_, assay = "SCT", features = unique(features), group.by = "Species")$SCT

## 保存数据
write.csv(avg, 'DEGs-xylem-heatmap.csv', quote=F)

# 作图
library(pheatmap)
library(RColorBrewer)
my_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(11)
my_palette <- colorRampPalette(c("white", "firebrick3"))(101)
my_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdPu")))(100)
my_palette <- colorRampPalette(brewer.pal(n = 7, name = "RdPu"))(100)

# avg <- read.csv("DEGs-xylem-heatmap.csv",header=T, row.names=1)

avg <- subset(as.data.frame(avg), subset = Ath<5 & Rice<5)

pheatmap(log10(avg+1), cluster_row = F, cluster_col = F,
         show_rownames = F, show_colnames = F,
         #scale = 'column',
         cellwidth = 40,
         color = my_palette, 
         border_color = NA,
         filename = 'DEGs-xylem-heatmap.pdf', width = 4, height = 5)






