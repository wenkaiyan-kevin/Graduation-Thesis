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

int.data <- readRDS('../02-Intergrated/atac-rna_integrated_liger.rds')

scRNA <- subset(x = int.data, subset = orig.ident == "rna")

# 数值结果
data <- as.data.frame(table(scRNA$cellType, scRNA$clusters))
write.csv(data, "cellType_clusters_value.csv", row.names=F, quote=F)

# 比例结果
data <- prop.table(table(scRNA$cellType, scRNA$clusters), 2)
write.csv(data, "cellType_clusters_propration.csv", row.names=T, quote=F)

library(reshape2)
library(dplyr)
library(ggpubr)
library(ggalluvial)

######
data <- read.csv('cellType_clusters_propration.csv',header = T, check.names = F)

melted <- melt(data) %>% rename(Count = value, Cluster = variable) %>% collect()

color_clusters <- c("#FC6E52","#E2C9B3","#86ADD2","#8596AA","#854DA8","#3B5EAB","#0B9BA4",
                    "#56B37F","#06E993","#CDC28B","#988B6F","#F9B1B1","#D94292")

p <- ggplot(melted, aes(x = Cluster, fill = Sample, stratum = Sample,
            alluvium = Sample, y = Count, label = Sample)) +
     ylab('Proportion') + 
     geom_flow() + geom_stratum() +
     scale_fill_manual(values = color_clusters) +
     theme_pubr(base_size = 16, legend = "right", border = T) + rotate_x_text(0)

ggsave("cellType_clusters_propration.pdf", p, width = 12, height = 4)

###########
data <- read.csv('cellType_clusters_value.csv',header = T, check.names = F)
p <- ggplot(data, aes(x=Var2, y=Freq, fill = Var1))+
            geom_bar(position = "stack", stat = "identity", width = 0.4) +
            theme_bw()+   
            scale_fill_manual(values=color_clusters) +
            scale_y_continuous(expand = c(0,0)) 

ggsave("cellType_clusters_value.pdf", p, width = 12, height = 4)

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


int.data <- readRDS('../02-Intergrated/atac-rna_integrated_liger.rds')
Idents(int.data) <- 'clusters'

1 <- xylem
2 <- epidermis/trichoblast
3 <- meristem
4 <- vascular cylinder
5 <- epidermis/trichoblast
6 <- epidermis/meristem
7 <- cortex
8 <- epidermis/atrichoblast
9 <- exodermis
10 <- pericycle
11 <- undefined
12 <- endodermis
13 <- endodermis
14 <- meristem
15 <- root cap
16 <- undefined
17 <- exodermis
18 <- endodermis
19 <- phloem
20 <- meristem

new.cluster.ids <- c("xylem","epidermis/trichoblast","meristem","vascular cylinder","epidermis/trichoblast",
                     "epidermis/meristem","cortex","epidermis/atrichoblast","exodermis","pericycle","undefined",
                     "endodermis","endodermis","meristem","root cap","undefined","exodermis","endodermis",
                     "phloem","meristem")

names(new.cluster.ids) <- levels(int.data)
int.data <- RenameIdents(int.data, new.cluster.ids)

mycolors <- c("#FC6E52","#E2C9B3","#86ADD2","#8596AA","#854DA8","#3B5EAB","#0B9BA4",
              "#56B37F","#06E993","#CDC28B","#988B6F","#F9B1B1","#D94292")

p <- DimPlot(int.data, reduction = "umap", label = F, cols = mycolors, pt.size=0.01)
ggsave("scRNA-RAM-atac-rna-integrated-annotation.pdf", p, width = 7, height = 4)


int.data$celltype <- Idents(int.data)
p <- DimPlot(int.data, group.by = "celltype", split.by = "orig.ident",
             cols = mycolors, pt.size = 0.01)
ggsave("atac-rna_integrated_liger-split-annotation.pdf", p, width = 9, height = 4)

# 保存数据
saveRDS(int.data, file = "atac-rna_integrated_liger-annotation.rds")