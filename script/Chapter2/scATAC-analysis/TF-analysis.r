### 特异性Peak数量作图

library(ggplot2)

data <- read.table('uniq-peaks-num.txt', header=T)
#data$Cluster <- factor(data$Cluster, levels=c('Cluster4','Cluster3','Cluster2','Cluster1','Cluster0'))

p <- ggplot(data, aes(x=Cluster, y=Number)) +
            geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
            coord_flip() +
            xlab("") +
            theme_classic()
ggsave("uniq-peak-cluster.pdf", plot = p, width = 5, height = 3)

### 提取区间
library(Seurat)
library(Signac)
library(rtracklayer)
library(tidyverse)
library(patchwork)
library(dplyr)
library(BSgenome.Osativa.IRGSP.IRGSP1)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 4)
set.seed(1234)

Peaks <- readRDS("../../02-DA-analysis/Peaks-cluster-special.rds")

endodermis <- as.data.frame(subset(Peaks, peak_called_in == "endodermis"))[,1:3]
write.table(endodermis, "endodermis.bed", sep='\t', row.names=F, col.names=F, quote=F)

epidermis_trichoblast <- as.data.frame(subset(Peaks, peak_called_in == "epidermis_trichoblast"))[,1:3]
write.table(epidermis_trichoblast, "epidermis_trichoblast.bed", sep='\t', row.names=F, col.names=F, quote=F)

exodermis <- as.data.frame(subset(Peaks, peak_called_in == "exodermis"))[,1:3]
write.table(exodermis, "exodermis.bed", sep='\t', row.names=F, col.names=F, quote=F)

meristem <- as.data.frame(subset(Peaks, peak_called_in == "meristem"))[,1:3]
write.table(meristem, "meristem.bed", sep='\t', row.names=F, col.names=F, quote=F)

xylem <- as.data.frame(subset(Peaks, peak_called_in == "xylem"))[,1:3]
write.table(xylem, "xylem.bed", sep='\t', row.names=F, col.names=F, quote=F)

### motif 分析作图
## 表皮
library(ggplot2)

data <- read.table('data.txt', header=T)

epidermis <- subset(data, cellType == "epidermis")

epidermis$Best_Match <- factor(epidermis$Best_Match, levels=rev(epidermis$Best_Match))

p <- ggplot(epidermis, aes(x=Best_Match, y=-log10(FDR))) +
            geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.6) +
            ylim(1, 2.2) +
            coord_flip() + scale_y_continuous(expand = c(0, 0)) +
            xlab("") +
            theme_classic()
ggsave("epidermis-motif.pdf", plot = p, width = 5, height = 3)
