library(UCell)
library(irGSEA)
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

# remotes::install_github("carmonalab/UCell", ref="v1.3")

# https://bioconductor.org/packages/release/bioc/src/contrib/AUCell_1.18.1.tar.gz
# 安装AUCell 最新版本
# install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/AUCell_1.18.1.tar.gz", repos=NULL, type="source")
scRNA <- readRDS("../../03-annotation/scRNA-RAM-annotation.rds")

data <- read.table('../01-gene-list/input-gene-list.txt', header=T)
# 构建需要输入的gene sets
geneSets <- list()
geneSets$Auxin_biosynthesis <- as.vector(CaseMatch(search=subset(data, Type == "Auxin_biosynthesis")[,2],
                                                   match=rownames(scRNA)))
geneSets$Response_to_auxin <- as.vector(CaseMatch(search=subset(data, Type == "Response_to_auxin")[,2],
                                                   match=rownames(scRNA)))
geneSets$Jasmonate_biosynthesis <- as.vector(CaseMatch(search=subset(data, Type == "Jasmonate_biosynthesis")[,2],
                                                   match=rownames(scRNA)))
geneSets$Response_to_jasmonate <- as.vector(CaseMatch(search=subset(data, Type == "Response_to_jasmonate")[,2],
                                                   match=rownames(scRNA)))
geneSets$Salicylic_acid_biosynthesis <- as.vector(CaseMatch(search=subset(data, Type == "Salicylic_acid_biosynthesis")[,2],
                                                   match=rownames(scRNA)))
geneSets$Response_to_salicylic_acid <- as.vector(CaseMatch(search=subset(data, Type == "Response_to_salicylic_acid")[,2],
                                                   match=rownames(scRNA)))
geneSets$ABA_biosynthesis <- as.vector(CaseMatch(search=subset(data, Type == "ABA_biosynthesis")[,2],
                                                   match=rownames(scRNA)))
geneSets$Response_to_ABA <- as.vector(CaseMatch(search=subset(data, Type == "Response_to_ABA")[,2],
                                                   match=rownames(scRNA)))
geneSets$Ethylene_biosynthesis <- as.vector(CaseMatch(search=subset(data, Type == "Ethylene_biosynthesis")[,2],
                                                   match=rownames(scRNA)))
geneSets$Response_to_ethylene <- as.vector(CaseMatch(search=subset(data, Type == "Response_to_ethylene")[,2],
                                                   match=rownames(scRNA)))
geneSets$GA_biosynthesis <- as.vector(CaseMatch(search=subset(data, Type == "GA_biosynthesis")[,2],
                                                   match=rownames(scRNA)))
geneSets$Response_to_GA <- as.vector(CaseMatch(search=subset(data, Type == "Response_to_GA")[,2],
                                                   match=rownames(scRNA)))
geneSets$BR_biosynthesis <- as.vector(CaseMatch(search=subset(data, Type == "BR_biosynthesis")[,2],
                                                   match=rownames(scRNA)))
geneSets$Response_to_BR <- as.vector(CaseMatch(search=subset(data, Type == "Response_to_BR")[,2],
                                                   match=rownames(scRNA)))

scRNA <- irGSEA.score(object = scRNA, assay = "SCT", 
                      slot = "data", seeds = 123, ncores = 10,
                      min.cells = 3, min.feature = 0,
                      custom = T, geneset = geneSets,
                      method = c("UCell", "ssgsea"),
                      aucell.MaxRank = NULL, ucell.MaxRank = NULL, kcdf = 'Gaussian')

# scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
#my_palette <- colorRampPalette(c("lightgray", "#FF7270", "#FF3838"))(101)
#my_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(101)
#scatterplot <- irGSEA.density.scatterplot(object = scRNA, method = "UCell",
#                             size = 0.5, pal = 'magma',
#                             show.geneset = "Auxin-biosynthesis",
#                             reduction = "umap") + 
#               scale_colour_gradientn(colours = my_palette)
#ggsave("test.pdf", scatterplot, width = 7, height = 6)
my_palette <- colorRampPalette(c("lightgray",brewer.pal(n = 9, name = "OrRd")))(101)
scatterplot <- irGSEA.density.scatterplot(object = scRNA, method = "ssgsea",
                             size = 0.5,
                             show.geneset = "Response-to-BR",
                             reduction = "umap") +
                             scale_colour_gradientn(colours = my_palette)
ggsave("Response-to-BR-scatterplot.pdf", scatterplot, width = 7, height = 6)

densityheatmap <- irGSEA.densityheatmap(object = scRNA,
                                        method = "ssgsea",
                                        show.geneset = "Response-to-BR")
ggsave("Response-to-BR-densityheatmap.pdf", densityheatmap, width = 7, height = 6)                                                   