

library(Seurat)

# 提取marker基因对象
scRNA <- readRDS("../02-SMC-Annotation/scRNA-SAM-SMC-annotation.rds")
Features <- read.table('scRNA-SAM-SMC-maker.txt',header=F)[,1]
scRNA <- subset(scRNA, features=Features)
# 计算平均表达
avg_exp <- as.data.frame(AverageExpression(scRNA)$RNA)
rowname <- read.table('scRNA-SAM-SMC-maker-SymId-TF.txt',header=F)[,2]
rownames(avg_exp) <- rowname
colnames(avg_exp) <- c('PZ','OC','CZ')
write.csv(avg_exp, 'avg_exp.csv', quote=F)

#-------------------------------- 作图 ----------------------------------------

library(ComplexHeatmap)
library(RColorBrewer)

col_fun <- circlize::colorRamp2(c(-4, 0, 4), c("navy", "white", "firebrick3"))

data <- read.csv('avg_exp.csv',header=T,row.names=1)

ha <- rowAnnotation(foo = anno_mark(at = c(7,16,131,136,150,180,189,193,199,219,220,239,244,250,252,261,271,272,304,323,341,456), 
                          labels = rownames(data)[c(7,16,131,136,150,180,189,193,199,219,220,239,244,250,252,261,271,272,304,323,341,456)],
                          extend = unit(0, "mm")),
                    gap = unit(2, "points"))

pdf(file="scRNA-SAM-SMC-marker-heatmap.pdf",onefile = FALSE,width=4, height=9)
p <- pheatmap(data, 
         cluster_row = T, cluster_col = F,
         show_rownames = F, show_colnames = T,
         scale = 'row',
         cutree_rows = 3,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "PRGn")))(100),
         right_annotation = ha,
         heatmap_legend_param = list(title='Z-Score'),
				 fontsize_row = 1,
         angle_col = '0')

draw(p, heatmap_legend_side = "left")

dev.off()

#-------------------------------- 调控网络分析 ----------------------------------------


library(Seurat)
library(GENIE3)
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

# 提取marker基因对象
scRNA <- readRDS("../../02-SMC-Annotation/scRNA-SAM-SMC-annotation.rds")
Features <- read.table('../scRNA-SAM-SMC-maker.txt',header=F)[,1]
scRNA <- subset(scRNA, features=Features)

exprMatr <- as.matrix(scRNA@assays$RNA@data)

rowname <- read.table('../scRNA-SAM-SMC-maker-SymId-TF.txt',header=F)[,2]
rownames(exprMatr) <- rowname


set.seed(1314)
regulators <- c(7,16,131,136,150,180,189,193,199,219,220,239,244,250,252,261,271,272,304,323,341,456)

weightMat <- GENIE3(exprMatr, regulators=regulators)

linkList <- getLinkList(weightMat)
> summary(linkList$weight)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.01518 0.03396 0.04072 0.04555 0.04926 0.29099

# 过滤
linkList <- getLinkList(weightMat， threshold=0.04926)