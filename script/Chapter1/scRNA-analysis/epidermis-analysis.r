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
p <- DimPlot(scRNA, reduction = "umap",
             cells.highlight = CellsByIdentities(object=scRNA, 
             idents=c("meristerm epidermis cell","leaf epidermal cell")))
mycolors <- c("grey50","#FC8D62","grey50","grey50","grey50","#E5C494","grey50")
p <- DimPlot(scRNA, reduction = "umap", cols=mycolors, pt.size=1)
ggsave('scRNA-EC-highlight.pdf', p, width = 8, height = 5)

####
scRNA_EC <- subset(scRNA, idents = c("meristerm epidermis cell", "leaf epidermal cell"))

DEGs <- FindMarkers(scRNA_EC, ident.1 = "meristerm epidermis cell", ident.2 = "leaf epidermal cell", min.pct = 0.5)
DEGs <- subset(DEGs, subset = p_val_adj < 0.05)

write.csv(DEGs,'scRNA-EC-DEGs.csv',quote=F)

# ------------------------------GO富集分析作图---------------------------------

library(ggplot2)

up <- read.csv('scRNA-EC-DEGs-UP-GO-Analysis.csv', header=T)
down <- read.csv('scRNA-EC-DEGs-DOWN-GO-Analysis.csv', header=T)

mycolors <- c("#FC8D62","#E5C494")

mytheme <-theme_bw() +  
          theme(plot.margin = unit(c(1,1,1,1),"cm"),
                plot.title = element_text(size = rel(7),hjust = 0.5),
                panel.grid.minor = element_blank(),
		            panel.grid.major = element_blank(),
		            panel.border = element_blank(),
		            panel.background = element_blank(),
                axis.line = element_line(color = 'black',size = 1),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.x = element_line(size=2),
								axis.ticks.length=unit(.55, "cm"),
								axis.text = element_text(size=rel(2)),
								axis.title = element_text(size = rel(5)))
#挑选感兴趣的iterm
p1 <- ggplot(up[c(1,4,5,14,21),], aes(y =-log10(fdr), x=reorder(Term,-log10(fdr)))) +
      geom_bar(stat="identity", fill='#FC8D62', width=0.9, alpha=0.7)+
      coord_flip() + 
      geom_text(aes(y=0,label = Term), hjust=0, vjust=0.5, size = 30) +
      ggtitle('DEGs UP GO Annotation') +
      ylab('-log10(padj)') + xlab('') +
      mytheme
p2 <- ggplot(down[c(1,2,7,9,11),], aes(y =-log10(fdr), x=reorder(Term,-log10(fdr)))) +
      geom_bar(stat="identity", fill='#E5C494', width=0.9, alpha=0.7)+
      coord_flip() + 
      geom_text(aes(y=0,label = Term), hjust=0, vjust=0.5, size = 30) +
      ggtitle('DEGs DOWN GO Annotation') +
      ylab('-log10(padj)') + xlab('') +
      mytheme


library(patchwork)

p <- p1 / p2

ggsave('scRNA-EC-DEGs-GO.pdf', p, width = 40, height = 40, limitsize = FALSE)