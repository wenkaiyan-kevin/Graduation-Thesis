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

scRNA <- readRDS("../../02-Clustering/02-Rice/03-SAM/02-Cluster/scRNA-SAM-Cluster.rds")
DefaultAssay(scRNA) <- 'RNA'

scRNA <- NormalizeData(scRNA, verbose = FALSE)

scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25,
                                test.use = "wilcox")


scRNA.markers <- subset(scRNA.markers, subset = pct.2 < 0.5 & p_val_adj < 0.05)

#scRNA.markers <- subset(scRNA.markers, subset = pct.1 > 0.5 & pct.2 < 0.5 & p_val_adj < 0.05)

> table(scRNA.markers$cluster)
	0   1   2   3   4   5   6
123 231 213 461  62 279   9

write.csv(scRNA.markers,'scRNA-SAM-marker.csv',quote=F)

# 作图------------------------------------------------------------------------------
library(RColorBrewer)
scRNA.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
p <- DotPlot(scRNA, features = unique(top5$gene)) + RotatedAxis() +
        scale_x_discrete("") + scale_y_discrete("") +
        scale_color_gradientn(colours = brewer.pal(9,"BuGn")) +
        #coord_flip() +
        theme_bw() +
        theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "inches"))+
        theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
ggsave("scRNA-SAM-top5-DotPlot.pdf", plot = p, width = 12, height = 5)


p <- DoHeatmap(scRNA, features = unique(scRNA.markers$gene), angle = 0, slot = "data") + NoLegend() + 
              theme(axis.text.y = element_blank())
ggsave("scRNA-SAM-marker-heatmap.pdf", plot = p, width = 16, height = 8)


# ------------------------------细胞数量与差异基因数量展示 ---------------------------------
library(ggplot2)
library(reshape2)

# 定义图形主题
#自定义图表主题，对图表做精细调整；
top.mar=0.2
right.mar=0.1
bottom.mar=0.2
left.mar=0.1
#隐藏纵轴，并对字体样式、坐标轴的粗细、颜色、刻度长度进行限定；
mytheme<-theme_classic()+
  theme(text=element_text(family = "sans",colour ="gray30",size = 16),
        axis.text.y = element_text(size = 12),
        axis.line.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.6,colour = "gray30"),
        axis.ticks = element_line(size = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))

# 读入数据
cell <- read.csv('cell-number.csv',header=T)
gene <- read.csv('gene-number.csv',header=T)

#指定绘图顺序
cell$Cluster <- factor(cell$Cluster,
                       levels = rev(cell$Cluster),
                       ordered = T)
gene$Cluster <- factor(gene$Cluster,
                       levels = rev(gene$Cluster),
                       ordered = T)

#绘制右侧的条形图
left <- ggplot(cell,aes(x=Cluster,y=Cell_Number))+
  coord_flip()+
  labs(x="",y="Cell Number",title="")+
  geom_bar(fill="orange",colour="black",size=0.5,width = 0.8, 
           stat = "identity", alpha=0.7)+
  scale_y_continuous(expand=expansion(add = c(0.1, 0.1)),
                     limits = c(-1800, 0),
                     breaks = c(-1500,-1000, -500, 0),
                     label = c("1500","1000", "500","0"))+mytheme

#### 绘制右边的图形

right <- ggplot(gene,aes(x=Cluster,y=Gene_Number))+
  coord_flip()+labs(x="",y="Marker Gene Number",title="")+
  geom_bar(fill="purple",stat = "identity",colour="black", alpha=0.7,
           size=0.5,width = 0.8)+
  scale_y_continuous(expand=expansion(add = c(0.1, 0.1)),
                     limits = c(0, 1000),breaks = c(0, 200, 400, 600, 800),
                     label = c("0","200","400","600","800"))+
  mytheme+theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())

#使用patchwork包将两个图合并成一个；
library(patchwork)
left+right

# ------------------------------差异基因分布统计---------------------------------
library(UpSetR)

mycolors <- c("#66C2A5","#FC8D62","#8DA0CB","#A6D854","#FFD92F","#E5C494","#B3B3B3")

DEGs_Table <- read.csv('scRNA-SAM-marker.csv', header=T, row.names=1)

listInput <- list()

for(i in seq(1:7)){
   listInput[[i]] <- subset(DEGs_Table, cluster==i-1)$gene
}

cluster_names <-  c("C-00","C-01","C-02","C-03","C-04","C-05","C-06")

names(listInput) <- cluster_names

pdf(file="scRNA-SAM-DEGs-upset.pdf",onefile = FALSE,width=13, height=8)

upset(fromList(listInput),
      nsets=7,
      sets=rev(cluster_names),
      sets.bar.color = rev(mycolors),
      mb.ratio = c(0.55, 0.45),
      keep.order = TRUE,
      text.scale = c(2, 1.5, 2, 1.5, 1.5, 2),
      sets.x.label = "Marker Genes Per Cluster")

dev.off()

# ------------------------------类群相关性分析---------------------------------
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

scRNA <- readRDS("../../02-Clustering/02-Rice/03-SAM/02-Cluster/scRNA-SAM-Cluster.rds")
DefaultAssay(scRNA) <- 'RNA'

scRNA <- NormalizeData(scRNA, verbose = FALSE)

DEGs_Genes <- read.csv('scRNA-SAM-marker.csv',row.names=1, header=T)

expr_matrix <- as.data.frame(AverageExpression(scRNA, assays='RNA', features = unique(DEGs_Genes$gene)))
write.csv(expr_matrix, 'scRNA-SAM-Cluster-expr.csv', quote=F)

expr_cor <- cor(expr_matrix)

cluster_names <-  c("C-00","C-01","C-02","C-03","C-04","C-05","C-06")

rownames(expr_cor) <- cluster_names
colnames(expr_cor) <- cluster_names

write.csv(expr_cor, 'scRNA-SAM-Cluster-Cor.csv', quote=F)

# ------------------------------GO富集分析作图---------------------------------
library(ggplot2)

cluster0 <- read.csv('Cluster0-GO-Analysis.csv', header=T)
cluster1 <- read.csv('Cluster1-GO-Analysis.csv', header=T)
cluster2 <- read.csv('Cluster2-GO-Analysis.csv', header=T)
cluster3 <- read.csv('Cluster3-GO-Analysis.csv', header=T)
cluster4 <- read.csv('Cluster4-GO-Analysis.csv', header=T)
cluster5 <- read.csv('Cluster5-GO-Analysis.csv', header=T)

mycolors <- c("#66C2A5","#FC8D62","#8DA0CB","#A6D854","#FFD92F","#E5C494","#B3B3B3")

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
								axis.text = element_text(size=rel(4)),
								axis.title = element_text(size = rel(5)))
#挑选感兴趣的iterm
p1 <- ggplot(cluster0[c(2,7,9,17,18),], aes(y =-log10(fdr), x=reorder(Term,-log10(fdr)))) +
      geom_bar(stat="identity", fill='#66C2A5', width=0.9, alpha=0.7)+
      coord_flip() + 
      geom_text(aes(y=0,label = Term), hjust=0, vjust=0.5, size = 30) +
      ggtitle('Cluster0 GO Annotation') +
      ylab('-log10(padj)') + xlab('') +
      mytheme
p2 <- ggplot(cluster1[c(1,3,11,15,24),], aes(y =-log10(fdr), x=reorder(Term,-log10(fdr)))) +
      geom_bar(stat="identity", fill='#FC8D62', width=0.9, alpha=0.7)+
      coord_flip() + 
      geom_text(aes(y=0,label = Term), hjust=0, vjust=0.5, size = 30) +
      ggtitle('Cluster1 GO Annotation') +
      ylab('-log10(padj)') + xlab('') +
      mytheme
p3 <- ggplot(cluster2[c(4,10,23,27,28),], aes(y =-log10(fdr), x=reorder(Term,-log10(fdr)))) +
      geom_bar(stat="identity", fill='#8DA0CB', width=0.9, alpha=0.7)+
      coord_flip() + 
      geom_text(aes(y=0,label = Term), hjust=0, vjust=0.5, size = 30) +
      ggtitle('Cluster2 GO Annotation') +
      ylab('-log10(padj)') + xlab('') +
      mytheme
p4 <- ggplot(cluster3[c(2,13,16,22,29),], aes(y =-log10(fdr), x=reorder(Term,-log10(fdr)))) +
      geom_bar(stat="identity", fill='#A6D854', width=0.9, alpha=0.7)+
      coord_flip() + 
      geom_text(aes(y=0,label = Term), hjust=0, vjust=0.5, size = 30) +
      ggtitle('Cluster3 GO Annotation') +
      ylab('-log10(padj)') + xlab('') +
      mytheme
p5 <- ggplot(cluster4[c(1,2,3,9,12),], aes(y =-log10(fdr), x=reorder(Term,-log10(fdr)))) +
      geom_bar(stat="identity", fill='#FFD92F', width=0.9, alpha=0.7)+
      coord_flip() + 
      geom_text(aes(y=0,label = Term), hjust=0, vjust=0.5, size = 30) +
      ggtitle('Cluster4 GO Annotation') +
      ylab('-log10(padj)') + xlab('') +
      mytheme
p6 <- ggplot(cluster5[c(1,2,3,4,5),], aes(y =-log10(fdr), x=reorder(Term,-log10(fdr)))) +
      geom_bar(stat="identity", fill='#E5C494', width=0.9, alpha=0.7)+
      coord_flip() + 
      geom_text(aes(y=0,label = Term), hjust=0, vjust=0.5, size = 30) +
      ggtitle('Cluster5 GO Annotation') +
      ylab('-log10(padj)') + xlab('') +
      mytheme

library(patchwork)

p <- (p1 | p2 | p3) / (p4 | p5 | p6)

ggsave('Cluster-GO.pdf', p, width = 80, height = 40, limitsize = FALSE)