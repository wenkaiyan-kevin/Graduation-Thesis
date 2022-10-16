# ------------------------------Create a gene activity matrix---------------------------------
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

scATAC <- readRDS("../04-ATAC-RNA-integrated/scATAC-SAM-Integrated.rds")
DefaultAssay(scATAC) <- 'peaks'

Peaks <- CallPeaks(object = scATAC,
                   group.by = "seurat_clusters",
                   extsize = 150,
                   effective.genome.size = 373245519,
                   additional.args = '-q 0.01',
                   macs2.path = "~/biosoftware/miniconda/envs/macs2/bin/macs2")
saveRDS(Peaks,'Peaks-cluster.rds')

# 针对每个类群中Unique的Peaks作图
library(ggplot2)

data <- read.table('uniq-peaks-num.txt', header=T)
data$Cluster <- factor(data$Cluster, levels=c('Cluster4','Cluster3','Cluster2','Cluster1','Cluster0'))


p <- ggplot(data, aes(x=Cluster, y=Number)) +
            geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
            coord_flip() +
            xlab("") +
            theme_classic()
ggsave("uniq-peak-cluster.pdf", plot = p, width = 4, height = 3)

# ------------------------------Peaks区间进行注释---------------------------------
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

scATAC <- readRDS("../04-ATAC-RNA-integrated/scATAC-SAM-Integrated.rds")
DefaultAssay(scATAC) <- 'peaks'

Peaks <- readRDS("Peaks-cluster.rds")

Cluster0 <- subset(Peaks, peak_called_in==0)
Cluster0_genes <- ClosestFeature(scATAC, regions = granges(Cluster0))
Cluster0_genes_pass <- subset(Cluster0_genes,
                              type=="gene" &
                              gene_biotype == "protein_coding" &
                              distance < 5000)
genes <- unique(Cluster0_genes_pass$gene_id)
write.table(genes, 'Cluster0-genes.txt', row.names=F, col.names=F, quote=F)
length(unique(Cluster0_genes_pass$gene_id))

#Cluster0 3321
#Cluster1 959
#Cluster2 388
#Cluster3 1209
#Cluster4 65
library(ggplot2)

data <- read.table('uniq-peaks-asso-gene-num.txt', header=T)
#data$Cluster <- factor(data$Cluster, levels=c('Cluster4','Cluster3','Cluster2','Cluster1','Cluster0'))


p <- ggplot(data, aes(x=Cluster, y=Number)) +
            geom_bar(stat="identity", fill="#20B2AA", alpha=.6, width=.4) +
            xlab("") +
            theme_classic()
ggsave("uniq-peaks-asso-gene.pdf", plot = p, width = 4, height = 2)


###类群差异基因之间的交集分析
library(UpSetR)

mycolors <- c("#E95C59", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3")

listInput <- list()

for(i in seq(1:5)){
   filename <- paste('Cluster',i-1,'-genes.txt',sep='')
   genes <- read.table(filename, header=F, stringsAsFactors=F)[,1]
   listInput[[i]] <- genes
}

cluster_names <-  c("C-00","C-01","C-02","C-03","C-04")

names(listInput) <- cluster_names

pdf(file="scATAC-SAM-DEGs-upset.pdf",onefile = FALSE,width=13, height=8)

upset(fromList(listInput),
      nsets=7,
      sets=rev(cluster_names),
      sets.bar.color = rev(mycolors),
      mb.ratio = c(0.55, 0.45),
      keep.order = TRUE,
      text.scale = c(2, 1.5, 2, 1.5, 1.5),
      sets.x.label = "Marker Genes Per Cluster")

dev.off()


# GO注释作图
library(ggplot2)
library(patchwork)

mycolors <- c("#E95C59", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3")

Cluseter0 <- read.csv('Cluster0-plot.csv',header=T)
p1 <- ggplot(Cluseter0, aes(x=reorder(ID,Number), y=Number)) +
             geom_bar(stat="identity", fill="#E95C59", alpha=.6, width=.4) +
             coord_flip() +
             xlab("") + ylab("") +
             theme_classic()

Cluseter1 <- read.csv('Cluster1-plot.csv',header=T)
p2 <- ggplot(Cluseter1, aes(x=reorder(ID,Number), y=Number)) +
             geom_bar(stat="identity", fill="#53A85F", alpha=.6, width=.4) +
             coord_flip() +
             xlab("") + ylab("") + 
             theme_classic()

Cluseter2 <- read.csv('Cluster2-plot.csv',header=T)
p3 <- ggplot(Cluseter2, aes(x=reorder(ID,Number), y=Number)) +
             geom_bar(stat="identity", fill="#F1BB72", alpha=.6, width=.4) +
             coord_flip() +
             xlab("") + ylab("") + 
             theme_classic()

Cluseter3 <- read.csv('Cluster3-plot.csv',header=T)
p4 <- ggplot(Cluseter3, aes(x=reorder(ID,Number), y=Number)) +
             geom_bar(stat="identity", fill="#F3B1A0", alpha=.6, width=.4) +
             coord_flip() +
             xlab("") + ylab("") + 
             theme_classic()

Cluseter4 <- read.csv('Cluster4-plot.csv',header=T)
p5 <- ggplot(Cluseter4, aes(x=reorder(ID,Number), y=Number)) +
             geom_bar(stat="identity", fill="#D6E7A3", alpha=.6, width=.4) +
             coord_flip() + 
             xlab("") + ylab("") +
             theme_classic()

p <- (p1 | p2 | p3) / (p4 | p5 | p3)

ggsave("uniq-peak-cluster-GO.pdf", plot = p, width = 4, height = 3)

# 每个类群中特异性表达的基因的气泡图
library(Seurat)
library(Signac)
library(rtracklayer)
library(tidyverse)
library(patchwork)
library(dplyr)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 4)
set.seed(1234)

scATAC <- readRDS("../04-ATAC-RNA-integrated/scATAC-SAM-Integrated.rds")

feature <- read.table('Cluster4-genes.txt', header=F)[,1]
head(FindMarkers(scATAC, ident.1="4", features=feature, only.pos = T,
                      logfc.threshold = 0.01, min.pct = 0.1))

# 0 Os11g0451051 Os04g0557500
# 1 Os07g0132600 Os03g0599400
# 2 Os12g0207600 Os07g0663200
# 3 Os03g0822000 Os09g0475400
# 4 Os03g0267800 Os04g0550900

library(RColorBrewer)

feature <- c("Os11g0451051","Os04g0557500","Os07g0132600","Os03g0599400",
             "Os12g0207600","Os07g0663200","Os03g0822000","Os09g0475400","Os03g0267800","Os04g0550900")
p <- DotPlot(scATAC, features = feature) + RotatedAxis() +
        scale_x_discrete("") + scale_y_discrete("") +
        scale_color_gradientn(colours = brewer.pal(9,"OrRd")) +
        #coord_flip() +
        theme_bw() +
        theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "inches"))+
        theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
ggsave("scATAC-SAM-top2-DotPlot.pdf", plot = p, width = 6, height = 3)


# ------------------------------特异性基因的作图---------------------------------

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

scATAC <- readRDS("../04-ATAC-RNA-integrated/scATAC-SAM-Integrated.rds")
DefaultAssay(scATAC) <- 'peaks'


colors <- c("#E95C59", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3")

# Os01g0600900
p <- CoveragePlot(object = scATAC, idents = "0",
                  region = '1-36813736-36815023',
                  features = "Os01g0854500",
                  extend.upstream = 2000,
                  extend.downstream = 2000,
                  window = 30)
p <- p & scale_fill_manual(values = colors)
ggsave("gene-regions.pdf", plot = p, width = 8, height = 6)



DefaultAssay(scATAC) <- 'RNA'
library(RColorBrewer)

feature <- c("Os01g0854500","Os01g0600900","Os02g0527200","Os02g0643200")
p <- DotPlot(scATAC, features = feature) + RotatedAxis() +
        scale_x_discrete("") + scale_y_discrete("") +
        #scale_color_gradientn(colours = brewer.pal(9,"OrRd")) +
        #coord_flip() +
        theme_bw() +
        theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "inches"))+
        theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
ggsave("test-DotPlot.pdf", plot = p, width = 6, height = 3)

