library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 10)
rm(list = ls())
set.seed(1234)

# 读取数据
scRNA <- readRDS("../01-integrated/scRNA-RAM-Integrate.rds")

DefaultAssay(scRNA) <- "SCT"
scRNA <- PrepSCTFindMarkers(scRNA, assay = "SCT")

scRNA.markers <- FindAllMarkers(scRNA, assay = "SCT", only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.35)
write.table(scRNA.markers,'scRNA-RAM-DEGs.txt',quote=F,row.names=F,col.names=T,sep='\t')

> table(scRNA.markers$cluster)
  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
240 287 376 220 282 288 412 311 514 426 274 124 334 393 348 419 270 400 374

scRNA.markers.specific <- subset(scRNA.markers, pct.2 <= 0.1)
write.table(scRNA.markers.specific,'scRNA-RAM-DEGs-specific.txt',quote=F,row.names=F,col.names=T,sep='\t')

> table(scRNA.markers.specific$cluster)
 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
44  51 126  22  18  68  34  97  84  60  61   9 152  48  91 129 102 105 177

# 保存 矫正后的数据
saveRDS(scRNA, file = "scRNA-RAM-Integrate-fixed.rds")

# 作图
library(UpSetR)

# 定义颜色
mycolors <- c('#D5EEBB','#7FC8A9','#5F7A61','#444941',
              '#F67280','#C06C84','#6C5B7B','#355C7D',
              '#F5E8C7','#DEBA9D','#9E7777','#6F4C5B',
              '#FBB448','#E3670C','#04009A','#77ACF1',
              '#3EDBF0','#867AE9','#FFF5AB')

DEGs_Table <- read.table('scRNA-RAM-DEGs.txt', header=T)

listInput <- list()

for(i in seq(1:21)){
   listInput[[i]] <- subset(DEGs_Table, cluster==i-1)$gene
}

cluster_names <-  c("C-00","C-01","C-02","C-03","C-04",
                    "C-05","C-06","C-07","C-08","C-09",
                    paste(rep("C-",9), seq(10,18), sep=''))

names(listInput) <- cluster_names

pdf(file="scRNA-RAM-DEGs-upset.pdf",onefile = FALSE,width=13, height=8)

upset(fromList(listInput),
      nsets=19,
      sets=rev(cluster_names),
      sets.bar.color = rev(mycolors),
      mb.ratio = c(0.55, 0.45),
      keep.order = TRUE,
      text.scale = c(2, 1.3, 2, 1.3, 1.2, 1.0),
      sets.x.label = "Marker Genes Per Cluster")

dev.off()

library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(RColorBrewer)
rm(list = ls())
set.seed(1234)

scRNA <- readRDS("scRNA-RAM-Integrate-fixed.rds")

Epidermi_genes <- c('Os01g0248900', 'Os10g0454200', 'Os10g0578200')
Exodermis_genes <- c("Os04g0125700","Os04g0452700","Os03g0107300")
cortex_genes <- c("Os01g0914100","Os03g0135700","Os03g0135700")
Endodermal_genes <- c("Os10g0155100","Os02g0745100","Os07g0134000","Os12g0122000")
Pericycle_genes <- c("Os07g0634400","Os10g0524300","Os03g0216700")
Phloem_genes <- c("Os10g0154700","Os06g0614000","Os01g0236300")
Xylem_genes <- c('Os07g0638500','Os05g0108600','Os01g0971400')
Vascular_genes <- c('Os01g0750300','Os01g0750300','Os01g0750300')
Root_cap_genes <- c('Os03g0247200','Os03g0247200','Os03g0247200')
Meristem_genes <- c('Os05g0438700','Os06g0159450','Os03g0146300','Os01g0854500')

p <- FeaturePlot(scRNA, features = Vascular_genes, order = T, ncol = 3, pt.size = 2, 
                 min.cutoff = 'q10', max.cutoff = 'q95', raster = T)

p <- FeaturePlot(scRNA, features = Meristem_genes, order = T, ncol = 3, pt.size = 2, raster = T)
ggsave("Meristem.pdf", p, width = 11, height = 6)

library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(RColorBrewer)
rm(list = ls())
set.seed(1234)

scRNA <- readRDS("scRNA-RAM-Integrate-fixed.rds")

CO <- c("Os03g0103200","Os02g0745100")
C1 <- c("Os01g0914300","Os03g0135700", "Os07g0556000")
C2 <- c("Os08g0561500","Os05g0438700", "Os06g0159450", "Os03g0146300", "Os03g0162200")
C3 <- c("Os04g0125700","Os10g0555900")
C4 <- c("Os08g0115800","Os05g0140100")
C5 <- c("Os01g0750300","Os08g0489300")
C6 <- c("Os02g0251900","Os01g0975000")
C7 <- c("Os01g0248900","Os03g0155900","Os02g0595900")
C8 <- c("Os09g0475800")
C9 <- c("Os08g0531000","Os08g0531000")
C1O <- c("Os10g0454200")
C12 <- c("Os10g0122600","Os10g0578200")
C14 <- c("Os03g0247200")
C15 <- c("Os10g0155100")
C16 <- c("Os07g0634400")
C17 <- c("Os01g0257300","Os10g0154700","Os01g0236300")
C18 <- c("Os07g0638500","Os05g0108600","Os01g0842500")

genes <- c(CO,C2,C3,C4,C5,C6,C7,C8,C9,C1O,C12,C14,C15,C16,C17,C18)
levels(scRNA) <- c('7', '10', '12', '3', '9', '1', '11', '0', '6', '8',
                   '15', '5', '2', '13', '4', '16', '18', '17', '14')

p <- DotPlot(scRNA, features = unique(genes)) + RotatedAxis() + scale_color_gradientn(colours = brewer.pal(9,"Reds"))
ggsave("scRNA-RAM-marker-Dotplot.pdf", p, width = 20, height =10)

# 细胞类群重命名
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

scRNA <- readRDS("../02-DEGs/scRNA-RAM-Integrate-fixed.rds")

new.cluster.ids <- c("endodermis (Osat)", "cortex (Osat)", "meristem (Osat)", "exodermis (Osat)",
                     "dividing cell (Osat)", "vascular cylinder (Osat)", "endodermis (Osat)",
                     "epidermis/atrichoblast (Osat)", "endodermis (Osat)", "exodermis (Osat)",
                     "epidermis/meristem (Osat)", "cortex (Osat)", "epidermis/trichoblast (Osat)",
                     "meristem (Osat)", "root cap (Osat)",
                     "endodermis (Osat)", "pericycle (Osat)", "phloem (Osat)", "xylem (Osat)")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA$cellType <- Idents(scRNA)

## save seurat object
saveRDS(scRNA, "scRNA-RAM-annotation.rds")

