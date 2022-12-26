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

new.cluster.ids <- c("maturation", "elongation", "meristem", "elongation", "elongation",
                     "maturation", "elongation", "maturation", "maturation", "maturation",
                     "elongation", "maturation", "maturation", "meristem", "maturation",
                     "elongation", "elongation", "elongation", "elongation")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)

colors <- c('#00A549','#8FDB96','#E1F6DE')
p <- DimPlot(scRNA, reduction = "umap", label = F, cols = colors, pt.size=0.001)
ggsave("scRNA_RAM_development.pdf", p, width = 7, height = 6)


## 模块化细胞发育轨迹分析
### 表皮
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(dplyr)
rm(list = ls())
set.seed(1234)

scRNA <- readRDS("../../../02-DEGs/scRNA-RAM-Integrate-fixed.rds")
scRNA_ <- scRNA[,scRNA$seurat_clusters %in% c("2", "7", "10", "12")]

scRNA.cds <- as.cell_data_set(scRNA_)
# Calculate size factors using built-in function in monocle3
scRNA.cds <- estimate_size_factors(scRNA.cds)
# Add gene names into CDS
rowData(scRNA.cds)$gene_name <- rownames(scRNA.cds)
rowData(scRNA.cds)$gene_short_name <- rowData(scRNA.cds)$gene_name
#
scRNA.cds <- cluster_cells(cds = scRNA.cds, reduction_method = "UMAP")
scRNA.cds <- learn_graph(scRNA.cds, use_partition = TRUE)

umap_ <- as.data.frame(scRNA_@reductions$umap@cell.embeddings)
start_cells <- row.names(subset(umap_, UMAP_1 < 3 & UMAP_1 > 2 & UMAP_2 < -4 & UMAP_2 > -5))

scRNA.cds <- order_cells(scRNA.cds, reduction_method = "UMAP", root_cells = start_cells)

# plot trajectories colored by pseudotime
p <- plot_cells(cds = scRNA.cds, color_cells_by = "pseudotime", show_trajectory_graph = F)
ggsave("monocle3-epidermis.pdf", p, width = 7, height = 6)

# 出现问题 'rBind' is defunct; 'base::rbind' handles S4 objects since R 3.2.0
# https://github.com/cole-trapnell-lab/monocle3/issues/509
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
Track_genes <- graph_test(scRNA.cds, neighbor_graph="principal_graph", cores=8)
Track_genes <- Track_genes[,c(5,7,3,4,1,6)] %>% filter(q_value < 1e-4)

write.csv(Track_genes,"Epidrmis-Trajectory-genes.csv", row.names = F)

## save seurat object
saveRDS(scRNA.cds, "scRNA-RAM-monocle3-epidermis.rds")

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(dplyr)
rm(list = ls())
set.seed(1234)

scRNA.cds <- readRDS("scRNA-RAM-monocle3-epidermis.rds")

genes <- c("Os08g0116500","Os02g0121300","Os02g0595900","Os01g0248900","Os10g0578200")

p <- plot_genes_in_pseudotime(scRNA.cds[genes,],color_cells_by="seurat_clusters",
                              min_expr=0.5, ncol= 2, cell_size = 0.1) +
     geom_smooth(method="loess", se = F) +
     scale_color_manual(values=c('#5F7A61', '#355C7D', '#9E7777', '#FBB448'))

# 去除散点图和原始折线图
p$layers[[1]] <- NULL
p$layers[[1]] <- NULL
p$layers[[1]] <- NULL

# 保存结果
ggsave("epidermis-1.pdf", p, width = 6, height = 5)
# 保存数据分布
p <- p + geom_rug(aes(color = seurat_clusters), sides='b', size = 0.001, position = "jitter")
ggsave("epidermis-2.pdf", p, width = 6, height = 5)

genes <- row.names(subset(Track_genes, q_value == 0 & morans_I > 0.25))
#Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
pt.matrix <- exprs(scRNA.cds)[match(genes,rownames(rowData(scRNA.cds))),order(pseudotime(scRNA.cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

pdf('test.pdf',width=5, height=5)
Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = FALSE,
  show_column_names            = FALSE,
  row_dend_reorder             = TRUE,
  #row_names_gp                = gpar(fontsize = 6),
  row_km                       = 3,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = TRUE,
  cluster_columns              = FALSE,
  column_split                 = NULL,
  use_raster                   = TRUE)
dev.off()

pdf('test1.pdf',width=5, height=5)
Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = FALSE,
  show_column_names            = FALSE,
  row_dend_reorder             = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
dev.off()

### 外皮层
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(dplyr)
rm(list = ls())
set.seed(1234)

scRNA <- readRDS("../../../02-DEGs/scRNA-RAM-Integrate-fixed.rds")

umap_ <- as.data.frame(scRNA@reductions$umap@cell.embeddings)
cells <- row.names(subset(umap_, UMAP_2 > 5))
scRNA <- scRNA[,cells]

scRNA_ <- scRNA[,scRNA$seurat_clusters %in% c("3", "9", "13")]

scRNA.cds <- as.cell_data_set(scRNA_)
# Calculate size factors using built-in function in monocle3
scRNA.cds <- estimate_size_factors(scRNA.cds)
# Add gene names into CDS
rowData(scRNA.cds)$gene_name <- rownames(scRNA.cds)
rowData(scRNA.cds)$gene_short_name <- rowData(scRNA.cds)$gene_name
#
scRNA.cds <- cluster_cells(cds = scRNA.cds, reduction_method = "UMAP")
scRNA.cds <- learn_graph(scRNA.cds, use_partition = TRUE)

umap_ <- as.data.frame(scRNA_@reductions$umap@cell.embeddings)
start_cells <- row.names(subset(umap_, UMAP_1 < -3 & UMAP_1 > -4 & UMAP_2 < 10 & UMAP_2 > 8))

scRNA.cds <- order_cells(scRNA.cds, reduction_method = "UMAP", root_cells = start_cells)

# plot trajectories colored by pseudotime
p <- plot_cells(cds = scRNA.cds, color_cells_by = "pseudotime", show_trajectory_graph = F)
ggsave("monocle3-exodermis.pdf", p, width = 7, height = 6)

# 出现问题 'rBind' is defunct; 'base::rbind' handles S4 objects since R 3.2.0
# https://github.com/cole-trapnell-lab/monocle3/issues/509
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
Track_genes <- graph_test(scRNA.cds, neighbor_graph="principal_graph", cores=8)
Track_genes <- Track_genes[,c(5,7,3,4,1,6)] %>% filter(q_value < 1e-4)

write.csv(Track_genes,"exodermis-Trajectory-genes.csv", row.names = F)

## save seurat object
saveRDS(scRNA.cds, "scRNA-RAM-monocle3-exodermis.rds")

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(dplyr)
rm(list = ls())
set.seed(1234)

scRNA.cds <- readRDS("scRNA-RAM-monocle3-exodermis.rds")

Track_genes <- read.csv('exodermis-Trajectory-genes.csv')
Track_genes <- Track_genes[order(-Track_genes$morans_I),]
genes <- subset(Track_genes, q_value == 0 & morans_I > 0.25)[,1]
#genes <- c("Os03g0177400","Os02g0653200","Os06g0143100","Os07g0174900","Os03g0368900")
p <- plot_genes_in_pseudotime(scRNA.cds[genes[1:6],],color_cells_by="seurat_clusters",
                              min_expr=0.5, ncol= 2, cell_size = 0.1) +
     geom_smooth(method="loess", se = F) +
     scale_color_manual(values=c('#444941', '#DEBA9D', '#E3670C'))

ggsave("test.pdf", p, width = 6, height = 3)

# 去除散点图和原始折线图
p$layers[[1]] <- NULL
p$layers[[1]] <- NULL
p$layers[[1]] <- NULL

# 保存结果
ggsave("exodermis-1.pdf", p, width = 6, height = 5)
# 保存数据分布
p <- p + geom_rug(aes(color = seurat_clusters), sides='b', size = 0.001, position = "jitter")
ggsave("exodermis-2.pdf", p, width = 6, height = 5)

### 皮层细胞
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(dplyr)
rm(list = ls())
set.seed(1234)

scRNA <- readRDS("../../../02-DEGs/scRNA-RAM-Integrate-fixed.rds")
scRNA_ <- scRNA[,scRNA$seurat_clusters %in% c("1", "2", "11", "13")]

scRNA.cds <- as.cell_data_set(scRNA_)
# Calculate size factors using built-in function in monocle3
scRNA.cds <- estimate_size_factors(scRNA.cds)
# Add gene names into CDS
rowData(scRNA.cds)$gene_name <- rownames(scRNA.cds)
rowData(scRNA.cds)$gene_short_name <- rowData(scRNA.cds)$gene_name
#
scRNA.cds <- cluster_cells(cds = scRNA.cds, reduction_method = "UMAP")
scRNA.cds <- learn_graph(scRNA.cds, use_partition = TRUE)

umap_ <- as.data.frame(scRNA_@reductions$umap@cell.embeddings)
start_cells <- row.names(subset(umap_, UMAP_1 < 3 & UMAP_1 > 2 & UMAP_2 < -4 & UMAP_2 > -5))

scRNA.cds <- order_cells(scRNA.cds, reduction_method = "UMAP", root_cells = start_cells)

# plot trajectories colored by pseudotime
p <- plot_cells(cds = scRNA.cds, color_cells_by = "pseudotime", show_trajectory_graph = F)
ggsave("monocle3-cortex.pdf", p, width = 7, height = 6)

# 出现问题 'rBind' is defunct; 'base::rbind' handles S4 objects since R 3.2.0
# https://github.com/cole-trapnell-lab/monocle3/issues/509
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
Track_genes <- graph_test(scRNA.cds, neighbor_graph="principal_graph", cores=8)
Track_genes <- Track_genes[,c(5,7,3,4,1,6)] %>% filter(q_value < 1e-4)

write.csv(Track_genes,"cortex-Trajectory-genes.csv", row.names = F)

## save seurat object
saveRDS(scRNA.cds, "scRNA-RAM-monocle3-cortex.rds")

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(dplyr)
rm(list = ls())
set.seed(1234)

scRNA.cds <- readRDS("scRNA-RAM-monocle3-cortex.rds")

Track_genes <- read.csv('cortex-Trajectory-genes.csv')
Track_genes <- Track_genes[order(-Track_genes$morans_I),]
genes <- subset(Track_genes, q_value == 0 & morans_I > 0.25)[,1]


#genes <- c("Os03g0424800","Os02g0662000","Os01g0914300","Os04g0554500","Os10g0552700")
p <- plot_genes_in_pseudotime(scRNA.cds[genes,],color_cells_by="seurat_clusters",
                              min_expr=0.5, ncol= 2, cell_size = 0.1) +
     geom_smooth(method="loess", se = F) +
     scale_color_manual(values=c('#7FC8A9','#5F7A61','#6F4C5B','#E3670C'))

ggsave("test.pdf", p, width = 6, height = 3)

# 去除散点图和原始折线图
p$layers[[1]] <- NULL
p$layers[[1]] <- NULL
p$layers[[1]] <- NULL

# 保存结果
ggsave("cortex-1.pdf", p, width = 6, height = 5)
# 保存数据分布
p <- p + geom_rug(aes(color = seurat_clusters), sides='b', size = 0.001, position = "jitter")
ggsave("cortex-2.pdf", p, width = 6, height = 5)

### 内皮层
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(dplyr)
rm(list = ls())
set.seed(1234)

scRNA <- readRDS("../../../02-DEGs/scRNA-RAM-Integrate-fixed.rds")
scRNA_ <- scRNA[,scRNA$seurat_clusters %in% c("0", "4", "6", "8", "13")]

scRNA.cds <- as.cell_data_set(scRNA_)
# Calculate size factors using built-in function in monocle3
scRNA.cds <- estimate_size_factors(scRNA.cds)
# Add gene names into CDS
rowData(scRNA.cds)$gene_name <- rownames(scRNA.cds)
rowData(scRNA.cds)$gene_short_name <- rowData(scRNA.cds)$gene_name
#
scRNA.cds <- cluster_cells(cds = scRNA.cds, reduction_method = "UMAP")
scRNA.cds <- learn_graph(scRNA.cds, use_partition = TRUE)

umap_ <- as.data.frame(scRNA_@reductions$umap@cell.embeddings)
start_cells <- row.names(subset(umap_, UMAP_1 < 4 & UMAP_1 > 2 & UMAP_2 < -7 & UMAP_2 > -9))

scRNA.cds <- order_cells(scRNA.cds, reduction_method = "UMAP", root_cells = start_cells)

# plot trajectories colored by pseudotime
p <- plot_cells(cds = scRNA.cds, color_cells_by = "pseudotime", show_trajectory_graph = F)
ggsave("monocle3-endodermis.pdf", p, width = 7, height = 6)

# 出现问题 'rBind' is defunct; 'base::rbind' handles S4 objects since R 3.2.0
# https://github.com/cole-trapnell-lab/monocle3/issues/509
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
Track_genes <- graph_test(scRNA.cds, neighbor_graph="principal_graph", cores=8)
Track_genes <- Track_genes[,c(5,7,3,4,1,6)] %>% filter(q_value < 1e-4)

write.csv(Track_genes,"endodermis-Trajectory-genes.csv", row.names = F)

## save seurat object
saveRDS(scRNA.cds, "scRNA-RAM-monocle3-endodermis.rds")

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(dplyr)
rm(list = ls())
set.seed(1234)

scRNA.cds <- readRDS("scRNA-RAM-monocle3-endodermis.rds")

Track_genes <- read.csv('endodermis-Trajectory-genes.csv')
Track_genes <- Track_genes[order(-Track_genes$morans_I),]
genes <- subset(Track_genes, q_value == 0 & morans_I > 0.25)[,1]

#genes <- c("Os03g0424800","Os02g0662000","Os01g0914300","Os04g0554500","Os10g0552700")
p <- plot_genes_in_pseudotime(scRNA.cds[genes[1:6],],color_cells_by="seurat_clusters",
                              min_expr=0.5, ncol= 2, cell_size = 0.1) +
     geom_smooth(method="loess", se = F) +
     scale_color_manual(values=c('#D5EEBB', '#F67280', '#6C5B7B', '#F5E8C7', '#E3670C'))

ggsave("test.pdf", p, width = 6, height = 3)

# 去除散点图和原始折线图
p$layers[[1]] <- NULL
p$layers[[1]] <- NULL
p$layers[[1]] <- NULL

# 保存结果
ggsave("endodermis-1.pdf", p, width = 6, height = 5)
# 保存数据分布
p <- p + geom_rug(aes(color = seurat_clusters), sides='b', size = 0.001, position = "jitter")
ggsave("endodermis-2.pdf", p, width = 6, height = 5)

### 维管
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(dplyr)
rm(list = ls())
set.seed(1234)

scRNA <- readRDS("../../../02-DEGs/scRNA-RAM-Integrate-fixed.rds")
scRNA_ <- scRNA[,scRNA$seurat_clusters %in% c("2", "4", "5", "13")]

scRNA.cds <- as.cell_data_set(scRNA_)
# Calculate size factors using built-in function in monocle3
scRNA.cds <- estimate_size_factors(scRNA.cds)
# Add gene names into CDS
rowData(scRNA.cds)$gene_name <- rownames(scRNA.cds)
rowData(scRNA.cds)$gene_short_name <- rowData(scRNA.cds)$gene_name
#
scRNA.cds <- cluster_cells(cds = scRNA.cds, reduction_method = "UMAP")
scRNA.cds <- learn_graph(scRNA.cds, use_partition = TRUE)

umap_ <- as.data.frame(scRNA_@reductions$umap@cell.embeddings)
start_cells <- row.names(subset(umap_, UMAP_1 < 3 & UMAP_1 > 2 & UMAP_2 < -4 & UMAP_2 > -5))

scRNA.cds <- order_cells(scRNA.cds, reduction_method = "UMAP", root_cells = start_cells)

# plot trajectories colored by pseudotime
p <- plot_cells(cds = scRNA.cds, color_cells_by = "pseudotime", show_trajectory_graph = F)
ggsave("monocle3-VC.pdf", p, width = 7, height = 6)

# 出现问题 'rBind' is defunct; 'base::rbind' handles S4 objects since R 3.2.0
# https://github.com/cole-trapnell-lab/monocle3/issues/509
# trace('calculateLW', edit = T, where = asNamespace("monocle3"))
Track_genes <- graph_test(scRNA.cds, neighbor_graph="principal_graph", cores=8)
Track_genes <- Track_genes[,c(5,7,3,4,1,6)] %>% filter(q_value < 1e-4)

write.csv(Track_genes,"VC-Trajectory-genes.csv", row.names = F)

## save seurat object
saveRDS(scRNA.cds, "scRNA-RAM-monocle3-VC.rds")

library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(dplyr)
rm(list = ls())
set.seed(1234)

scRNA.cds <- readRDS("scRNA-RAM-monocle3-VC.rds")

Track_genes <- read.csv('VC-Trajectory-genes.csv')
Track_genes <- Track_genes[order(-Track_genes$morans_I),]
genes <- subset(Track_genes, q_value == 0 & morans_I > 0.25)[,1]

genes <- c("Os05g0160200","Os08g0559200","Os11g0629400","Os11g0167800","Os03g0416200")
p <- plot_genes_in_pseudotime(scRNA.cds[genes,],color_cells_by="seurat_clusters",
                              min_expr=0.5, ncol= 2, cell_size = 0.1) +
     geom_smooth(method="loess", se = F) +
     scale_color_manual(values=c('#5F7A61','#F67280','#C06C84','#E3670C'))

ggsave("test.pdf", p, width = 6, height = 3)

# 去除散点图和原始折线图
p$layers[[1]] <- NULL
p$layers[[1]] <- NULL
p$layers[[1]] <- NULL

# 保存结果
ggsave("VC-1.pdf", p, width = 6, height = 5)
# 保存数据分布
p <- p + geom_rug(aes(color = seurat_clusters), sides='b', size = 0.001, position = "jitter")
ggsave("VC-2.pdf", p, width = 6, height = 5)