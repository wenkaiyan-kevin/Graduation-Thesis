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

scRNA <- readRDS("../../01-integrated/scRNA-RAM-Integrate.rds")

umap_ <- as.data.frame(scRNA@reductions$umap@cell.embeddings)
cells <- row.names(subset(umap_, UMAP_1 > 5 & UMAP_2 < 5))

scRNA_ <- scRNA[,scRNA$seurat_clusters %in% c("7", "10", "12")]

## save seurat object
saveRDS(scRNA_, "scRNA-RAM-epidermis.rds")
## 作图
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

scRNA <- readRDS("../../01-integrated/scRNA-RAM-Integrate.rds")

colors <- c(rep("lightgrey",7),"#355C7D","lightgrey","lightgrey","#9E7777",
            "lightgrey","#FBB448",rep("lightgrey",6))
p <- DimPlot(scRNA, reduction = "umap", label = F, cols = colors, pt.size=0.001)
ggsave("scRNA-RAM-epidermis.pdf", p, width = 7, height = 6)

## monocle2分析
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(monocle)
library(RColorBrewer)
library(SeuratWrappers)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 10)
rm(list = ls())
set.seed(1234)

scRNA <- readRDS("../01-subset-cluster/scRNA-RAM-epidermis.rds")
scRNA$celltype <- Idents(scRNA)
## source这个关键词在monocle2中貌似有冲突
#pData(cds)$sample <- pData(cds)$source
DefaultAssay(scRNA) <- 'RNA'

expr_matrix <- as(as.matrix(scRNA@assays$RNA@counts), 'sparseMatrix')
p_data <- scRNA@meta.data
f_data <- data.frame(gene_short_name = row.names(scRNA),row.names = row.names(scRNA))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.1,
                      expressionFamily = negbinomial.size())

pData(cds)$sample <- pData(cds)$source

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
dim(fData(cds)) #此时有28444个基因
#过滤掉在小于10个细胞中表达的基因，还剩19595个基因。
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10)) 

disp_table <- dispersionTable(cds)
## 1148个基因
ordergene <- subset(disp_table, mean_expression >= 0.1 & mean_expression <= 1 & 
                    dispersion_empirical >= dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, ordergene)

cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')

cds <- orderCells(cds)
## 根据生物学只是重新指定root
# cds <- orderCells(cds, root_state = '1')

# 按照时间发育作图
p <- plot_complex_cell_trajectory(cds, color_by="Pseudotime", cell_size = 0.5) + 
                                  scale_color_gradient2(low="#feebe2", mid="#fa9fb5", high="#c51b8a")
ggsave("scRNA-epidermis-Pseudotime.pdf",p,width = 7,height = 7)
# 按照数据来源作图
p <- plot_complex_cell_trajectory(cds, color_by="seurat_clusters", cell_size = 0.5) + 
                                  scale_color_manual(values = c("#355C7D","#9E7777","#FBB448"))
ggsave("scRNA-epidermis-source.pdf",p,width = 7,height = 7)

saveRDS(cds, "scRNA-RAM-epidermis-CDS.rds")

# BEAM分析
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(monocle)
library(RColorBrewer)
library(SeuratWrappers)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 10)
rm(list = ls())
set.seed(1234)

cds <- readRDS("scRNA-RAM-epidermis-CDS.rds")

disp_table <- dispersionTable(cds)
## 1148个基因
ordergene <- subset(disp_table, mean_expression >= 0.1 & mean_expression <= 1 & 
                    dispersion_empirical >= dispersion_fit)$gene_id

BEAM_res <- BEAM(cds[ordergene,], branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

#1148个基因可视化
BEAM_genes <- top_n(BEAM_res, n = 1148, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = 1,
                                 branch_labels = c("fate 1(Trichoblasts)", "fate 2(Atrichoblasts)"),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"),
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(60),
                                 num_clusters = 4, show_rownames = F, return_heatmap = T)
ggsave("scRNA_BEAM_epidermis.pdf", p$ph_res, width = 6, height = 5)

# 保存聚类结果
data <- p$annotation_row
write.csv(data, 'scRNA_BEAM_epidermis.csv', quote=F)
cluster1 <- rownames(subset(data, Cluster==1))

#cluster5 <- gene[which(p$annotation_row$Cluster==5)]
#write.table(cluster5,'cluster5-genes.txt',quote=F,row.names=F,col.names=F)

## GO富集分析
library(topGO)

# 设置输入文件，后面直接在这个地方修改文件名称就可以直接运行了
input="Cluster4-gene.txt"  #差异基因名称的列表
mapfile="/gss1/home/yanwk/seqlib/GO-Annotation-list/rice/IRGSP-Gene2Go.txt"    #所有基因GO map结果

# 读取数据
gene2GO <- readMappings(file = mapfile) #如果是读取其他文件格式，后面参数还需要修改
gene_names <- names(gene2GO)
my_genes <- read.table(input, header=F)[,1]

genelist <- factor(as.integer(gene_names %in% my_genes))
names(genelist) <- gene_names

# 开始分析
top_diff_genes <- function(allScore){return(allScore<0.05)}

GO_data <- new("topGOdata",ontology='BP',allGenes=genelist,
               annot=annFUN.gene2GO, gene2GO=gene2GO, geneSel=top_diff_genes)

result <- runTest(GO_data, algorithm = "classic", statistic = "fisher")

#提取基因 table
allGO <- usedGO(object = GO_data)
allres <- GenTable(GO_data,classicFisher=result, orderBy="classic",ranksOf="classicFisher", topNodes=50)

# 统计检验矫正
fdr <- p.adjust(p=allres[,"classicFisher"], method="fdr")
final_res <- cbind(allres,fdr)

# 过滤得到显著
final_res <- subset(final_res, fdr < 0.1)

#输出结果
write.csv(final_res, "Cluster4_GO_Analysis.csv", quote=F, row.names=F)

# 作图
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggsci)
library(monocle)
library(RColorBrewer)
library(SeuratWrappers)
library(future)
options(future.globals.maxSize = 50*1024^3)
options(future.rng.onMisuse="ignore")
plan("multicore", workers = 10)
rm(list = ls())
set.seed(1234)

cds <- readRDS("../02-monocle2-analysis/scRNA-RAM-epidermis-CDS.rds")

genes <- c("Os01g0248900","Os09g0417600","Os10g0578200","Os10g0454300")

p <- plot_genes_branched_pseudotime(cds[genes], branch_point = 1, color_by = "celltype",
                               cell_size=1, ncol = 2) + 
                               scale_color_manual(values = c("#355C7D","#9E7777","#FBB448"))
ggsave("scRNA-epidermis-gene-Pseudotime.pdf", p, width = 10, height = 5)












