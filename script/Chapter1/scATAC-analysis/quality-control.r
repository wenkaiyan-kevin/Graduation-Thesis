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

# 读取peaks
Peaks <- readRDS('Peaks.rds')

# 读取metadata
meta_file <- "/gss1/home/yanwk/ywk_Graduation_Project/01-Preprocessing/02-Rice/02-Rice-Shoot/02-Mydata-ATAC/outs/singlecell.csv"
metadata <- read.csv(file = meta_file, header = TRUE, row.names = 1)
metadata <- subset(metadata, is__cell_barcode == 1)

# 构建counts矩阵
fragment_file <- '/gss1/home/yanwk/ywk_Graduation_Project/01-Preprocessing/02-Rice/02-Rice-Shoot/02-Mydata-ATAC/outs/fragments.tsv.gz'
fragment <- CreateFragmentObject(path = fragment_file, cells = rownames(metadata))

counts <- FeatureMatrix(fragments = fragment,
                        features = Peaks,
                        cells = rownames(metadata))

# 构建seurat对象
chrom_assay <- CreateChromatinAssay(counts = counts,
                                    sep = c("-", "-"),
                                    genome = seqinfo(BSgenome.Osativa.IRGSP.IRGSP1),
                                    fragments = fragment_file,
                                    min.cells = 10,
                                    min.features = 200)

scATAC <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)
scATAC$lib <- 'scATAC'

# An object of class Seurat
# 72098 features across 8001 samples within 1 assay
# Active assay: peaks (72098 features, 0 variable features)


# 添加注释
gtf <- import('/gss1/home/yanwk/seqlib/cellRanger_genome/rice/Oryza_sativa.IRGSP-1.0.51.gtf')
gtf$gene_name <- gtf$gene_id
gtf$tx_id <- gtf$gene_id
Annotation(scATAC) <- gtf

# 过滤细胞
# compute nucleosome signal score per cell
scATAC <- NucleosomeSignal(object = scATAC)
# compute TSS enrichment score per cell
scATAC <- TSSEnrichment(object = scATAC, fast = TRUE)
# fraction of reads in peaks
scATAC$pct_reads_in_peaks <- scATAC$peak_region_fragments / scATAC$passed_filters * 100

scATAC <- subset(scATAC, subset = peak_region_fragments > 3000 & peak_region_fragments < 20000)

# An object of class Seurat
# 72098 features across 8000 samples within 1 assay
# Active assay: peaks (72098 features, 0 variable features)

# 针对指标作图
p <- VlnPlot(object = scATAC, ncol = 4, pt.size = 0,
             features = c('pct_reads_in_peaks', 'peak_region_fragments',
                          'TSS.enrichment', 'nucleosome_signal'))

ggsave("scATAC-QC.pdf", plot = p, width = 12, height = 6)

# 保存数据
saveRDS(scATAC, file = "scATAC-SAM-QC.rds")