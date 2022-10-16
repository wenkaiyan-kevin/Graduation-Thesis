## 1、构建水稻BSgenome

工作目录：/gss1/home/yanwk/seqlib/cellRanger_genome/rice/BSgenome

### 1.1 首先，压缩基因组文件


faToTwoBit Oryza_sativa.IRGSP-1.0.dna.toplevel.chr.fa IRGSP-1.0_genome.2bit


### 1.2 创建`IRGSP-1.0_genome.seed`文件


Package: BSgenome.Osativa.IRGSP.IRGSP1
Title: Full genome sequences for Oryza sativa (IRGSP version IRGSP-1.0)
Description: Full genome sequences for Oryza sativa(rice) as provided by RAP-DB (IRGSP-1.0, Jun. 2020) and stored in Biostrings objects.
Version: 1.0.0
License: GPL (>= 2)
organism: Oryza sativa
common_name: Rice
provider: MUS
provider_version: MUS7
release_date: Jun. 2020
release_name: Rice Annotation Project Database
organism_biocview: Oryza_sativa
BSgenomeObjname: Osativa
seqnames: c("1","2","3","4","5","6","7","8","9","10","11","12","Mt","Pt")
circ_seqs: c("Mt","Pt")
SrcDataFiles: IRGSP-1.0_genome.fasta.gz from ftp://ftp.ensemblgenomes.org/pub/plants/release-48/
seqs_srcdir: /gss1/home/yanwk/seqlib/cellRanger_genome/rice/BSgenome
seqfile_name: IRGSP-1.0_genome.2bit


#然后在R中执行

library(BSgenome)
forgeBSgenomeDataPkg("IRGSP-1.0_genome.seed")


### 1.3 创建R包和安装


R CMD build BSgenome.Osativa.IRGSP.IRGSP1

R CMD check BSgenome.Osativa.IRGSP.IRGSP1_1.0.0.tar.gz

# 安装
R CMD INSTALL BSgenome.Osativa.IRGSP.IRGSP1_1.0.0.tar.gz
