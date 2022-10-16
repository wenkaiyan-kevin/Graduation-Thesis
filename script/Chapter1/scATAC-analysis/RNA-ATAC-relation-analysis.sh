工作目录：/gss1/home/yanwk/ywk_work/5_other_Analysis/1_Shoot_ATAC_RNA_relation


#我首先需要把scRNA-seq数据中每个基因的表达量给计算出来
##使用R中的seurat包读取数据
library(Seurat)

counts <- Read10X_h5(filename = "/gss1/home/yanwk/ywk_work/4_scRNA_Analysis/1_MP_Analysis/Root/outs/filtered_feature_bc_matrix.h5")

counts_mat <- as.matrix(counts)

row_sum <- as.data.frame(rowSums(counts_mat))

write.table(row_sum, 'sc-RNA-gene-count.txt', row.names=T, col.names=F, sep='\t', quote=F)

##然后把输出的结果中关于线粒体和叶绿体的基因删除

grep -v -E '^LOC_Osp|LOC_Osm' sc-RNA-gene-count.txt > temp

rm -rf sc-RNA-gene-count.txt

mv temp sc-RNA-gene-count.txt

##现在需要把基因的长度加入到我们的文件中去

paste sc-RNA-gene-count.txt gene.length.txt | cut -f 1,2,7 > sc-RNA-gene-count-length.txt

#给文件加入表头
gene	count	length

##使用python的bioinfokit把数据
from bioinfokit.analys import norm, get_data
import pandas as pd

sc_exp = pd.read_table("sc-RNA-gene-count-length.txt", index_col = 'gene')
nm = norm()
nm.tpm(df=sc_exp, gl='length')
tpm_df = nm.tpm_norm
tpm_df.to_csv("sc-RNA-gene-count-TPM.txt", index_label="gene", sep='\t')


Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    0.000     0.088     1.007    35.261     7.641 18788.432



#################
#接下来是对bam文件进行处理，过滤后进行使用
工作目录:/gss1/home/yanwk/ywk_work/3_scATAC_Analysis/1_Shoot_Analysis/1_cellranger_result/Shoot/outs/filtered_bam
samtools view -@ 24 -ShuF 4 -f 2 -q 30 ../possorted_bam.bam | samtools sort -@ 24 - -o Shoot.aln_f2q30.sort.bam

samtools index -@ 24 Shoot.aln_f2q30.sort.bam


##去重复
picard MarkDuplicates -INPUT Shoot.aln_f2q30.sort.bam -OUTPUT Shoot.aln_f2q30.sort.pmd.bam -REMOVE_DUPLICATES true -ASSUME_SORTED true -METRICS_FILE Shoot.aln_f2q30_pmd.out

samtools index -@ 24 Shoot.aln_f2q30.sort.pmd.bam

#################
bamCoverage -p 24 -b Shoot.aln_f2q30.sort.pmd.bam \
            --binSize 50 \
            --smoothLength 150 \
            --normalizeUsing BPM \
            --maxFragmentLength 150 \
            -o Shoot.aln_f2q30.sort.pmd.bw

computeMatrix scale-regions -p 24 --maxThreshold 50 -R gene.zone.bed -S ../1_BamToBw/GL-input.bw ../1_BamToBw/GL-1.bw ../1_BamToBw/GL-2.bw --samplesLabel GL_input GL_1 GL_2 -b 2000 -a 2000 --regionBodyLength 4000 --skipZeros -o GL_matrix.mat_Region.gz


computeMatrix scale-regions -p 24 -S Shoot.aln_f2q30.sort.pmd.bw \
                            -R gene_exp_level_0-0.bed gene_exp_level_0-0.bed gene_exp_level_25-50.bed gene_exp_level_50-75.bed gene_exp_level_75-100.bed \
                            --beforeRegionStartLength 2000 \
                            --regionBodyLength 3000 \
                            --afterRegionStartLength 2000 \
                            --skipZeros --maxThreshold 10 -o matrix.mat.gz

plotProfile -m matrix.mat.gz \
            -out RNA_ATAC_relation.pdf \
            --samplesLabel "Relationship of RNA and ATAC"