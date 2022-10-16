
#-------------------------------TSS富集分数作图-------------------------------

samtools view -@ 24 -ShuF 4 -f 2 -q 30 /gss1/home/yanwk/ywk_Graduation_Project/01-Preprocessing/02-Rice/02-Rice-Shoot/02-Mydata-ATAC/outs/possorted_bam.bam \
| samtools sort -@ 24 - -o Shoot.aln_f2q30.sort.bam
samtools index -@ 24 Shoot.aln_f2q30.sort.bam

##去重复
picard MarkDuplicates -INPUT Shoot.aln_f2q30.sort.bam -OUTPUT Shoot.aln_f2q30.sort.pmd.bam -REMOVE_DUPLICATES true -ASSUME_SORTED true -METRICS_FILE Shoot.aln_f2q30_pmd.out
samtools index -@ 24 Shoot.aln_f2q30.sort.pmd.bam
rm -rf Shoot.aln_f2q30.sort.bam

#samtools view -h Shoot.aln_f2q30.sort.pmd.bam | grep -v 'Mt' | grep -v "Pt" | samtools view -bS -o Shoot.final.bam
#samtools index -@ 24 Shoot.final.bam
#rm -rf Shoot.aln_f2q30.sort.pmd.bam


bamCoverage -p 24 -b Shoot.aln_f2q30.sort.pmd.bam \
            --binSize 50 \
            --smoothLength 150 \
            --normalizeUsing BPM \
            --maxFragmentLength 150 \
            -o Shoot.aln_f2q30.sort.pmd.bw

computeMatrix scale-regions -p 24 -S Shoot.aln_f2q30.sort.pmd.bw \
                            -R ~/seqlib/cellRanger_genome/rice/Genome-Features/gene.bed \
                            --beforeRegionStartLength 2000 \
                            --regionBodyLength 3000 \
                            --afterRegionStartLength 2000 \
                            --missingDataAsZero \
                            --skipZeros --maxThreshold 10 -o matrix.mat.gz

plotHeatmap -m matrix.mat.gz --heatmapHeight 16 --heatmapWidth 6 \
            --colorMap BuGn \
            -out TSS-enrichment.pdf
#-------------------------------Peak区间注释-------------------------------

bedtools intersect -a ../Peaks.bed \
-b ~/seqlib/cellRanger_genome/rice/Genome-Features/dwonstream-200bp.bed -wa | \
sort -n -k1 -k2 | uniq | wc -l
#3553

bedtools intersect -a ../Peaks.bed \
-b ~/seqlib/cellRanger_genome/rice/Genome-Features/five_prime_utr.bed -wa | \
sort -n -k1 -k2 | uniq | wc -l
#14192

bedtools intersect -a ../Peaks.bed \
-b ~/seqlib/cellRanger_genome/rice/Genome-Features/intergenic.bed -wa | \
sort -n -k1 -k2 | uniq | wc -l
#40741

bedtools intersect -a ../Peaks.bed \
-b ~/seqlib/cellRanger_genome/rice/Genome-Features/three_prime_utr.bed -wa | \
sort -n -k1 -k2 | uniq | wc -l
#5285

bedtools intersect -a ../Peaks.bed \
-b ~/seqlib/cellRanger_genome/rice/Genome-Features/exon.bed -wa | \
sort -n -k1 -k2 | uniq | wc -l
#27080

bedtools intersect -a ../Peaks.bed \
-b ~/seqlib/cellRanger_genome/rice/Genome-Features/gene.bed -wa | \
sort -n -k1 -k2 | uniq | wc -l
#28923

bedtools intersect -a ../Peaks.bed \
-b ~/seqlib/cellRanger_genome/rice/Genome-Features/promoter-1kB.bed -wa | \
sort -n -k1 -k2 | uniq | wc -l
#18552


# 作图
library(ggplot2)

data <- read.table('annotation.txt', header=T)

colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f")
          
p <- ggplot(data, aes_string(x = 1, fill = "Type", y = "Fraction"))
p <- p + geom_bar(stat="identity") + coord_flip() + theme_bw()
p <- p + ylab("Percentage(%)") + xlab('') + ggtitle("Feature Distribution")
p <- p + scale_x_continuous(breaks=NULL)
p <- p + scale_fill_manual(values=colors, guide=guide_legend(reverse=TRUE))

ggsave("Peaks-annotation.pdf", plot = p, width = 8, height = 6)