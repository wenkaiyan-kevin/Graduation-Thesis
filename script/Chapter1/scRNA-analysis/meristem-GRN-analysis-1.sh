# 基因symbol ID的添加以及转录因子的注释
vim scRNA-SAM-SMC-maker.txt

# 添加symbol ID
awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($1 in a){print $0"\t"a[$1]}else{print $0}}' \
    ~/seqlib/cellRanger_genome/rice/gene.name.txt scRNA-SAM-SMC-maker.txt \
    > scRNA-SAM-SMC-maker-SymId.txt

# 添加转录因子信息
# 从数据库中下载水稻转录因子信息：http://planttfdb.gao-lab.org/download.php
# 由于下载的数据基因是LOC基因，需要转换成Os基因, 然后注释我们的基因
cut -f 2,3 Osj_TF_list.txt | sed  '1d' | sort | uniq > Osj_TF_list_filter.txt

awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($2 in a){print $0"\t"a[$2]}}' \
    Osj_TF_list_filter.txt ~/seqlib/cellRanger_genome/rice/RAP-MSU_2021-05-10-v3.txt |  \
    cut -f 1,2,4 > TF-info.txt 

awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($1 in a){print $0"\t"a[$2]}else{print $0"\tNone\tNone\tNone"}}' \
    TF-info.txt scRNA-SAM-SMC-maker-SymId.txt | cut -f 1,3,6 > scRNA-SAM-SMC-maker-SymId-TF.txt