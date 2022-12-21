# 水稻顶端分生组织单细胞转录组与染色质可及性研究
## 论文内容  
本论文的主要研究内容为

## 仓库信息 
在该仓库中可以查询到论文中各个章节中所涉及的主要分析程序代码与分析结果

## 分析环境的搭建  
推荐使用conda搭建论文的生信分析环境，关于conda的安装与使用请[点击](https://docs.conda.io/en/latest/)。  
**创建环境**
```
conda create -n sc-env
```

**安装R分析环境所需软件**
```
conda install -c bioconda -c conda-forge -c bu_cnio -c bioturing r-seuratdata r-seuratwrappers \
        bioconductor-cicero bioconductor-rsamtools bioconductor-motifmatchr bioconductor-jaspar2020 \
        bioconductor-tfbstools bioconductor-monocle r-ggplot2 r-seurat r-seuratdisk  \
        r-seuratObject r-velocyto.r r-signac r-monocle3 r-wgcna r-igraph r-devtools r-base r-cairo
```
**安装python分析环境所需软件**
```
conda install -c bioconda -c conda-forge python scanpy scvelo numpy pandas matplotlib loompy velocyto.py
```


## 一、水稻茎端分生组织单细胞转录组与染色质可及性图谱绘制
### 1.1 单细胞转录组分析结果
**1.1.1 数据预处理**  
本研究使用10x Genomics官方软件Cell Ranger（v.6.0.2）构建单细胞基因表达矩阵。
```
cellranger count --id=scRNA_Rice_SAM --transcriptome=/gss1/home/yanwk/seqlib/cellRanger_genome/rice/cellranger-RNA-IRGSP-genome/Oryza_sativa_IRGSP \
                 --fastqs=/gss1/home/yanwk/ywk_Graduation_Project/00-SingleCell-data/02-Rice/02-Rice-Shoot/01-Mydata \
                 --force-cells=10000   
```
**1.1.2 细胞的聚类与分群**  
基于单细胞基因表达矩阵，我们将测序所得细胞进行过滤、分群与聚类，并对每个类群中的`marker`基因进行鉴定。  
- [质控流程](script/Chapter1/scRNA-analysis/quality-control.r)
- [聚类分群流程](script/Chapter1/scRNA-analysis/cell-clustering.r)  

**1.1.3 marker基因鉴定**  
- [差异基因分析流程](script/Chapter1/scRNA-analysis/DEGs-analysis.r)  

**1.1.4 RNA速率分析**  
- [loom文件的生成](script/Chapter1/scRNA-analysis/RNA-velocity-1.sh)
- [转录剪切矩阵的生成](script/Chapter1/scRNA-analysis/RNA-velocity-2.r)
- [RNA速率的计算](script/Chapter1/scRNA-analysis/RNA-velocity-3.py)  

**1.1.5 茎端表皮细胞转录组分析**  
- [分析内容](script/Chapter1/scRNA-analysis/epidermis-analysis.r)  

**1.1.6 初始分生细胞转录组分析**  
- [聚类分析](script/Chapter1/scRNA-analysis/meristem-analysis-1.r)
- [拟时序分析](script/Chapter1/scRNA-analysis/meristem-analysis-2.py)
- [调控网络分析](script/Chapter1/scRNA-analysis/meristem-GRN-analysis-2.r)

### 1.2 单细胞染色质可及性分析结果
**1.2.1 数据预处理** 
```
cellranger-atac count --id SAM --reference /home/wkyan/ywk_lab/04-scATAC-analysis/00-seqlib/IRGSP \
                      --fastqs /home/wkyan/ywk_lab/04-scATAC-analysis/01-raw-data/03-Shoot/data \
                      --sample S_20201215NA --force-cells 8000
```
**1.2.2 细胞的聚类与分群**  
基于单细胞染色质可及性矩阵，我们将测序所得细胞进行过滤、分群与聚类，并对每个类群中的`marker peak`进行鉴定。
- [BSgenome的构建](script/Chapter1/scATAC-analysis/BSgenome-create.sh)
- [基本信息统计](script/Chapter1/scATAC-analysis/data-info-analysis.sh)
- [质控流程](script/Chapter1/scATAC-analysis/quality-control.r)
- [聚类分群流程](script/Chapter1/scATAC-analysis/cell-clustering.r)
- [染色质可及性与基因表达关系分析](script/Chapter1/scATAC-analysis/RNA-ATAC-relation-analysis.sh)


**1.2.3 scRNA数据对scATAC细胞类群注释**  
分析流程如下：  
- [注释流程](script/Chapter1/scATAC-analysis/scRNA-scATAC-annotation.r)


## 二、水稻根端分生组织单细胞转录组与染色质可及性图谱绘制
### 2.1 单细胞转录组分析结果


### 2.2 单细胞染色质可及性分析结果

## 三、单细胞精度下茎端与根端分生组织差异性研究

