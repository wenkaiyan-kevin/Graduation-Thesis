# 水稻顶端分生组织单细胞转录组与染色质可及性研究
## 论文内容简介  
&emsp;在植物的整个生命周期中，顶端分生组织是维持其不断生长和发育的重要组织。植物的顶端分生组织分为茎端分生组织（Shoot Apical Meristem，SAM）和根端分生组织（Root Apical Meristem，RAM），它们主要分布于植物的茎尖和根尖。茎端分生组织和根端分生组织都具有自我更新和无限分化的潜能，能够通过其组织内干细胞的分裂和分化实现自身的增殖与其他组织器官的形成。水稻是世界主要的粮食作物，在国家粮食安全保障体系中占有重要地位。随着世界人口数量的不断增加，耕地面积的不断减少，培育优质高产水稻品种是提高粮食产量的一个重要途径。在水稻中，顶端分生组织的生长发育影响着众多的农业性状（如穗粒数，穗粒重，植株高度，分蘖数等）。因此，深入研究水稻顶端分生组织发生及发育的分子机制，将为进一步提升水稻产量和品质提供坚实的理论基础。  
&emsp;单细胞测序是一种新的测序方法，该方法以单个细胞为对象，通过对细胞内基因组、转录组等遗传物质进行扩增测序，进而完成对单个细胞中多种组学信息的检测。为更加精确地解析茎端或根端分生组织中不同类型的细胞在时间和空间梯度上的差异，构建组织内细胞的分化发育轨迹，揭示顶端分生组织不断生长背后所涉及的转录调控机制，本研究以茎端和根端分生组织为对象，采用单细胞转录组和单细胞染色质可及性测序技术开展了多方面的研究。

## 仓库信息 
在该仓库中可以查询到与论文各个章节主要结果相关的生信分析流程与分析程序。由于仓库仍处于开发中，如遇问题可发送邮件至610498675@qq.com。

## 生信分析环境的搭建  
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
### 本章节主要内容概览
```mermaid
        graph LR
                A[scRNA-seq data] --> |cellRanger| C1(数据预处理)
                B[scATAC-seq data] --> |cellRanger ATAC| C2(数据预处理)
                C1(数据预处理) --> |质控| D1(茎尖细胞异质性分析)
                C1(数据预处理) --> |质控| E1(细胞簇标记基因功能分析)
                C1(数据预处理) --> |质控| F1(细胞簇分化轨迹分析) 
                D1(茎尖细胞异质性分析) --> G1(茎尖单细胞转录组图谱的构建)
                E1(细胞簇标记基因功能分析) --> G1(茎尖单细胞转录组图谱的构建)
                F1(细胞簇分化轨迹分析) --> G1(茎尖单细胞转录组图谱的构建)
                G1(茎尖单细胞转录组图谱的构建) --> H1(茎尖表皮细胞亚群转录组分析)
                G1(茎尖单细胞转录组图谱的构建) --> I1(茎尖干细胞亚群异质性分析)
                G1(茎尖单细胞转录组图谱的构建) --> J1(茎尖干细胞发育轨迹分析)
                G1(茎尖单细胞转录组图谱的构建) --> K1(茎尖干细胞基因调控网络的构建)

                C2(数据预处理) --> |质控| D2(茎尖细胞异质性分析)
                D2(茎尖细胞异质性分析) --> E2{细胞类型注释}
                E2{细胞类型注释} --> F2(多组学数据整合)
                G1(茎尖单细胞转录组图谱的构建) --> F2(多组学数据整合)
                E2{细胞类型注释} --> G2(类群特异性开放区间的鉴定)
                E2{细胞类型注释} --> H2(类群特异性表达基因的鉴定)
```
___
### 1.1 单细胞转录组分析结果
#### 1.1.1 数据预处理
本研究使用10x Genomics官方软件Cell Ranger（v.6.0.2）构建单细胞基因表达矩阵。
```
cellranger count --id=scRNA_Rice_SAM --transcriptome=/gss1/home/yanwk/seqlib/cellRanger_genome/rice/cellranger-RNA-IRGSP-genome/Oryza_sativa_IRGSP \
                 --fastqs=/gss1/home/yanwk/ywk_Graduation_Project/00-SingleCell-data/02-Rice/02-Rice-Shoot/01-Mydata \
                 --force-cells=10000   
```
#### 1.1.2 茎尖细胞异质性分析
基于`cellranger`生成的单细胞基因表达矩阵，我们将对其再次进行质量控制与过滤。随后，完成对细胞的聚类和分群。  

[数据质控流程](script/Chapter1/scRNA-analysis/quality-control.r)：**图2-4**(不同质控条件下的数据分布)  
[聚类分群流程](script/Chapter1/scRNA-analysis/cell-clustering.r)：**图2-5**(水稻茎尖单细胞转录组图谱的构建)  

#### 1.1.3 类群标记基因的鉴定
[细胞簇标记基因分析流程](script/Chapter1/scRNA-analysis/DEGs-analysis.r)：**图2-5**(水稻茎尖单细胞转录组图谱的构建)；**图2-6**(水稻茎尖单细胞图谱中不同细胞类群在其功能上的差异性)

#### 1.1.4 细胞簇分化轨迹分析
本小节分析流程主要由三个步骤组成，分别是：（1）[loom文件的生成](script/Chapter1/scRNA-analysis/RNA-velocity-1.sh)；（2）[转录剪切矩阵的生成](script/Chapter1/scRNA-analysis/RNA-velocity-2.r)；（3）[RNA速率的计算](script/Chapter1/scRNA-analysis/RNA-velocity-3.py)。主要分析结果参见**图2-7**(水稻茎端细胞动态发育轨迹的构建)

#### 1.1.5 水稻茎尖表皮细胞亚群转录组分析
本小节的分析流程中包括单细胞亚群的提取、差异表达基因分析、火山图、柱状图的绘制。  
[主要分析流程](script/Chapter1/scRNA-analysis/epidermis-analysis.r)：**图2-9**(表皮细胞亚群间的差异表达分析)

#### 1.1.6 水稻茎尖干细胞亚群异质性分析
本小节分析主要涉及单细胞亚群的提取、再次聚类分群、基因集富集表达分析。  
主要分析流程请点击：[干胞亚群聚类分群流程](script/Chapter1/scRNA-analysis/meristem-analysis-1.r)。  
主要分析结果参见**图2-10**(分生细胞类群的异质性)。

#### 1.1.7 水稻茎尖干细胞亚群发育轨迹和基因调控网络的构建
**Palantir**的使用、调控网路的构建请点击[水稻茎尖干细胞亚群发育轨迹分析](script/Chapter1/scRNA-analysis/meristem-analysis-2.py)和[水稻茎尖干细胞亚群发育轨迹调控网络分析](script/Chapter1/scRNA-analysis/meristem-GRN-analysis-2.r)。主要分析结果参见**图2-11**(茎尖分生细胞分化轨迹)和**图2-12**(茎尖分生细胞分化过程中转录调控网络的构建)。
___
### 1.2 单细胞染色质可及性分析结果
#### 1.2.1 测序数据预处理
```
cellranger-atac count --id SAM --reference /home/wkyan/ywk_lab/04-scATAC-analysis/00-seqlib/IRGSP \
                      --fastqs /home/wkyan/ywk_lab/04-scATAC-analysis/01-raw-data/03-Shoot/data \
                      --sample S_20201215NA --force-cells 8000
```
#### 1.2.2 基于染色质可及性的水稻茎尖细胞异质性分析
本论文中的单细胞染色质可及性数据主要使用Signac软件进行分析，由于该软件的分析流程中需要使用基因组信息，因此需先在R中构建水稻的`BSgenome`包，方便后续的加载与使用。构建流程请点击：[BSgenome的构建](script/Chapter1/scATAC-analysis/BSgenome-create.sh)  
本节分析主要由以下三方面内容组成，详细分析流程如下所示：  
[数据信息统计](script/Chapter1/scATAC-analysis/data-info-analysis.sh)：**图2-13**(水稻茎尖单细胞染色质可及性数据质量)  
[数据质控流程](script/Chapter1/scATAC-analysis/quality-control.r)：**图2-14**(不同质控条件下的数据分布)  
[聚类分群流程](script/Chapter1/scATAC-analysis/cell-clustering.r)：**图2-15**(水稻茎尖单细胞染色质可及性图谱的构建)  

#### 1.2.3 细胞类型的注释
本节分析首先从染色质可及性与基因表达水平之间的相关性入手，通过整合RNA与ATAC两类单细胞数据对，细胞类型进行注释。  
[染色质可及性与基因表达关系分析](script/Chapter1/scATAC-analysis/RNA-ATAC-relation-analysis.sh)：**图2-16A**(水稻茎尖单细胞染色质可及性图谱中细胞类型的注释)  
[注释工作的分析流程](script/Chapter1/scATAC-analysis/scRNA-scATAC-annotation.r)：**图2-16**(水稻茎尖单细胞染色质可及性图谱中细胞类型的注释)；**图2-17**(五个细胞类群在功能上的差异性)  


## 二、水稻根端分生组织单细胞转录组与染色质可及性图谱绘制
### 本章节主要内容概览
```mermaid
        graph LR
                A[两组scRNA-seq数据集] --> |Seurat| C1(数据整合)
                B[scATAC-seq data] --> |cellRanger ATAC| C2(数据预处理)
                C1(数据整合) --> |去除批次效应| D1(水稻根尖细胞异质性研究)
                D1(水稻根尖细胞异质性研究) --> E1(水稻根尖细胞发育图谱的构建)
                D1(水稻根尖细胞异质性研究) --> F1(根尖表皮细胞发育轨迹)
                D1(水稻根尖细胞异质性研究) --> G1(植物激素合成及响应基因在根尖细胞中的表达与分布)

                C2(数据预处理) --> |质控| D2(单细胞多组学数据整合方法间的比较)
                D1(水稻根尖细胞异质性研究) --> D2(单细胞多组学数据整合方法间的比较)
                D2(单细胞多组学数据整合方法间的比较) --> E2(细胞类型的注释)
                E2(细胞类型的注释) --> F2(细胞类型中特异性染色质开放区间的分析)
```
___

### 2.1 单细胞转录组分析结果
#### 2.1.1 水稻根尖细胞异质性研究
通过整合已发表的水稻根尖单细胞转录组数据，构建出更为详细的的根尖单细胞图谱。主要分析内容和流程如下所示：  
[两组水稻根尖单细胞转录组数据的整合流程](script/Chapter2/scRNA-analysis/scRNA-integrated.r)：**图3-1**(水稻根尖单细胞图谱的构建)  
[细胞类群的注释流程](script/Chapter2/scRNA-analysis/scRNA-integrated.r)：**图3-1**(水稻根尖单细胞图谱的构建)

#### 2.1.2 水稻根尖细胞发育图谱的构建
通过使用Monocle3软件构建出表皮，外皮层，皮层，内皮层以及维管细胞的分化发育图谱。详细的分析流程请点击：[根尖细胞发育图谱的构建流程](script/Chapter2/scRNA-analysis/scRNA-Development-analysis.r)。详细分析结果如**图3-2**(水稻根尖细胞的动态分化过程)所示。

#### 2.1.3 根尖表皮细胞发育轨迹
提取表皮细胞亚群，刻画出根毛和非根毛细胞的分化路径。详细分析流程如下所示：  
[根尖表皮细胞发育轨迹分析流程](script/Chapter2/scRNA-analysis/scRNA-epidermal-analysis.r)：**图3-3**(水稻根尖表皮细胞分化轨迹分析)

#### 2.1.4 植物激素合成及响应基因在根尖细胞中的表达与分布
本节主要是观察生长素、细胞分裂素等植物激素合成和响应基因的表达模式。详细分析流程如下所示：  
[植物激素合成及响应基因分析流程](script/Chapter2/scRNA-analysis/scRNA-Phytohormones-analysis.r)：**图3-4**(植物激素合成及响应基因表达模式的UMAP可视化)
___

### 2.2 单细胞染色质可及性分析结果

#### 2.2.1 单细胞多组学数据整合方法间的比较
本节内容主要是使用`Seurat`、`liger`、`scglue`三种方法完成对scRNA-seq与scATAC-seq数据的整合。三种方法的分析流程如下所示：  
- [Seurat分析流程](script/Chapter2/scATAC-analysis/Seurat-analysis.r)：**图3-8**(不同数据整合方法下水稻根尖scRNA与scATAC数据的联合分析)
- [liger分析流程](script/Chapter2/scATAC-analysis/liger-analysis.r)：**图3-8**(不同数据整合方法下水稻根尖scRNA与scATAC数据的联合分析)
- [scglue分析流程](script/Chapter2/scATAC-analysis/glue-analysis.py)：**图3-8**(不同数据整合方法下水稻根尖scRNA与scATAC数据的联合分析)

#### 2.2.2 细胞类型的注释
依据liger的分析结果进一步对细胞的类型进行注释，详细的分析流程请点击：[细胞注释流程](script/Chapter2/scATAC-analysis/cell-annotation.r)。分析结果如**图3-9**(单细胞转录组与单细胞染色质可及性数据联合分析以鉴定细胞类型)。

#### 2.2.3 细胞类型中特异性染色质开放区间的分析
本节通过提取每个细胞簇特异性peak区间，然后对其结合的转录因子进行富集分析，部分分析流程请点击：[转录因子富集分析流程](script/Chapter2/scATAC-analysis/TF-analysis.r)。分析结果如**图3-10**(细胞类型特异性Peaks区间内转录因子结合motifs的富集分析)。

## 三、单细胞精度下茎端与根端分生组织差异性研究
### 本章节主要内容概览
```mermaid
        graph TD
                A1[五组拟南芥根尖scRNA-seq数据整合] --> C1(水稻与拟南芥根尖单细胞转录组比较分析)
                A2[两组水稻根尖scRNA-seq数据整合] --> C1(水稻与拟南芥根尖单细胞转录组比较分析)

                A3[拟南芥茎尖scRNA-seq数据] --> C2(水稻与拟南芥茎尖单细胞转录组比较分析)
                A4[水稻茎尖scRNA-seq数据] --> C2(水稻与拟南芥茎尖单细胞转录组比较分析)

                A5[水稻茎尖scRNA-seq数据] --> C3(水稻茎尖与根尖单细胞转录组比较分析)
                A6[水稻根尖scRNA-seq数据] --> C3(水稻茎尖与根尖单细胞转录组比较分析)
```
___
### 3.1 拟南芥根尖细胞异质性分析
通过收集整理与发表的五组拟南芥根尖单细胞转录组数据，构建出更为细致的细胞图谱。整合结果如**图4-1**(拟南芥根尖单细胞图谱的构建)所示，详细的分析流程请点击：[数据整合流程](script/Chapter3/Ath-integrated-analysis.r)。

### 3.2 水稻与拟南芥根尖单细胞比较转录组分析
分析结果如**图4-2**(单细胞精度下水稻与拟南芥根尖比较分析)所示。  
[物种间根尖单细胞比较转录组分析流程](script/Chapter3/root-analysis.r)

### 3.3 水稻与拟南芥茎尖单细胞比较转录组分析
分析结果如**图4-3**(单细胞精度下水稻与拟南芥茎尖比较分析)所示。  
[物种间茎尖单细胞比较转录组分析流程](script/Chapter3/shoot-analysis.r)

### 3.4 水稻茎尖与根尖单细胞比较转录组分析
分析结果如**图4-4**(单细胞精度下水稻根尖与茎尖比较分析)所示。  
[物种内根尖与茎尖单细胞比较转录组分析流程](script/Chapter3/rice-root-shoot-analysis.r)

## 四、基于深度学习的植物转录因子结合位点预测
模型详细代码信息请访问：https://github.com/wenkaiyan-kevin/PlantBind