### 首先是在R中提取PCA降维
"""
library(Seurat)
scRNA <- readRDS("../scRNA-SMC-Cluster.rds")
DefaultAssay(scRNA) <- 'RNA'

# 输出该细胞类型的count矩阵 cell X gene
write.table(t(as.matrix(GetAssayData(object = scRNA, slot = "counts"))), 
            'SMC_counts.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

# 保存降维的数据
embed_pca <- Embeddings(scRNA,"pca")
write.table(embed_pca, "SMC_embed_pca.csv", sep = ',', row.names = T, col.names = T, quote = F)

"""
### 开始分析内容
# 加载所需的python包
import palantir
import scanpy as sc
import numpy as np
import os
import pandas as pd

# Plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# 加载数据集 ------------------------------------------------------------------
counts = palantir.io.from_csv('SMC_counts.csv')
counts
# 数据的归一化和标准化 ---------------------------------------------------------
norm_df = palantir.preprocess.normalize_counts(counts)
norm_df = palantir.preprocess.log_transform(norm_df)

# 数据的PCA降维 --------------------------------------------------------------
# 💡 PCA reduction, 这里注意需要导入之前Seurat中Harmony的降维数据
# pca_projections, _ = palantir.utils.run_pca(norm_df, use_hvg=False)
pca_projections = pd.read_csv("SMC_embed_pca.csv", header=0, index_col=0)
# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=50)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
# MAGIC imputation：使用MAGIC算法对单细胞的表达数据进行imputation处理
imp_df = palantir.utils.run_magic_imputation(norm_df, dm_res)

# tSNE降维与可视化 ------------------------------------------------------------
tsne = palantir.utils.run_tsne(ms_data)

fig, ax = palantir.plot.plot_tsne(tsne)
plt.savefig('scRNA-SMC-Palantir-TSNE.pdf')

# 绘制特征基因的表达谱 --------------------------------------------------------
# Os03g0155500:OsEXPA25  Os01g0248900:EXPA8  Os10g0578200:CSLD1 
fig = palantir.plot.plot_gene_expression(imp_df, tsne, ['Os08g0532800','Os03g0181500','Os07g0194300','Os02g0814200'])
plt.savefig('scRNA-SMC-CZ-Gene.pdf')

fig = palantir.plot.plot_gene_expression(imp_df, tsne, ['Os12g0274700','Os01g0760900','Os01g0914300','Os10g0554800','Os11g0150400',
              'Os10g0159700','Os01g0638000','Os03g0845000','Os02g0729700'])
plt.savefig('scRNA-SMC-PZ-Gene.pdf')

# 对降维后的数据进行聚类分群 ----------------------------------------------------
# 数据聚类
clusters = palantir.utils.determine_cell_clusters(pca_projections, k=200)
# 聚类结果可视化
fig = palantir.plot.plot_cell_clusters(tsne, clusters)
plt.savefig('scRNA-SMC-Palantir-Cluster.pdf')

# 运行Palantir进行分化发育轨迹推断 -----------------------------------------------
## 可以指定一个近似的最早的起始细胞来运行Palantir。Palantir可以自动确定终末分化状态的细胞，
## 也可以使用termine_states参数指定它们。
## 💡这里可以把之前的TSNE结果给输出出来，然后在EXCEL中进行作图，这样就可以很容易把最前面的cell的名字给找出来
tsne.to_csv("scRNA-SMC-Palantir-TSNE.csv",index=True,sep=',')
# 运行Palantir
start_cell = 'scRNA-SAM_TCAGTCCCAGTAACGG-1'
terminal_states = pd.Series(['CZ', 'PZ'],index=['scRNA-SAM_TCAAGCATCATAAGGA-1','scRNA-SAM_TGAGGTTAGGCTAGCA-1'])

pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500, 
                                    terminal_states=terminal_states.index, use_early_cell_as_start=True)

# 根据已有的生物学知识对终末分化的细胞进行重命名, 这里需要根据自己的具体情况修改
mapping = pd.Series(index=['CZ', 'PZ'])
mapping['CZ'] = "scRNA-SAM_TCAAGCATCATAAGGA-1"
mapping['PZ'] = "scRNA-SAM_TGAGGTTAGGCTAGCA-1"
mapping = pd.Series(mapping.index, index=mapping)
pr_res.branch_probs.columns = mapping[pr_res.branch_probs.columns]
pr_res.branch_probs = pr_res.branch_probs.loc[:, ['CZ', 'PZ']]

# 可视化Palantir的结果
fig = palantir.plot.plot_palantir_results(pr_res, tsne)
plt.savefig('scRNA-SMC-Palantir-Result.pdf')

# 基因表达趋势分析 ------------------------------------------------------------------
# Palantir使用Generalized Additive Models(GAMs)模型计算基因在不同分化细胞中的表达趋势
# 选择特征基因
genes = ['Os08g0532800','Os03g0181500','Os07g0194300','Os02g0814200']
# 计算基因的表达趋势
gene_trends = palantir.presults.compute_gene_trends(pr_res, imp_df.loc[:, genes])
# 基因表达趋势的可视化
fig = palantir.plot.plot_gene_trends(gene_trends)
plt.savefig('scRNA-SMC-Palantir-CZ-GeneTrends.pdf')
# 绘制基因表达趋势热图
fig = palantir.plot.plot_gene_trend_heatmaps(gene_trends)
plt.savefig('palantir_gene_trends_heatmaps.pdf')