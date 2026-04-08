
## pandas处理表格数据（excel级别）

处理表达矩阵/DEG结果

```python
import pandas as pd
df = pd.read_csv("expr.tsv", sep="\t")   # df一个表格
print(df.head())    #看数据

df[df["logFC"]>2]   #自动处理表格，之前写的是for+if
df[(df["logFC"] > 2) & (df["pvalue"] < 0.05)]
df[(df["logFC"] > 2) & (df["pvalue"] < 0.05)] ["gene"]  #只输出gene列
df[abs(df["logFC"])>2]
df["gene"]   #取某一列
df.to_csv("result.tsv",sep="\t",index=false)
```

#### 应用：处理DEG结果

```python
df = pd.read_csv("DEG.tsv", sep="\t")
df_up = de[(df["logfc"]>1) & (df["FDR"] < 0.05)]
```

## numpy 底层计算

	处理矩阵
	数值运算
	被pandas/sklearn调用

```python
import numpy as np

arr = np.array([1,2,3])
print(arr * 2)

[2,4,6]#输出
```

## scipy 统计工具箱

	t-test
	相关性
	分布计算

## biopython 生信

	处理生物序列

FASTA/FASTQ/序列分析

## pysam 处理测序数据 

	BAM
	VCF
	SAM

## seaborn 统计图可视化

## scikit-learn 机器学习

机器学习+降维+聚类

	PCA
	聚类（K-means)
	分类