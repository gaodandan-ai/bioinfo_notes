## 1.Preprocessing of Raw Data

确保原始数据的测序质量

### 评估读数质量

FastQC 是一种针对高通量测序数据（例如 FASTQ 文件）的质量控制工具，可用于生成有关原始测序读数整体质量的报告。

```python
fastqc sample.fastq
```
### 修剪低质量读数

通常使用Trimmomatic 或者 cutadapt等工具删除低质量碱基和接头序列，确保下游中仅适用高质量读数

```python
cutadapt -q 20 -o trimmed_output.fastq sample.fastq
```
### 将读数映射到参考基因组

将原始读数与参考基因组或者转录组进行比对，将原始核苷酸序列映射到已知的基因位置

#### 第一步：索引参考基因组
#### 第二步：比对读数

比对过程中生成SAM或者BAM文件，其中储存的是映射的读数与其比对的位置

#### 第三步：生成基因表达矩阵

## 用于预处理和比对的python工具

### snapy

```python
import scanpy as sc  
# Read 10x Genomics data 
adata = sc.read_10x_mtx('path_to_data_directory/', var_  
names='gene_symbols', cache=True, prefix='your_file_name_prefix')
```
### anndata
### HTSeq

## Quality Control and Filtering of RNA-seq Data

### Normalization and Scaling of Gene Expression Data

标准化调整细胞间的基因表达值，以考虑测序深度和捕获的差异
使用酶百万技术（CPM）或者对数标准化来调整这些差异

## Dimensionality

目的：保留重要信息，同时简化分析
常用方法：PCA/t-SNE/UMAP

### PCA
将高维基因表达数据转换为较少数量的主成分，每个主成分捕获方差最大的几个方向，前几个成分能解释大部分方差

