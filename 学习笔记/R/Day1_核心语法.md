### 赋值 <-
```R
x <- 10
gene <- "geneA"
#把右边的值赋给左边变量
```

```R
gene_matrix <- read.csv("gene_matrix.csv")
#读取文件，赋值给左边变量
```

### function()

结构
```
函数名（参数1，参数2）
```
常见函数
```R
read.csv()#读取文件
write.csv()#写文件
cpm()#标准化
glmFit()#差异分析
```

### 数据框data.frame
```R
df <- dataframe(
gene = c("dnaA","dnaB","dnaC"),
logFC = c(-2,4,4)
)
```
df
	gene logFC
	dnaA -2
	dnaB 4
	dnaC 4

### 取数据
```R
df$gene#取df里面gene这一列
```

```R
gene_list <-df$gene
#将df里边的gene列赋给gene_list
```
gene_list
	dnaA -2
	dnaB 4
	dnaC 4



学会用 ggplot2 画散点图。

---

## 代码

```r
library(ggplot2)

ggplot(df, aes(x=logFC, y=-log10(pvalue))) +
  geom_point()
```

---

## 解释

```
ggplot() 创建图
aes() 设置坐标轴
geom_point() 画散点
```

---

## 练习

```r
ggplot(df, aes(x=gene, y=count)) +
  geom_bar(stat="identity")
```