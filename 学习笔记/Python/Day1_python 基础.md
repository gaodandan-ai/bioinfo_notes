
## 学习目标

	变量
	字符串
	列表
	for循环
	文件读写
	spilt
	数据筛选
	字典

---

### 变量

```python
gene = "danK"
fc = 2.3
```

---

### 字符串

字符串=一串文本（序列、名字等）

```python
seq = "ACTTATGAG"

print(seq)
print(len(seq))
print(seq(0))   #第一个字符
print(seq[1:4])  #第2~4个字符
```

---

### 列表

```python
genes = ["danA","groEL","clpB"]

print(genes[0])   # 第一个
print(len(genes))
```

```
danA
3
```

---

## for循环

```python
genes = ["dnaK", "groEL", "clpB"]

for gene in genes:
    print("当前基因是:", gene)
```


```python
for i in range(3):
    print("A")

print("B")
```
输出
```
A
A
A
B
```

```python
for i in range(3):
    if i > 1:
        print("A")
    print("B")
```

```
B  #i = 0  不打印A
B
A  #i = 2  打印A
B
```
### （）和 [ ]

[ ] = 索引/查找/取值
（）= 调用函数

```python
genes[0]   #genes不是函数
len(genes)  #len是函数
```
### 换行strip

```python
line.strip()
```

### 拆列spilt

```python
line.spilt("\t")
```

## 文件读写

### 读文件

DEG分析代码

```python
with open("expr.tsv", "r") as f:  #打开文件，叫它f，以后用f操作文件，with:代码结束自动关闭
    header = next(f)   # 跳过表头

    for line in f:                   #对剩下的每一行依次处理
        line = line.strip()          #去掉行末换行符
        fields = line.split("\t")    #按照tab拆成几列

        gene = fields[0]             #第一列：gene名
        logfc = float(fields[1])     #第二列：logFC,转成小数
        pvalue = float(fields[2])    #第三列：pvalue,转成小数

        if logfc > 2 and pvalue < 0.05:#同时满足筛选条件
            print(gene, logfc, pvalue)
```

### 写输出代码

```python
with open("result.tsv", "w") as out:
    out.write("gene\tlogFC\tpvalue\n")

    with open("expr.tsv", "r") as f:
        header = next(f)

        for line in f:
            line = line.strip()
            fields = line.split("\t")

            gene = fields[0]
            logfc = float(fields[1])
            pvalue = float(fields[2])

            if logfc > 2 and pvalue < 0.05:
                out.write(f"{gene}\t{logfc}\t{pvalue}\n")
```


## 练习

打印expr.tsv的每一行gene，只打印第一列
```python 
with open("expr.tsv", "r") as f: 
	header = next(f)
	
	for line in f:                  
        line = line.strip()          
        fields = line.split("\t") 
        
        gene = fields[0]
	    print (gene)
```

修改后
```python
with open("expr.tsv", "r") as f:
    next(f)  # 不用存变量

    for line in f:
        fields = line.strip().split("\t")
        print(fields[0])
```

从表达矩阵中筛选上调基因，然后写到一个新的文件里

expr.tsv

```
gene	logFC	pvalue
dnaK	2.1	0.01
groEL	1.5	0.08
clpB	3.2	0.001
```

筛选条件
logFC >2
pvalue < 0.05

符合条件的基因输出到新文件 up_genes.tsv

```python
with open("expr.tsv", "r") as infile, open("up_genes.tsv", "w") as outfile:
    header = next(f)    #跳过原文件表头  gene logFC  pvalue
    out.write("gene\tlogFC\tpvalue\n")

    for line in f:    #对原文件剩下的每一行进行轮流处理
        line = line.strip()    #去掉每一行后的\n
        fields = line.split("\t")  #按照tab拆列

        gene = fields[0]
        logfc = float(fields[1])
        pvalue = float(fields[2])

        if logfc > 2 and pvalue < 0.05:
            out.write(f"{gene}\t{logfc}\t{pvalue}\n")
```

```python
 out.write(f"{gene}\t{logfc}\t{pvalue}\n")
 
f"..."  #将变量的值塞进字符串中
```

```python
print (f"gene is {gene}")  #推荐

print("gene is" + gene)
```

### 字典 dict

字典=对应关系

```python
gene2fc = {
    "dnaK": 2.1,
    "groEL": 1.5,
    "clpB": 3.2
}
```

```
dnaK  → 2.1
groEL → 1.5
clpB  → 3.2
```

dict=用 key 找到 value

#### 怎么用字典

查值
```python
print(gene2fc["dnaK"])
```

遍历字典

```python
for gene , fc in gene2fc.items():
	print（gene，fc）
```

添加元素

```python
gene2fc["newGene"] = 4.5
```

判断是否存在

```python 
if "dnaK" in gene2fc:
	print("存在")
```

练习：dict，输出logfc >2的基因

```python
gene2fc = {
    "dnaK": 2.1,
    "groEL": 1.5,
    "clpB": 3.2
}
```

```python 
for gene,logfc in gene2fc.items():
	if logfc >2:
		print(gene)

```