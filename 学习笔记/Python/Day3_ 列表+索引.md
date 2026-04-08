```python
genes = ["dnaK","groEL","clpB"]
#  0 1 2
```

```python
print(genes[0])

dnaK
```

###  列表长度
```python
print len(genes())

3
```

### for+列表
```python
genes = ["dnaK","groEL","clpB"]
	for gene in genes;
	print(gene)

```
### 列表++=

```python
genes = []

genes += ["dnak"]     #= genes.append("dnaK")
genes += ["groEL"]

["dnaK“，"groEL"]
```


## 生信中的用法

### 存筛选到的基因

```python
selected = []

genes = ["dnaK", "groEL", "clpB"]
logfc = [2.1, 1.5, 3.2]

for i in range(3):    #i=0 1 2, 同时访问两个列表genes[i],logfc[i],让数据对齐
    if logfc[i] > 2:
        selected.append(genes[i])

print(selected)
```

输出
```
['dnaK', 'clpB']
```

### 练习1
```python
genes = ["dnaK", "groEL", "clpB"]
```

打印第二个基因
```python
genes = ["dnaK", "groEL", "clpB"]
print (genes[1])
```

用for打印所有的基因

```python 
genes = ["dnaK", "groEL", "clpB"]
for gene in genes:
	print (gene)   #for要定格，print要缩进
```

只打印logFC >2的基因

```python
genes = ["dnaK", "groEL", "clpB"]
logfc = [2.1, 1.5, 3.2]
```

```python
genes = ["dnaK", "groEL", "clpB"]
logfc = [2.1, 1.5, 3.2]

selected=[]

for i in range(3):
	
	if logfc[i] > 2:
        selected.append(genes[i])
        
print(selected)

```

更高级一点的写法

```python
genes = ["dnaK", "groEL", "clpB"]
logfc = [2.1, 1.5, 3.2]
for gene, fc in zip(genes,logfc):
	if fc>2:
		print(gene)
```

