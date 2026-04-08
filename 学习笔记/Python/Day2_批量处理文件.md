# 主题

	1.input
	2.print
	
---

## 现实问题

```python
sample1.fastq
sample2.fastq
sample3.fastq
```

使用python 处理所有的文件
```python
files = ["sample1.fastq", "sample2.fastq", "sample3.fastq"]

for file in files:
    print("正在处理:", file)
```

### 自动找文件

```python
from pathlib import Path

for file in Path(".").glob("*.fastq"):    #path(.)当前文件夹   glob("*.fastq")找所有的fastq文件
    print(file)
```

### 统计每个fastq的reads数

```python
from pathlib import Path   #从一个叫pathlib的文件夹里拿出path的工具

for file in Path(".").glob("*.fastq"):
    lines = 0    #每个文件重新计数

    with open(file, "r") as f:   #打开当前文件
        for line in f:
            lines += 1

    reads = lines // 4

    print(file, "reads数:", reads)
```


### 统计每个fastq的平均reads长度

```python
from pathlib import Path
for file in Path(".").glob("*.fastq"):
	total_length = 0
	count = 0     #累加所有的序列长度，统计reads数量
    
    with open(file,"r") as f:
	    for i, line in enumerate(f,start = 1):  #给每一行编号
	    line =line.strip()
	   
	   if i%4 == 2:    #每隔4行的第二行=序列
	   total_length += len(line)  
	   count +=1
	
	if count >0:
	avg_length = total_length /count
	 
	 print(file,"平均长度“，avg_length)
```