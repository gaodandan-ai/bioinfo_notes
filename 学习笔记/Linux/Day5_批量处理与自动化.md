

# 一、学习目标

- for循环
- xargs
- 批量简单处理 
---

### for 循环

基础语法
```bash
for 变量 in 列表
do 
   命令
done
```

例子：批量统计FASTQ reads
```bash
for file in *.fastq
do  
  wc -l $file  #$file 当前文件
done
```

### xargs

作用：将一堆结果变成命令参数

```bash
ls *.txt | xargs wc -l  #统计所有文件行数
```

等价于

```bash
wc -l file1.txt file2.txt file3.txt
```


## 实例

```bash
for sample in *.fastq
do 
   fastp -i $sample -o clean_$sample  #fastp：测序数据质控工具（QC）
   # -i 输入文件   -o输出文件
done
#对当前目录的所有fastq文件，进行质量控制处理
```
in
	sample1.fastq
	sample2.fastq
	sample3.fastq

out
	clean_sample1.fastq
	clean_sample2.fastq
	clean_sample3.fastq