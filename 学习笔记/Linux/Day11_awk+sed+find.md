# 学习重点

	1. awk
	2. sed 批量修改
	3. find 找文件

## awk

```bash
awk '$2>2' DEG.txt
```

### 生信应用
很多生信文件都是

	tsv: tab
	csv: 逗号

```bash
awk -F "\t" '$3 > 2' DEG.tsv
awk -f "," '$3 > 2' DEG.csv
```


筛选差异基因表达表

```bash
awk -F "\t" '$2>2 && $3<0.05' DEG.tsv
```

打印指定列

```bash
awk '{print $1, $3}' file.txt
awk -F "\t" '{print $1"\t"$3}' file.tsv  # tsv写法
```

改列
```bash
awk '{print $1, $2*2}' file.txt
```

计算reads
```bash
awk '{print $1/4}' file.txt
```

添加表头

```bash
awk 'BEGIN{print "Gene\tlogFC"} {print $1"\t"$2}' file.txt
```

输出文件

```bash
awk '$2 >2 {print $1}' file.txt > up_genes.txt
```

统计行数

```bash
awk 'END{print NR}' file.txt
```

求和
```bash
awk '{sum+=$2} END{print sum}' file.txt
```

平均值
```bash
awk '{sum+=$2} END{print sum/NR}' file.txt
```

分组统计
```bash
awk '{count[$1]++} END{for (i in count) print i, count[i]}' file.txt
```
	统计某基因出现的次数
	统计突变位点

## sed批量修改

```bash
sed 's/old/new/g' file.txt
sed -i 's/old/new/g' file.txt #直接改文件
```

删除某一行
```bash
sed '1d' file.txt  #删掉第一行
```


## find

```bash
find  -name "*.fastq.gz"  #找所有的fastq文件
find . -size +1G  #找大于1G的文件
```

```bash
find -name "*.fastq.gz" -exec fastq {}\; #找到文件后直接跑fastqc
#-exec 对找到的每个文件都执行这个命令
#\;表示结束命令
# factqc file1  fastqc file2  fastqc file3
find . -name "*.fastq.gz" -exec fastqc {} +  #更快,一次把多个文件传递给命令
#fastqc file1 file2 file3
```

## \转义符

意思：让后边的符号不被bash当做特殊符号

```bash
\;# 转义；
My\ file.txt  #等价于
"My file.txt"

#换行
command\
	arg1\
	arg2
#表示下一行还是同一条命令

```