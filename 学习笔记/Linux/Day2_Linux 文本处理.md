"# 主题

例如：

Day1 Linux基础命令  
Python 文件读取  
R ggplot绘图  

---

# 一、学习目标

	grep
	cut
	sort
	uniq
	awk


---

# 二、知识点总结


## 1. 命令 / 函数

### grep  查找文本

在文件中查找包含某个字符串的行

```bash
grep "关键词"文件
```

示例输出：

```bash
grep -i "dnaK" genes.txt  #-i忽略大小写
grep -n "dnaK"  genes.txt  #-n显示行号
grep "Genea" Variants.vcf
```

### wc (word count)用于统计文本

```bash
wc -l #统计行数
wc -w #统计单词数
wc -c #统计字符数
```

### cut 提取列

```bash
cut -f列号 文件

cut -f1,3 varent.vcf #选择多列
```

如果是csv文件

```bash
cut -d "," -f1 file.csv
#-d指定分隔符
#-f指定列
```

---

### sort 排序

```bash
sort 文件
```

常用参数

```bash
sort -n file.txt
#数值排序
sort -ke file.txt
#按照第二例排序
```

---

### uniq 去重

uniq只统计相邻重复，需要搭配sort

```bash
uniq file.txt
```

示例
```
geneA
geneA
geneB
geneC
```

```bash
uniq gene.txt
uniq -c gene.txt#-c统计出现次数
```

```
geneA
geneB
geneC
2 geneA
1 geneB
1 geneB
```


### awk

作用：
提取列
筛选数据
统计数据
简单计算

基本语法
```bash
awk '{print $1}' file
#$1第一列  $2第二列
#核心原则：awk '程序' 文件
#''将awk的程序整体包起来，''是shell用的
#{}定义执行的动作，{}是awk用的
```

{}什么时候要
	写动作语句的时候，如print $1
	条件语句就可以不要， 如$2>10
awk的默认条件式print $0，打印整行

实例
	geneA 12
	geneB 20
	geneC 15

```bash
awk "{print $1}" genes.txt
```
	geneA
	geneB
	geneC

```bash
awk '{$2>15}' genes.txt | cut -f1
```
	geneB 20

生信例子
```bash
awk '{sum += $2} END {print sum}' counts.txt
```

计数下调的差异蛋白
```bash
awk -F ',' 'NR>1 && $2<0 {count++} END{print "logFC < 0 下调数量：", count}' 4h_DEG.csv
```
- `-F ','`：按**逗号**分割 csv 文件
- `NR>1`：**跳过第一行表头**（不统计标题行）
- `$2<0`：**第 2 列数值小于 0**（就是你要的 logFC < 0）
- `{count++}`：符合条件就计数 +1
- `END{print ...}`：最后输出总数

## Linux命令组合使用

Linux很强的地方是管道（|），把前一个命令的结果传输给后一个

实例

```bash
cut -f2 gene.txt | sort | uniq
#提取第2列，排序，去重
```

```bash
cut -f1 genes.txt | sort | uniq -c
```
	2 geneA
	1 geneB
	1 geneC



---

# 三、组合示例

```bash
grep -v "^#" file.vcf | wc -l
#grep -v "^#"去掉注释行
#wc -l 统计突变数量
```

---

# 总结

用自己的话写：

```
grep  → 搜索
cut   → 提取列
sort  → 排序
uniq  → 去重
awk   → 文本处理
```

---

