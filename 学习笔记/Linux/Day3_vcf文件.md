

# 一、学习目标

- 管道 | 用法
- head
- tail
- vcf 文件
---

### 管道 |

基础语法
```bash
cut -f1 genes.txt | sort | uniq -c
#提取列，排序，统计
```
	3 geneA
	2 geneB
	2 geneC


### head

```bash
head file.txt#默认显示前10行
head -5 file.txt#指定前5行
```

默认显示前10行

```bash
/data/gaodd
```

### tail 同理 head

## VCF文件

全名：variant call format
突变信息文件
用于记录基因组中检测到的突变（SNP INDEL)

基本结构(注释信息+突变数据)
```
##fileformat=VCFv4.2        #文件版本
##source=bcftools           #工具来源
#CHROM POS ID REF ALT QUAL FILTER INFO   #染色体 突变位置 突变编号 参考碱基 突变碱基 质量值 过滤状态 附加信息

chr1  12345  .  A  T  60  PASS  DP=20#测序深度20
chr1  23456  .  G  C  50  PASS  DP=18
chr2  34567  .  T  TA 40  PASS  DP=12
```

### 重测序流程

```
FASTQ（测序数据）
        ↓
比对到参考基因组
        ↓
BAM（比对结果）
        ↓
variant calling
        ↓
VCF（突变结果）#最终得到 varient.vcf
```

### linux处理vcf文件

```bash
head varient.vcf
grep -v "^#" varients.vcf | wc -l #去除开头注释信息,统计突变数量
grep -v "^#" varients.vcf | cut -f1 | sort | uniq -c #统计每条染色体的突变数量

```

# 三、代码示例

```bash
mkdir bioinfo_learning
cd bioinfo_learning
mkdir python
mkdir r
mkdir data
ls
```

---

# 四、我的理解

用自己的话写：

```
pwd 就是查看我现在在哪个目录
ls 就是查看当前目录有什么文件
cd 是进入目录
mkdir 是创建文件夹
```

---

# 五、实践练习

今天实际做过的操作：

```bash
mkdir test_linux
cd test_linux
touch test.txt
ls
```

---

# 六、易错点

```
cd 后面必须是目录
rm 删除文件后不能恢复
Linux 区分大小写
```

---

# 七、扩展知识

例如：

```
ls -lh 可以更容易看文件大小
```

---

# 八、以后可能用到的场景

例如：

```
生信分析时经常需要在服务器创建目录
存放 fastq / bam / vcf 数据
```