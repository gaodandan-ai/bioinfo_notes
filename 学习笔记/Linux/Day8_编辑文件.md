# 主题

	1.sed
	2.awk
	3.vim
	4.
	重定向
---

## sed 流编辑器

### 在命令行中“查找+替换+编辑”

基本语法

```bash
sed 's/旧名称/新名称' 文件
```
	geneA
	geneB

```bash
sed 's/geneA/geneC' genes.txt
```
	geneC
	geneB

不会修改原文件，只是修改输出结果

保存

```bash
sed -i 's/geneA/geneC' genes.txt   #-i 直接修改文件
```

全部替换
```bash
sed 's/geneA/geneC/g' genes.txt
```

### 删除内容

```bash
sed 'sd' file.txt   #删除第二行
```

删除包含某个词的行

```bash
sed '/geneA/d' file.txt
```

## awk

修改某一列
```bash
awk '{print $1, $2*2}' file.txt #第2列×1
awk '{print $1,$2*2}' file.txt > new file.txt
```

上调或者下调的基因写成新表

```bash
awk -F ',' 'NR==1 || $2>0' G:/thermal_protemic/24h/24h_DEG.csv > G:/thermal_protemic/24h/24h_up.csv
```
- `-F ','`：按逗号分割 CSV
- `NR==1`：保留表头
- `$2>0`：第 2 列是 logFC，大于 0 → 上调
- `$2<0`：第 2 列小于 0 → 下调
- `> 24h_up.csv`：输出新文件

```bash
awk -F ',' 'BEGIN{OFS=","} NR==1 || $2<0 {print $1,$2,$3,$10,$15}' G:/thermal_protemic/24h/24h_DEG.csv > G:/thermal_protemic/24h/24h_down_selected.csv
```
- `-F ','` 按逗号分隔
- `OFS=","` 输出也是逗号分隔
- `print $1,$2,$3,$10,$15` **只输出这 5 列**
- `NR==1` 保留表头
- `$2>0` 只留上调
- `$2<0` 只留下调
## 重定向

```
> 覆盖写入
>>追加写入
```

```bash
echo "hello" > file.txt
echo "hello" >>file.txt
```

## 生信应用

#### 去掉VCF注释并保存

```bash
grep -v "^#"varients.vcf > clean.vcf
```
#### 修改FASTQ.header
```bash
sed "s/>/>chr_" genome.fa > new.fa
```
#### 提取修改表达矩阵
```bash
awk '{print $1, $2/100}' expr.txt > norm.txt
```

## vim

操作流程
	打开文件 -->i -->输入 -->esc -->:wq or :q

进入vim
```bash
vim file.txt
```
输入i 开始编辑
esc退出输入模式
：wq保存并退出
：q不保存退出
：w只保存

```
i      → 输入
Esc    → 退出输入
:wq    → 保存退出
:q!    → 强制退出
dd     → 删除一行
yy     → 复制
p      → 粘贴
/xxx   → 搜索
```

## nano简单文件编辑器

```bash
nano file.txt
```
```
Ctrl + O → 保存
Ctrl + X → 退出
Ctrl + W → 搜索
Ctrl + K → 删除一行
```