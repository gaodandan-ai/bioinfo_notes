文件名：nano genes.txt
```
geneA 12
geneB 20
geneC 15
geneA 18
geneB 8
geneD 30
geneC 5
geneA 25
```

# 查看文件内容

```bash
cat nano genes.txt
```

# 提取第一列

```bash
cut -f1 nano genes.txt
```

# 找出所有geneA，统计次数

```bash
grep "geneA" nano gene.txt
grep "geneA" nano gene.txt | wc -l
```

# 打印第一列

```bash
awk '{print $1}' nano gene.txt
```


# 筛选表达量大于15，打印基因名

```bash
awk '$2 >15 [print $1]' nono genes.txt
```

# 统计每个基因出现的次数

```bash
cut -f1 nanom genes.txt | sort | uniq -c
```

# 寻找表达量最高的基因

```bash
sort -k2 -n nano genes.txt |tail -1 
#-k2 按照第二列进行排序
# -n 按照数字进行排序，不然会按照字符串进行排序
#tail 查看文件末尾，tail会显示最后十行，tail -1显示最后一行
#同理，head与tail相反
```

