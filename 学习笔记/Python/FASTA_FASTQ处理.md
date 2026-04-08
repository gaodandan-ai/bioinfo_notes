
FASTA文件格式

test.fa
```
>seq1
ATGCGTAC
>seq2
GGGTTTAA
```

### 统计有多少条序列

思路：每次遇到一个>就是一条序列

```python
count = 0

with open("test.fa", "r") as f:
    for line in f:
        line = line.strip()

        if line.startswith(">"):
            count += 1

print("序列数:", count)
```

### FASTA统计长度

```python
seq_id = ""   #当前序列名字
seq = ""      #当前序列碱基

with open("test.fa", "r") as f:
    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if seq_id != "":
                print(seq_id, len(seq))

            seq_id = line[1:]
            seq = ""

        else:
            seq += line

if seq_id != "":
    print(seq_id, len(seq))
```

## FASTQ文件

基本结构
```
@read1       #第一行  read名
ATGCGTAC     #第二行  序列
+            #第三行  +
IIIIIIII     #第四行  质量值
@read2
GGGTTTAA
+
HHHHHHHH
```

### 统计reads数量
思路：每四行=1条read

```python
lines = 0

with open("test.fastq", "r") as f:
    for line in f:
        lines += 1

print("reads数:", lines // 4)
```

### 统计每条read长度

```python
with open("test.fastq", "r") as f:
    for i, line in enumerate(f, start=1):
        line = line.strip()

        if i % 4 == 2:
            print("长度:", len(line))
```