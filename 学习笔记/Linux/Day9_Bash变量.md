# 核心点
	
	bash变量
	if判断
	批量样本处理
	mini reseq pipeline

## bash变量

```bash
sample="T32"  #将T32存入sample
echo $sample  #打印sample
```

```bash
sample="T32"  #等号两边不能有空格
echo "${sample}_R1.fastq.gz"
echo "${sample}_R2.fastq.gz"   #{sample}比sample更稳
```

常见写法

```bash
ref="genome.fna"
raw_dir="/data/gaodd/cg_sequencing/01_raw_data"
out_dir="/data/gaodd/cg_sequencing/03_qc"
threads=8
```

用变量拼命令

```bash
sample="T32"
raw_dir="/data/raw/"
echo "${raw_dir}/{$sample}_R1.fastq.gz"
```

	/data/raw/T32_R1.fastq.gz

## if

```bash
if [ 条件 ]; then
    命令
fi
```