	1.原始FASTQ做FastQC
	2.Trim Galore 修剪接头和低质量序列
	3.修剪后的FASTQ再做FastQC，确认接头是否去掉，低质量区域是否改善以及修剪后的数据是否更加干净
	4.用BWA比对到参考基因组
	5.用samtools排序并建立索引
	6.用bcftools做初步变异检测

```bash
#!/bin/bash   #脚本用bash执行
set -euo pipefail

# =========================
# 1. 基本参数
# =========================
raw_dir="/data/gaodd/cg_sequencing/01_raw_data"
trim_dir="/data/gaodd/cg_sequencing/04_clean_data"     #修剪后的clean reads 、、、、 
qc_raw_dir="/data/gaodd/cg_sequencing/03_qc/raw"       #原始数据FastQC后的数据
qc_trim_dir="/data/gaodd/cg_sequencing/03_qc/trimmed"  #修剪后的数据FastQC后的数据
bam_dir="/data/gaodd/cg_sequencing/05_alignment/sorted_bam"  #比对后的文件目录
vcf_dir="/data/gaodd/cg_sequencing/06_variant_calling/raw"   #变异检测后的文件目录

ref="/data/gaodd/cg_sequencing/00_reference/index/genome"
threads=4

# =========================
# 2. 创建输出目录
# =========================
mkdir -p "$trim_dir" "$qc_raw_dir" "$qc_trim_dir" "$bam_dir" "$vcf_dir"

# =========================
# 3. 进入原始数据目录
# =========================
cd "$raw_dir" || exit 1   #||或者，否则，失败就取消

# =========================
# 4. 遍历所有 R1 文件
# =========================
for r1 in *_R1.fastq.gz
do
    sample=${r1%_R1.fastq.gz}    #将后缀_R1.fastq.gz丢掉，获得样本名，bash的字符截断语法，需要用{}括起来
    r2="${sample}_R2.fastq.gz"

    echo "=================================="
    echo "[INFO] 开始处理样本: $sample"

    # 检查 R2 是否存在
    if [ ! -f "$r2" ]; then      #if判断需要加上引号，不然文件名存在一些特殊字符会炸
        echo "[WARN] 跳过 $sample: 找不到配对文件 $r2"
        continue
    fi

    # 定义中间文件和输出文件
    trim_r1="${trim_dir}/${sample}_R1_val_1.fq.gz"
    trim_r2="${trim_dir}/${sample}_R2_val_2.fq.gz"
    bam="${bam_dir}/${sample}.sorted.bam"
    vcf="${vcf_dir}/${sample}.raw.vcf.gz"

    # =========================
    # Step 1. 原始数据质控
    # =========================
    echo "[INFO] Step 1: FastQC 原始数据"
    fastqc -t "$threads" -o "$qc_raw_dir" "$r1" "$r2"

    # =========================
    # Step 2. 质控修剪
    # =========================
    echo "[INFO] Step 2: Trim Galore 修剪"
    trim_galore --paired -o "$trim_dir" "$r1" "$r2"

    # 检查修剪后文件是否生成
    if [ ! -f "$trim_r1" ] || [ ! -f "$trim_r2" ]; then
        echo "[ERROR] $sample 修剪后文件不存在，停止处理该样本"
        continue
    fi

    # =========================
    # Step 3. 修剪后质控
    # =========================
    echo "[INFO] Step 3: FastQC 修剪后数据"
    fastqc -t "$threads" -o "$qc_trim_dir" "$trim_r1" "$trim_r2"

    # =========================
    # Step 4. 比对 + 排序   核心步骤
    # =========================
    echo "[INFO] Step 4: BWA MEM 比对 + samtools sort"
    bwa mem -t "$threads" "$ref" "$trim_r1" "$trim_r2" | \   #第一步，将修剪后的双端reads比对到参考基因组
        samtools sort -@ "$threads" -o "$bam"  #将上游传过来的比对结果排序，输出为BAM文件

    # 建 BAM 索引
    echo "[INFO] Step 5: samtools index"
    samtools index "$bam"

    # =========================
    # Step 5. 变异检测   将比对结果变成变异结果
    # =========================
    echo "[INFO] Step 6: bcftools call"
    bcftools mpileup -Ou -f "$ref" "$bam" | \
        bcftools call -mv -Oz -o "$vcf"

    # 建 VCF 索引
    bcftools index "$vcf"

    echo "[INFO] 样本 $sample 处理完成"
done

echo "=================================="
echo "[INFO] 全部样本处理完成"
```