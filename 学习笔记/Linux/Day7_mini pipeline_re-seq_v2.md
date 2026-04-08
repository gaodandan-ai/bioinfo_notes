```bash
#!/bin/bash
set -euo pipefail
shopt -s nullglob   #如果通配符没有匹配到文件，就让它变成空，而不是保留原样

# =========================================================
# Corynebacterium glutamicum resequencing pipeline v4
# 功能：
# 1) raw FASTQ 质控
# 2) Trim Galore 修剪
# 3) trimmed FASTQ 质控
# 4) BWA MEM 比对
# 5) samtools sort + index
# 6) Picard 去重
# 7) bcftools mpileup + call
# 8) 输出 summary 表
# =========================================================

# -----------------------------
# 0. 基本参数
# -----------------------------
RAW_DIR="/data/gaodd/cg_sequencing/01_raw_data"
SCRIPT_DIR="/data/gaodd/cg_sequencing/02_scripts"
REF_FA="/data/gaodd/cg_sequencing/00_reference/genome.fna"
BWA_PREFIX="/data/gaodd/cg_sequencing/00_reference/index/genome"

QC_RAW_DIR="/data/gaodd/cg_sequencing/03_qc/raw"
QC_TRIM_DIR="/data/gaodd/cg_sequencing/03_qc/trimmed"
CLEAN_DIR="/data/gaodd/cg_sequencing/04_clean_data"

ALIGN_DIR="/data/gaodd/cg_sequencing/05_alignment"
SAM_DIR="${ALIGN_DIR}/sam"
SORTED_BAM_DIR="${ALIGN_DIR}/sorted_bam"
DEDUP_BAM_DIR="${ALIGN_DIR}/dedup_bam"

VCF_DIR="/data/gaodd/cg_sequencing/06_variant_calling"
RAW_VCF_DIR="${VCF_DIR}/raw"

LOG_DIR="/data/gaodd/cg_sequencing/08_logs"
SUMMARY_DIR="/data/gaodd/cg_sequencing/09_summary"

THREADS=8

# 工具路径（按你的环境写）
FASTQC="/data/gaodd/miniconda3/envs/cg_seq/bin/fastqc"
TRIM_GALORE="/data/gaodd/miniconda3/envs/cg_seq/bin/trim_galore"
BWA="/data/gaodd/miniconda3/envs/cg_seq/bin/bwa"
SAMTOOLS="/data/gaodd/miniconda3/envs/cg_seq/bin/samtools"
BCFTOOLS="/data/gaodd/miniconda3/envs/cg_seq/bin/bcftools"
PICARD_JAR="/data/gaodd/miniconda3/envs/cg_seq/share/picard-3.4.0-0/picard.jar"

# bcftools 参数
MIN_MAPQ=20
MIN_BASEQ=20
MAX_DEPTH=1000

# summary 文件
SUMMARY_TSV="${SUMMARY_DIR}/pipeline_summary.tsv"

# -----------------------------
# 1. 创建目录
# -----------------------------
mkdir -p \
  "$QC_RAW_DIR" "$QC_TRIM_DIR" "$CLEAN_DIR" \
  "$SAM_DIR" "$SORTED_BAM_DIR" "$DEDUP_BAM_DIR" \
  "$RAW_VCF_DIR" "$LOG_DIR" "$SUMMARY_DIR"

# -----------------------------
# 2. 检查必要文件
# -----------------------------
if [ ! -f "$REF_FA" ]; then
    echo "[ERROR] 参考基因组不存在: $REF_FA" >&2
    exit 1
fi

# 检查 BWA index
if [ ! -f "${BWA_PREFIX}.bwt" ]; then
    echo "[INFO] 未发现 BWA index，开始构建..."
    mkdir -p "$(dirname "$BWA_PREFIX")"
    "$BWA" index -p "$BWA_PREFIX" "$REF_FA"
    echo "[INFO] BWA index 构建完成"
fi

# 检查 summary 文件是否存在，不存在则写表头
if [ ! -f "$SUMMARY_TSV" ]; then
    echo -e "sample\traw_reads\tclean_reads\tbam_reads\tvcf_sites\tstatus" > "$SUMMARY_TSV"
fi

# -----------------------------
# 3. 进入原始数据目录
# -----------------------------
cd "$RAW_DIR" || {
    echo "[ERROR] 无法进入原始数据目录: $RAW_DIR" >&2
    exit 1
}

# -----------------------------
# 4. 遍历样本
# -----------------------------
for R1 in *_R1.fastq.gz *_R1_001.fastq.gz *_1.fastq.gz *_R1.fq.gz *_R1_001.fq.gz *_1.fq.gz
do
    [ -e "$R1" ] || continue

    SAMPLE=""

    case "$R1" in
        *_R1.fastq.gz) SAMPLE="${R1%_R1.fastq.gz}"; R2="${SAMPLE}_R2.fastq.gz" ;;
        *_R1_001.fastq.gz) SAMPLE="${R1%_R1_001.fastq.gz}"; R2="${SAMPLE}_R2_001.fastq.gz" ;;
        *_1.fastq.gz) SAMPLE="${R1%_1.fastq.gz}"; R2="${SAMPLE}_2.fastq.gz" ;;
        *_R1.fq.gz) SAMPLE="${R1%_R1.fq.gz}"; R2="${SAMPLE}_R2.fq.gz" ;;
        *_R1_001.fq.gz) SAMPLE="${R1%_R1_001.fq.gz}"; R2="${SAMPLE}_R2_001.fq.gz" ;;
        *_1.fq.gz) SAMPLE="${R1%_1.fq.gz}"; R2="${SAMPLE}_2.fq.gz" ;;
        *) 
            echo "[WARN] 无法识别命名格式: $R1"
            continue
            ;;
    esac

    SAMPLE_LOG="${LOG_DIR}/${SAMPLE}.log"

    # 输出文件
    TRIM_R1="${CLEAN_DIR}/${SAMPLE}_R1_val_1.fq.gz"
    TRIM_R2="${CLEAN_DIR}/${SAMPLE}_R2_val_2.fq.gz"

    SAM_OUT="${SAM_DIR}/${SAMPLE}.sam"
    SORTED_BAM="${SORTED_BAM_DIR}/${SAMPLE}.sorted.bam"
    DEDUP_BAM="${DEDUP_BAM_DIR}/${SAMPLE}.dedup.bam"
    DEDUP_METRICS="${DEDUP_BAM_DIR}/${SAMPLE}.dedup.metrics.txt"
    VCF_OUT="${RAW_VCF_DIR}/${SAMPLE}.raw.vcf.gz"

    {
        echo "=================================================="
        echo "[INFO] 开始处理样本: $SAMPLE"
        echo "[INFO] R1: $R1"
        echo "[INFO] R2: $R2"
        echo "[INFO] 时间: $(date '+%F %T')"

        # 检查 R2
        if [ ! -f "$R2" ]; then
            echo "[WARN] 缺少配对文件，跳过样本: $SAMPLE"
            echo -e "${SAMPLE}\tNA\tNA\tNA\tNA\tmissing_R2" >> "$SUMMARY_TSV"
            continue
        fi

        # 如果最终结果已存在，直接跳过
        if [ -f "$VCF_OUT" ] && [ -f "${VCF_OUT}.csi" ]; then
            echo "[INFO] VCF 已完成，跳过样本: $SAMPLE"

            RAW_READS=$(( $("${SAMTOOLS}" view -c "$DEDUP_BAM" 2>/dev/null || echo 0) ))
            VCF_SITES=$(( $("${BCFTOOLS}" view -H "$VCF_OUT" 2>/dev/null | wc -l || echo 0) ))

            echo -e "${SAMPLE}\tNA\tNA\t${RAW_READS}\t${VCF_SITES}\tskipped_existing_vcf" >> "$SUMMARY_TSV"
            continue
        fi

        # -----------------------------
        # Step 1. 原始数据 FastQC
        # -----------------------------
        echo "[INFO] Step 1: FastQC(raw)"
        "$FASTQC" -t "$THREADS" -o "$QC_RAW_DIR" "$R1" "$R2"

        # 原始 reads 数量（FASTQ 四行一条）
        RAW_LINES_R1=$(zcat "$R1" | wc -l)
        RAW_READS=$(( RAW_LINES_R1 / 4 ))
        echo "[INFO] raw_reads(R1): $RAW_READS"

        # -----------------------------
        # Step 2. Trim Galore
        # -----------------------------
        echo "[INFO] Step 2: Trim Galore"
        "$TRIM_GALORE" \
            --paired \
            --cores "$THREADS" \
            -o "$CLEAN_DIR" \
            "$R1" "$R2"

        if [ ! -f "$TRIM_R1" ] || [ ! -f "$TRIM_R2" ]; then
            echo "[ERROR] Trim 后文件不存在，样本失败: $SAMPLE"
            echo -e "${SAMPLE}\t${RAW_READS}\tNA\tNA\tNA\ttrim_failed" >> "$SUMMARY_TSV"
            continue
        fi

        CLEAN_LINES_R1=$(zcat "$TRIM_R1" | wc -l)
        CLEAN_READS=$(( CLEAN_LINES_R1 / 4 ))
        echo "[INFO] clean_reads(R1): $CLEAN_READS"

        # -----------------------------
        # Step 3. 修剪后 FastQC
        # -----------------------------
        echo "[INFO] Step 3: FastQC(trimmed)"
        "$FASTQC" -t "$THREADS" -o "$QC_TRIM_DIR" "$TRIM_R1" "$TRIM_R2"

        # -----------------------------
        # Step 4. BWA MEM 比对
        # -----------------------------
        echo "[INFO] Step 4: BWA MEM"
        "$BWA" mem \
            -t "$THREADS" \
            -M \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
            "$BWA_PREFIX" \
            "$TRIM_R1" "$TRIM_R2" \
            > "$SAM_OUT"

        # -----------------------------
        # Step 5. samtools sort + index
        # -----------------------------
        echo "[INFO] Step 5: samtools sort/index"
        "$SAMTOOLS" sort -@ "$THREADS" -o "$SORTED_BAM" "$SAM_OUT"
        "$SAMTOOLS" index "$SORTED_BAM"

        # 这里把大的 sam 删掉，节省空间
        rm -f "$SAM_OUT"

        # -----------------------------
        # Step 6. Picard MarkDuplicates
        # -----------------------------
        echo "[INFO] Step 6: Picard MarkDuplicates"
        java -jar "$PICARD_JAR" MarkDuplicates \
            I="$SORTED_BAM" \
            O="$DEDUP_BAM" \
            M="$DEDUP_METRICS" \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=SILENT

        "$SAMTOOLS" index "$DEDUP_BAM"

        BAM_READS=$("$SAMTOOLS" view -c "$DEDUP_BAM")
        echo "[INFO] dedup_bam_reads: $BAM_READS"

        # -----------------------------
        # Step 7. bcftools mpileup + call
        # -----------------------------
        echo "[INFO] Step 7: bcftools mpileup/call"
        "$BCFTOOLS" mpileup \
            -Ou \
            -f "$REF_FA" \
            -a FORMAT/DP,FORMAT/AD \
            -d "$MAX_DEPTH" \
            -q "$MIN_MAPQ" \
            -Q "$MIN_BASEQ" \
            -C 50 \
            "$DEDUP_BAM" | \
        "$BCFTOOLS" call \
            -mv \
            --ploidy 1 \
            -Oz \
            -o "$VCF_OUT"

        "$BCFTOOLS" index "$VCF_OUT"

        VCF_SITES=$("$BCFTOOLS" view -H "$VCF_OUT" | wc -l)
        echo "[INFO] vcf_sites: $VCF_SITES"

        # -----------------------------
        # Step 8. summary
        # -----------------------------
        echo -e "${SAMPLE}\t${RAW_READS}\t${CLEAN_READS}\t${BAM_READS}\t${VCF_SITES}\tsuccess" >> "$SUMMARY_TSV"
        echo "[INFO] 样本完成: $SAMPLE"

    } > "$SAMPLE_LOG" 2>&1

    echo "[INFO] $SAMPLE 处理结束，日志: $SAMPLE_LOG"
done

echo "[INFO] 全部样本处理完成"
echo "[INFO] Summary: $SUMMARY_TSV"
```