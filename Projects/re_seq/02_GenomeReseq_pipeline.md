输入：原始双端测序数据
输出：
- 质控报告
- 修剪后的 clean reads
- 比对后的 BAM
- 去重后的 BAM
- 检测到的突变 VCF
- 按 AF 分组的突变结果
- snpEff 注释后的突变文件
- 最后统计报告

流程
```
原始FASTQ
→ FastQC质控
→ Trim Galore修剪
→ 修剪后再质控
→ BWA比对到参考基因组
→ samtools排序
→ Picard去重
→ bcftools做变异检测
→ 基础质量过滤
→ 按AF拆分为多态突变 / 高频突变
→ snpEff注释
→ bcftools stats统计
```

```bash
#!/bin/bash
# ============================================================================
# 谷氨酸棒杆菌重测序分析脚本（双输出版，升级稳健版）
# 支持原始双端命名：
#   1) xxx_1.fastq.gz / xxx_2.fastq.gz
#   2) xxx_R1.fastq.gz / xxx_R2.fastq.gz
#   3) xxx_R1_001.fastq.gz / xxx_R2_001.fastq.gz
#   4) xxx_1.fq.gz / xxx_2.fq.gz
#   5) xxx_R1.fq.gz / xxx_R2.fq.gz
#   6) xxx_R1_001.fq.gz / xxx_R2_001.fq.gz
#
# 流程：
# 原始质控 → 修剪 → 修剪后质控 → 比对 → 排序 → 去重 → 变异检测
# → 基础质量过滤（不按AF硬过滤）
# → 按AF拆分为 polymorphic / highfreq
# → 两套结果分别注释与统计
# ============================================================================

set -euo pipefail   
#set -e只要某一步命令报错，脚本立刻停止 set -u只要用了一个没被定义的变量，就报错
#set -o pipefail 对管道命名更加严格，A|B只要有一个报错就停止
shopt -s nullglob

# -------------------------- 1. 环境与工具配置 --------------------------
export JAVA_OPTS="-Xmx16g -Xms4g"   #snpeff是JAVA程序，给他设置内存
export SNPEFF_CONFIG="/data/gaodd/miniconda3/envs/cg_seq/snpEff/snpEff.config"

bwa="/data/gaodd/miniconda3/envs/cg_seq/bin/bwa"
trim_galore="/data/gaodd/miniconda3/envs/cg_seq/bin/trim_galore"
fastqc="/data/gaodd/miniconda3/envs/cg_seq/bin/fastqc"
multiqc="/data/gaodd/miniconda3/envs/cg_seq/bin/multiqc"
samtools="/data/gaodd/miniconda3/envs/cg_seq/bin/samtools"
bcftools="/data/gaodd/miniconda3/envs/cg_seq/bin/bcftools"
cutadapt="/data/gaodd/miniconda3/envs/cg_seq/bin/cutadapt"
snpeff="/data/gaodd/miniconda3/envs/cg_seq/bin/snpEff"
picard="/data/gaodd/miniconda3/envs/cg_seq/bin/picard"
python3_bin="/data/gaodd/miniconda3/envs/cg_seq/bin/python3"

# -------------------------- 2. 参数设置 --------------------------
reference="/data/gaodd/cg_sequencing/00_reference/genome.fna"
bwa_index="/data/gaodd/cg_sequencing/00_reference/index/genome"

SNPEFF_DATA_DIR="/data/gaodd/cg_sequencing/07_annotation/"
CG_DB_NAME="Corynebacterium_glutamicum_ATCC13032"

fastq_dir="/data/gaodd/cg_sequencing/01_raw_data/"
fastqc_raw="/data/gaodd/cg_sequencing/03_qc/raw/"
fastqc_clean="/data/gaodd/cg_sequencing/03_qc/trimmed/"
trimgalore_dir="/data/gaodd/cg_sequencing/04_clean_data/"
alignment="/data/gaodd/cg_sequencing/05_alignment/sam/"
sorted="/data/gaodd/cg_sequencing/05_alignment/sorted_bam/"
dedup="/data/gaodd/cg_sequencing/05_alignment/dedup_bam/"

variant_raw="/data/gaodd/cg_sequencing/06_variant_calling/raw/"
variant_qc="/data/gaodd/cg_sequencing/06_variant_calling/quality_filtered/"
variant_poly="/data/gaodd/cg_sequencing/06_variant_calling/polymorphic/"
variant_high="/data/gaodd/cg_sequencing/06_variant_calling/highfreq/"

variant_poly_annot="/data/gaodd/cg_sequencing/06_variant_calling/annotated/polymorphic/"
variant_high_annot="/data/gaodd/cg_sequencing/06_variant_calling/annotated/highfreq/"

variant_tables="/data/gaodd/cg_sequencing/06_variant_calling/tables/"
helper_dir="/data/gaodd/cg_sequencing/02_scripts/helpers/"
log_dir="/data/gaodd/cg_sequencing/08_logs/sequencing_analysis/"

# 阈值参数
MIN_DP_CALL=10
MIN_DP_AF=20
MIN_AF_POLY=0.01
MAX_AF_POLY=0.99
MIN_AF_HIGH=0.80
MIN_QUAL=30

# -------------------------- 3. 初始化 --------------------------
echo "=== 初始化：创建目录 & 检查依赖 ==="

mkdir -p \
  "$fastqc_raw" "$fastqc_clean" "$trimgalore_dir" \
  "$alignment" "$sorted" "$dedup" \
  "$variant_raw" "$variant_qc" "$variant_poly" "$variant_high" \
  "$variant_poly_annot" "$variant_high_annot" \
  "$variant_tables" "$helper_dir" "$log_dir"

tools=("$bwa" "$trim_galore" "$fastqc" "$multiqc" "$samtools" "$bcftools" "$cutadapt" "$snpeff" "$picard" "$python3_bin")
for tool in "${tools[@]}"; do
    if [ ! -e "$tool" ]; then
        echo "❌ 错误：工具不存在 → $tool"
        exit 1
    fi
done

if [ ! -f "$reference" ]; then
    echo "❌ 错误：参考基因组不存在 → $reference"
    exit 1
fi

if [ ! -f "${bwa_index}.bwt" ]; then
    echo "⚠️ BWA索引不存在，自动构建..."
    "$bwa" index -p "$bwa_index" "$reference"
fi

if [ ! -d "${SNPEFF_DATA_DIR}/${CG_DB_NAME}" ] || [ ! -f "${SNPEFF_DATA_DIR}/${CG_DB_NAME}/snpEffectPredictor.bin" ]; then
    echo "❌ 错误：snpEff自定义数据库不存在！"
    echo "请先构建数据库："
    echo "$snpeff build -Xmx16g -gff3 -v -noCheckCds -noCheckProtein -dataDir $SNPEFF_DATA_DIR $CG_DB_NAME"
    exit 1
else
    echo "✅ snpEff数据库检查通过：$CG_DB_NAME"
fi

if find "$fastq_dir" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | grep -q .; then
    fastq_n=$(find "$fastq_dir" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | wc -l)
    echo "✅ 原始数据检查通过：找到 ${fastq_n} 个 FASTQ 文件"
else
    echo "❌ 错误：$fastq_dir 下没有 .fastq.gz 或 .fq.gz 文件"
    exit 1
fi

# -------------------------- 4. 写入辅助Python脚本 --------------------------

# 4.1 按 AF 拆分 VCF
split_py="${helper_dir}/split_vcf_by_af.py"

cat > "$split_py" << 'PYEOF'
#!/usr/bin/env python3
from __future__ import annotations
import argparse
import gzip
from pathlib import Path

def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def parse_info(info_str):
    d = {}
    if info_str == "." or not info_str:
        return d
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = "1"
    return d

def safe_int(x):
    if x is None:
        return None
    try:
        return int(float(x))
    except Exception:
        return None

def pick_ao_dp(info, fmt_keys, sample_fields, alt_index=0):
    fmt = dict(zip(fmt_keys, sample_fields)) if fmt_keys and sample_fields else {}
    ao = None
    dp = None

    # FORMAT/DP
    if "DP" in fmt:
        dp = safe_int(fmt.get("DP"))

    # FORMAT/AO
    if "AO" in fmt:
        vals = fmt.get("AO", "").split(",")
        if alt_index < len(vals):
            ao = safe_int(vals[alt_index])

    # FORMAT/AD = ref,alt1,alt2...
    if (ao is None or dp is None) and "AD" in fmt:
        vals = [safe_int(x) for x in fmt.get("AD", "").split(",")]
        if len(vals) > 1 + alt_index:
            ao = vals[1 + alt_index]
        if dp is None:
            valid = [x for x in vals if x is not None]
            if valid:
                dp = sum(valid)

    # INFO fallback
    if dp is None and "DP" in info:
        dp = safe_int(info.get("DP"))
    if ao is None and "AO" in info:
        vals = info.get("AO", "").split(",")
        if alt_index < len(vals):
            ao = safe_int(vals[alt_index])

    return ao, dp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--poly_sites", required=True)
    ap.add_argument("--high_sites", required=True)
    ap.add_argument("--poly_tsv", required=True)
    ap.add_argument("--high_tsv", required=True)
    ap.add_argument("--min_dp", type=int, default=20)
    ap.add_argument("--min_af_poly", type=float, default=0.01)
    ap.add_argument("--max_af_poly", type=float, default=0.99)
    ap.add_argument("--min_af_high", type=float, default=0.8)
    args = ap.parse_args()

    poly_sites = []
    high_sites = []
    poly_rows = []
    high_rows = []

    with open_text(Path(args.vcf)) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, _id, ref, alts, qual, flt, info_str = cols[:8]
            fmt_keys = []
            sample_fields = []

            if len(cols) >= 10:
                fmt_keys = cols[8].split(":")
                sample_fields = cols[9].split(":")

            info = parse_info(info_str)

            for alt_index, alt in enumerate(alts.split(",")):
                ao, dp = pick_ao_dp(info, fmt_keys, sample_fields, alt_index)

                if dp is None or dp < args.min_dp:
                    continue

                if ao is None:
                    continue

                af = ao / dp if dp > 0 else 0.0
                site = f"{chrom}\t{pos}\n"
                row = [chrom, pos, ref, alt, dp, ao, f"{float(af):.6f}"]

                if args.min_af_poly <= af < args.max_af_poly:
                    poly_sites.append(site)
                    poly_rows.append(row)

                if af >= args.min_af_high:
                    high_sites.append(site)
                    high_rows.append(row)

    def write_unique_lines(path, lines):
        seen = set()
        with open(path, "w") as out:
            for x in lines:
                if x not in seen:
                    out.write(x)
                    seen.add(x)

    write_unique_lines(args.poly_sites, poly_sites)
    write_unique_lines(args.high_sites, high_sites)

    header = "CHROM\tPOS\tREF\tALT\tDP\tAO\tAF\n"
    with open(args.poly_tsv, "w") as out:
        out.write(header)
        for r in poly_rows:
            out.write("\t".join(map(str, r)) + "\n")

    with open(args.high_tsv, "w") as out:
        out.write(header)
        for r in high_rows:
            out.write("\t".join(map(str, r)) + "\n")

if __name__ == "__main__":
    main()
PYEOF

chmod +x "$split_py"

# 4.2 识别多种 FASTQ 命名并输出样本配对表
pair_py="${helper_dir}/detect_fastq_pairs.py"

cat > "$pair_py" << 'PYEOF'
#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import re
import sys

PATTERNS = [
    (re.compile(r"^(?P<sample>.+)_1\.(fastq|fq)\.gz$"), "R1"),
    (re.compile(r"^(?P<sample>.+)_2\.(fastq|fq)\.gz$"), "R2"),
    (re.compile(r"^(?P<sample>.+)_R1\.(fastq|fq)\.gz$"), "R1"),
    (re.compile(r"^(?P<sample>.+)_R2\.(fastq|fq)\.gz$"), "R2"),
    (re.compile(r"^(?P<sample>.+)_R1_001\.(fastq|fq)\.gz$"), "R1"),
    (re.compile(r"^(?P<sample>.+)_R2_001\.(fastq|fq)\.gz$"), "R2"),
]

def detect_pairs(indir: Path):
    samples = {}
    unmatched = []

    files = sorted(set(list(indir.glob("*.fastq.gz")) + list(indir.glob("*.fq.gz"))))

    for p in files:
        name = p.name
        matched = False
        for rx, read_tag in PATTERNS:
            m = rx.match(name)
            if m:
                sample = m.group("sample")
                samples.setdefault(sample, {})
                if read_tag in samples[sample]:
                    print(f"ERROR: duplicated {read_tag} for sample {sample}: {name}", file=sys.stderr)
                    sys.exit(1)
                samples[sample][read_tag] = str(p.resolve())
                matched = True
                break
        if not matched:
            unmatched.append(name)

    good = []
    bad = []

    for sample, d in sorted(samples.items()):
        r1 = d.get("R1")
        r2 = d.get("R2")
        if r1 and r2:
            good.append((sample, r1, r2))
        else:
            bad.append((sample, r1, r2))

    return good, bad, unmatched

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    good, bad, unmatched = detect_pairs(Path(args.indir))

    if unmatched:
        print("WARNING: unmatched FASTQ files:", file=sys.stderr)
        for x in unmatched:
            print(f"  {x}", file=sys.stderr)

    if bad:
        print("ERROR: some samples do not have complete R1/R2 pairs:", file=sys.stderr)
        for sample, r1, r2 in bad:
            print(f"  {sample}\tR1={r1}\tR2={r2}", file=sys.stderr)
        sys.exit(1)

    if not good:
        print("ERROR: no paired FASTQ files detected.", file=sys.stderr)
        sys.exit(1)

    with open(args.out, "w") as out:
        out.write("sample\tfq1\tfq2\n")
        for sample, fq1, fq2 in good:
            out.write(f"{sample}\t{fq1}\t{fq2}\n")

if __name__ == "__main__":
    main()
PYEOF

chmod +x "$pair_py"

# -------------------------- 5. 原始数据FastQC --------------------------
echo -e "\n=== Step 1: 原始数据FastQC质控 ==="
mapfile -t raw_fastqs < <(find "$fastq_dir" -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)

if [ "${#raw_fastqs[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到原始FASTQ文件"
    exit 1
fi

"$fastqc" -t 4 -o "$fastqc_raw" "${raw_fastqs[@]}"
"$multiqc" -o "$fastqc_raw" "$fastqc_raw"
echo "✅ 原始数据质控完成"

# -------------------------- 6. 识别双端配对 --------------------------
echo -e "\n=== Step 2: 识别原始FASTQ双端配对 ==="
pair_tsv="${trimgalore_dir}/fastq_pairs.tsv"
"$python3_bin" "$pair_py" --indir "$fastq_dir" --out "$pair_tsv"

sample_count=$(tail -n +2 "$pair_tsv" | wc -l)
echo "✅ 已识别到 $sample_count 个双端样本"

# -------------------------- 7. Trim Galore --------------------------
echo -e "\n=== Step 3: Trim Galore! 修剪 ==="

tail -n +2 "$pair_tsv" | while IFS=$'\t' read -r sample fq1 fq2; do
    echo "正在修剪样本：$sample"

    "$trim_galore" -q 25 --phred33 --length 50 \
        --paired --cores 4 --path_to_cutadapt "$cutadapt" \
        -o "$trimgalore_dir" "$fq1" "$fq2"
done

echo "✅ 修剪完成：共处理 $sample_count 个样本"

# -------------------------- 8. 修剪后FastQC --------------------------
echo -e "\n=== Step 4: 修剪后FastQC质控 ==="
trimmed_fastqs=( "$trimgalore_dir"/*val*.fq.gz )

if [ "${#trimmed_fastqs[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到修剪后的 fq.gz 文件"
    exit 1
fi

"$fastqc" -t 4 -o "$fastqc_clean" "${trimmed_fastqs[@]}"
"$multiqc" -o "$fastqc_clean" "$fastqc_clean"
echo "✅ 修剪后质控完成"

# -------------------------- 9. BWA比对 --------------------------
echo -e "\n=== Step 5: BWA-MEM比对 ==="

tail -n +2 "$pair_tsv" | while IFS=$'\t' read -r sample fq1 fq2; do
    trimmed_fq1=""
    trimmed_fq2=""

    for candidate in \
        "${trimgalore_dir}/${sample}_1_val_1.fq.gz" \
        "${trimgalore_dir}/${sample}_R1_val_1.fq.gz" \
        "${trimgalore_dir}/${sample}_R1_001_val_1.fq.gz"
    do
        if [ -f "$candidate" ]; then
            trimmed_fq1="$candidate"
            break
        fi
    done

    for candidate in \
        "${trimgalore_dir}/${sample}_2_val_2.fq.gz" \
        "${trimgalore_dir}/${sample}_R2_val_2.fq.gz" \
        "${trimgalore_dir}/${sample}_R2_001_val_2.fq.gz"
    do
        if [ -f "$candidate" ]; then
            trimmed_fq2="$candidate"
            break
        fi
    done

    if [ -z "$trimmed_fq1" ] || [ -z "$trimmed_fq2" ]; then
        echo "❌ 错误：未找到修剪后的双端文件：$sample"
        echo "R1=$trimmed_fq1"
        echo "R2=$trimmed_fq2"
        exit 1
    fi

    sam_out="${alignment}${sample}.sam"
    log_out="${alignment}${sample}_bwa.log"

    echo "正在比对样本：$sample"

    "$bwa" mem -t 4 -M -T 30 \
        -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}_lib" \
        "$bwa_index" "$trimmed_fq1" "$trimmed_fq2" > "$sam_out" 2> "$log_out"
done

bwa_logs=( "$alignment"/*.log )
if [ "${#bwa_logs[@]}" -gt 0 ]; then
    "$multiqc" -o "$alignment" "$alignment"/*.log || true
fi
echo "✅ 比对完成"

# -------------------------- 10. SAM转BAM并排序 --------------------------
echo -e "\n=== Step 6: SAM转BAM并排序 ==="
cd "$alignment"

sam_files=( *.sam )
if [ "${#sam_files[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到 SAM 文件"
    exit 1
fi

for sam in "${sam_files[@]}"; do
    [ -f "$sam" ] || continue
    sample=$(basename "$sam" ".sam")
    sorted_bam="${sorted}${sample}.sorted.bam"

    "$samtools" sort -O bam -@ 4 -m 1G -o "$sorted_bam" "$sam"
    "$samtools" index "$sorted_bam"
    rm -f "$sam"
done

echo "✅ 排序完成"

# -------------------------- 11. 去重 --------------------------
echo -e "\n=== Step 7: Picard去重 ==="
cd "$sorted"

sorted_bams=( *.sorted.bam )
if [ "${#sorted_bams[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到 sorted BAM 文件"
    exit 1
fi

for bam in "${sorted_bams[@]}"; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" ".sorted.bam")
    dedup_bam="${dedup}${sample}.dedup.bam"
    dedup_metrics="${dedup}${sample}_dup_metrics.txt"

    "$picard" MarkDuplicates \
        I="$bam" \
        O="$dedup_bam" \
        M="$dedup_metrics" \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

    "$samtools" index "$dedup_bam"
done

echo "✅ 去重完成"

# -------------------------- 12. 变异检测 --------------------------
echo -e "\n=== Step 8: 变异检测 ==="
cd "$dedup"

dedup_bams=( *.dedup.bam )
if [ "${#dedup_bams[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到 dedup BAM 文件"
    exit 1
fi

for bam in "${dedup_bams[@]}"; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" ".dedup.bam")
    raw_vcf_gz="${variant_raw}${sample}.raw.vcf.gz"
    mpileup_log="${variant_raw}${sample}_mpileup.log"

    echo "正在检测变异：$sample"

    "$bcftools" mpileup \
        -a FORMAT/DP,FORMAT/AD \
        -f "$reference" \
        -d 1000 -q 20 -Q 20 -C 50 \
        "$bam" 2> "$mpileup_log" | \
    "$bcftools" call -mv -O z --ploidy 1 -o "$raw_vcf_gz"

    "$bcftools" index -f "$raw_vcf_gz"
done

echo "✅ 变异检测完成：$variant_raw"

# -------------------------- 13. 基础质量过滤 --------------------------
echo -e "\n=== Step 9: 基础质量过滤（保留多态突变） ==="
cd "$variant_raw"

raw_vcfs=( *.raw.vcf.gz )
if [ "${#raw_vcfs[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到 raw VCF 文件"
    exit 1
fi

for vcf_gz in "${raw_vcfs[@]}"; do
    [ -f "$vcf_gz" ] || continue
    sample=$(basename "$vcf_gz" ".raw.vcf.gz")
    qc_vcf_gz="${variant_qc}${sample}.qc.vcf.gz"

    echo "正在基础过滤：$sample"

    "$bcftools" filter \
        -e "QUAL<${MIN_QUAL} || FORMAT/DP<${MIN_DP_CALL}" \
        -O z -o "$qc_vcf_gz" "$vcf_gz"

    if [ -s "$qc_vcf_gz" ]; then
        "$bcftools" index -f "$qc_vcf_gz"
    else
        echo "⚠️ 过滤后为空：$sample"
    fi
done

echo "✅ 基础质量过滤完成：$variant_qc"

# -------------------------- 14. 按AF拆分：polymorphic / highfreq --------------------------
echo -e "\n=== Step 10: 按AF拆分为多态突变和高频突变 ==="
cd "$variant_qc"

qc_vcfs=( *.qc.vcf.gz )
if [ "${#qc_vcfs[@]}" -eq 0 ]; then
    echo "❌ 错误：未找到 qc VCF 文件"
    exit 1
fi

for vcf_gz in "${qc_vcfs[@]}"; do
    [ -f "$vcf_gz" ] || continue
    sample=$(basename "$vcf_gz" ".qc.vcf.gz")

    poly_sites="${variant_tables}${sample}.polymorphic.sites.tsv"
    high_sites="${variant_tables}${sample}.highfreq.sites.tsv"
    poly_tsv="${variant_tables}${sample}.polymorphic.table.tsv"
    high_tsv="${variant_tables}${sample}.highfreq.table.tsv"

    poly_vcf_gz="${variant_poly}${sample}.polymorphic.vcf.gz"
    high_vcf_gz="${variant_high}${sample}.highfreq.vcf.gz"

    echo "正在AF拆分：$sample"

    "$python3_bin" "$split_py" \
        --vcf "$vcf_gz" \
        --poly_sites "$poly_sites" \
        --high_sites "$high_sites" \
        --poly_tsv "$poly_tsv" \
        --high_tsv "$high_tsv" \
        --min_dp "$MIN_DP_AF" \
        --min_af_poly "$MIN_AF_POLY" \
        --max_af_poly "$MAX_AF_POLY" \
        --min_af_high "$MIN_AF_HIGH"

    if [ -s "$poly_sites" ]; then
        "$bcftools" view -R "$poly_sites" -O z -o "$poly_vcf_gz" "$vcf_gz"
        "$bcftools" index -f "$poly_vcf_gz"
        echo "✅ polymorphic 结果已输出：$sample"
    else
        echo "⚠️ polymorphic 结果为空：$sample"
    fi

    if [ -s "$high_sites" ]; then
        "$bcftools" view -R "$high_sites" -O z -o "$high_vcf_gz" "$vcf_gz"
        "$bcftools" index -f "$high_vcf_gz"
        echo "✅ highfreq 结果已输出：$sample"
    else
        echo "⚠️ highfreq 结果为空：$sample"
    fi
done

echo "✅ AF拆分完成"

# -------------------------- 15. 注释函数 --------------------------
annotate_dir () {
    local input_dir="$1"
    local output_dir="$2"
    local label="$3"

    echo -e "\n=== Step 11 (${label}): snpEff注释 ==="
    cd "$input_dir"

    local vcfs=( *.vcf.gz )
    if [ "${#vcfs[@]}" -eq 0 ]; then
        echo "⚠️ ${label} 目录下没有 VCF 文件，跳过注释"
        return 0
    fi

    for vcf_gz in "${vcfs[@]}"; do
        [ -f "$vcf_gz" ] || continue
        sample=$(basename "$vcf_gz" ".vcf.gz")

        temp_vcf="${output_dir}${sample}.tmp.vcf"
        ann_vcf_gz="${output_dir}${sample}.annotated.vcf.gz"
        ann_log="${output_dir}${sample}_snpEff.log"
        stats_html="${output_dir}${sample}_snpEff_stats.html"

        echo "正在注释：$sample"

        "$bcftools" view -O v "$vcf_gz" | \
        "$snpeff" ann -Xmx16g \
            -config "$SNPEFF_CONFIG" \
            -dataDir "$SNPEFF_DATA_DIR" \
            -stats "$stats_html" \
            -v "$CG_DB_NAME" - > "$temp_vcf" 2> "$ann_log"

        if [ -s "$temp_vcf" ]; then
            "$bcftools" view -O z -o "$ann_vcf_gz" "$temp_vcf"
            "$bcftools" index -f "$ann_vcf_gz"
            rm -f "$temp_vcf"
            echo "✅ 注释完成：$sample"
        else
            echo "⚠️ 注释结果为空：$sample"
            rm -f "$temp_vcf"
        fi
    done

    "$multiqc" -o "$output_dir" "$output_dir" || true
}

annotate_dir "$variant_poly" "$variant_poly_annot" "polymorphic"
annotate_dir "$variant_high" "$variant_high_annot" "highfreq"

# -------------------------- 16. 统计 --------------------------
echo -e "\n=== Step 12: 变异统计 ==="

for ann_dir in "$variant_poly_annot" "$variant_high_annot"; do
    cd "$ann_dir"

    ann_vcfs=( *.annotated.vcf.gz )
    if [ "${#ann_vcfs[@]}" -eq 0 ]; then
        echo "⚠️ $ann_dir 下没有 annotated VCF，跳过统计"
        continue
    fi

    for vcf_gz in "${ann_vcfs[@]}"; do
        [ -f "$vcf_gz" ] || continue
        sample=$(basename "$vcf_gz" ".annotated.vcf.gz")
        "$bcftools" stats "$vcf_gz" > "${ann_dir}${sample}.variant_stats.txt"
    done

    "$multiqc" -o "$ann_dir" "$ann_dir" || true
done

# -------------------------- 17. 清理临时文件 --------------------------
echo -e "\n=== Step 13: 清理临时文件 ==="
rm -f "$pair_tsv"

# -------------------------- 18. 完成提示 --------------------------
echo -e "\n=============================================="
echo "✅ 谷氨酸棒杆菌重测序分析全流程完成（双输出版，多命名兼容，升级稳健版）！"
echo "=============================================="
echo "结果路径汇总："
echo "1. 原始质控：$fastqc_raw"
echo "2. 修剪后质控：$fastqc_clean"
echo "3. 清洁数据：$trimgalore_dir"
echo "4. 去重BAM：$dedup"
echo "5. 原始VCF：$variant_raw"
echo "6. 基础质量过滤VCF：$variant_qc"
echo "7. 多态突变VCF：$variant_poly"
echo "8. 高频突变VCF：$variant_high"
echo "9. 多态突变注释：$variant_poly_annot"
echo "10. 高频突变注释：$variant_high_annot"
echo "11. AF拆分表格：$variant_tables"
echo "12. 日志目录：$log_dir"
echo "=============================================="

exit 0


```