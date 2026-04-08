```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CRISPRi双基因筛选数据分析流程（正式生产版）

特点：
1. 仅保留当前数据中验证有效的 direct 3-repeat 提取逻辑
2. 使用已验证的方向规则：
   - 后一个 repeat 是 repeat1 -> 前面的 spacer reverse complement
   - 后一个 repeat 是 repeat2 -> 前面的 spacer 保持原样
3. 修复 FASTA header / spacer1-spacer2 配对问题
4. gene level 与 gRNA level 统计按原始 read 配对后计数
5. 适合正式全量运行
"""

import os
import sys
import re
import gzip
import time
import logging
import argparse
import subprocess
from typing import Dict, List, Optional, Tuple

import pandas as pd
import pysam
from Bio import SeqIO
from fuzzysearch import find_near_matches


# =========================
# 日志设置
# =========================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("crispri_dual_gene_analysis.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class CRISPRiDualGeneAnalyzer:
    def __init__(
        self,
        rawdata_dir: str,
        reference_fasta: str,
        anno_file: str,
        output_dir: str,
        threads: int = 10,
        max_reads: Optional[int] = None,
        use_compressed: bool = True,
        repeat_max_errors: int = 3,
    ):
        self.rawdata_dir = rawdata_dir
        self.library_file = reference_fasta
        self.anno_file = anno_file
        self.result_dir = output_dir
        self.threads = threads
        self.max_reads = max_reads
        self.use_compressed = use_compressed
        self.repeat_max_errors = repeat_max_errors

        self.repeat1 = "ATCTACAACAGTAGAAATTC"
        self.repeat2 = "GAATTTCTACTGTTGTAGAT"
        self.repeat_len = 20
        self.ideal_spacer_len = 23
        self.ideal_repeat_gap = self.repeat_len + self.ideal_spacer_len  # 43

        self.samples = self.get_samples()
        logger.info(f"发现 {len(self.samples)} 个样本: {', '.join(self.samples)}")

        self.create_result_directories()
        self.load_annotation()

        self.index_prefix = None

    # =========================
    # 基础函数
    # =========================
    @staticmethod
    def strip_fastq_suffix(filename: str) -> str:
        for suf in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
            if filename.endswith(suf):
                return filename[:-len(suf)]
        return os.path.splitext(filename)[0]

    @staticmethod
    def reverse_complement(seq: str) -> str:
        complement = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
        return seq.translate(complement)[::-1]

    @staticmethod
    def hamming_distance(a: str, b: str) -> int:
        if len(a) != len(b):
            return max(len(a), len(b))
        return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

    @staticmethod
    def safe_int(x, default=0):
        try:
            return int(x)
        except Exception:
            return default

    @staticmethod
    def parse_query_name(query_name: str) -> Tuple[Optional[str], Optional[str]]:
        m = re.match(r"(.+)_spacer([12])$", query_name)
        if not m:
            return None, None
        return m.group(1), m.group(2)

    def sanitize_read_id(self, read_name: str) -> str:
        read_name = str(read_name).split()[0]
        read_name = read_name.replace("/", "_")
        return read_name

    # =========================
    # 目录与样本
    # =========================
    def create_result_directories(self):
        dirs = [
            self.result_dir,
            os.path.join(self.result_dir, "01_spacer"),
            os.path.join(self.result_dir, "02_bowtie2"),
            os.path.join(self.result_dir, "03_sam_extract"),
            os.path.join(self.result_dir, "04_count"),
            os.path.join(self.result_dir, "05_dual_gRNA_matrix"),
            os.path.join(self.result_dir, "05_dual_gene_matrix"),
        ]

        for sample in self.samples:
            dirs.extend([
                os.path.join(self.result_dir, "01_spacer", sample),
                os.path.join(self.result_dir, "02_bowtie2", sample),
                os.path.join(self.result_dir, "03_sam_extract", sample),
                os.path.join(self.result_dir, "04_count", sample),
            ])

        for d in dirs:
            os.makedirs(d, exist_ok=True)

    def get_samples(self) -> List[str]:
        samples = []
        if not os.path.isdir(self.rawdata_dir):
            raise FileNotFoundError(f"原始数据目录不存在: {self.rawdata_dir}")

        for item in os.listdir(self.rawdata_dir):
            item_path = os.path.join(self.rawdata_dir, item)
            if os.path.isdir(item_path):
                fastq_files = [
                    f for f in os.listdir(item_path)
                    if f.endswith((".fastq.gz", ".fq.gz", ".fq", ".fastq"))
                ]
                if fastq_files:
                    samples.append(item)
        return sorted(samples)

    # =========================
    # 注释
    # =========================
    def load_annotation(self):
        logger.info("加载注释文件...")
        anno_df = pd.read_csv(self.anno_file)

        self.gene_annotation = {}
        for _, row in anno_df.iterrows():
            gene_id = row.get("cgl-id", "")
            if pd.isna(gene_id) or gene_id == "":
                continue
            self.gene_annotation[gene_id] = {
                "gene_name": row.get("gene-name", ""),
                "product": row.get("product", ""),
                "start": row.get("cgl_start", 0),
                "end": row.get("cgl_end", 0),
            }

        logger.info(f"成功加载 {len(self.gene_annotation)} 个基因注释")

    # =========================
    # repeat 识别
    # =========================
    def deduplicate_overlapping_matches(self, all_matches: List[Dict]) -> List[Dict]:
        all_matches.sort(key=lambda x: (x["start"], x["edit_distance"]))
        filtered_matches = []

        for match in all_matches:
            overlapping_idx = []
            for i, old in enumerate(filtered_matches):
                if not (match["end"] <= old["start"] or match["start"] >= old["end"]):
                    overlapping_idx.append(i)

            if not overlapping_idx:
                filtered_matches.append(match)
            else:
                overlapping = [filtered_matches[i] for i in overlapping_idx]
                if all(match["edit_distance"] < old["edit_distance"] for old in overlapping):
                    new_filtered = []
                    for i, old in enumerate(filtered_matches):
                        if i not in overlapping_idx:
                            new_filtered.append(old)
                    new_filtered.append(match)
                    filtered_matches = new_filtered

        filtered_matches.sort(key=lambda x: x["start"])
        return filtered_matches

    def collect_repeat_matches(self, seq_str: str, max_errors: int) -> List[Dict]:
        all_matches = []
        for repeat_seq in (self.repeat1, self.repeat2):
            matches = find_near_matches(repeat_seq, seq_str, max_l_dist=max_errors)
            for match in matches:
                all_matches.append({
                    "start": match.start,
                    "end": match.end,
                    "repeat_seq": repeat_seq,
                    "edit_distance": match.dist
                })

        if not all_matches:
            return []

        return self.deduplicate_overlapping_matches(all_matches)

    def score_triplet(self, m1: Dict, m2: Dict, m3: Dict) -> float:
        spacer1_len = m2["start"] - m1["end"]
        spacer2_len = m3["start"] - m2["end"]
        gap12 = m2["start"] - m1["start"]
        gap23 = m3["start"] - m2["start"]

        score = 0.0
        score += m1["edit_distance"] + m2["edit_distance"] + m3["edit_distance"]
        score += abs(spacer1_len - self.ideal_spacer_len) * 2.0
        score += abs(spacer2_len - self.ideal_spacer_len) * 2.0
        score += abs(gap12 - self.ideal_repeat_gap) * 1.0
        score += abs(gap23 - self.ideal_repeat_gap) * 1.0
        return score

    def extract_from_triplet(self, seq_str: str, m1: Dict, m2: Dict, m3: Dict) -> Optional[Dict]:
        spacer1 = seq_str[m1["end"]:m2["start"]]
        spacer2 = seq_str[m2["end"]:m3["start"]]

        if not (20 <= len(spacer1) <= 26 and 20 <= len(spacer2) <= 26):
            return None

        # 已验证规则
        if m2["repeat_seq"] == self.repeat1:
            spacer1 = self.reverse_complement(spacer1)
        if m3["repeat_seq"] == self.repeat1:
            spacer2 = self.reverse_complement(spacer2)

        return {
            "spacer1": spacer1,
            "spacer2": spacer2,
            "repeat_types": [
                "repeat1" if m1["repeat_seq"] == self.repeat1 else "repeat2",
                "repeat1" if m2["repeat_seq"] == self.repeat1 else "repeat2",
                "repeat1" if m3["repeat_seq"] == self.repeat1 else "repeat2",
            ]
        }

    def find_best_triplet(self, seq_str: str, matches: List[Dict]) -> Optional[Dict]:
        if len(matches) < 3:
            return None

        best = None
        best_score = None
        n = len(matches)

        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    m1, m2, m3 = matches[i], matches[j], matches[k]

                    if not (m1["end"] <= m2["start"] and m2["end"] <= m3["start"]):
                        continue

                    extracted = self.extract_from_triplet(seq_str, m1, m2, m3)
                    if extracted is None:
                        continue

                    score = self.score_triplet(m1, m2, m3)
                    if best_score is None or score < best_score:
                        best_score = score
                        best = {
                            "spacer1": extracted["spacer1"],
                            "spacer2": extracted["spacer2"],
                            "repeat_types": extracted["repeat_types"],
                            "score": score,
                        }

        return best

    # =========================
    # spacer 提取
    # =========================
    def extract_dual_spacer(self, fastq_file: str, output_file: str) -> Dict[str, int]:
        logger.info(f"开始提取双基因spacer序列: {os.path.basename(fastq_file)}")
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        dual_spacers = {}
        processed_count = 0
        matched_triplet_count = 0
        t0 = time.time()

        try:
            if fastq_file.endswith(".gz"):
                handle = gzip.open(fastq_file, "rt")
            else:
                handle = open(fastq_file, "r")

            with handle:
                for seq_record in SeqIO.parse(handle, "fastq"):
                    if self.max_reads and processed_count >= self.max_reads:
                        break

                    seq_str = str(seq_record.seq).upper()
                    read_id = self.sanitize_read_id(seq_record.id)

                    matches = self.collect_repeat_matches(seq_str, max_errors=self.repeat_max_errors)
                    result = self.find_best_triplet(seq_str, matches) if matches else None

                    if result is not None:
                        dual_spacers[read_id] = {
                            "spacer1": result["spacer1"],
                            "spacer2": result["spacer2"],
                            "repeat_types": result["repeat_types"],
                            "full_seq": seq_str
                        }
                        matched_triplet_count += 1

                    processed_count += 1

                    if processed_count % 100000 == 0:
                        logger.info(
                            f"{os.path.basename(fastq_file)}: 已处理 {processed_count:,} reads, "
                            f"提取到 {len(dual_spacers):,} 个双基因reads"
                        )

        except Exception as e:
            logger.error(f"处理FASTQ文件失败: {e}")
            return {
                "processed_reads": processed_count,
                "dual_reads": 0,
                "matched_triplet_reads": 0,
                "elapsed_sec": int(time.time() - t0),
            }

        with open(output_file, "w") as f:
            for key, data in dual_spacers.items():
                header1 = (
                    f">{key}_spacer1 I "
                    f"spacer1={data['spacer1']} spacer2={data['spacer2']} "
                    f"repeat_types={data['repeat_types']}\n"
                )
                f.write(header1)
                f.write(data["spacer1"] + "\n")

                header2 = (
                    f">{key}_spacer2 I "
                    f"spacer1={data['spacer1']} spacer2={data['spacer2']} "
                    f"repeat_types={data['repeat_types']}\n"
                )
                f.write(header2)
                f.write(data["spacer2"] + "\n")

        elapsed = int(time.time() - t0)
        logger.info(
            f"提取完成: {len(dual_spacers)} 个双基因reads -> {output_file}; "
            f"processed={processed_count}, direct={matched_triplet_count}, elapsed={elapsed}s"
        )

        return {
            "processed_reads": processed_count,
            "dual_reads": len(dual_spacers),
            "matched_triplet_reads": matched_triplet_count,
            "elapsed_sec": elapsed,
        }

    # =========================
    # Bowtie2
    # =========================
    def build_bowtie2_index(self) -> str:
        index_dir = os.path.join(self.result_dir, "02_bowtie2")
        index_prefix = os.path.join(index_dir, "reference_index")

        if os.path.exists(f"{index_prefix}.1.bt2"):
            logger.info("bowtie2索引已存在，跳过构建")
            return index_prefix

        cmd = f"bowtie2-build {self.library_file} {index_prefix}"
        logger.info(f"构建bowtie2索引: {cmd}")

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            logger.info("bowtie2索引构建成功")
            return index_prefix
        raise RuntimeError(f"bowtie2索引构建失败: {result.stderr}")

    def run_bowtie2_alignment(self, fasta_file: str, sam_file: str, index_prefix: str) -> bool:
        cmd = (
            f"bowtie2 -f -a -p {self.threads} --local "
            f"-x {index_prefix} -U {fasta_file} -S {sam_file} "
            f"-N 1 -L 22 --quiet"
        )
        logger.info(f"运行bowtie2比对: {os.path.basename(fasta_file)}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            logger.info(f"bowtie2比对完成: {sam_file}")
            return True
        logger.error(f"bowtie2比对失败: {result.stderr}")
        return False

    # =========================
    # SAM提取
    # =========================
    def extract_from_sam(self, sam_file: str, output_file: str) -> Dict[str, int]:
        logger.info(f"从SAM文件提取有效比对: {os.path.basename(sam_file)}")
        logger.info("过滤条件: mapped, query_alignment_length=20-26bp")

        valid_alignments = []
        total_reads = 0
        valid_reads = 0

        try:
            with pysam.AlignmentFile(sam_file, "r") as samfile:
                for read in samfile:
                    total_reads += 1

                    if read.is_unmapped:
                        continue
                    if read.cigarstring is None:
                        continue

                    qalen = read.query_alignment_length
                    if qalen is None:
                        continue

                    if 20 <= qalen <= 26:
                        valid_reads += 1
                        valid_alignments.append({
                            "query_name": read.query_name,
                            "reference_name": read.reference_name,
                            "mapping_quality": read.mapping_quality,
                            "query_alignment_length": qalen,
                            "query_length": read.query_length if read.query_length else 0,
                            "cigar": read.cigarstring
                        })

        except Exception as e:
            logger.error(f"处理SAM文件失败: {e}")
            return {"total_records": total_reads, "valid_records": 0}

        with open(output_file, "w") as f:
            f.write("query_name\treference_name\tmapping_quality\tquery_alignment_length\tquery_length\tcigar\n")
            for aln in valid_alignments:
                f.write(
                    f"{aln['query_name']}\t{aln['reference_name']}\t{aln['mapping_quality']}\t"
                    f"{aln['query_alignment_length']}\t{aln['query_length']}\t{aln['cigar']}\n"
                )

        percent = (valid_reads / total_reads * 100) if total_reads > 0 else 0.0
        logger.info(f"总记录数: {total_reads}, 通过过滤: {valid_reads} ({percent:.1f}%)")
        logger.info(f"提取完成: {len(valid_alignments)} 个有效比对 -> {output_file}")

        return {"total_records": total_reads, "valid_records": valid_reads}

    # =========================
    # reference -> gene
    # =========================
    def extract_gene_id_from_reference(self, reference_name: str) -> str:
        if reference_name.startswith("cgb_"):
            m = re.match(r"cgb_([^-]+)", reference_name)
            if m:
                try:
                    return f"Cgl{int(m.group(1))}"
                except Exception:
                    return f"unknown_{reference_name[:20]}"
        elif reference_name.startswith("Cgl"):
            m = re.match(r"(Cgl\d+)", reference_name)
            if m:
                return m.group(1)

        return f"unknown_{reference_name[:20]}"

    # =========================
    # 读取 spacer header
    # =========================
    def load_spacer_pair_mapping(self, spacer_file: str) -> Dict[str, Dict[str, str]]:
        pair_map = {}

        with open(spacer_file, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    continue

                header = line.strip()[1:]
                header_first = header.split()[0]
                parent_id, spacer_no = self.parse_query_name(header_first)
                if parent_id is None:
                    continue

                m1 = re.search(r"spacer1=([ATCGNatcgn]+)", line)
                m2 = re.search(r"spacer2=([ATCGNatcgn]+)", line)
                if not (m1 and m2):
                    continue

                if parent_id not in pair_map:
                    pair_map[parent_id] = {}

                pair_map[parent_id]["spacer1"] = m1.group(1).upper()
                pair_map[parent_id]["spacer2"] = m2.group(1).upper()

        return pair_map

    # =========================
    # 唯一配对命中
    # =========================
    def collect_paired_unique_hits(self, extract_file: str) -> Dict[str, Dict[str, str]]:
        if not os.path.exists(extract_file) or os.path.getsize(extract_file) == 0:
            return {}

        df = pd.read_csv(extract_file, sep="\t")
        if df.empty:
            return {}

        tmp = {}
        for _, row in df.iterrows():
            query_name = str(row["query_name"])
            reference_name = str(row["reference_name"])

            parent_id, spacer_no = self.parse_query_name(query_name)
            if parent_id is None or spacer_no is None:
                continue

            if parent_id not in tmp:
                tmp[parent_id] = {"1": set(), "2": set()}

            tmp[parent_id][spacer_no].add(reference_name)

        paired_hits = {}
        for parent_id, hits in tmp.items():
            if len(hits["1"]) == 1 and len(hits["2"]) == 1:
                paired_hits[parent_id] = {
                    "1": list(hits["1"])[0],
                    "2": list(hits["2"])[0]
                }

        return paired_hits

    # =========================
    # 计数
    # =========================
    def count_dual_gene_abundance(self, extract_file: str, output_file: str) -> int:
        logger.info(f"统计双基因组合丰度: {os.path.basename(extract_file)}")
        paired_hits = self.collect_paired_unique_hits(extract_file)
        gene_combinations = {}

        for _, hit_pair in paired_hits.items():
            gene1 = self.extract_gene_id_from_reference(hit_pair["1"])
            gene2 = self.extract_gene_id_from_reference(hit_pair["2"])
            genes = sorted([gene1, gene2])
            combo = f"{genes[0]}|{genes[1]}"
            gene_combinations[combo] = gene_combinations.get(combo, 0) + 1

        with open(output_file, "w") as f:
            f.write("dual_gene_combination\tcount\tgene1\tgene2\n")
            for combo, count in sorted(gene_combinations.items()):
                gene1, gene2 = combo.split("|", 1)
                f.write(f"{combo}\t{count}\t{gene1}\t{gene2}\n")

        logger.info(f"双基因组合丰度统计完成: {len(gene_combinations)} 个组合")
        return len(gene_combinations)

    def count_dual_gRNA_abundance(self, extract_file: str, spacer_file: str, output_file: str) -> int:
        logger.info(f"统计双基因gRNA组合丰度: {os.path.basename(extract_file)}")
        spacer_pair_map = self.load_spacer_pair_mapping(spacer_file)
        paired_hits = self.collect_paired_unique_hits(extract_file)
        gRNA_combinations = {}

        for parent_id in paired_hits:
            if parent_id not in spacer_pair_map:
                continue

            spacer1 = spacer_pair_map[parent_id].get("spacer1", "")
            spacer2 = spacer_pair_map[parent_id].get("spacer2", "")
            if not spacer1 or not spacer2:
                continue

            combo = f"{spacer1}|{spacer2}"
            gRNA_combinations[combo] = gRNA_combinations.get(combo, 0) + 1

        with open(output_file, "w") as f:
            f.write("dual_gRNA_combination\tcount\tspacer1\tspacer2\n")
            for combo, count in sorted(gRNA_combinations.items()):
                spacer1, spacer2 = combo.split("|", 1)
                f.write(f"{combo}\t{count}\t{spacer1}\t{spacer2}\n")

        logger.info(f"双基因gRNA组合丰度统计完成: {len(gRNA_combinations)} 个组合")
        return len(gRNA_combinations)

    # =========================
    # 单样本处理
    # =========================
    def process_sample(self, sample: str):
        logger.info(f"开始处理样本: {sample}")
        sample_dir = os.path.join(self.rawdata_dir, sample)

        fastq_files = []
        for file in os.listdir(sample_dir):
            if self.use_compressed and file.endswith((".fq.gz", ".fastq.gz")):
                fastq_files.append(os.path.join(sample_dir, file))
            elif (not self.use_compressed) and file.endswith((".fq", ".fastq")):
                fastq_files.append(os.path.join(sample_dir, file))

        fastq_files = sorted(fastq_files)
        if not fastq_files:
            logger.warning(f"样本 {sample} 中未找到合适的fastq文件")
            return

        for fastq_file in fastq_files:
            file_basename = self.strip_fastq_suffix(os.path.basename(fastq_file))

            spacer_file = os.path.join(self.result_dir, "01_spacer", sample, f"{file_basename}_spacer.fasta")
            spacer_stats = self.extract_dual_spacer(fastq_file, spacer_file)
            if spacer_stats["dual_reads"] == 0:
                logger.warning(f"样本 {sample} 文件 {file_basename} 未提取到有效spacer")
                continue

            sam_file = os.path.join(self.result_dir, "02_bowtie2", sample, f"{file_basename}.sam")
            if not self.run_bowtie2_alignment(spacer_file, sam_file, self.index_prefix):
                logger.warning(f"样本 {sample} 文件 {file_basename} bowtie2比对失败")
                continue

            extract_file = os.path.join(self.result_dir, "03_sam_extract", sample, f"{file_basename}_extract.txt")
            sam_stats = self.extract_from_sam(sam_file, extract_file)
            if sam_stats["valid_records"] == 0:
                logger.warning(f"样本 {sample} 文件 {file_basename} 未提取到有效比对")
                continue

            dual_gene_count_file = os.path.join(self.result_dir, "04_count", sample, f"{file_basename}_dual_gene_count.txt")
            dual_gene_count = self.count_dual_gene_abundance(extract_file, dual_gene_count_file)

            dual_gRNA_count_file = os.path.join(self.result_dir, "04_count", sample, f"{file_basename}_dual_gRNA_count.txt")
            dual_gRNA_count = self.count_dual_gRNA_abundance(extract_file, spacer_file, dual_gRNA_count_file)

            logger.info(
                f"样本 {sample} 文件 {file_basename} 完成："
                f"processed_reads={spacer_stats['processed_reads']}, "
                f"dual_reads={spacer_stats['dual_reads']}, "
                f"valid_alignments={sam_stats['valid_records']}, "
                f"dual_gene_combos={dual_gene_count}, "
                f"dual_gRNA_combos={dual_gRNA_count}"
            )

    # =========================
    # 样本内合并
    # =========================
    def merge_dual_gene_counts(self, sample: str) -> str:
        logger.info(f"合并样本双基因组合计数: {sample}")
        count_dir = os.path.join(self.result_dir, "04_count", sample)
        merged_file = os.path.join(count_dir, f"{sample}_merged_dual_gene_count.txt")

        gene_combinations = {}
        for file in os.listdir(count_dir):
            if file.endswith("_dual_gene_count.txt") and not file.startswith(f"{sample}_merged"):
                file_path = os.path.join(count_dir, file)
                with open(file_path, "r") as f:
                    for i, line in enumerate(f):
                        if i == 0 or not line.strip():
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            combo, count = parts[0], self.safe_int(parts[1], 0)
                            gene_combinations[combo] = gene_combinations.get(combo, 0) + count

        with open(merged_file, "w") as f:
            f.write("dual_gene_combination\tcount\tgene1\tgene2\n")
            for combo, count in sorted(gene_combinations.items()):
                gene1, gene2 = combo.split("|", 1)
                f.write(f"{combo}\t{count}\t{gene1}\t{gene2}\n")

        logger.info(f"样本 {sample} 双基因组合合并完成，共 {len(gene_combinations)} 个组合")
        return merged_file

    def merge_dual_gRNA_counts(self, sample: str) -> str:
        logger.info(f"合并样本双基因gRNA组合计数: {sample}")
        count_dir = os.path.join(self.result_dir, "04_count", sample)
        merged_file = os.path.join(count_dir, f"{sample}_merged_dual_gRNA_count.txt")

        gRNA_combinations = {}
        for file in os.listdir(count_dir):
            if file.endswith("_dual_gRNA_count.txt") and not file.startswith(f"{sample}_merged"):
                file_path = os.path.join(count_dir, file)
                with open(file_path, "r") as f:
                    for i, line in enumerate(f):
                        if i == 0 or not line.strip():
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            combo, count = parts[0], self.safe_int(parts[1], 0)
                            gRNA_combinations[combo] = gRNA_combinations.get(combo, 0) + count

        with open(merged_file, "w") as f:
            f.write("dual_gRNA_combination\tcount\tspacer1\tspacer2\n")
            for combo, count in sorted(gRNA_combinations.items()):
                spacer1, spacer2 = combo.split("|", 1)
                f.write(f"{combo}\t{count}\t{spacer1}\t{spacer2}\n")

        logger.info(f"样本 {sample} 双基因gRNA组合合并完成，共 {len(gRNA_combinations)} 个组合")
        return merged_file

    # =========================
    # 多样本矩阵
    # =========================
    def merge_multi_sample_dual_gRNA_matrix(self) -> str:
        logger.info("合并多个样本的双基因gRNA计数，生成汇总表达矩阵")

        all_combos = set()
        sample_counts = {}

        for sample in self.samples:
            count_dir = os.path.join(self.result_dir, "04_count", sample)
            merged_file = os.path.join(count_dir, f"{sample}_merged_dual_gRNA_count.txt")

            sample_data = {}
            if os.path.exists(merged_file):
                with open(merged_file, "r") as f:
                    for i, line in enumerate(f):
                        if i == 0 or not line.strip():
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            combo = parts[0]
                            count = self.safe_int(parts[1], 0)
                            sample_data[combo] = count
                            all_combos.add(combo)
            sample_counts[sample] = sample_data

        output_dir = os.path.join(self.result_dir, "05_dual_gRNA_matrix")
        os.makedirs(output_dir, exist_ok=True)
        merged_matrix_file = os.path.join(output_dir, "all_samples_dual_gRNA_merged_matrix.txt")

        with open(merged_matrix_file, "w") as f:
            header = ["dual_gRNA_combination"] + self.samples
            f.write("\t".join(header) + "\n")
            for combo in sorted(all_combos):
                row = [combo] + [str(sample_counts[sample].get(combo, 0)) for sample in self.samples]
                f.write("\t".join(row) + "\n")

        logger.info(f"多样本双gRNA矩阵已保存至: {merged_matrix_file}")
        return merged_matrix_file

    def merge_multi_sample_dual_gene_matrix(self) -> str:
        logger.info("合并多个样本的双基因计数，生成汇总表达矩阵")

        all_combos = set()
        sample_counts = {}

        for sample in self.samples:
            count_dir = os.path.join(self.result_dir, "04_count", sample)
            merged_file = os.path.join(count_dir, f"{sample}_merged_dual_gene_count.txt")

            sample_data = {}
            if os.path.exists(merged_file):
                with open(merged_file, "r") as f:
                    for i, line in enumerate(f):
                        if i == 0 or not line.strip():
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            combo = parts[0]
                            count = self.safe_int(parts[1], 0)
                            sample_data[combo] = count
                            all_combos.add(combo)
            sample_counts[sample] = sample_data

        output_dir = os.path.join(self.result_dir, "05_dual_gene_matrix")
        os.makedirs(output_dir, exist_ok=True)
        merged_matrix_file = os.path.join(output_dir, "all_samples_dual_gene_merged_matrix.txt")

        with open(merged_matrix_file, "w") as f:
            header = ["dual_gene_combination"] + self.samples
            f.write("\t".join(header) + "\n")
            for combo in sorted(all_combos):
                row = [combo] + [str(sample_counts[sample].get(combo, 0)) for sample in self.samples]
                f.write("\t".join(row) + "\n")

        logger.info(f"多样本双基因矩阵已保存至: {merged_matrix_file}")
        return merged_matrix_file

    # =========================
    # 总流程
    # =========================
    def run_complete_analysis(self):
        logger.info("开始CRISPRi双基因完整分析流程")
        self.index_prefix = self.build_bowtie2_index()

        for sample in self.samples:
            self.process_sample(sample)
            self.merge_dual_gene_counts(sample)
            self.merge_dual_gRNA_counts(sample)

        spacer_level_matrix = self.merge_multi_sample_dual_gRNA_matrix()
        gene_level_matrix = self.merge_multi_sample_dual_gene_matrix()

        logger.info("CRISPRi双基因分析流程完成！")
        logger.info(f"双gRNA组合表达矩阵: {spacer_level_matrix}")
        logger.info(f"双基因组合表达矩阵: {gene_level_matrix}")

        return {
            "spacer_level_matrix": spacer_level_matrix,
            "gene_level_matrix": gene_level_matrix
        }


def main():
    parser = argparse.ArgumentParser(description="CRISPRi双基因筛选数据分析（正式生产版）")
    parser.add_argument("--rawdata-dir", required=True, help="原始数据目录")
    parser.add_argument("--reference-fasta", required=True, help="参考序列文件")
    parser.add_argument("--anno-file", required=True, help="基因注释文件")
    parser.add_argument("--output-dir", required=True, help="输出结果目录")
    parser.add_argument("--threads", type=int, default=10, help="线程数")
    parser.add_argument("--repeat-max-errors", type=int, default=3, help="repeat匹配最大编辑距离，默认3")
    parser.add_argument("--test", action="store_true", help="测试模式（限制处理reads数量）")
    parser.add_argument("--max-reads", type=int, default=None, help="最大处理reads数量")
    parser.add_argument("--use-uncompressed", action="store_true", help="使用未压缩fastq文件")
    args = parser.parse_args()

    use_compressed = not args.use_uncompressed
    max_reads = args.max_reads

    if args.test:
        max_reads = 50000
        logger.info("运行测试模式，限制处理 50,000 reads")

    try:
        analyzer = CRISPRiDualGeneAnalyzer(
            rawdata_dir=args.rawdata_dir,
            reference_fasta=args.reference_fasta,
            anno_file=args.anno_file,
            output_dir=args.output_dir,
            threads=args.threads,
            max_reads=max_reads,
            use_compressed=use_compressed,
            repeat_max_errors=args.repeat_max_errors,
        )

        result = analyzer.run_complete_analysis()

        print("\n分析完成！")
        print(f"双基因组合表达矩阵（spacer级别）: {result['spacer_level_matrix']}")
        print(f"双基因组合表达矩阵（基因级别）: {result['gene_level_matrix']}")

    except Exception as e:
        print(f"分析失败: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
```