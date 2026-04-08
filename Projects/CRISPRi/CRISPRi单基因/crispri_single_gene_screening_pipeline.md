```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CRISPRi单基因筛选数据分析流程

"""

import os
import sys
import subprocess
import pandas as pd
import pysam
import re
import gzip
from Bio import SeqIO
from collections import Counter
import bisect
import logging
from tqdm import tqdm

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CRISPRiAnalyzerOptimized:
    def __init__(self, rawdata_dir, reference_fasta, anno_file, output_dir, threads=10, max_reads=None, use_compressed=True):
        self.rawdata_dir = rawdata_dir
        self.result_dir = output_dir
        self.library_file = reference_fasta
        self.anno_file = anno_file
        self.threads = threads
        self.use_compressed = use_compressed  # 是否使用压缩文件
        self.max_reads = max_reads  # 限制处理的reads数量（用于测试）
        
        # 获取样本列表
        self.samples = self.get_samples()
        logger.info(f"发现样本: {self.samples}")
        
        # 创建结果目录
        self.create_directories()
        
        # 加载注释数据
        self.load_annotation()
        
    def create_directories(self):
        """创建分析结果目录"""
        dirs = [
            os.path.join(self.result_dir, '01_spacer'),
            os.path.join(self.result_dir, '02_bowtie2'),
            os.path.join(self.result_dir, '03_sam_extract'),
            os.path.join(self.result_dir, '04_count'),
            os.path.join(self.result_dir, '05_gene_matrix'),
            os.path.join(self.result_dir, '05_gRNA_matrix')
        ]
        
        # 为所有样本创建目录
        for sample in self.samples:
            dirs.extend([
                os.path.join(self.result_dir, '01_spacer', sample),
                os.path.join(self.result_dir, '02_bowtie2', sample),
                os.path.join(self.result_dir, '03_sam_extract', sample),
                os.path.join(self.result_dir, '04_count', sample)
            ])
        
        for dir_path in dirs:
            os.makedirs(dir_path, exist_ok=True)
            
    def get_samples(self):
        """获取样本列表"""
        samples = []
        for item in os.listdir(self.rawdata_dir):
            item_path = os.path.join(self.rawdata_dir, item)
            # 识别所有包含时间点信息的样本目录
            if os.path.isdir(item_path):
                samples.append(item)
        return sorted(samples)
    
    def load_annotation(self):
        """加载CGL注释文件"""
        logger.info("加载CGL注释文件...")
        try:
            anno_df = pd.read_csv(self.anno_file)
            self.anno_cgl_id = anno_df['cgl-id'].tolist()
            self.anno_start = anno_df['cgl_start'].tolist()
            self.anno_end = anno_df['cgl_end'].tolist()
            
            self.anno_start_sorted = sorted(self.anno_start)
            self.anno_values_sorted = sorted(self.anno_end)
            self.start_cgl_id = dict(zip(self.anno_start, self.anno_cgl_id))
            self.end_cgl_id = dict(zip(self.anno_end, self.anno_cgl_id))
            
            # 创建基因ID到注释信息的映射
            self.gene_annotation = {}
            for _, row in anno_df.iterrows():
                gene_id = row['cgl-id']
                self.gene_annotation[gene_id] = {
                    'gene_name': row.get('gene-name', ''),
                    'product': row.get('product', ''),
                    'start': row.get('cgl_start', 0),
                    'end': row.get('cgl_end', 0)
                }
            logger.info(f"成功加载 {len(self.anno_cgl_id)} 个基因注释，建立 {len(self.gene_annotation)} 个基因映射")
        except Exception as e:
            logger.error(f"加载注释文件失败: {e}")
            self.gene_annotation = {}
            raise
    
    def count_fastq_reads(self, fastq_file):
        """统计FASTQ文件中的reads数量"""
        count = 0
        try:
            if fastq_file.endswith('.gz'):
                with gzip.open(fastq_file, 'rt') as handle:
                    for record in SeqIO.parse(handle, 'fastq'):
                        count += 1
                        if self.max_reads and count >= self.max_reads:
                            break
            else:
                for record in SeqIO.parse(fastq_file, 'fastq'):
                    count += 1
                    if self.max_reads and count >= self.max_reads:
                        break
        except Exception as e:
            logger.error(f"统计reads失败: {e}")
        return count
    
    def extract_spacer(self, fastq_file, output_file, index='I'):
        """提取spacer序列"""
        logger.info(f"开始提取spacer序列: {os.path.basename(fastq_file)}")
        
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # 定义重复序列模式
        repeat1 = 'ATCTACAACAGTAGAAATTC'
        repeat2 = 'ATCTACAACAG' #具体测序时重复序列
        repeat3 = 'GAATTTCTACTGTTGTAGAT'
        
        valid_spacers = {}
        
        try:
            # 直接处理fq.gz文件
            if fastq_file.endswith('.gz'):
                file_handle = gzip.open(fastq_file, 'rt')
            else:
                file_handle = open(fastq_file, 'r')
            
            with file_handle as handle:
                processed_count = 0
                for seq_record in SeqIO.parse(handle, 'fastq'):
                    if self.max_reads and processed_count >= self.max_reads:
                        break
                    
                    sequence = str(seq_record.seq)
                    spacer_info = None
                    
                    # 模式1: repeat1 -> repeat2（严格顺序）
                    repeat1_positions = [m.start() for m in re.finditer(repeat1, sequence)]
                    if repeat1_positions:
                        for repeat1_pos in repeat1_positions:
                            repeat1_end = repeat1_pos + len(repeat1)
                            remaining_sequence = sequence[repeat1_end:]
                            repeat2_matches = list(re.finditer(repeat2, remaining_sequence))
                            if repeat2_matches:
                                # 取第一个匹配的 repeat2
                                repeat2_relative_pos = repeat2_matches[0].start()
                                repeat2_absolute_pos = repeat1_end + repeat2_relative_pos
                                spacer_seq = sequence[repeat1_end:repeat2_absolute_pos]
                                if len(spacer_seq) > 0:
                                    spacer_info = {
                                        'spacer': spacer_seq,
                                        'spacer_length': len(spacer_seq),
                                        'type': 'repeat12',
                                        'repeat1_pos': repeat1_pos,
                                        'repeat2_pos': repeat2_absolute_pos
                                    }
                                    break  # 找到第一个有效组合即停止
                    
                    # 模式2: 两个 repeat3 之间的序列（若模式1未命中再尝试）
                    if spacer_info is None:
                        repeat3_positions = [m.start() for m in re.finditer(repeat3, sequence)]
                        if len(repeat3_positions) >= 2:
                            pos1 = repeat3_positions[0]
                            pos2 = repeat3_positions[1]
                            start = pos1 + len(repeat3)
                            end = pos2
                            if end > start:
                                spacer_seq = sequence[start:end]
                                if len(spacer_seq) > 0:
                                    spacer_info = {
                                        'spacer': spacer_seq,
                                        'spacer_length': len(spacer_seq),
                                        'type': 'repeat3pair',
                                        'repeat3_pos1': pos1,
                                        'repeat3_pos2': pos2
                                    }
                    
                    # 记录命中
                    if spacer_info is not None:
                        valid_spacers[seq_record.description] = spacer_info
                    
                    processed_count += 1
                    # 每处理10万条记录输出一次进度
                    if processed_count % 100000 == 0:
                        logger.info(f"已处理 {processed_count:,} reads，找到 {len(valid_spacers)} 个有效spacer")
        
        except Exception as e:
            logger.error(f"处理FASTQ文件失败: {e}")
            return 0
        
        # 保存结果（根据命中类型写不同的元信息）
        with open(output_file, 'w') as f:
            for seq_id, info in valid_spacers.items():
                if info['type'] == 'repeat12':
                    id_line = f">{seq_id} {index} spacer_len:{info['spacer_length']} repeat1_pos:{info['repeat1_pos']} repeat2_pos:{info['repeat2_pos']}\n"
                else:
                    id_line = f">{seq_id} {index} spacer_len:{info['spacer_length']} repeat3_pos1:{info['repeat3_pos1']} repeat3_pos2:{info['repeat3_pos2']}\n"
                f.write(id_line)
                f.write(str(info['spacer']) + '\n')
        
        logger.info(f"提取完成: {len(valid_spacers)} 个spacer序列 -> {output_file}")
        return len(valid_spacers)

    def build_bowtie2_index(self):
        """构建bowtie2索引"""
        index_dir = os.path.join(self.result_dir, '02_bowtie2')
        index_prefix = os.path.join(index_dir, 'reference_index')
        
        # 检查索引是否已存在
        if os.path.exists(f"{index_prefix}.1.bt2"):
            logger.info("bowtie2索引已存在，跳过构建")
            return index_prefix
        
        cmd = f"bowtie2-build {self.library_file} {index_prefix}"
        logger.info(f"构建bowtie2索引: {cmd}")
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                logger.info("bowtie2索引构建成功")
                return index_prefix
            else:
                logger.error(f"bowtie2索引构建失败: {result.stderr}")
                raise Exception(f"bowtie2索引构建失败: {result.stderr}")
        except Exception as e:
            logger.error(f"执行bowtie2-build命令失败: {e}")
            raise
    
    def run_bowtie2_alignment(self, fasta_file, sam_file, index_prefix):
        """运行bowtie2比对"""
        # 完整的bowtie2参数设置
        # -f: 输入为FASTA格式
        # -a: 报告所有比对结果
        # -p: 使用指定线程数
        # --local: 局部比对模式
        # -N 1: 允许1个错配
        # -L 22: seed长度为22
        # --quiet: 静默模式
        cmd = f"bowtie2 -f -a -p {self.threads} --local -x {index_prefix} -U {fasta_file} -S {sam_file} -N 1 -L 22 --quiet"
        logger.info(f"运行bowtie2比对: {os.path.basename(fasta_file)} (参数: -f -a -p {self.threads} --local -N 1 -L 22)")
        
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                logger.info(f"bowtie2比对完成: {sam_file}")
                return True
            else:
                logger.error(f"bowtie2比对失败: {result.stderr}")
                return False
        except Exception as e:
            logger.error(f"执行bowtie2命令失败: {e}")
            return False
    
    def separate_and_sum(self, cigar_string):
        """解析CIGAR字符串"""
        separated = re.findall(r'\d+[MD]', cigar_string)
        sum_M = sum_D = 0
        for item in separated:
            if item.endswith('M'):
                sum_M += int(item[:-1])
            elif item.endswith('D'):
                sum_D += int(item[:-1])
        return sum_M + sum_D
    
    def extract_from_sam(self, sam_file, output_file):
        """从SAM文件提取有效比对"""
        logger.info(f"从SAM文件提取有效比对: {os.path.basename(sam_file)}")
        logger.info("过滤条件: mapping_quality=255, CIGAR长度22-24bp")
        sam_data = []
        
        try:
            samfile = pysam.AlignmentFile(sam_file, 'r')
            total_reads = 0
            valid_reads = 0
            cgl_reads = 0
            cgb_reads = 0
            
            for read in samfile:
                total_reads += 1
                # 严格的过滤条件：mapping_quality=255 且 CIGAR长度在22-24之间
                if (read.mapping_quality == 255 and 
                    22 <= self.separate_and_sum(read.cigarstring) <= 24):
                    
                    ref_name = read.reference_name if read.reference_name is not None else 'N'
                    
                    # 统计不同类型的基因
                    if ref_name.startswith('Cgl'):
                        cgl_reads += 1
                    elif ref_name.startswith('cgb_'):
                        cgb_reads += 1
                    
                    sam_data.append([
                        read.query_name if read.query_name else 'N',
                        read.flag if read.flag is not None else 'N',
                        ref_name,
                        read.pos if read.pos is not None else 'N',
                        read.mapping_quality if read.mapping_quality is not None else 'N',
                        read.cigarstring if read.cigarstring else 'N',
                        read.query_sequence if read.query_sequence else 'N'
                    ])
                    valid_reads += 1
                
                # 每处理10万条记录输出一次进度
                if total_reads % 100000 == 0:
                    logger.info(f"已处理 {total_reads:,} 比对，找到 {valid_reads} 个有效比对 (Cgl: {cgl_reads}, cgb: {cgb_reads})")
            
            samfile.close()
            
            # 保存结果
            with open(output_file, 'w') as f:
                for data in sam_data:
                    line = ', '.join(str(x) for x in data) + '\n'
                    f.write(line)
            
            logger.info(f"提取完成: {len(sam_data)} 个有效比对 -> {output_file}")
            logger.info(f"基因类型分布: Cgl基因 {cgl_reads} 个, cgb基因 {cgb_reads} 个")
            return len(sam_data)
            
        except Exception as e:
            logger.error(f"处理SAM文件失败: {e}")
            return 0
    
    def count_gene_abundance(self, extract_file, output_file):
        """统计基因丰度"""
        logger.info(f"统计基因丰度: {os.path.basename(extract_file)}")
        
        try:
            # 读取基因ID信息（第3列，索引为2）
            df = pd.read_csv(extract_file, sep=', ', usecols=[2], header=None)
            gene_refs = df.iloc[:, 0].tolist()
            
            cgl_ids = []
            cgb_count = 0
            unmapped_count = 0
            
            for gene_ref in tqdm(gene_refs, desc="映射基因"):
                if pd.isna(gene_ref) or gene_ref == 'N':
                    cgl_ids.append('')
                    unmapped_count += 1
                    continue
                
                try:
                    # 从基因引用中提取基因ID
                    # 格式: Cgl1953-822-5***86***+ 或 cgb_00105-133-133-1***26***+
                    gene_ref_str = str(gene_ref)
                    
                    if gene_ref_str.startswith('Cgl'):
                        # 提取Cgl基因ID（如Cgl1953）
                        gene_id = gene_ref_str.split('-')[0]
                        # 验证基因ID是否在注释文件中
                        if hasattr(self, 'gene_annotation') and gene_id in self.gene_annotation:
                            cgl_ids.append(gene_id)
                        else:
                            # 即使不在注释文件中，也保留Cgl基因
                            cgl_ids.append(gene_id)
                    elif gene_ref_str.startswith('cgb_'):
                        # cgb基因被过滤掉，不计入最终统计
                        cgb_count += 1
                        cgl_ids.append('')
                    else:
                        cgl_ids.append('')
                        unmapped_count += 1
                        
                except Exception as e:
                    cgl_ids.append('')
                    unmapped_count += 1
            
            # 统计基因丰度
            gene_counts = Counter(cgl_ids)
            if '' in gene_counts:
                del gene_counts['']
            
            # 保存结果（包含注释信息）
            with open(output_file, 'w') as f:
                # 写入表头
                f.write("gene_id\tcount\tgene_name\tproduct\n")
                for gene_id, count in sorted(gene_counts.items()):
                    # 获取注释信息
                    if hasattr(self, 'gene_annotation') and gene_id in self.gene_annotation:
                        gene_name = self.gene_annotation[gene_id].get('gene_name', '')
                        product = self.gene_annotation[gene_id].get('product', '')
                    else:
                        gene_name = ''
                        product = ''
                    f.write(f"{gene_id}\t{count}\t{gene_name}\t{product}\n")
            
            logger.info(f"基因丰度统计完成: {len(gene_counts)} 个Cgl基因，{cgb_count} 个cgb基因被过滤，{unmapped_count} 个未映射位点")
            return len(gene_counts)
            
        except Exception as e:
            logger.error(f"统计基因丰度失败: {e}")
            return 0
    
    def count_gRNA_abundance(self, extract_file, output_file):
        """统计gRNA丰度(完整的基因引用ID)"""
        logger.info(f"统计gRNA丰度: {os.path.basename(extract_file)}")
        
        try:
            # 读取gRNA信息（第3列，索引为2）
            df = pd.read_csv(extract_file, sep=', ', usecols=[2], header=None)
            gRNA_refs = df.iloc[:, 0].tolist()
            
            cgl_gRNAs = []
            cgb_count = 0
            unmapped_count = 0
            
            for gRNA_ref in tqdm(gRNA_refs, desc="统计gRNA"):
                if pd.isna(gRNA_ref) or gRNA_ref == 'N':
                    cgl_gRNAs.append('')
                    unmapped_count += 1
                    continue
                
                try:
                    # 保留完整的gRNA引用ID
                    # 格式: Cgl1953-822-5***86***+ 或 cgb_00105-133-133-1***26***+
                    gRNA_ref_str = str(gRNA_ref)
                    
                    if gRNA_ref_str.startswith('Cgl'):
                        # 保留完整的Cgl gRNA ID
                        cgl_gRNAs.append(gRNA_ref_str)
                    elif gRNA_ref_str.startswith('cgb_'):
                        # cgb基因被过滤掉，不计入最终统计
                        cgb_count += 1
                        cgl_gRNAs.append('')
                    else:
                        cgl_gRNAs.append('')
                        unmapped_count += 1
                        
                except Exception as e:
                    cgl_gRNAs.append('')
                    unmapped_count += 1
            
            # 统计gRNA丰度
            gRNA_counts = Counter(cgl_gRNAs)
            if '' in gRNA_counts:
                del gRNA_counts['']
            
            # 保存结果
            with open(output_file, 'w') as f:
                # 写入表头
                f.write("gRNA_id\tcount\tgene_id\n")
                for gRNA_id, count in sorted(gRNA_counts.items()):
                    # 从gRNA_id中提取基因ID
                    gene_id = gRNA_id.split('-')[0] if '-' in gRNA_id else gRNA_id
                    f.write(f"{gRNA_id}\t{count}\t{gene_id}\n")
            
            logger.info(f"gRNA丰度统计完成: {len(gRNA_counts)} 个Cgl gRNA，{cgb_count} 个cgb gRNA被过滤，{unmapped_count} 个未映射位点")
            return len(gRNA_counts)
            
        except Exception as e:
            logger.error(f"统计gRNA丰度失败: {e}")
            return 0
    
    def process_sample(self, sample):
        """处理单个样本"""
        logger.info(f"开始处理样本: {sample}")
        
        sample_dir = os.path.join(self.rawdata_dir, sample)
        
        # 处理R1和R2文件
        for read_num in ['1', '2']:
            # 选择输入文件（优先使用压缩文件）
            if self.use_compressed:
                fastq_file = os.path.join(sample_dir, f"{sample}_{read_num}.fq.gz")
            else:
                fastq_file = os.path.join(sample_dir, f"{sample}_{read_num}.fq")
            
            if not os.path.exists(fastq_file):
                # 尝试另一种格式
                alt_file = os.path.join(sample_dir, f"{sample}_{read_num}.fq.gz" if not self.use_compressed else f"{sample}_{read_num}.fq")
                if os.path.exists(alt_file):
                    fastq_file = alt_file
                else:
                    logger.warning(f"未找到文件: {fastq_file}")
                    continue
            
            # 输出文件路径
            spacer_file = os.path.join(self.result_dir, '01_spacer', sample, f"{sample}_{read_num}_spacer.fasta")
            sam_file = os.path.join(self.result_dir, '02_bowtie2', sample, f"{sample}_{read_num}.sam")
            extract_file = os.path.join(self.result_dir, '03_sam_extract', sample, f"{sample}_{read_num}_extract.txt")
            gene_count_file = os.path.join(self.result_dir, '04_count', sample, f"{sample}_{read_num}_gene_count.txt")
            gRNA_count_file = os.path.join(self.result_dir, '04_count', sample, f"{sample}_{read_num}_gRNA_count.txt")
            
            # 1. 提取spacer
            spacer_count = self.extract_spacer(fastq_file, spacer_file)
            if spacer_count == 0:
                logger.warning(f"样本 {sample}_{read_num} 未提取到spacer序列")
                continue
            
            # 2. bowtie2比对
            if not self.run_bowtie2_alignment(spacer_file, sam_file, self.index_prefix):
                logger.warning(f"样本 {sample}_{read_num} bowtie2比对失败")
                continue
            
            # 3. 从SAM文件提取有效比对
            extract_count = self.extract_from_sam(sam_file, extract_file)
            if extract_count == 0:
                logger.warning(f"样本 {sample}_{read_num} 未提取到有效比对")
                continue
            
            # 4. 统计基因丰度
            gene_count = self.count_gene_abundance(extract_file, gene_count_file)
            
            # 5. 统计gRNA丰度
            gRNA_count = self.count_gRNA_abundance(extract_file, gRNA_count_file)
            
            logger.info(f"样本 {sample}_{read_num} 处理完成，统计到 {gene_count} 个基因，{gRNA_count} 个gRNA")
    
    def merge_sample_counts(self, sample):
        """合并样本的R1和R2计数"""
        logger.info(f"合并样本计数: {sample}")
        
        count_dir = os.path.join(self.result_dir, '04_count', sample)
        r1_file = os.path.join(count_dir, f"{sample}_1_gene_count.txt")
        r2_file = os.path.join(count_dir, f"{sample}_2_gene_count.txt")
        merged_file = os.path.join(count_dir, f"{sample}_merged_gene_count.txt")
        
        gene_counts = {}
        gene_annotations = {}
        
        # 读取R1计数（跳过表头）
        if os.path.exists(r1_file):
            with open(r1_file, 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if i == 0 and line.startswith('gene_id'):
                        continue  # 跳过表头
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            gene_id, count = parts[0], parts[1]
                            gene_counts[gene_id] = gene_counts.get(gene_id, 0) + int(count)
                            # 保存注释信息
                            if len(parts) >= 4:
                                gene_annotations[gene_id] = {
                                    'gene_name': parts[2],
                                    'product': parts[3]
                                }
        
        # 读取R2计数（跳过表头）
        if os.path.exists(r2_file):
            with open(r2_file, 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if i == 0 and line.startswith('gene_id'):
                        continue  # 跳过表头
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            gene_id, count = parts[0], parts[1]
                            gene_counts[gene_id] = gene_counts.get(gene_id, 0) + int(count)
                            # 保存注释信息
                            if len(parts) >= 4 and gene_id not in gene_annotations:
                                gene_annotations[gene_id] = {
                                    'gene_name': parts[2],
                                    'product': parts[3]
                                }
        
        # 保存合并结果（包含注释信息）
        with open(merged_file, 'w') as f:
            f.write("gene_id\tcount\tgene_name\tproduct\n")
            for gene_id, count in sorted(gene_counts.items()):
                if gene_id in gene_annotations:
                    gene_name = gene_annotations[gene_id]['gene_name']
                    product = gene_annotations[gene_id]['product']
                else:
                    gene_name = ''
                    product = ''
                f.write(f"{gene_id}\t{count}\t{gene_name}\t{product}\n")
        
        logger.info(f"样本 {sample} 合并完成，共 {len(gene_counts)} 个基因")
        return merged_file
    
    def merge_gRNA_counts(self, sample):
        """合并样本的R1和R2 gRNA计数"""
        logger.info(f"合并样本gRNA计数: {sample}")
        
        count_dir = os.path.join(self.result_dir, '04_count', sample)
        r1_file = os.path.join(count_dir, f"{sample}_1_gRNA_count.txt")
        r2_file = os.path.join(count_dir, f"{sample}_2_gRNA_count.txt")
        merged_file = os.path.join(count_dir, f"{sample}_merged_gRNA_count.txt")
        
        gRNA_counts = {}
        
        # 读取R1计数（跳过表头）
        if os.path.exists(r1_file):
            with open(r1_file, 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if i == 0 and line.startswith('gRNA_id'):
                        continue  # 跳过表头
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            gRNA_id, count = parts[0], parts[1]
                            gRNA_counts[gRNA_id] = gRNA_counts.get(gRNA_id, 0) + int(count)
        
        # 读取R2计数（跳过表头）
        if os.path.exists(r2_file):
            with open(r2_file, 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if i == 0 and line.startswith('gRNA_id'):
                        continue  # 跳过表头
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            gRNA_id, count = parts[0], parts[1]
                            gRNA_counts[gRNA_id] = gRNA_counts.get(gRNA_id, 0) + int(count)
        
        # 保存合并结果
        with open(merged_file, 'w') as f:
            f.write("gRNA_id\tcount\tgene_id\n")
            for gRNA_id, count in sorted(gRNA_counts.items()):
                # 从gRNA_id中提取基因ID
                gene_id = gRNA_id.split('-')[0] if '-' in gRNA_id else gRNA_id
                f.write(f"{gRNA_id}\t{count}\t{gene_id}\n")
        
        logger.info(f"样本 {sample} gRNA合并完成，共 {len(gRNA_counts)} 个gRNA")
        return merged_file
    
    def create_gene_expression_matrix(self):
        """创建基因表达矩阵（不包含注释信息）"""
        logger.info("创建基因表达矩阵...")
        
        # 收集所有基因ID
        all_genes = set()
        sample_data = {}
        
        for sample in self.samples:
            merged_file = os.path.join(self.result_dir, '04_count', sample, f"{sample}_merged_gene_count.txt")
            if os.path.exists(merged_file):
                sample_counts = {}
                with open(merged_file, 'r') as f:
                    lines = f.readlines()
                    for i, line in enumerate(lines):
                        if i == 0 and line.startswith('gene_id'):
                            continue  # 跳过表头
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= 2:
                                gene_id, count = parts[0], parts[1]
                                sample_counts[gene_id] = int(count)
                                all_genes.add(gene_id)
                sample_data[sample] = sample_counts
            else:
                logger.warning(f"未找到样本 {sample} 的基因合并计数文件")
                sample_data[sample] = {}
        
        # 创建基因表达矩阵
        all_genes = sorted(list(all_genes))
        matrix_data = []
        
        for gene_id in all_genes:
            row = [gene_id]
            # 添加样本计数
            for sample in sorted(self.samples):
                count = sample_data[sample].get(gene_id, 0)
                row.append(count)
            matrix_data.append(row)
        
        # 保存为CSV文件
        output_file = os.path.join(self.result_dir, '05_gene_matrix', 'gene_expression_matrix.csv')
        
        columns = ['gene_id'] + sorted(self.samples)
        df = pd.DataFrame(matrix_data, columns=columns)
        df.to_csv(output_file, index=False)
        
        logger.info(f"基因表达矩阵已保存: {output_file}")
        logger.info(f"矩阵维度: {len(all_genes)} 基因 × {len(self.samples)} 样本")
        
        # 显示矩阵摘要
        print("\n=== 基因表达矩阵摘要 ===")
        print(f"基因数量: {len(all_genes)}")
        print(f"样本数量: {len(self.samples)}")
        print(f"样本列表: {', '.join(sorted(self.samples))}")
        print(f"结果文件: {output_file}")
        
        return output_file
    
    def create_gRNA_expression_matrix(self):
        """创建gRNA表达矩阵（不包含注释信息）"""
        logger.info("创建gRNA表达矩阵...")
        
        # 收集所有gRNA ID
        all_gRNAs = set()
        sample_data = {}
        
        for sample in self.samples:
            merged_file = os.path.join(self.result_dir, '04_count', sample, f"{sample}_merged_gRNA_count.txt")
            if os.path.exists(merged_file):
                sample_counts = {}
                with open(merged_file, 'r') as f:
                    lines = f.readlines()
                    for i, line in enumerate(lines):
                        if i == 0 and line.startswith('gRNA_id'):
                            continue  # 跳过表头
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= 2:
                                gRNA_id, count = parts[0], parts[1]
                                sample_counts[gRNA_id] = int(count)
                                all_gRNAs.add(gRNA_id)
                sample_data[sample] = sample_counts
            else:
                logger.warning(f"未找到样本 {sample} 的gRNA合并计数文件")
                sample_data[sample] = {}
        
        # 创建gRNA表达矩阵
        all_gRNAs = sorted(list(all_gRNAs))
        matrix_data = []
        
        for gRNA_id in all_gRNAs:
            row = [gRNA_id]
            # 添加样本计数
            for sample in sorted(self.samples):
                count = sample_data[sample].get(gRNA_id, 0)
                row.append(count)
            matrix_data.append(row)
        
        # 保存为CSV文件
        output_file = os.path.join(self.result_dir, '05_gRNA_matrix', 'gRNA_expression_matrix.csv')
        
        columns = ['gRNA_id'] + sorted(self.samples)
        df = pd.DataFrame(matrix_data, columns=columns)
        df.to_csv(output_file, index=False)
        
        logger.info(f"gRNA表达矩阵已保存: {output_file}")
        logger.info(f"矩阵维度: {len(all_gRNAs)} gRNA × {len(self.samples)} 样本")
        
        # 显示矩阵摘要
        print("\n=== gRNA表达矩阵摘要 ===")
        print(f"gRNA数量: {len(all_gRNAs)}")
        print(f"样本数量: {len(self.samples)}")
        print(f"样本列表: {', '.join(sorted(self.samples))}")
        print(f"结果文件: {output_file}")
        
        return output_file
    
    def run_analysis(self, samples):
        """运行指定样本的分析流程"""
        logger.info(f"开始CRISPRi分析流程，处理 {len(samples)} 个样本")
        
        try:
            # 1. 构建bowtie2索引
            self.index_prefix = self.build_bowtie2_index()
            
            # 2. 处理指定样本
            for sample in samples:
                if sample in self.samples:
                    self.process_sample(sample)
                    self.merge_sample_counts(sample)
                    self.merge_gRNA_counts(sample)
                else:
                    logger.warning(f"样本 {sample} 不在可用样本列表中，跳过处理")
            
            # 3. 创建输出目录
            os.makedirs(os.path.join(self.result_dir, '05_gene_matrix'), exist_ok=True)
            os.makedirs(os.path.join(self.result_dir, '05_gRNA_matrix'), exist_ok=True)
            
            # 4. 创建基因和gRNA表达矩阵
            gene_matrix_file = self.create_gene_expression_matrix()
            gRNA_matrix_file = self.create_gRNA_expression_matrix()
            
            logger.info("CRISPRi分析流程完成！")
            logger.info(f"基因表达矩阵: {gene_matrix_file}")
            logger.info(f"gRNA表达矩阵: {gRNA_matrix_file}")
            
            return gene_matrix_file, gRNA_matrix_file
            
        except Exception as e:
            logger.error(f"分析流程失败: {e}")
            raise
    
    def run_complete_analysis(self):
        """运行完整分析流程"""
        logger.info("开始CRISPRi完整分析流程（优化版本）")
        
        try:
            # 1. 构建bowtie2索引
            self.index_prefix = self.build_bowtie2_index()
            
            # 2. 处理所有样本
            for sample in self.samples:
                self.process_sample(sample)
                self.merge_sample_counts(sample)
                self.merge_gRNA_counts(sample)
            
            # 3. 创建输出目录
            os.makedirs(os.path.join(self.result_dir, '05_gene_matrix'), exist_ok=True)
            os.makedirs(os.path.join(self.result_dir, '05_gRNA_matrix'), exist_ok=True)
            
            # 4. 创建基因和gRNA表达矩阵
            gene_matrix_file = self.create_gene_expression_matrix()
            gRNA_matrix_file = self.create_gRNA_expression_matrix()
            
            logger.info("CRISPRi分析流程完成！")
            logger.info(f"基因表达矩阵: {gene_matrix_file}")
            logger.info(f"gRNA表达矩阵: {gRNA_matrix_file}")
            
            return gene_matrix_file, gRNA_matrix_file
            
        except Exception as e:
            logger.error(f"分析流程失败: {e}")
            raise

def main():
    """主函数"""
    # 解析命令行参数
    import argparse
    parser = argparse.ArgumentParser(description='CRISPRi筛选数据分析')
    parser.add_argument('--rawdata-dir', required=True, help='原始数据目录')
    parser.add_argument('--reference-fasta', required=True, help='参考序列文件')
    parser.add_argument('--anno-file', required=True, help='基因注释文件')
    parser.add_argument('--output-dir', required=True, help='输出结果目录')
    parser.add_argument('--threads', type=int, default=10, help='并行线程数')
    parser.add_argument('--test', action='store_true', help='测试模式（限制处理reads数量）')
    parser.add_argument('--max-reads', type=int, default=None, help='最大处理reads数量')
    parser.add_argument('--use-uncompressed', action='store_true', help='使用未压缩的fastq文件')
    
    args = parser.parse_args()
    
    # 设置参数
    use_compressed = not args.use_uncompressed
    max_reads = args.max_reads
    
    if args.test:
        max_reads = 50000  # 测试模式限制5万reads
        logger.info("运行测试模式，限制处理50,000 reads")
    
    try:
        analyzer = CRISPRiAnalyzerOptimized(
            rawdata_dir=args.rawdata_dir,
            reference_fasta=args.reference_fasta,
            anno_file=args.anno_file,
            output_dir=args.output_dir,
            threads=args.threads,
            max_reads=max_reads,
            use_compressed=use_compressed
        )
        gene_matrix_file, gRNA_matrix_file = analyzer.run_complete_analysis()
        print(f"\n分析完成！")
        print(f"基因表达矩阵: {gene_matrix_file}")
        print(f"gRNA表达矩阵: {gRNA_matrix_file}")
        
    except Exception as e:
        print(f"分析失败: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
```