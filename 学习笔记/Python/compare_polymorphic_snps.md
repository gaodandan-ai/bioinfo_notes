脚本任务：把多个样品的 VCF 读进来，提取每个位点的测序深度和突变支持数，算出 AF，然后把不同样品放在一起比较
核心逻辑：读文件 → 拆字段 → 过滤 → 合并 → 输出
	DP：这个位点的总覆盖深度
	AO：支持ALT的reads数量
	AF：ALT allele frequency，计算方式是AF=AO/AP
	AF=20/100=0.2，大约有20%的reads支持突变

该脚本
自动分样品为
	对照组：C*
	实验组：T*
	baseline：C0 and T0

读取每个文件的VCF，对每个位点计算DP/AO/AF
合并样品所有点
输出表：
	long.tsv
	wide.tsv
	events.tsv
	treat_specific.tsv
	control_specific.tsv


```python
#!/usr/bin/env python3
from __future__ import annotations

import argparse    #用来接受命令行参数
import gzip        #用来读取.gz文件
import re          #正则表达块，主要用来识别样品名
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import csv


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")
#自动判断文件是不是。gz文件，因为普通的open()读取不了。gz

def parse_info(info_str: str) -> Dict[str, str]:
    d = {}
    if info_str == "." or not info_str:
        return d
    # 如果 INFO 列是 "." 或空，就返回空字典
    for item in info_str.split(";"):  
    # INFO 里通常是 ; 分隔的多个键值对  
	# 例如: DP=100;AO=20;TYPE=snp
	
        if "=" in item: #如果有等号，说明是K=value的形式
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = "1"
    return d
#parse——info(info_str)把vcf的info列解析成字典
#vcf文件里的格式可能是DP=84;AO=12;RO=72;TYPE=snp
#解析后变成{
#    "DP": "84",
#    "AO": "12",
#    "RO": "72",
#    "TYPE": "snp"

def safe_int(x: Optional[str]) -> Optional[int]:
    if x is None or x in {"", "."}:
        return None  # 如果是 None、空字符串、"."，直接认为没值
    try:
        return int(float(x))  # 先转 float 再转 int，兼容 "20.0" 这种情况
    except Exception:
        return None  #转换失败就返回none,而不是报错
#把字符串转成整数

def pick_alt_depth_and_dp(   
    info: Dict[str, str],
    fmt_keys: List[str],
    sample_fields: List[str],
    alt_index: int = 0,
) -> Tuple[Optional[int], Optional[int]]:

    fmt = dict(zip(fmt_keys, sample_fields)) if fmt_keys and sample_fields else {}
# 把 FORMAT 字段名和样品列值一一对应，做成字典  
# 比如:  
# fmt_keys = ["GT", "DP", "AO"]  
# sample_fields = ["0/1", "100", "20"]  
# 结果:  
# fmt = {"GT": "0/1", "DP": "100", "AO": "20"}
    ao = None
    dp = None
#尽可能从不同格式里的VCF文件中啥选出AO和DP


    # FORMAT: DP
    if "DP" in fmt:
        dp = safe_int(fmt.get("DP"))

    # FORMAT: AO
    if "AO" in fmt:
        ao_vals = fmt.get("AO", "")
        if ao_vals:
            parts = ao_vals.split(",")
            if alt_index < len(parts):
                ao = safe_int(parts[alt_index])

    # FORMAT: AD = ref,alt1,alt2...
    if (ao is None or dp is None) and "AD" in fmt:
        ad_vals = fmt.get("AD", "")
        if ad_vals:
            parts = [safe_int(p) for p in ad_vals.split(",")]
            if len(parts) > 1 + alt_index:
                ao = parts[1 + alt_index]
            if dp is None:
                valid_parts = [p for p in parts if p is not None]
                dp = sum(valid_parts) if valid_parts else None

    # INFO fallback
    if dp is None and "DP" in info:
        dp = safe_int(info.get("DP"))

    if ao is None and "AO" in info:
        ao_vals = info.get("AO", "")
        if ao_vals:
            parts = ao_vals.split(",")
            if alt_index < len(parts):
                ao = safe_int(parts[alt_index])

    return ao, dp


def iter_vcf_records(vcf_path: Path, snp_only: bool = False):  #逐条读取vcf记录
    with open_text(vcf_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue

            cols = line.split("\t")
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
                if snp_only and (len(ref) != 1 or len(alt) != 1):
                    continue

                yield {
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "alt_index": alt_index,
                    "qual": qual,
                    "filter": flt,
                    "info": info,
                    "fmt_keys": fmt_keys,
                    "sample_fields": sample_fields,
                }


def load_sample(                #load_smple()把一个样品的vcf文件全部整理完，整理成一个字典
    vcf_path: Path,
    min_dp: int,
    snp_only: bool = False,
    pass_only: bool = False,
) -> Dict[Tuple[str, int, str, str], Dict[str, float]]:
    out = {}
    for rec in iter_vcf_records(vcf_path, snp_only=snp_only):
        if pass_only and rec["filter"] not in {"PASS", "."}:  #如果你指定了 `--pass_only`，就只保留 `PASS` 或 `.` 的位点
            continue

        ao, dp = pick_alt_depth_and_dp(
            rec["info"], rec["fmt_keys"], rec["sample_fields"], rec["alt_index"]
        )

        if dp is None or dp < min_dp:         #太低深度不保留
            continue
        if ao is None:
            continue
        if dp <= 0:
            continue
#深度过滤或者确实过滤
        af = ao / dp
        key = (rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
        out[key] = {"AO": float(ao), "DP": float(dp), "AF": float(af)}

    return out


def write_tsv(path: Path, header: List[str], rows: List[List[object]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
#将表格写成tsv文件

def detect_vcfs(indir: Path, pattern: str = "*.raw.vcf*") -> Dict[str, Path]:
    samples = {}
    for subdir in sorted(indir.iterdir()):
        if not subdir.is_dir():
            continue

        matched = sorted(subdir.glob(pattern))
        if len(matched) == 1:
            samples[subdir.name] = matched[0]
        elif len(matched) > 1:
            raise SystemExit(f"Multiple matched VCFs in {subdir}: {[str(x) for x in matched]}")
    return samples


def classify_samples(sample_names: List[str]) -> Tuple[str, str, List[str], List[str]]:         #按照命名规则自动分组
    """
    Returns:
      c0, t0, control_samples, treatment_samples
    """
    if "C0" not in sample_names:
        raise SystemExit("Missing required control baseline sample: C0")
    if "T0" not in sample_names:
        raise SystemExit("Missing required treatment baseline sample: T0")

    control_samples = sorted(
        [s for s in sample_names if re.fullmatch(r"C\d+", s)],
        key=lambda x: int(x[1:])
    )
    treatment_samples = sorted(
        [s for s in sample_names if re.fullmatch(r"T\d+", s)],
        key=lambda x: int(x[1:])     #按照数字排序，防止按照字符串排序
    )

    return "C0", "T0", control_samples, treatment_samples


def build_long_table(         #生成长表
    data: Dict[str, Dict[Tuple[str, int, str, str], Dict[str, float]]],
    sample_order: List[str],
    min_af: float,
    max_af: float,
) -> List[List[object]]:
    keys = sorted(set().union(*[set(d.keys()) for d in data.values()]))
    rows = []

    for key in keys:
        chrom, pos, ref, alt = key
        for sid in sample_order:
            rec = data[sid].get(key)
            if rec is None:
                continue
            af = rec["AF"]
            if af < min_af or af > max_af:
                continue
            rows.append([chrom, pos, ref, alt, sid, int(rec["DP"]), int(rec["AO"]), f"{af:.6f}"])
    return rows


def build_wide_table(      #生成宽表
    data: Dict[str, Dict[Tuple[str, int, str, str], Dict[str, float]]],
    sample_order: List[str],
    min_af: float,
    max_af: float,
) -> Tuple[List[str], List[List[object]]]:
    keys = sorted(set().union(*[set(d.keys()) for d in data.values()]))

    header = ["CHROM", "POS", "REF", "ALT"]
    for sid in sample_order:
        header += [f"{sid}_DP", f"{sid}_AO", f"{sid}_AF"]

    rows = []
    for key in keys:
        chrom, pos, ref, alt = key
        row = [chrom, pos, ref, alt]
        afs = []

        for sid in sample_order:
            rec = data[sid].get(key)
            if rec is None:
                row.extend([".", ".", "."])
                afs.append(None)
            else:
                af = rec["AF"]
                row.extend([str(int(rec["DP"])), str(int(rec["AO"])), f"{af:.6f}"])
                afs.append(af)

        if any(a is not None and min_af <= a <= max_af for a in afs):
            rows.append(row)

    return header, rows


def build_events_table(      #给每个位点打上标签，根据C0/T0以及后续C/T样品，给每个位点添加解释标签
    data: Dict[str, Dict[Tuple[str, int, str, str], Dict[str, float]]],
    control_samples: List[str],
    treatment_samples: List[str],
    min_af: float,
) -> Tuple[List[str], List[List[object]]]:
    all_samples = control_samples + treatment_samples
    keys = sorted(set().union(*[set(d.keys()) for d in data.values()]))

    header = ["CHROM", "POS", "REF", "ALT"]
    for sid in all_samples:
        header.append(f"{sid}_AF")
    header.append("LABELS")

    rows = []

    for key in keys:
        chrom, pos, ref, alt = key

        afs = {}
        for sid in all_samples:
            rec = data[sid].get(key)
            afs[sid] = None if rec is None else float(rec["AF"])

        def present(af: Optional[float]) -> bool:
            return af is not None and af >= min_af

        labels = []

        c0_af = afs.get("C0")
        t0_af = afs.get("T0")

        # 在后续样品中是否出现
        present_controls = [s for s in control_samples if s != "C0" and present(afs[s])]
        present_treatments = [s for s in treatment_samples if s != "T0" and present(afs[s])]

        if not present(c0_af) and present_controls:
            labels.append("appears_after_C0")
        if not present(t0_af) and present_treatments:
            labels.append("appears_after_T0")

        if present_controls and not present_treatments:      #是否组特异
            labels.append("control_specific")
        if present_treatments and not present_controls:
            labels.append("treatment_specific")
        if present_controls and present_treatments:
            labels.append("shared_between_groups")

        # 相对各自 baseline 的 AF 升高
        for sid in control_samples:
            if sid == "C0":
                continue
            if c0_af is not None and afs[sid] is not None and afs[sid] - c0_af >= 0.20:      #是否相对baselineAF升高
                labels.append(f"AF_increase_vs_C0:{sid}")

        for sid in treatment_samples:
            if sid == "T0":
                continue
            if t0_af is not None and afs[sid] is not None and afs[sid] - t0_af >= 0.20:
                labels.append(f"AF_increase_vs_T0:{sid}")

        # 哪些样品里存在
        present_samples = [sid for sid in all_samples if present(afs[sid])]
        if present_samples:
            labels.append("present_in:" + ",".join(present_samples))

        row = [chrom, pos, ref, alt]
        for sid in all_samples:
            row.append("." if afs[sid] is None else f"{afs[sid]:.6f}")
        row.append(",".join(labels) if labels else "no_label")
        rows.append(row)

    return header, rows


def build_group_specific_tables(    #挑出特异性位点
    data: Dict[str, Dict[Tuple[str, int, str, str], Dict[str, float]]],
    control_samples: List[str],
    treatment_samples: List[str],
    min_af: float,
) -> Tuple[List[List[object]], List[List[object]]]:
    """
    Returns:
      treatment_specific_rows, control_specific_rows
    """
    all_samples = control_samples + treatment_samples
    keys = sorted(set().union(*[set(d.keys()) for d in data.values()]))

    treatment_specific_rows = []
    control_specific_rows = []

    def present(x: Optional[float]) -> bool:
        return x is not None and x >= min_af

    for key in keys:
        chrom, pos, ref, alt = key

        afs = {}
        for sid in all_samples:
            rec = data[sid].get(key)
            afs[sid] = None if rec is None else float(rec["AF"])

        present_controls = any(present(afs[s]) for s in control_samples if s != "C0")
        present_treatments = any(present(afs[s]) for s in treatment_samples if s != "T0")

        row = [chrom, pos, ref, alt]
        for sid in all_samples:
            row.append("." if afs[sid] is None else f"{afs[sid]:.6f}")

        if present_treatments and not present_controls:
            treatment_specific_rows.append(row)
        if present_controls and not present_treatments:
            control_specific_rows.append(row)

    return treatment_specific_rows, control_specific_rows


def main():
    ap = argparse.ArgumentParser(
        description="Compare multi-sample VCFs using fixed naming rules: C0/T0 baselines, C* control group, T* treatment group."
    )
    ap.add_argument("--indir", required=True, help="Directory containing sample subdirectories")
    ap.add_argument("--pattern", default="*.raw.vcf*", help="VCF glob pattern inside each sample subdirectory")
    ap.add_argument("--outdir", default="snp_compare", help="Output directory")
    ap.add_argument("--min_dp", type=int, default=20, help="Minimum DP")
    ap.add_argument("--min_af", type=float, default=0.01, help="Minimum AF")
    ap.add_argument("--max_af", type=float, default=0.99, help="Maximum AF for long/wide tables")
    ap.add_argument("--snp_only", action="store_true", help="Keep SNPs only")
    ap.add_argument("--pass_only", action="store_true", help="Keep PASS/. variants only")
    args = ap.parse_args()

    samples = detect_vcfs(Path(args.indir), args.pattern)
    if not samples:
        raise SystemExit("No VCF files found.")

    sample_names = list(samples.keys())
    c0, t0, control_samples, treatment_samples = classify_samples(sample_names)

    sample_order = control_samples + treatment_samples
    outdir = Path(args.outdir)

    print("Detected samples:")
    for sid in sample_order:
        print(f"  {sid}\t{samples[sid]}")

    data = {}
    for sid in sample_order:
        print(f"Loading {sid}: {samples[sid]}")
        data[sid] = load_sample(
            samples[sid],
            min_dp=args.min_dp,
            snp_only=args.snp_only,
            pass_only=args.pass_only,
        )

    # 1. long
    long_rows = build_long_table(data, sample_order, args.min_af, args.max_af)
    write_tsv(
        outdir / "polymorphic_snps.long.tsv",
        ["CHROM", "POS", "REF", "ALT", "SAMPLE", "DP", "AO", "AF"],
        long_rows,
    )

    # 2. wide
    wide_header, wide_rows = build_wide_table(data, sample_order, args.min_af, args.max_af)
    write_tsv(outdir / "polymorphic_snps.wide.tsv", wide_header, wide_rows)

    # 3. events
    events_header, events_rows = build_events_table(
        data=data,
        control_samples=control_samples,
        treatment_samples=treatment_samples,
        min_af=args.min_af,
    )
    write_tsv(outdir / "polymorphic_snps.events.tsv", events_header, events_rows)

    # 4. group specific
    group_header = ["CHROM", "POS", "REF", "ALT"] + [f"{sid}_AF" for sid in sample_order]
    treatment_specific_rows, control_specific_rows = build_group_specific_tables(
        data=data,
        control_samples=control_samples,
        treatment_samples=treatment_samples,
        min_af=args.min_af,
    )
    write_tsv(outdir / "treatment_specific.tsv", group_header, treatment_specific_rows)
    write_tsv(outdir / "control_specific.tsv", group_header, control_specific_rows)

    print("Done. Wrote:")
    print(f"  {outdir / 'polymorphic_snps.long.tsv'}")
    print(f"  {outdir / 'polymorphic_snps.wide.tsv'}")
    print(f"  {outdir / 'polymorphic_snps.events.tsv'}")
    print(f"  {outdir / 'treatment_specific.tsv'}")
    print(f"  {outdir / 'control_specific.tsv'}")


if __name__ == "__main__":
    main()
```