目录结构：
```
/data/gaodd/cg_seq_run_vcf/raw_vcf_data/
├── T32/
│   └── T32.raw.vcf.gz
├── T33/
│   └── T33.raw.vcf.gz
├── T34/
│   └── T34.raw.vcf.gz
...

```

run
```bash
python compare_snps.py \
  --indir /data/gaodd/cg_seq_run_vcf/raw_vcf_data \
  -o /data/gaodd/cg_seq_run_vcf/results/snp_compare.tsv

```

scripts

```python
#!/usr/bin/env python3

import argparse
from pathlib import Path
import csv
import gzip


def open_vcf(vcf_file):
    """支持普通 vcf 和 vcf.gz"""
    if str(vcf_file).endswith(".gz"):
        return gzip.open(vcf_file, "rt")
    else:
        return open(vcf_file, "r")


def parse_vcf(vcf, min_dp=0):
    data = {}

    with open_vcf(vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            chrom, pos, _, ref, alt, _, _, info, fmt, sample = fields[:10]

            fmt_keys = fmt.split(":")
            sample_vals = sample.split(":")
            fmt_dict = dict(zip(fmt_keys, sample_vals))

            dp = None
            ao = None

            if "DP" in fmt_dict and fmt_dict["DP"] not in {".", ""}:
                try:
                    dp = int(fmt_dict["DP"])
                except ValueError:
                    dp = None

            if "AO" in fmt_dict and fmt_dict["AO"] not in {".", ""}:
                try:
                    ao = int(fmt_dict["AO"].split(",")[0])
                except ValueError:
                    ao = None

            if dp is not None and ao is not None and dp >= min_dp and dp > 0:
                af = ao / dp
                key = (chrom, pos, ref, alt)

                data[key] = {
                    "DP": dp,
                    "AO": ao,
                    "AF": af
                }

    return data


def find_samples(indir):
    """
    自动识别样品：
    /path/to/indir/T32/T32.raw.vcf.gz
    /path/to/indir/T33/T33.raw.vcf.gz
    """
    samples = {}

    for subdir in sorted(indir.iterdir()):
        if not subdir.is_dir():
            continue

        sample_name = subdir.name
        vcf_path = subdir / f"{sample_name}.raw.vcf.gz"

        if vcf_path.exists():
            samples[sample_name] = vcf_path
        else:
            # 兼容非压缩版本
            vcf_path_unzipped = subdir / f"{sample_name}.raw.vcf"
            if vcf_path_unzipped.exists():
                samples[sample_name] = vcf_path_unzipped

    return samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--indir",
        default="/data/gaodd/cg_seq_run_vcf/raw_vcf_data",
        help="包含各样品子目录的总目录"
    )
    parser.add_argument("--min_af", type=float, default=0.01)
    parser.add_argument("--max_af", type=float, default=0.99)
    parser.add_argument("--min_dp", type=int, default=20)
    parser.add_argument("--out", "-o", default="snp_compare.tsv")

    args = parser.parse_args()

    indir = Path(args.indir)

    if not indir.exists():
        raise FileNotFoundError(f"Input directory not found: {indir}")

    samples = find_samples(indir)

    if not samples:
        raise FileNotFoundError(
            f"No sample VCF files found in {indir}"
        )

    print("Detected samples:")
    for s, v in samples.items():
        print(f"  {s}\t{v}")

    data = {}
    for s, v in samples.items():
        print("Loading", v)
        data[s] = parse_vcf(v, min_dp=args.min_dp)

    all_sites = set()
    for s in data:
        all_sites.update(data[s].keys())

    rows = []

    for site in sorted(all_sites, key=lambda x: (x[0], int(x[1]), x[2], x[3])):
        chrom, pos, ref, alt = site
        row = [chrom, pos, ref, alt]
        afs = {}

        for s in samples:
            rec = data[s].get(site)

            if rec:
                dp = rec["DP"]
                ao = rec["AO"]
                af = rec["AF"]
                row.extend([dp, ao, round(af, 4)])
                afs[s] = af
            else:
                row.extend([".", ".", "."])
                afs[s] = None

        if any(a is not None and args.min_af <= a <= args.max_af for a in afs.values()):
            rows.append(row)

    header = ["CHROM", "POS", "REF", "ALT"]
    for s in samples:
        header.extend([f"{s}_DP", f"{s}_AO", f"{s}_AF"])

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)

    print("Finished. Output:", args.out)


if __name__ == "__main__":
    main()
    
```