
# CRISPRi Screening Analysis Pipeline

This project provides an automated pipeline for analyzing CRISPRi screening data, supporting both Single Gene and Dual Gene screening workflows.

##  Features

- **Single Gene Screening Analysis**: Process data from single sgRNA libraries.
- **Dual Gene Screening Analysis**: Process data from dual sgRNA libraries (combinatorial screening).
- **Automated Workflow**: Automated generation of count matrices from raw sequencing data (FastQ).
- **Multi-threading Support**: Parallel processing to improve analysis efficiency.

##  Installation

This project uses Conda for environment management. Please ensure Anaconda or Miniconda is installed.

1. **Clone the repository (if not already downloaded)**
   ```bash
   cd /data/zuoll/1.project/00.CRISPRi/crispri-screening-pipeline
   ```

2. **Create Conda Environment**
   Create the runtime environment using the `environment.yml` file in the project root:
   ```bash
   conda env create -f environment.yml
   ```

3. **Activate Environment**
   ```bash
   conda activate crispri
   # Note: Check the 'name' field in environment.yml for the exact environment name, usually 'crispri' or similar.
   ```

##  Directory Structure

```
.
├── docs/               # Documentation and flowchart (workflow.png)
├── scripts/            # Analysis scripts
│   ├── crispri_single_gene_screening_pipeline.py  # Single gene pipeline script
│   └── crispri_dual_gene_screening_pipeline.py       # Dual gene pipeline script
├── environment.yml     # Environment configuration file
└── README.md           # Documentation
```

##  Usage

All analysis scripts are located in the `scripts/` directory.

### 1. Single Gene Screening Analysis

For analyzing single sgRNA library data.

**Basic Usage:**
```bash
python scripts/crispri_single_gene_screening_pipeline.py \
    --rawdata-dir <RAW_DATA_DIR> \
    --reference-fasta <REFERENCE_FASTA> \
    --anno-file <ANNOTATION_FILE> \
    --output-dir <OUTPUT_DIR> \
    --threads 8
```

**Arguments:**

| Argument             | Required | Description                                                               |
| :------------------- | :------- | :------------------------------------------------------------------------ |
| `--rawdata-dir`      | ✅        | Directory path containing raw FastQ files                                 |
| `--reference-fasta`  | ✅        | Reference sequence FASTA file (usually sgRNA library sequences)           |
| `--anno-file`        | ✅        | Gene annotation file (mapping sgRNAs to genes)                            |
| `--output-dir`       | ✅        | Directory for output results                                              |
| `--threads`          | ❌        | Number of parallel threads (Default: 1)                                   |
| `--test`             | ❌        | Enable test mode (process only a small number of reads for quick testing) |
| `--max-reads`        | ❌        | Maximum number of reads to process                                        |
| `--use-uncompressed` | ❌        | Add this flag if input fastq files are uncompressed                       |

### 2. Dual Gene Screening Analysis

For analyzing dual sgRNA library (combinatorial screening) data.

**Basic Usage:**
```bash
python scripts/crispri_dual_gene_screening_pipeline.py \
    --rawdata-dir <RAW_DATA_DIR> \
    --reference-fasta <REFERENCE_FASTA> \
    --anno-file <ANNOTATION_FILE> \
    --output-dir <OUTPUT_DIR> \
    --threads 8
```

**Arguments:**
(Same arguments as the single gene pipeline)

##  Workflow

The workflow diagram is located at `docs/workflow.png`, detailing the data processing steps.

##  Examples

Assuming your data is in `data/raw`, reference sequences in `refs/library.fa`, annotation file in `refs/library.csv`, and output to `results/`.

**Run Single Gene Pipeline:**
```bash
conda activate crispri
cd /data/gaodd/CRISPRi_thermal/crispri_screening_pipeline

nohup python3 scripts/single_gene/crispri_single_gene_screening_pipeline.py \
--rawdata-dir /data/gaodd/CRISPRi_thermal/raw_data/single_gene \
--reference-fasta /data/gaodd/CRISPRi_thermal/library/library3_slt.fasta \
--anno-file /data/gaodd/CRISPRi_thermal/library/cgl-anno.csv \
--output-dir /data/gaodd/CRISPRi_thermal/results/single_gene_results2 \
--threads 8 \
> run.log 2>&1 &

tail -f run.log
```

**Run in Test Mode:**
```bash
python scripts/crispri_single_gene_screening_pipeline.py \
    --rawdata-dir ./data/raw \
    --reference-fasta ./refs/library.fa \
    --anno-file ./refs/library.csv \
    --output-dir ./results/test_run \
    --test
```