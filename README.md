# DELeGANce Screening Pipeline

DELeGANce is an end-to-end pipeline for DNA-encoded library (DEL) sequencing analysis: FASTQ preprocessing, barcode decoding, hit calling (GLM-based), and interactive reporting. This repository provides a ready-to-run workflow you can use as a template for new targets. The target is fully user-defined via your input FASTQs and sample columns.

---

## Quick Start

1) **Create a clean environment**
```bash
conda create -n delegance python=3.10
conda activate delegance
conda install -c conda-forge -c bioconda fastp
conda install -c conda-forge numpy pandas scipy matplotlib bokeh
# Optional (recommended for speed + structure images)
conda install -c pytorch pytorch torchvision torchaudio
conda install -c conda-forge rdkit
```

2) **Run the pipeline**
```bash
python3 run_delegance_pipeline.py \
  --fastq-dir 00_input_fastq \
  --bbinfo 00_BB_information.txt \
  --output-dir DELeGANce_out/my_run \
  --threads 6 \
  --mismatch hp_op_cp \
  --r1 R1C1 R1C2 R1C3 R1C4 \
  --r2 R2C1 R2C2 R2C3 R2C4 \
  --neg NEG_R1 NEG_R2 \
  --del2 DEL2
```

3) **Open the results**
- `DELeGANce_out/my_run/index.html`
- `DELeGANce_out/my_run/03_normalized/.../interactive_hits.html`

---

## Repository Layout

- `run_delegance_pipeline.py` - Orchestrator for preprocess -> decode -> hit calling.
- `01_preprocess_reads.pl` - FASTQ cleaning + barcode table reconciliation (fastp wrapper).
- `02_decode_reads.pl` - Maps merged reads to tags and builds raw/scaled count matrices.
- `03_call_hits.py` - GLM hit calling + plots + static HTML report.
- `04_build_interactive_report.py` - Interactive Bokeh report (Top-N explorer).
- `export_beginner_qc_report.py` - Beginner-friendly QC HTML/TSV.
- `export_final_excel.py` - Excel export with a guide tab (for wet-lab review).
- `run_*_autopilot.sh` - Optional target-specific helper script (if present).
- `00_*_input_fastq/` - Input FASTQs (expected format: `<SAMPLE>_1.fastq.gz`, `<SAMPLE>_2.fastq.gz`).
- `00_*_BB_information*.txt` - BB metadata tables.
- `DELeGANce_out/` - Outputs (TSV/plots/HTML/logs).

---

## Running the Pipeline (Step by Step)

### 1) Verify toolchain
```bash
fastp --version
perl -v
python3 --version
```

### 2) Make sure your FASTQs are named consistently
Each sample should have paired reads:
```
<SAMPLE>_1.fastq.gz
<SAMPLE>_2.fastq.gz
```
The `<SAMPLE>` names must match the columns you pass to the runner (`--r1`, `--r2`, `--neg`, `--del2`).
If you do not know the column names yet, run preprocess+decode first and inspect the header of
`DELeGANce_out/<run_name>/02_decoded/raw_counts_matrix.tsv`, then run hit calling only.

### 3) Launch a full run
```bash
python3 run_delegance_pipeline.py \
  --fastq-dir 00_input_fastq \
  --bbinfo 00_BB_information.txt \
  --output-dir DELeGANce_out/my_run \
  --threads 6 \
  --mismatch hp_op_cp \
  --r1 R1C1 R1C2 R1C3 R1C4 \
  --r2 R2C1 R2C2 R2C3 R2C4 \
  --neg NEG_R1 NEG_R2 \
  --del2 DEL2
```

---

## Outputs (Where to Look)

After a run, open `DELeGANce_out/<run_name>/index.html` for a summary page.
Key files:

- `03_normalized/<preset>/05_hybrid_annot.tsv` - main results table
- `03_normalized/<preset>/report.html` - static summary report
- `03_normalized/<preset>/interactive_hits.html` - interactive Top-N explorer
- `Beginner_QC_Report.html` - quick QC overview (includes Tier/PickGroup recommendations)
- `Beginner_QC_TopHits.tsv` - QC table with Tier/PickGroup columns
- `DELeGANce_final_results.xlsx` - Excel export for wet-lab review

---

## Generic Setup Template (Target-Agnostic)

### 1) Folder layout
```
project_root/
  00_input_fastq/          # paired FASTQs
    SAMPLE_A_1.fastq.gz
    SAMPLE_A_2.fastq.gz
    SAMPLE_B_1.fastq.gz
    SAMPLE_B_2.fastq.gz
  00_BB_information.txt     # BB info (tab-delimited)
  DELeGANce_out/
```

### 2) BB information file format (tab-delimited)
Minimum 7 columns in this exact order (header optional):
```
type    seq     bb_id   cycle   tag_id  lib_id  smiles
```
Notes:
- `type` includes CODON/HP/OP/CP (CODON rows are used for BB mapping)
- `cycle` is 1/2/3/4 for each BB
- `lib_id` is required if multiple libraries are present

### 3) Sample column names
After decode, column names come from FASTQ sample names.
Check header of:
```
DELeGANce_out/<run_name>/02_decoded/raw_counts_matrix.tsv
```
Use those exact names in `--r1`, `--r2`, `--neg`, `--del2`.

### 4) Minimal run command (example template)
```bash
python3 run_delegance_pipeline.py \
  --fastq-dir 00_input_fastq \
  --bbinfo 00_BB_information.txt \
  --output-dir DELeGANce_out/my_run \
  --threads 6 \
  --mismatch hp_op_cp \
  --r1 R1C1 R1C2 R1C3 R1C4 \
  --r2 R2C1 R2C2 R2C3 R2C4 \
  --neg NEG_R1 NEG_R2 \
  --del2 DEL2
```

### 5) Quick validation
- `02_decoded/raw_counts_matrix.tsv` has `lib_id` and your sample columns
- `03_normalized/.../05_hybrid_annot.tsv` is produced
- `Beginner_QC_Report.html` opens without errors

---

## Useful Stand-Alone Commands

### Hit calling only (after decode)
```bash
python3 03_call_hits.py \
  --run_root DELeGANce_out/my_run \
  --r1_cols R1C1 R1C2 R1C3 R1C4 \
  --r2_cols R2C1 R2C2 R2C3 R2C4 \
  --neg_r1_col NEG_R1 \
  --neg_r2_col NEG_R2 \
  --del2_col DEL2 \
  --preset balanced
```

### Regenerate the interactive report
```bash
python3 04_build_interactive_report.py \
  --master_tsv DELeGANce_out/my_run/03_normalized/<preset>/05_hybrid_annot.tsv \
  --bbinfo DELeGANce_out/my_run/BB_information_fixed.tsv \
  --out DELeGANce_out/my_run/03_normalized/<preset>/interactive_hits.html
```

### Export beginner QC + Excel
```bash
python3 export_beginner_qc_report.py \
  --run_root DELeGANce_out/my_run \
  --out_html DELeGANce_out/Beginner_QC_Report.html \
  --out_tsv DELeGANce_out/Beginner_QC_TopHits.tsv \
  --neg_high_quantile 0.90 \
  --recommend_a 50 \
  --recommend_b 20 \
  --recommend_diverse 50 \
  --diverse_key BB1_BB2_BB3

python3 export_final_excel.py \
  --run_root DELeGANce_out/my_run \
  --out DELeGANce_out/DELeGANce_final_results.xlsx \
  --top_n 1000
```

### Run hit calling + postprocess for multiple runs
```bash
bash run_hits_then_postprocess.sh \
  --run-root DELeGANce_out/run_A \
  --run-root DELeGANce_out/run_B \
  --r1 "R1C1 R1C2 R1C3 R1C4" \
  --r2 "R2C1 R2C2 R2C3 R2C4" \
  --neg "NEG_R1 NEG_R2" \
  --del2 DEL2
```

### Subsample paired FASTQs (for fast testing)
```bash
python3 subsample_fastq_pairs.py \
  --input-dir 00_input_fastq \
  --output-dir 00_input_fastq_subsampled \
  --n-pairs 50000 \
  --mode random \
  --seed 42
```

---

## Troubleshooting

- **fastp not found**: install it or add `--skip-fastp` (not recommended for fresh data).
- **Torch/CUDA issues**: the pipeline falls back to CPU automatically. You can force CPU by setting `DELEGANCE_DISABLE_TORCH=1`.
- **RDKit missing**: only affects molecule thumbnails in the interactive report.
- **Large matrices**: the runner auto-selects GLM top mode for huge datasets; tune with `--glm_mode` and `--glm_top_*`.

---

## GitHub Hygiene

- Outputs and FASTQs should not be committed. Keep them locally.
- A `.gitignore` is included to prevent large data from being tracked.

---

## Citation / Notes

If you use this pipeline in a manuscript, cite the DEL screening methodology and note the following:
- R2 negative controls are treated as the primary NEG for final hit gating.
- R1 NEG is reported as QC only (to show round-to-round distribution shifts).

---

## License

This project is released under a non-commercial academic research license.
See `LICENSE.txt` for the full terms.

---

If you need a tailored README for another target, tell me the folder names and sample columns and I will generate it.
