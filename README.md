# DELeGANce KRAS Screening Pipeline

DELeGANce is an end-to-end pipeline for DNA-encoded library (DEL) sequencing analysis: FASTQ preprocessing, barcode decoding, hit calling (GLM-based), and interactive reporting. This repository is a ready-to-run KRAS workflow with example input folders and scripts you can use as a template for new targets.

---

## Quick Start (KRAS demo)

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

2) **Run the pipeline (KRAS_both)**
```bash
python3 run_delegance_pipeline.py \
  --fastq-dir 00_KRAS_input_fastq \
  --bbinfo 00_DELeGANce_KRASMAT2A_BB_information.txt \
  --output-dir DELeGANce_out/KRAS_both \
  --threads 6 \
  --mismatch hp_op_cp \
  --r1 K_R2C1 K_R2C2 K_R2C3 K_R2C4 \
  --r2 K_R3C1 K_R3C2 K_R3C3 K_R3C4 \
  --neg K_R2C5 K_R3C5 \
  --del2 DEL234
```

3) **Open the results**
- `DELeGANce_out/KRAS_both/index.html`
- `DELeGANce_out/KRAS_both/03_normalized/.../interactive_hits.html`

---

## Repository Layout

- `run_delegance_pipeline.py` - Orchestrator for preprocess -> decode -> hit calling.
- `01_preprocess_reads.pl` - FASTQ cleaning + barcode table reconciliation (fastp wrapper).
- `02_decode_reads.pl` - Maps merged reads to tags and builds raw/scaled count matrices.
- `03_call_hits.py` - GLM hit calling + plots + static HTML report.
- `04_build_interactive_report.py` - Interactive Bokeh report (Top-N explorer).
- `export_beginner_qc_report.py` - Beginner-friendly QC HTML/TSV.
- `export_final_excel.py` - Excel export with a guide tab (for wet-lab review).
- `run_kras_6del_full_autopilot.sh` - Convenience wrapper for the KRAS_6DEL_full run.
- `00_*_input_fastq/` - Input FASTQs (expected format: `<SAMPLE>_1.fastq.gz`, `<SAMPLE>_2.fastq.gz`).
- `00_DELeGANce_KRASMAT2A_BB_information.txt`, `00_6DEL_BB_information_20241013.txt` - BB metadata tables.
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
  --fastq-dir 00_KRAS_input_fastq \
  --bbinfo 00_DELeGANce_KRASMAT2A_BB_information.txt \
  --output-dir DELeGANce_out/KRAS_test_run \
  --threads 6 \
  --mismatch hp_op_cp \
  --r1 K_R2C1 K_R2C2 K_R2C3 K_R2C4 \
  --r2 K_R3C1 K_R3C2 K_R3C3 K_R3C4 \
  --neg K_R2C5 K_R3C5 \
  --del2 DEL234
```

### 4) Optional: KRAS_6DEL_full shortcut
This script reuses merged FASTQs from `KRAS_both` and runs the 6-DEL BB table.
```bash
bash run_kras_6del_full_autopilot.sh run
```
Environment overrides are supported:
```bash
RUN_NAME=KRAS_6DEL_full \
FASTQ_DIR=00_KRAS_input_fastq \
BBINFO=00_6DEL_BB_information_20241013.txt \
bash run_kras_6del_full_autopilot.sh run
```

---

## Outputs (Where to Look)

After a run, open `DELeGANce_out/<run_name>/index.html` for a summary page.
Key files:

- `03_normalized/<preset>/05_hybrid_annot.tsv` - main results table
- `03_normalized/<preset>/report.html` - static summary report
- `03_normalized/<preset>/interactive_hits.html` - interactive Top-N explorer
- `Beginner_QC_Report.html` - quick QC overview
- `DELeGANce_final_results.xlsx` - Excel export for wet-lab review

---

## Useful Stand-Alone Commands

### Hit calling only (after decode)
```bash
python3 03_call_hits.py \
  --run_root DELeGANce_out/my_run \
  --r1_cols K_R2C1 K_R2C2 K_R2C3 K_R2C4 \
  --r2_cols K_R3C1 K_R3C2 K_R3C3 K_R3C4 \
  --neg_r1_col K_R2C5 \
  --neg_r2_col K_R3C5 \
  --del2_col DEL234 \
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
  --run_root DELeGANce_out/KRAS_both \
  --run_root DELeGANce_out/KRAS_6DEL_full \
  --out_html DELeGANce_out/Beginner_QC_Report.html \
  --out_tsv DELeGANce_out/Beginner_QC_TopHits.tsv

python3 export_final_excel.py \
  --run_root DELeGANce_out/KRAS_both \
  --run_root DELeGANce_out/KRAS_6DEL_full \
  --out DELeGANce_out/DELeGANce_final_results.xlsx \
  --top_n 1000
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

If you need a tailored README for another target (e.g., GPCR, MAT2A), tell me the desired folder names and sample columns and I will generate it.
