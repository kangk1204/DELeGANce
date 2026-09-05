#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RUN_NAME="${RUN_NAME:-full_run}"
RUN_ROOT="${RUN_ROOT:-$ROOT_DIR/DELeGANce_out/$RUN_NAME}"

# Reuse already-merged FASTQs + fastp JSONs from an existing run to avoid re-running fastp.
SRC_FASTP_OUT="${SRC_FASTP_OUT:-$ROOT_DIR/DELeGANce_out/base_run/01_fastp_out}"

DEFAULT_FASTQ_DIR=""
for cand in \
  "$ROOT_DIR/00_input_fastq" \
  "$ROOT_DIR/00_input_fastq_set1" \
  "$ROOT_DIR/00_input_fastq_set2" \
  "$ROOT_DIR/00_input_fastq_set3" \
  "$ROOT_DIR/00_input_fastq_set4" \
  "$ROOT_DIR/00_input_fastq_set5"
do
  if [[ -d "$cand" ]]; then
    DEFAULT_FASTQ_DIR="$cand"
    break
  fi
done
FASTQ_DIR="${FASTQ_DIR:-$DEFAULT_FASTQ_DIR}"
BBINFO="${BBINFO:-$ROOT_DIR/00_BB_information.txt}"
THREADS="${THREADS:-6}"
MISMATCH_MODE="${MISMATCH_MODE:-hp_op_cp}"

# NOTE: Update these sample column names to match your own decoded matrix,
# or override via env (space-separated): R1_COLS="R1C1 R1C2" R2_COLS="" NEG_COLS="NEG_R1"
IFS=' ' read -r -a R1_COLS <<< "${R1_COLS:-R1C1 R1C2 R1C3 R1C4}"
IFS=' ' read -r -a R2_COLS <<< "${R2_COLS-R2C1 R2C2 R2C3 R2C4}"   # `-` (not `:-`): R2_COLS="" means R1-only
IFS=' ' read -r -a NEG_COLS <<< "${NEG_COLS:-NEG_R1 NEG_R2}"
DEL2_COL="${DEL2_COL:-DEL2}"

usage() {
  cat <<EOF
Usage:
  $(basename "$0") run        # link merged FASTQs/fastp json + run full pipeline
  $(basename "$0") link       # only create run_root + symlinks into 01_fastp_out
  $(basename "$0") tail       # tail pipeline log

Env overrides:
  RUN_NAME, RUN_ROOT, SRC_FASTP_OUT, FASTQ_DIR, BBINFO, THREADS, MISMATCH_MODE, DEL2_COL
  R1_COLS, R2_COLS, NEG_COLS   (space-separated column lists; R2_COLS="" for R1-only)
  DELEGANCE_DISABLE_TORCH      (default 1: hit-caller uses the lightweight CPU GLM fallback,
                                which 03_call_hits.py reports as reduced accuracy; set 0 to use torch)
EOF
}

link_fastp_outputs() {
  if [[ ! -d "$SRC_FASTP_OUT" ]]; then
    echo "[ERROR] SRC_FASTP_OUT not found: $SRC_FASTP_OUT" >&2
    exit 2
  fi
  mkdir -p "$RUN_ROOT/01_fastp_out"
  # Absolute source path: a relative SRC_FASTP_OUT would otherwise produce dangling symlinks
  # (link targets are resolved relative to $RUN_ROOT/01_fastp_out, not to the cwd).
  SRC_FASTP_OUT="$(cd "$SRC_FASTP_OUT" && pwd -P)"

  shopt -s nullglob
  local src
  for src in "$SRC_FASTP_OUT"/*; do
    local base dest
    base="$(basename "$src")"
    dest="$RUN_ROOT/01_fastp_out/$base"
    if [[ -e "$dest" ]]; then
      continue
    fi
    ln -s "$src" "$dest"
  done

  local merged_files json_files
  # Same set of merged-FASTQ names the pipeline/decoder accept (MERGED_FASTQ_RE / AUTO_FASTQ_REGEX):
  # *_merged.fq[.gz], *_merged.fastq[.gz], *.fpmerged.fq[.gz]
  shopt -s nocaseglob
  merged_files=("$RUN_ROOT/01_fastp_out"/*_merged.fq "$RUN_ROOT/01_fastp_out"/*_merged.fq.gz
                "$RUN_ROOT/01_fastp_out"/*_merged.fastq "$RUN_ROOT/01_fastp_out"/*_merged.fastq.gz
                "$RUN_ROOT/01_fastp_out"/*.fpmerged.fq "$RUN_ROOT/01_fastp_out"/*.fpmerged.fq.gz)
  json_files=("$RUN_ROOT/01_fastp_out"/*.fastp.json)
  shopt -u nocaseglob
  local n_merged="${#merged_files[@]}"
  local n_json="${#json_files[@]}"
  shopt -u nullglob
  echo "[INFO] Linked outputs into: $RUN_ROOT/01_fastp_out  (merged=$n_merged, fastp_json=$n_json)"
  if [[ "$n_merged" -lt 1 ]]; then
    echo "[ERROR] No merged FASTQ (*_merged.fq[.gz], *_merged.fastq[.gz], *.fpmerged.fq[.gz]) found under $RUN_ROOT/01_fastp_out" >&2
    exit 2
  fi
}

run_pipeline() {
  mkdir -p "$RUN_ROOT"

  if [[ ! -f "$BBINFO" ]]; then
    echo "[ERROR] BBINFO not found: $BBINFO" >&2
    exit 2
  fi
  if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "[ERROR] FASTQ_DIR not found: $FASTQ_DIR" >&2
    exit 2
  fi

  link_fastp_outputs

  export PYTHONUNBUFFERED=1
  # Default 1 keeps the historical behaviour (CPU GLM fallback in 03_call_hits.py). The pipeline's
  # auto-opt now honours this variable too, so the output dir label matches the device actually used.
  export DELEGANCE_DISABLE_TORCH="${DELEGANCE_DISABLE_TORCH:-1}"

  python3 "$ROOT_DIR/run_delegance_pipeline.py" \
    --fastq-dir "$FASTQ_DIR" \
    --bbinfo "$BBINFO" \
    --output-dir "$RUN_ROOT" \
    --threads "$THREADS" \
    --mismatch "$MISMATCH_MODE" \
    --skip-fastp \
    --r1 "${R1_COLS[@]}" \
    ${R2_COLS[@]:+--r2 "${R2_COLS[@]}"} \
    --neg "${NEG_COLS[@]}" \
    --del2 "$DEL2_COL"
}

cmd="${1:-run}"
case "$cmd" in
  run)  run_pipeline ;;
  link) mkdir -p "$RUN_ROOT"; link_fastp_outputs ;;
  tail)
    if [[ ! -f "$RUN_ROOT/00_pipeline.log" ]]; then
      echo "[ERROR] log not found: $RUN_ROOT/00_pipeline.log" >&2
      exit 2
    fi
    tail -n 200 -f "$RUN_ROOT/00_pipeline.log"
    ;;
  -h|--help|help) usage ;;
  *)
    echo "[ERROR] Unknown command: $cmd" >&2
    usage
    exit 2
    ;;
esac
