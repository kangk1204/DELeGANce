#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RUN_NAME="${RUN_NAME:-KRAS_6DEL_full}"
RUN_ROOT="${RUN_ROOT:-$ROOT_DIR/DELeGANce_out/$RUN_NAME}"

# Reuse already-merged FASTQs + fastp JSONs from an existing run to avoid re-running fastp.
SRC_FASTP_OUT="${SRC_FASTP_OUT:-$ROOT_DIR/DELeGANce_out/KRAS_both/01_fastp_out}"

FASTQ_DIR="${FASTQ_DIR:-$ROOT_DIR/00_KRAS_input_fastq}"
BBINFO="${BBINFO:-$ROOT_DIR/00_6DEL_BB_information_20241013.txt}"
THREADS="${THREADS:-6}"
MISMATCH_MODE="${MISMATCH_MODE:-hp_op_cp}"

R1_COLS=(K_R2C1 K_R2C2 K_R2C3 K_R2C4)
R2_COLS=(K_R3C1 K_R3C2 K_R3C3 K_R3C4)
NEG_COLS=(K_R2C5 K_R3C5)
DEL2_COL="${DEL2_COL:-DEL234}"

usage() {
  cat <<EOF
Usage:
  $(basename "$0") run        # link merged FASTQs/fastp json + run full pipeline
  $(basename "$0") link       # only create run_root + symlinks into 01_fastp_out
  $(basename "$0") tail       # tail pipeline log

Env overrides:
  RUN_NAME, RUN_ROOT, SRC_FASTP_OUT, FASTQ_DIR, BBINFO, THREADS, MISMATCH_MODE, DEL2_COL
EOF
}

link_fastp_outputs() {
  if [[ ! -d "$SRC_FASTP_OUT" ]]; then
    echo "[ERROR] SRC_FASTP_OUT not found: $SRC_FASTP_OUT" >&2
    exit 2
  fi
  mkdir -p "$RUN_ROOT/01_fastp_out"

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
  merged_files=("$RUN_ROOT/01_fastp_out"/*_merged.fq.gz)
  json_files=("$RUN_ROOT/01_fastp_out"/*.fastp.json)
  local n_merged="${#merged_files[@]}"
  local n_json="${#json_files[@]}"
  shopt -u nullglob
  echo "[INFO] Linked outputs into: $RUN_ROOT/01_fastp_out  (merged=$n_merged, fastp_json=$n_json)"
  if [[ "$n_merged" -lt 1 ]]; then
    echo "[ERROR] No *_merged.fq.gz found under $RUN_ROOT/01_fastp_out" >&2
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
  export DELEGANCE_DISABLE_TORCH="${DELEGANCE_DISABLE_TORCH:-1}"

  python3 "$ROOT_DIR/run_delegance_pipeline.py" \
    --fastq-dir "$FASTQ_DIR" \
    --bbinfo "$BBINFO" \
    --output-dir "$RUN_ROOT" \
    --threads "$THREADS" \
    --mismatch "$MISMATCH_MODE" \
    --skip-fastp \
    --r1 "${R1_COLS[@]}" \
    --r2 "${R2_COLS[@]}" \
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
