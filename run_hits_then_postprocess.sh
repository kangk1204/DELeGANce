#!/usr/bin/env bash
set -euo pipefail

# Resolve sibling scripts relative to this file so the script works from any cwd.
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export PYTHONUNBUFFERED=1
export TORCHDYNAMO_DISABLE=1
export TORCH_COMPILE_DISABLE=1

usage() {
  cat <<EOF
Usage:
  $(basename "$0") --run-root <RUN_ROOT> [--run-root <RUN_ROOT> ...] \\
    --r1 "R1C1 R1C2" [--r2 "R2C1 R2C2"] --neg "NEG_R1 NEG_R2" --del2 DEL2

Env overrides:
  R1_COLS, R2_COLS, NEG_COLS, DEL2_COL
EOF
}

RUN_ROOTS=()
R1_COLS_STR="${R1_COLS:-}"
R2_COLS_STR="${R2_COLS:-}"
NEG_COLS_STR="${NEG_COLS:-}"
DEL2_COL="${DEL2_COL:-}"

while [[ $# -gt 0 ]]; do
  # Options that take a value must have one (otherwise set -u reports a bare "$2: unbound variable").
  case "$1" in
    --run-root|--r1|--r2|--neg|--del2)
      if [[ $# -lt 2 ]]; then
        echo "[ERROR] $1 requires a value." >&2
        usage; exit 2
      fi ;;
  esac
  case "$1" in
    --run-root)
      RUN_ROOTS+=("$2"); shift 2 ;;
    --r1)
      R1_COLS_STR="$2"; shift 2 ;;
    --r2)
      R2_COLS_STR="$2"; shift 2 ;;
    --neg)
      NEG_COLS_STR="$2"; shift 2 ;;
    --del2)
      DEL2_COL="$2"; shift 2 ;;
    -h|--help|help)
      usage; exit 0 ;;
    *)
      echo "[ERROR] Unknown arg: $1" >&2
      usage; exit 2 ;;
  esac
done

if [[ ${#RUN_ROOTS[@]} -eq 0 ]]; then
  echo "[ERROR] --run-root is required." >&2
  usage; exit 2
fi
if [[ -z "$R1_COLS_STR" || -z "$NEG_COLS_STR" || -z "$DEL2_COL" ]]; then
  echo "[ERROR] --r1, --neg, --del2 are required." >&2
  usage; exit 2
fi

IFS=' ' read -r -a R1_COLS <<< "$R1_COLS_STR"
IFS=' ' read -r -a R2_COLS <<< "$R2_COLS_STR"
IFS=' ' read -r -a NEG_COLS <<< "$NEG_COLS_STR"

run_one() {
  local run_root="$1"
  # The fixed folder label promises "GLM full / CPU / float64"; pin those knobs so auto-opt
  # cannot silently switch to glm top / CUDA for large matrices while keeping this name.
  local outdir="${run_root}/03_normalized/glm_full_dev_cpu_fp64"
  local cmd=(python3 "$ROOT_DIR/run_delegance_pipeline.py"
    --only hit
    --force-hit
    --auto-opt 0
    --glm-mode full
    --device cpu
    --dtype float64
    --output-dir "${run_root}"
    --hit-out "${outdir}"
    --r1 "${R1_COLS[@]}"
    --neg "${NEG_COLS[@]}"
    --del2 "${DEL2_COL}"
  )
  if [[ ${#R2_COLS[@]} -gt 0 && -n "${R2_COLS_STR}" ]]; then
    cmd+=(--r2 "${R2_COLS[@]}")
  fi
  "${cmd[@]}"
}

for run_root in "${RUN_ROOTS[@]}"; do
  run_one "${run_root}"
done

# DEL2_COL must reach the child via the environment: postprocess_after_hits.sh reads it from env only.
DEL2_COL="${DEL2_COL}" bash "$ROOT_DIR/postprocess_after_hits.sh" "${RUN_ROOTS[@]}"
