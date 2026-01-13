#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 RUN_ROOT [RUN_ROOT ...]" >&2
  exit 2
fi

DEL2_COL="${DEL2_COL:-DEL234}"
TOP_HITSCORE="${TOP_HITSCORE:-10000}"
PLOT_HEIGHT="${PLOT_HEIGHT:-260}"

runs=("$@")
patterns=()
for run_root in "${runs[@]}"; do
  patterns+=("03_call_hits.py --run_root ${run_root}")
done

while :; do
  any=0
  for pat in "${patterns[@]}"; do
    if pgrep -f "${pat}" >/dev/null; then
      any=1
      break
    fi
  done
  if [ "${any}" -eq 0 ]; then
    break
  fi
  echo "[$(date '+%F %T')] Waiting for hit-caller runs to finish..."
  sleep 60
done

for run_root in "${runs[@]}"; do
  norm="${run_root}/03_normalized/glm_full_dev_cpu_fp64"
  python3 make_display_hybrid_tsv.py --run_root "${run_root}"
  python3 export_beginner_qc_report.py \
    --run_root "${run_root}" \
    --del2_col "${DEL2_COL}" \
    --out_html "${norm}/report.html" \
    --out_tsv "${norm}/report.tsv"
  python3 export_final_excel.py --run_root "${run_root}" --out "${run_root}/final_hits.xlsx"
  python3 04_build_interactive_report.py \
    --master_tsv "${norm}/05_hybrid_annot.tsv" \
    --bbinfo "${run_root}/BB_information_fixed.tsv" \
    --out "${norm}/interactive_hits.html" \
    --top_hitscore "${TOP_HITSCORE}" \
    --plot_height "${PLOT_HEIGHT}"
done

echo "[$(date '+%F %T')] Postprocess complete."
