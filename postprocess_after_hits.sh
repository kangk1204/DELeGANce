#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 RUN_ROOT [RUN_ROOT ...]" >&2
  exit 2
fi

# Repository-wide default DEL2 column is "DEL2" (README, 03_call_hits.py --del2_col).
DEL2_COL="${DEL2_COL:-DEL2}"
TOP_HITSCORE="${TOP_HITSCORE:-10000}"
PLOT_HEIGHT="${PLOT_HEIGHT:-260}"

runs=("$@")
patterns=()
for run_root in "${runs[@]}"; do
  # run_delegance_pipeline.py passes --run_root as an ABSOLUTE path; match that, anchored at the
  # end of the argument so "run1" does not also wait for "run10". Regex metacharacters escaped.
  abs_root="$(cd "${run_root}" 2>/dev/null && pwd -P || printf '%s' "${run_root}")"
  esc_root="$(printf '%s' "${abs_root}" | sed -e 's/[][\.*^$+?(){}|/]/\\&/g')"
  patterns+=("03_call_hits\.py --run_root ${esc_root}( |\$)")
done

while :; do
  any=0
  for pat in "${patterns[@]}"; do
    if pgrep -f -- "${pat}" >/dev/null; then
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
  norm_base="${run_root}/03_normalized"
  if [[ ! -d "$norm_base" ]]; then
    echo "[ERROR] Missing: $norm_base" >&2
    exit 2
  fi
  annot="$(python3 - <<'PY' "$norm_base"
import os, sys
root = sys.argv[1]
cands = []
for dirpath, dirnames, filenames in os.walk(root):
    # Skip archived copies moved aside by run_delegance_pipeline.py's legacy migration.
    dirnames[:] = [d for d in dirnames if not d.startswith("legacy_")]
    if "05_hybrid_annot.tsv" in filenames:
        p = os.path.join(dirpath, "05_hybrid_annot.tsv")
        cands.append(p)
if not cands:
    sys.exit(2)
cands.sort(key=lambda p: os.path.getmtime(p), reverse=True)
print(cands[0])
PY
  )" || { echo "[ERROR] 05_hybrid_annot.tsv not found under $norm_base" >&2; exit 2; }
  norm="$(dirname "$annot")"
  python3 "$ROOT_DIR/make_display_hybrid_tsv.py" --in_tsv "$annot"
  # beginner_qc_report.html: must NOT be report.html, which is 03_call_hits.py's all-in-one report
  # in the same directory (and the target of index.html's "All-in-one report" link).
  python3 "$ROOT_DIR/export_beginner_qc_report.py" \
    --annot_tsv "$annot" \
    --del2_col "${DEL2_COL}" \
    --out_html "${norm}/beginner_qc_report.html" \
    --out_tsv "${norm}/beginner_qc_tophits.tsv"
  # Optional exporters: a missing openpyxl/xlsxwriter or bokeh must not abort the remaining runs.
  python3 "$ROOT_DIR/export_final_excel.py" --annot_tsv "$annot" --out "${run_root}/final_hits.xlsx" \
    || echo "[WARN] export_final_excel.py failed for ${run_root} (continuing)" >&2
  python3 "$ROOT_DIR/04_build_interactive_report.py" \
    --master_tsv "$annot" \
    --bbinfo "${run_root}/BB_information_fixed.tsv" \
    --out "${norm}/interactive_hits.html" \
    --top_hitscore "${TOP_HITSCORE}" \
    --plot_height "${PLOT_HEIGHT}" \
    || echo "[WARN] 04_build_interactive_report.py failed for ${run_root} (continuing)" >&2
done

echo "[$(date '+%F %T')] Postprocess complete."
