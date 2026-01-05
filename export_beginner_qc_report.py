#!/usr/bin/env python3
import argparse
import html
import re
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd


def _resolve_annot(run_root: Path, prefer_dir: str) -> Path:
    if prefer_dir:
        preferred = run_root / prefer_dir / "05_hybrid_annot.tsv"
        if preferred.exists():
            return preferred
    fallback = run_root / "03_normalized" / "05_hybrid_annot.tsv"
    if fallback.exists():
        return fallback

    candidates = list(run_root.rglob("05_hybrid_annot.tsv"))
    if not candidates:
        raise FileNotFoundError(f"[ERROR] 05_hybrid_annot.tsv not found under {run_root}")
    if len(candidates) == 1:
        return candidates[0]
    msg = "\n".join(str(c) for c in sorted(candidates))
    raise FileNotFoundError(
        "[ERROR] Multiple 05_hybrid_annot.tsv files found. "
        "Specify --annot_tsv or --prefer_dir.\n" + msg
    )


def _guess_run_name(annot_path: Path) -> str:
    parts = list(annot_path.parts)
    if "03_normalized" in parts:
        idx = parts.index("03_normalized")
        if idx > 0:
            return parts[idx - 1]
    return annot_path.parent.name


def _pick_col(cols: List[str], preferred: List[str]) -> Optional[str]:
    for c in preferred:
        if c in cols:
            return c
    return None


def _add_hit_score(df: pd.DataFrame) -> pd.DataFrame:
    if "HitScore" in df.columns:
        df["HitScore"] = pd.to_numeric(df["HitScore"], errors="coerce")
        return df
    if "HitScore_GLM" in df.columns:
        df["HitScore"] = pd.to_numeric(df["HitScore_GLM"], errors="coerce")
        return df
    if "HitScore_RS" in df.columns:
        df["HitScore"] = pd.to_numeric(df["HitScore_RS"], errors="coerce")
        return df
    raise ValueError("[ERROR] No HitScore/HitScore_GLM/HitScore_RS column found.")


def _normalize_bool(series: pd.Series) -> pd.Series:
    if series.dtype == object:
        s = series.astype(str).str.lower()
        return s.isin(["true", "1", "yes", "y", "t"])
    return series.astype(bool)


def _auto_del2_col(cols: List[str]) -> Optional[str]:
    # Prefer non-CPM DEL columns
    for c in cols:
        if c.startswith("DEL") and not c.endswith("_CPM"):
            return c
    return None


_LIBDEL_RE = re.compile(r"_LIB[\w\.-]+", re.IGNORECASE)


def _strip_libdel(value) -> str:
    if pd.isna(value):
        return ""
    return _LIBDEL_RE.sub("", str(value))


def _pick_bb_value(row: pd.Series, base: str) -> str:
    for c in (base, f"{base}_x", f"{base}_y"):
        if c in row.index:
            return _strip_libdel(row.get(c, ""))
    return ""


def _make_display_id(lib: str, bb1: str, bb2: str, bb3: str, bb4: str, fallback_id: str) -> str:
    lib = (lib or "").strip()
    if lib:
        return f"{lib}_{bb1}_{bb2}_{bb3}_{bb4}"
    return _strip_libdel(fallback_id)


def _qc_flags(row: pd.Series, del2_col: Optional[str]) -> str:
    flags = []

    if del2_col and pd.notna(row.get(del2_col)):
        if float(row.get(del2_col, 0)) < 10:
            flags.append("LOW_DEL2")

    lfc_neg = row.get("LFC_NEG_centered")
    if pd.notna(lfc_neg) and float(lfc_neg) >= 0:
        flags.append("NEG_HIGH")

    lfc_neg_r1 = row.get("LFC_NEG_centered_R1")
    if pd.notna(lfc_neg_r1) and float(lfc_neg_r1) >= 0:
        flags.append("NEG_R1_HIGH")

    lfc_r1 = row.get("LFC_R1_vs_DEL2")
    lfc_r2 = row.get("LFC_R2_vs_DEL2")
    if pd.notna(lfc_r1) and pd.notna(lfc_r2):
        if float(lfc_r2) < float(lfc_r1) - 1.0:
            flags.append("R2_DROP")

    for base in ["BB1", "BB2", "BB3"]:
        val = _pick_bb_value(row, base).strip().upper()
        if val in ["", "NA", "NAN", "NONE"]:
            flags.append("BB_MISSING")
            break

    neg_fail = row.get("NEG_hard_fail")
    if pd.notna(neg_fail) and bool(neg_fail):
        flags.append("NEG_HARD_FAIL")

    return "OK" if not flags else ";".join(flags)


def _build_top_table(df: pd.DataFrame, top_n: int) -> Tuple[pd.DataFrame, str]:
    df = _add_hit_score(df.copy())
    df = df[np.isfinite(df["HitScore"])].sort_values("HitScore", ascending=False)
    df = df.head(top_n).copy()

    cols = list(df.columns)
    id_col = _pick_col(cols, ["ID", "id", "ID_x", "id_x", "ID_y", "id_y"])
    lib_col = _pick_col(cols, ["LIB_ID", "LIB_ID_x", "LIB_ID_y"])

    bb_cols = None
    if all(c in cols for c in ["BB1_x", "BB2_x", "BB3_x", "BB4_x"]):
        bb_cols = ["BB1_x", "BB2_x", "BB3_x", "BB4_x"]
    elif all(c in cols for c in ["BB1", "BB2", "BB3", "BB4"]):
        bb_cols = ["BB1", "BB2", "BB3", "BB4"]

    del2_col = _auto_del2_col(cols)

    display_cols = []
    if lib_col:
        display_cols.append(lib_col)
    if id_col:
        display_cols.append(id_col)
    if "cycles" in cols:
        display_cols.append("cycles")
    if bb_cols:
        display_cols.extend(bb_cols)

    if del2_col:
        display_cols.append(del2_col)

    metrics = [
        "HitScore",
        "LFC_R1_vs_DEL2",
        "LFC_R2_vs_DEL2",
        "LFC_NEG_centered",
        "LFC_NEG_centered_R1",
        "LFC_NEG_centered_R2",
        "GLM_hit",
        "RS_pass",
        "Consensus_hit",
        "NEG_hard_fail",
    ]
    for m in metrics:
        if m in cols:
            display_cols.append(m)

    table = df[display_cols].copy()
    table.rename(columns={lib_col: "LibID", id_col: "ID"}, inplace=True)
    if bb_cols:
        table.rename(columns={c: c.replace("_x", "") for c in bb_cols}, inplace=True)
    if del2_col:
        table.rename(columns={del2_col: "DEL2_raw"}, inplace=True)

    table.insert(0, "Rank", range(1, len(table) + 1))
    # Clean LIB suffix in BB fields and build display ID with LibID prefix
    raw_id = table["ID"].copy() if "ID" in table.columns else pd.Series([""] * len(table))
    for c in ["BB1", "BB2", "BB3", "BB4"]:
        if c in table.columns:
            table[c] = table[c].apply(_strip_libdel)
            table[c] = table[c].replace("", "NA")
    if "ID" in table.columns and "LibID" in table.columns:
        table["ID"] = table.apply(
            lambda r: _make_display_id(
                r.get("LibID", ""),
                r.get("BB1", "NA"),
                r.get("BB2", "NA"),
                r.get("BB3", "NA"),
                r.get("BB4", "NA"),
                raw_id.loc[r.name],
            ),
            axis=1,
        )
    table["QC_Flags"] = df.apply(_qc_flags, axis=1, del2_col=del2_col)

    return table, del2_col or ""


def _summary_counts(df: pd.DataFrame) -> dict:
    out = {"rows": int(len(df))}
    for c in ["GLM_hit", "RS_pass", "NEG_hard_fail", "Consensus_hit"]:
        if c in df.columns:
            out[c] = int(_normalize_bool(df[c]).sum())
    return out


def _html_table(df: pd.DataFrame) -> str:
    return df.to_html(index=False, escape=True, border=0, classes="data-table")


def main() -> None:
    ap = argparse.ArgumentParser(description="Beginner-friendly QC report for DELeGANce outputs")
    ap.add_argument("--run_root", action="append", default=[], help="Run root (e.g., DELeGANce_out/my_run)")
    ap.add_argument("--annot_tsv", action="append", default=[], help="Explicit 05_hybrid_annot.tsv path (repeatable)")
    ap.add_argument("--prefer_dir", default="03_normalized/glm_full_dev_cpu_fp64",
                    help="Preferred subdir under run_root for 05_hybrid_annot.tsv")
    ap.add_argument("--out_html", default="DELeGANce_out/Beginner_QC_Report.html")
    ap.add_argument("--out_tsv", default="DELeGANce_out/Beginner_QC_TopHits.tsv")
    ap.add_argument("--top_n", type=int, default=200, help="Top N hits per run (HitScore ranking)")
    args = ap.parse_args()

    annot_paths = [Path(p) for p in args.annot_tsv]
    run_roots = [Path(p) for p in args.run_root]
    if not annot_paths and not run_roots:
        raise SystemExit("[ERROR] --run_root or --annot_tsv is required.")
    out_html = Path(args.out_html)
    out_tsv = Path(args.out_tsv)

    summaries = []
    top_tables = []
    per_run_html = []

    for annot_path in annot_paths:
        if not annot_path.exists():
            raise FileNotFoundError(f"[ERROR] annot_tsv not found: {annot_path}")
        df = pd.read_csv(annot_path, sep="\t", low_memory=False)
        run_name = _guess_run_name(annot_path)
        summary = _summary_counts(df)
        summary["run"] = run_name
        summaries.append(summary)

        top_table, del2_col = _build_top_table(df, args.top_n)
        top_table.insert(1, "Run", run_name)
        top_tables.append(top_table)

        per_run_html.append(
            f"<h2>{html.escape(run_name)}</h2>"
            f"<p>Top {args.top_n} by HitScore (HitScore_GLM/HitScore_RS fallback).</p>"
            f"{_html_table(top_table)}"
        )

    for run_root in run_roots:
        annot_path = _resolve_annot(run_root, args.prefer_dir)
        df = pd.read_csv(annot_path, sep="\t", low_memory=False)

        run_name = run_root.name
        summary = _summary_counts(df)
        summary["run"] = run_name
        summaries.append(summary)

        top_table, del2_col = _build_top_table(df, args.top_n)
        top_table.insert(1, "Run", run_name)
        top_tables.append(top_table)

        per_run_html.append(
            f"<h2>{html.escape(run_name)}</h2>"
            f"<p>Top {args.top_n} by HitScore (HitScore_GLM/HitScore_RS fallback).</p>"
            f"{_html_table(top_table)}"
        )

    # Summary table
    summary_df = pd.DataFrame(summaries)[["run", "rows", "GLM_hit", "RS_pass", "NEG_hard_fail", "Consensus_hit"]]
    summary_html = _html_table(summary_df.fillna(""))

    # TSV output
    all_top = pd.concat(top_tables, ignore_index=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    all_top.to_csv(out_tsv, sep="\t", index=False, na_rep="NA")

    # HTML output
    out_html.parent.mkdir(parents=True, exist_ok=True)
    html_text = f"""<!doctype html>
<html lang="ko">
<head>
  <meta charset="utf-8" />
  <title>DELeGANce QC 리포트 (초보자용)</title>
  <style>
    body {{ font-family: Arial, sans-serif; line-height: 1.5; padding: 20px; color: #222; }}
    h1, h2 {{ color: #1f4e79; }}
    .note {{ background: #f7f9fc; border: 1px solid #e2e6ee; padding: 12px; }}
    .data-table {{ border-collapse: collapse; width: 100%; margin: 10px 0 30px 0; }}
    .data-table th, .data-table td {{ border: 1px solid #dfe3ea; padding: 6px 8px; font-size: 12px; }}
    .data-table th {{ background: #f1f4f9; text-align: left; }}
    .small {{ font-size: 12px; color: #555; }}
  </style>
</head>
<body>
  <h1>DELeGANce QC 리포트 (초보자용)</h1>
  <div class="note">
    <p><b>이 리포트의 목적</b>: 최종 hit 결과를 초보자도 빠르게 이해하고, 위험 신호(NEG 높음, R2 급감 등)를 한눈에 확인하도록 돕습니다.</p>
    <p class="small">주의: QC 플래그는 “의심 지점” 표시입니다. 최종 판단은 실험 맥락과 함께 해석해 주세요.</p>
  </div>

  <h2>핵심 지표 설명 (간단 버전)</h2>
  <ul class="small">
    <li><b>HitScore</b>: 여러 지표를 종합한 점수 (높을수록 우수)</li>
    <li><b>LFC_R1_vs_DEL2</b>: 1차 결합(R1)이 baseline(DEL2) 대비 얼마나 증가했는지 (log2)</li>
    <li><b>LFC_R2_vs_DEL2</b>: 재결합/정제(R2)에서의 증가 정도 (log2)</li>
    <li><b>LFC_NEG_centered</b>: NEG 대비 특이성 지표 (0 이상이면 NEG가 높을 가능성)</li>
    <li><b>GLM_hit / RS_pass / Consensus_hit</b>: 통계/규칙 기반 필터 통과 여부</li>
    <li><b>NEG_hard_fail</b>: NEG가 너무 높아 탈락한 항목</li>
  </ul>

  <h2>전체 요약</h2>
  {summary_html}

  <h2>QC 플래그 기준</h2>
  <ul class="small">
    <li><b>LOW_DEL2</b>: DEL2 raw count &lt; 10 (기저 카운트가 낮아 변동성 큼)</li>
    <li><b>NEG_HIGH</b>: LFC_NEG_centered ≥ 0 (NEG가 충분히 낮지 않음)</li>
    <li><b>R2_DROP</b>: LFC_R2_vs_DEL2가 R1보다 1.0 이상 낮음 (재결합 단계에서 급감)</li>
    <li><b>BB_MISSING</b>: BB1~BB3 중 누락</li>
    <li><b>NEG_HARD_FAIL</b>: NEG_hard_fail = TRUE</li>
  </ul>

  {''.join(per_run_html)}
</body>
</html>
"""
    out_html.write_text(html_text, encoding="utf-8")
    print(f"[OK] HTML: {out_html}")
    print(f"[OK] TSV: {out_tsv}")


if __name__ == "__main__":
    main()
