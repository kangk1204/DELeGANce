#!/usr/bin/env python3
import argparse
import html
import json
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
    if series is None:
        return pd.Series([], dtype=bool)
    if series.dtype == bool:
        return series.fillna(False)
    if series.dtype == object:
        s = series.astype(str).str.strip().str.lower()
        return s.isin(["true", "1", "yes", "y", "t"])
    s = pd.to_numeric(series, errors="coerce").fillna(0.0)
    return s > 0

def _truthy_value(value) -> bool:
    if value is None:
        return False
    try:
        if isinstance(value, float) and np.isnan(value):
            return False
    except Exception:
        pass
    if isinstance(value, bool):
        return value
    s = str(value).strip().lower()
    return s in ["true", "1", "yes", "y", "t"]


# Derived DEL2 columns produced by 03_call_hits.py (never raw counts)
_DEL2_DERIVED = {"DEL2_norm", "DEL2_sum"}


def _load_del2_from_params(annot_path: Path) -> Optional[str]:
    """Read normalized_columns.del2 from hit_params.json next to the annot file (written by the orchestrator)."""
    hp = annot_path.parent / "hit_params.json"
    if not hp.exists():
        return None
    try:
        params = json.loads(hp.read_text(encoding="utf-8"))
    except Exception:
        return None
    del2 = (params.get("normalized_columns") or {}).get("del2")
    return str(del2) if del2 else None


def _auto_del2_col(cols: List[str]) -> Optional[str]:
    if "DEL2_OVERRIDE" in cols:
        return "DEL2_OVERRIDE"
    def _is_derived(c: str) -> bool:
        return c.endswith("_CPM") or c.endswith("_norm") or c.endswith("_sum") or c in _DEL2_DERIVED

    # 1) raw sample column named DEL* (skip CPM and *_norm/*_sum derived columns from 03_call_hits)
    for c in cols:
        if c.startswith("DEL") and not _is_derived(c):
            return c
    # 2) raw sample column (has a <c>_CPM counterpart) whose name contains "DEL" (e.g. K_DEL234)
    colset = set(cols)
    for c in cols:
        if "DEL" in c.upper() and not _is_derived(c) and f"{c}_CPM" in colset:
            return c
    return None


# Anchored: strips only a trailing "_LIB<lib>" namespace token from a single BB value.
_LIBDEL_RE = re.compile(r"_LIB[\w.-]+$")
# Fallback for full tag IDs (cycles_BB1_BB2_BB3[_BB4]) when no lib_id is known: removes each "_LIB<lib>"
# token without consuming the following "_BB" tokens. It is ambiguous for lib_ids that contain "_"
# (01_preprocess allows [A-Za-z0-9_.-]), so callers pass the run's known lib_ids whenever available.
_LIBDEL_TOKEN_RE = re.compile(r"_LIB[^_]+(?=_|$)")
_NA_LIB_VALUES = {"", "na", "nan", "none", "<na>"}


def _known_libs(series) -> Tuple[str, ...]:
    """Distinct real lib_ids from a LIB_ID column (NA-like values dropped), longest first."""
    if series is None:
        return ()
    vals = {str(v).strip() for v in series.dropna().unique()}
    vals = {v for v in vals if v.lower() not in _NA_LIB_VALUES}
    return tuple(sorted(vals, key=len, reverse=True))


def _strip_libdel(value) -> str:
    if pd.isna(value):
        return ""
    return _LIBDEL_RE.sub("", str(value))


def _strip_libdel_anywhere(value, known_libs: Tuple[str, ...] = ()) -> str:
    """Remove "_LIB<lib>" tokens from a full tag ID. With known lib_ids only the exact tokens
    "_LIB<lib_id>" (followed by "_" or end) are removed; the generic regex is used only when no
    lib_id is known for the run."""
    if pd.isna(value):
        return ""
    s = str(value)
    if known_libs:
        for lib in known_libs:
            s = re.sub(r"_LIB" + re.escape(lib) + r"(?=_|$)", "", s)
        return s
    return _LIBDEL_TOKEN_RE.sub("", s)


def _pick_bb_value(row: pd.Series, base: str) -> str:
    for c in (base, f"{base}_x", f"{base}_y"):
        if c in row.index:
            return _strip_libdel(row.get(c, ""))
    return ""


def _make_display_id(lib: str, bb1: str, bb2: str, bb3: str, bb4: str, fallback_id: str,
                     known_libs: Tuple[str, ...] = ()) -> str:
    # NaN (pandas reads a literal "NA" lib_id as NaN) is truthy, so test with isna before str()
    lib = "" if (lib is None or (isinstance(lib, float) and np.isnan(lib))) else str(lib).strip()
    if lib and lib.lower() not in _NA_LIB_VALUES:
        return f"{lib}_{bb1}_{bb2}_{bb3}_{bb4}"
    return _strip_libdel_anywhere(fallback_id, known_libs)


def _neg_quantile_threshold(series: pd.Series, quantile: float) -> Optional[float]:
    if series is None:
        return None
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return None
    q = float(s.quantile(quantile))
    if not np.isfinite(q):
        return None
    return max(0.0, q)


def _select_diverse_indices(df: pd.DataFrame, key_cols: List[str], max_n: int) -> List[int]:
    seen = set()
    picked = []
    for idx, row in df.iterrows():
        key = "|".join(str(row.get(c, "NA")) for c in key_cols)
        if key in seen:
            continue
        seen.add(key)
        picked.append(idx)
        if max_n is not None and len(picked) >= max_n:
            break
    return picked


def _qc_flags(row: pd.Series, del2_col: Optional[str],
              neg_high_thr: Optional[float],
              neg_r1_high_thr: Optional[float]) -> str:
    flags = []

    if del2_col:
        del2_val = pd.to_numeric(pd.Series([row.get(del2_col)]), errors="coerce").iloc[0]
        if pd.notna(del2_val) and float(del2_val) < 10:
            flags.append("LOW_DEL2")

    lfc_neg = row.get("LFC_NEG_centered")
    neg_thr = 0.0 if neg_high_thr is None else float(neg_high_thr)
    if pd.notna(lfc_neg) and float(lfc_neg) >= neg_thr:
        flags.append("NEG_HIGH")

    lfc_neg_r1 = row.get("LFC_NEG_centered_R1")
    neg_r1_thr = 0.0 if neg_r1_high_thr is None else float(neg_r1_high_thr)
    if pd.notna(lfc_neg_r1) and float(lfc_neg_r1) >= neg_r1_thr:
        flags.append("NEG_R1_HIGH")

    lfc_r1 = row.get("LFC_R1_vs_DEL2_used", row.get("LFC_R1_vs_DEL2"))
    lfc_r2 = row.get("LFC_R2_vs_DEL2_used", row.get("LFC_R2_vs_DEL2"))
    if pd.notna(lfc_r1) and pd.notna(lfc_r2):
        if float(lfc_r2) < float(lfc_r1) - 1.0:
            flags.append("R2_DROP")

    for base in ["BB1", "BB2", "BB3"]:
        val = _pick_bb_value(row, base).strip().upper()
        if val in ["", "NA", "NAN", "NONE"]:
            flags.append("BB_MISSING")
            break

    neg_fail = row.get("NEG_hard_fail")
    if _truthy_value(neg_fail):
        flags.append("NEG_HARD_FAIL")

    return "OK" if not flags else ";".join(flags)


def _build_top_table(df: pd.DataFrame, top_n: int, neg_high_quantile: float,
                     recommend_a: int, recommend_b: int, recommend_diverse: int,
                     diverse_key: str) -> Tuple[pd.DataFrame, str, dict, dict]:
    df_all = _add_hit_score(df.copy())
    df_all = df_all[np.isfinite(df_all["HitScore"])].copy()
    # Coerce the metric columns used by _qc_flags once (non-numeric strings become NaN instead of raising)
    for c in ("LFC_NEG_centered", "LFC_NEG_centered_R1", "LFC_R1_vs_DEL2_used", "LFC_R2_vs_DEL2_used",
              "LFC_R1_vs_DEL2", "LFC_R2_vs_DEL2"):
        if c in df_all.columns:
            df_all[c] = pd.to_numeric(df_all[c], errors="coerce")
    lib_col_all = _pick_col(list(df_all.columns), ["LIB_ID", "LIB_ID_x", "LIB_ID_y"])
    known_libs = _known_libs(df_all[lib_col_all]) if lib_col_all else ()
    neg_thr = _neg_quantile_threshold(df_all.get("LFC_NEG_centered"), neg_high_quantile)
    neg_r1_thr = _neg_quantile_threshold(df_all.get("LFC_NEG_centered_R1"), neg_high_quantile)

    # Stable sort with a deterministic tie-breaker so Rank/PickGroup are reproducible across runs
    tie_col = _pick_col(list(df_all.columns), ["ID", "id", "ID_x", "id_x", "ID_y", "id_y"])
    sort_keys = ["HitScore"] + ([tie_col] if tie_col else [])
    df = df_all.sort_values(sort_keys, ascending=[False] + [True] * (len(sort_keys) - 1),
                            kind="mergesort").head(top_n).copy()
    df = df.reset_index(drop=True)

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
                known_libs,
            ),
            axis=1,
        )
    table["QC_Flags"] = df.apply(
        _qc_flags,
        axis=1,
        del2_col=del2_col,
        neg_high_thr=neg_thr,
        neg_r1_high_thr=neg_r1_thr,
    )

    # Tiering (A: consensus + no QC flags; B: consensus + has QC flags)
    consensus = _normalize_bool(df.get("Consensus_hit", pd.Series([False] * len(df))))
    neg_hard = _normalize_bool(df.get("NEG_hard_fail", pd.Series([False] * len(df))))
    no_flags = table["QC_Flags"].astype(str).str.upper().eq("OK")
    tier = np.where(consensus & (~neg_hard) & no_flags, "A",
                    np.where(consensus & (~neg_hard) & (~no_flags), "B", "Other"))
    table["Tier"] = tier

    # Pick groups
    pick_group = [""] * len(table)
    tier_a_idx = table.index[table["Tier"] == "A"].tolist()
    tier_b_idx = table.index[table["Tier"] == "B"].tolist()

    pick_a = tier_a_idx[: max(0, int(recommend_a))]
    pick_b = tier_b_idx[: max(0, int(recommend_b))]

    key_cols = ["BB1", "BB2", "BB3"] if diverse_key == "BB1_BB2_BB3" else \
               ["BB1", "BB2"] if diverse_key == "BB1_BB2" else \
               ["BB1", "BB2", "BB3", "BB4"] if diverse_key == "BB1_BB2_BB3_BB4" else \
               ["BB1"]
    diverse_candidates = table[table["Tier"] == "A"].copy()
    diverse_idx = _select_diverse_indices(diverse_candidates, key_cols, max(0, int(recommend_diverse)))

    def _add_tag(idx_list, tag):
        for i in idx_list:
            if pick_group[i]:
                pick_group[i] += ";" + tag
            else:
                pick_group[i] = tag

    _add_tag(pick_a, "TierA_Top")
    _add_tag(pick_b, "TierB_Control")
    _add_tag(diverse_idx, "TierA_Diverse")
    table["PickGroup"] = [v if v else "NA" for v in pick_group]

    # Reorder so Tier/PickGroup are near front
    front = ["Rank", "Tier", "PickGroup"]
    cols_final = front + [c for c in table.columns if c not in front]
    table = table[cols_final]

    thresholds = {
        "neg_high_thr": neg_thr,
        "neg_r1_high_thr": neg_r1_thr,
        "neg_high_quantile": neg_high_quantile,
    }
    rec_tables = {
        "tier_a_top": table.loc[pick_a].copy(),
        "tier_b_control": table.loc[pick_b].copy(),
        "tier_a_diverse": table.loc[diverse_idx].copy(),
    }
    return table, del2_col or "", thresholds, rec_tables


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
    ap.add_argument("--del2_col", default="", help="DEL2 column name override (optional)")
    ap.add_argument("--out_html", default="DELeGANce_out/Beginner_QC_Report.html")
    ap.add_argument("--out_tsv", default="DELeGANce_out/Beginner_QC_TopHits.tsv")
    ap.add_argument("--top_n", type=int, default=200, help="Top N hits per run (HitScore ranking)")
    ap.add_argument("--neg_high_quantile", type=float, default=0.90,
                    help="Quantile for NEG_HIGH/NEG_R1_HIGH flags (per-run, fallback to 0 if lower)")
    ap.add_argument("--recommend_a", type=int, default=50, help="Tier A top picks per run")
    ap.add_argument("--recommend_b", type=int, default=20, help="Tier B control picks per run")
    ap.add_argument("--recommend_diverse", type=int, default=50, help="Tier A diversity picks per run")
    ap.add_argument("--diverse_key", choices=["BB1", "BB1_BB2", "BB1_BB2_BB3", "BB1_BB2_BB3_BB4"],
                    default="BB1_BB2_BB3", help="Key for diversity picks")
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

    # Resolve (annot_path, run_name) jobs once; the same processing applies to both input styles
    annot_jobs: List[Tuple[Path, str]] = []
    for annot_path in annot_paths:
        if not annot_path.exists():
            raise FileNotFoundError(f"[ERROR] annot_tsv not found: {annot_path}")
        annot_jobs.append((annot_path, _guess_run_name(annot_path)))
    for run_root in run_roots:
        annot_jobs.append((_resolve_annot(run_root, args.prefer_dir), run_root.name))

    for annot_path, run_name in annot_jobs:
        df = pd.read_csv(annot_path, sep="\t", low_memory=False)
        summary = _summary_counts(df)
        summary["run"] = run_name
        summaries.append(summary)

        # DEL2 column: explicit --del2_col > hit_params.json (orchestrator/03) > name heuristic.
        # An explicit name that is not a column falls through to the next source (with a warning)
        # instead of silently skipping hit_params.json.
        del2_from_params = _load_del2_from_params(annot_path)
        del2_override = ""
        for cand in (args.del2_col, del2_from_params):
            if cand and cand in df.columns:
                del2_override = cand
                break
        if args.del2_col and args.del2_col not in df.columns:
            nxt = f"hit_params.json ({del2_from_params})" if del2_override else "name heuristic"
            print(f"[WARN] --del2_col {args.del2_col!r} is not a column of {annot_path}; using {nxt}")
        if del2_override:
            df = df.rename(columns={del2_override: "DEL2_OVERRIDE"})
        top_table, del2_col, th, rec = _build_top_table(
            df, args.top_n, args.neg_high_quantile,
            args.recommend_a, args.recommend_b, args.recommend_diverse, args.diverse_key
        )
        top_table.insert(1, "Run", run_name)
        top_tables.append(top_table)

        def _fmt_thr(v) -> str:
            return f"≥ {v:.3f}" if v is not None else "n/a (column missing)"

        thr_note = (
            f"<p class='small'>QC thresholds: NEG_HIGH {_fmt_thr(th['neg_high_thr'])}, "
            f"NEG_R1_HIGH {_fmt_thr(th['neg_r1_high_thr'])} "
            f"(quantile={th['neg_high_quantile']:.2f}, min=0)</p>"
        )
        rec_html = (
            f"<h3>Recommended Tier A (Top {args.recommend_a})</h3>"
            f"{_html_table(rec['tier_a_top'])}"
            f"<h3>Recommended Tier A Diversity (Top {args.recommend_diverse}, key={args.diverse_key})</h3>"
            f"{_html_table(rec['tier_a_diverse'])}"
            f"<h3>Recommended Tier B Controls (Top {args.recommend_b})</h3>"
            f"{_html_table(rec['tier_b_control'])}"
        )
        per_run_html.append(
            f"<h2>{html.escape(run_name)}</h2>"
            f"<p>Top {args.top_n} by HitScore (HitScore_GLM/HitScore_RS fallback).</p>"
            f"{thr_note}"
            f"{_html_table(top_table)}"
            f"{rec_html}"
        )

    # Summary table
    # reindex (not []) so runs lacking some flag columns do not raise KeyError
    summary_df = pd.DataFrame(summaries).reindex(columns=["run", "rows", "GLM_hit", "RS_pass", "NEG_hard_fail", "Consensus_hit"])
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
    <li><b>Tier</b>: A=합성 우선(Consensus_hit + QC 플래그 없음), B=컨트롤 후보(Consensus_hit + QC 플래그 존재, NEG_hard_fail 제외 → Other)</li>
    <li><b>PickGroup</b>: TierA_Top, TierA_Diverse, TierB_Control (아래 추천 리스트와 동일)</li>
  </ul>

  <h2>전체 요약</h2>
  {summary_html}

  <h2>QC 플래그 기준</h2>
  <ul class="small">
    <li><b>LOW_DEL2</b>: DEL2 raw count &lt; 10 (기저 카운트가 낮아 변동성 큼)</li>
    <li><b>NEG_HIGH</b>: LFC_NEG_centered ≥ Q{int(args.neg_high_quantile * 100)} (per-run, 최소 0)</li>
    <li><b>NEG_R1_HIGH</b>: LFC_NEG_centered_R1 ≥ Q{int(args.neg_high_quantile * 100)} (per-run, 최소 0)</li>
    <li><b>R2_DROP</b>: LFC_R2_vs_DEL2가 R1보다 1.0 이상 낮음 (재결합 단계에서 급감)</li>
    <li><b>BB_MISSING</b>: BB1~BB3 중 누락</li>
    <li><b>NEG_HARD_FAIL</b>: NEG_hard_fail = TRUE</li>
  </ul>

  <h2>추천 선택 가이드 (초보자용)</h2>
  <ul class="small">
    <li><b>Tier A Top</b>: 합성 우선 후보 (Consensus_hit + QC 플래그 없음)</li>
    <li><b>Tier A Diversity</b>: BB 조합을 다양화한 추천 (key={args.diverse_key})</li>
    <li><b>Tier B Control</b>: 리스크/컨트롤 후보 (Consensus_hit이지만 QC 플래그 존재)</li>
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
