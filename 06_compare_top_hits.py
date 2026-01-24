#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare top hits across multiple runs and generate an interactive summary.
"""

import argparse
import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

_HAS_BOKEH = True
try:
    from bokeh.io import output_file, save
    from bokeh.layouts import column
    from bokeh.models import ColumnDataSource, DataTable, Div, HoverTool, NumberFormatter, Panel, TableColumn, Tabs, TabPanel
    from bokeh.plotting import figure
except Exception:
    _HAS_BOKEH = False

_HAS_TABPANEL = False
try:
    from bokeh.models import TabPanel as _BokehTabPanel  # type: ignore
    _HAS_TABPANEL = True
except Exception:
    _HAS_TABPANEL = False


def _strip_lib_suffix(val: str) -> str:
    if val is None:
        return "NA"
    s = str(val).strip()
    if s in ("", "NA", "nan", "None"):
        return "NA"
    return re.sub(r"_LIB[\w\.-]+$", "", s)


def _pick_latest(paths: List[str]) -> Optional[str]:
    if not paths:
        return None
    return max(paths, key=lambda p: os.path.getmtime(p))


def resolve_hybrid_path(base: str, preset: Optional[str]) -> str:
    if os.path.isfile(base):
        return base
    if os.path.isdir(base):
        direct = os.path.join(base, "05_hybrid_annot.tsv")
        if os.path.isfile(direct):
            return direct
        norm = os.path.join(base, "03_normalized")
        if os.path.isdir(norm):
            if preset:
                cand = os.path.join(norm, preset, "05_hybrid_annot.tsv")
                if os.path.isfile(cand):
                    return cand
            paths = []
            for root, _, files in os.walk(norm):
                if "05_hybrid_annot.tsv" in files:
                    paths.append(os.path.join(root, "05_hybrid_annot.tsv"))
            cand = _pick_latest(paths)
            if cand:
                return cand
    raise FileNotFoundError(f"05_hybrid_annot.tsv not found under: {base}")


def infer_score_col(paths: List[str], preset: Optional[str], override: Optional[str]) -> str:
    if override:
        return override
    has_glm = True
    has_rs = True
    for p in paths:
        hp = resolve_hybrid_path(p, preset)
        cols = pd.read_csv(hp, sep="\t", nrows=0).columns.tolist()
        has_glm = has_glm and ("HitScore_GLM" in cols)
        has_rs = has_rs and ("HitScore_RS" in cols)
    if has_glm:
        return "HitScore_GLM"
    if has_rs:
        return "HitScore_RS"
    raise SystemExit("[ERROR] No common score column across runs (need HitScore_GLM or HitScore_RS).")


def _to_int_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").fillna(0).astype(int)


def _make_compound_key(df: pd.DataFrame) -> pd.Series:
    bb_cols = [c for c in ["BB1_x", "BB2_x", "BB3_x", "BB4_x"] if c in df.columns]
    if not bb_cols:
        return df.get("ID_x", pd.Series(["NA"] * len(df))).astype(str)
    for c in bb_cols:
        df[c] = df[c].fillna("NA").astype(str).map(_strip_lib_suffix)
    return df[bb_cols].agg("|".join, axis=1)


def _apply_recommend_filter(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "Consensus_hit" in out.columns:
        out = out[_to_int_series(out["Consensus_hit"]) == 1]
    if "GLM_hit" in out.columns:
        out = out[_to_int_series(out["GLM_hit"]) == 1]
    if "RS_pass" in out.columns:
        out = out[_to_int_series(out["RS_pass"]) == 1]
    if "NEG_hard_fail" in out.columns:
        out = out[_to_int_series(out["NEG_hard_fail"]) == 0]
    return out


def _sanitize_label(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_") or "run"


def _table_columns(df: pd.DataFrame, cols: List[str]) -> List[TableColumn]:
    out = []
    for c in cols:
        if c not in df.columns:
            continue
        if df[c].dtype.kind in "if":
            out.append(TableColumn(field=c, title=c, formatter=NumberFormatter(format="0.0000")))
        else:
            out.append(TableColumn(field=c, title=c))
    return out


def _make_panel(child, title: str):
    if _HAS_TABPANEL:
        return _BokehTabPanel(child=child, title=title)
    return Panel(child=child, title=title)


def build_interactive_html(summary_html: str, pattern_df: pd.DataFrame,
                           overlap_df: pd.DataFrame, run_tabs: List[Panel],
                           out_html: str) -> None:
    if not _HAS_BOKEH:
        print("[WARN] bokeh is not available; skipping interactive HTML output.")
        return

    items = [Div(text=summary_html)]

    if not pattern_df.empty:
        src = ColumnDataSource(pattern_df)
        p = figure(
            title="Top-N overlap patterns",
            x_range=pattern_df["pattern"].astype(str).tolist(),
            height=320,
            tools="pan,reset,save",
        )
        p.vbar(x="pattern", top="count", width=0.9, source=src, color="#74a9cf")
        p.xaxis.major_label_orientation = 1.0
        p.add_tools(HoverTool(tooltips=[("pattern", "@pattern"), ("count", "@count"), ("pct", "@pct{0.0}%")]))
        items.append(p)

    if not overlap_df.empty:
        cols = [c for c in overlap_df.columns if c not in ("compound_key",)]
        table_cols = _table_columns(overlap_df, ["compound_key"] + cols)
        table = DataTable(source=ColumnDataSource(overlap_df), columns=table_cols, height=360, width=1200, index_position=None)
        items.append(Div(text="<h3>Overlapping compounds (>=2 runs)</h3>"))
        items.append(table)

    comp_panel = _make_panel(column(*items, sizing_mode="stretch_width"), "Comparison")
    tabs = Tabs(tabs=[comp_panel] + run_tabs)
    output_file(out_html, title="Top-hit comparison")
    save(tabs)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True, help="Run roots or 05_hybrid_annot.tsv paths")
    ap.add_argument("--labels", nargs="*", default=None, help="Optional labels for runs (same order)")
    ap.add_argument("--preset", default=None, help="Preset under 03_normalized (optional)")
    ap.add_argument("--top-n", type=int, default=100, help="Top N per run (default)")
    ap.add_argument("--recommend-n", type=int, default=10, help="Recommended hits per run (default)")
    ap.add_argument("--top-n-list", default=None,
                    help="Comma/space-separated list of top-N per run (e.g., 100,100,200)")
    ap.add_argument("--recommend-n-list", default=None,
                    help="Comma/space-separated list of recommended-N per run (e.g., 10,10,20)")
    ap.add_argument("--score-col", default=None, help="Override score column")
    ap.add_argument("--out-dir", default=".", help="Output directory (default: current dir)")
    ap.add_argument("--out-prefix", default=None, help="Output prefix (default: compare_topN)")
    ap.add_argument("--min-overlap", type=int, default=2, help="Min runs to consider overlap")
    ap.add_argument("--no-html", action="store_true", help="Skip interactive HTML output")
    args = ap.parse_args()

    run_paths = args.runs
    labels = args.labels or []
    if labels and len(labels) != len(run_paths):
        raise SystemExit("[ERROR] --labels must match number of --runs")
    if not labels:
        labels = [os.path.basename(os.path.normpath(p)) or f"run{i+1}" for i, p in enumerate(run_paths)]
    labels = [_sanitize_label(l) for l in labels]

    score_col = infer_score_col(run_paths, args.preset, args.score_col)
    def _parse_int_list(s: Optional[str]) -> Optional[List[int]]:
        if not s:
            return None
        parts = [p for p in re.split(r"[,\s]+", str(s).strip()) if p]
        return [max(1, int(p)) for p in parts]

    top_n = max(1, int(args.top_n))
    rec_n = max(1, int(args.recommend_n))
    top_n_list = _parse_int_list(args.top_n_list)
    rec_n_list = _parse_int_list(args.recommend_n_list)
    if top_n_list and len(top_n_list) != len(run_paths):
        raise SystemExit("[ERROR] --top-n-list must match number of runs")
    if rec_n_list and len(rec_n_list) != len(run_paths):
        raise SystemExit("[ERROR] --recommend-n-list must match number of runs")
    min_overlap = max(2, int(args.min_overlap))

    out_dir = args.out_dir or "."
    os.makedirs(out_dir, exist_ok=True)
    name_top = max(top_n_list) if top_n_list else top_n
    prefix = args.out_prefix or os.path.join(out_dir, f"compare_top{name_top}")
    html_path = f"{prefix}_interactive.html"

    run_info = []
    for idx, (path, label) in enumerate(zip(run_paths, labels)):
        hybrid = resolve_hybrid_path(path, args.preset)
        header = pd.read_csv(hybrid, sep="\t", nrows=0).columns.tolist()
        want = [
            "LIB_ID_x", "ID_x", "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x",
            "HitScore_GLM", "HitScore_RS",
            "GLM_hit", "RS_pass", "Consensus_hit", "NEG_hard_fail",
        ]
        usecols = [c for c in want if c in header] + [score_col]
        usecols = list(dict.fromkeys(usecols))
        df = pd.read_csv(hybrid, sep="\t", usecols=usecols)
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(float("-inf"))
        df["compound_key"] = _make_compound_key(df)
        run_top_n = top_n_list[idx] if top_n_list else top_n
        run_rec_n = rec_n_list[idx] if rec_n_list else rec_n
        df_top = df.sort_values(score_col, ascending=False).head(run_top_n).copy()
        df_top["compound_key"] = _make_compound_key(df_top)
        df_rec = _apply_recommend_filter(df_top)
        df_rec = df_rec.sort_values(score_col, ascending=False)
        df_rec = df_rec.drop_duplicates("compound_key").head(run_rec_n).copy()
        run_info.append({
            "label": label,
            "hybrid": hybrid,
            "top": df_top,
            "rec": df_rec,
            "top_n": run_top_n,
            "rec_n": run_rec_n,
        })

    # Save per-run tables
    for info in run_info:
        label = info["label"]
        top_path = f"{prefix}_{label}_top{info['top_n']}.tsv"
        rec_path = f"{prefix}_{label}_recommended{info['rec_n']}.tsv"
        info["top"].to_csv(top_path, sep="\t", index=False)
        info["rec"].to_csv(rec_path, sep="\t", index=False)
        info["top_path"] = top_path
        info["rec_path"] = rec_path

    # Build comparison tables
    top_sets = {r["label"]: set(r["top"]["compound_key"].astype(str)) for r in run_info}
    union_keys = sorted(set().union(*top_sets.values()))

    pattern_rows = []
    overlap_rows = []
    union_rows = []
    for key in union_keys:
        present = [lbl for lbl in labels if key in top_sets.get(lbl, set())]
        pattern = "+".join(present)
        pattern_rows.append(pattern)
        row = {"compound_key": key, "pattern": pattern, "n_runs": len(present)}
        # Attach per-run scores
        for info in run_info:
            lbl = info["label"]
            scores = info["top"].loc[info["top"]["compound_key"] == key, score_col]
            row[f"{lbl}_score"] = float(scores.max()) if not scores.empty else np.nan
        union_rows.append(row)
        if len(present) >= min_overlap:
            overlap_rows.append(row.copy())

    pattern_df = pd.DataFrame(pattern_rows, columns=["pattern"])
    if not pattern_df.empty:
        pattern_df = (
            pattern_df.groupby("pattern", as_index=False)
            .size().rename(columns={"size": "count"})
            .sort_values("count", ascending=False)
        )
        pattern_df["pct"] = (pattern_df["count"] / pattern_df["count"].sum()) * 100.0

    union_df = pd.DataFrame(union_rows)
    overlap_df = pd.DataFrame(overlap_rows)

    union_path = f"{prefix}_union_top{name_top}.tsv"
    overlap_path = f"{prefix}_overlap_top{name_top}.tsv"
    pattern_path = f"{prefix}_pattern_counts.tsv"
    union_df.to_csv(union_path, sep="\t", index=False)
    overlap_df.to_csv(overlap_path, sep="\t", index=False)
    pattern_df.to_csv(pattern_path, sep="\t", index=False)

    # Build HTML
    if not args.no_html:
        summary_lines = ["<h2>Top-hit comparison</h2>", "<ul>"]
        summary_lines.append(f"<li>Score column: <b>{score_col}</b></li>")
        summary_lines.append(f"<li>Top-N default: <b>{top_n}</b></li>")
        summary_lines.append(f"<li>Min overlap: <b>{min_overlap}</b> runs</li>")
        for info in run_info:
            summary_lines.append(
                f"<li>{info['label']}: top={len(info['top'])} (target {info['top_n']}), "
                f"recommended={len(info['rec'])} (target {info['rec_n']})</li>"
            )
        summary_lines.append(f"<li>Union size: <b>{len(union_df)}</b></li>")
        summary_lines.append(f"<li>Overlap (>= {min_overlap}): <b>{len(overlap_df)}</b></li>")
        summary_lines.append("</ul>")

        run_tabs = []
        for info in run_info:
            top_df = info["top"].copy()
            rec_df = info["rec"].copy()
            cols_top = ["compound_key", score_col, "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x"]
            cols_rec = cols_top
            cols_top = [c for c in cols_top if c in top_df.columns]
            cols_rec = [c for c in cols_rec if c in rec_df.columns]

            top_table = DataTable(
                source=ColumnDataSource(top_df[cols_top]),
                columns=_table_columns(top_df, cols_top),
                height=320,
                width=1200,
                index_position=None,
            )
            rec_table = DataTable(
                source=ColumnDataSource(rec_df[cols_rec]),
                columns=_table_columns(rec_df, cols_rec),
                height=260,
                width=1200,
                index_position=None,
            )
            panel = _make_panel(
                column(
                    Div(text=f"<h3>{info['label']}</h3>"),
                    Div(text=f"<p>Top-N table (n={len(top_df)})</p>"),
                    top_table,
                    Div(text=f"<p>Recommended (n={len(rec_df)})</p>"),
                    rec_table,
                    sizing_mode="stretch_width",
                ),
                info["label"],
            )
            run_tabs.append(panel)

        build_interactive_html("\n".join(summary_lines), pattern_df, overlap_df, run_tabs, html_path)

    print(f"[INFO] score_col={score_col}")
    print(f"[INFO] outputs: {union_path}, {overlap_path}, {pattern_path}")
    for info in run_info:
        print(f"[INFO] {info['label']}: {info['top_path']}, {info['rec_path']}")
    if not args.no_html:
        print(f"[INFO] html: {html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
