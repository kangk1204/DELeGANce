#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Summarize top hits from 05_hybrid_annot.tsv.

- Ranks by HitScore_GLM (fallback: HitScore_RS)
- Shows most frequent BBs in top N
- Recommends top X compounds after simple pass filters
"""

import argparse
import base64
import os
import re
import sys
from io import BytesIO
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

_HAS_BOKEH = True
try:
    from bokeh.io import output_file, save
    from bokeh.layouts import column
    from bokeh.models import ColumnDataSource, DataTable, Div, HoverTool, NumberFormatter, TableColumn
    from bokeh.plotting import figure
except Exception:
    _HAS_BOKEH = False

_HAS_RDKIT = True
_RDKIT_IMPORT_ERR = ""
try:
    from rdkit import Chem  # type: ignore
    from rdkit.Chem import Draw  # type: ignore
except Exception as e:
    _HAS_RDKIT = False
    _RDKIT_IMPORT_ERR = str(e)

_TRANSPARENT_PNG = (
    "data:image/png;base64,"
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR4nGNgYAAAAAMA"
    "ASsJTYQAAAAASUVORK5CYII="
)


def _pick_latest(paths: List[str]) -> Optional[str]:
    if not paths:
        return None
    return max(paths, key=lambda p: os.path.getmtime(p))


def resolve_hybrid_path(base: str, preset: Optional[str]) -> str:
    # base can be:
    # - 05_hybrid_annot.tsv
    # - 03_normalized/<preset>
    # - run root (contains 03_normalized)
    if os.path.isfile(base):
        return base
    if os.path.isdir(base):
        direct = os.path.join(base, "05_hybrid_annot.tsv")
        if os.path.isfile(direct):
            return direct
        # base is run root
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


def _to_int_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").fillna(0).astype(int)


def _strip_lib_suffix(val: str) -> str:
    if val is None:
        return "NA"
    s = str(val).strip()
    if s in ("", "NA", "nan", "None"):
        return "NA"
    return re.sub(r"_LIB[\w\.-]+$", "", s)


def _pick_y_column(cols: List[str], score_col: str) -> str:
    for cand in ["LFC_R2_vs_DEL2", "LFC_R1_vs_DEL2", "mean_R2_norm", "mean_R1_norm", "DEL2_norm"]:
        if cand in cols:
            return cand
    return score_col


def _safe_float_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def smiles_to_base64(smiles: str, size=(140, 140)) -> str:
    if not isinstance(smiles, str) or not smiles.strip():
        return _TRANSPARENT_PNG
    if not _HAS_RDKIT:
        return _TRANSPARENT_PNG
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return _TRANSPARENT_PNG
        img = Draw.MolToImage(mol, size=size)
        buf = BytesIO()
        img.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{b64}"
    except Exception:
        return _TRANSPARENT_PNG


def build_bb_smiles_map(df_top: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    bb_smiles_cols = {
        "BB1_x": "bb1_smiles",
        "BB2_x": "bb2_smiles",
        "BB3_x": "bb3_smiles",
        "BB4_x": "bb4_smiles",
    }
    out: Dict[str, Dict[str, str]] = {}
    for bb_col, smi_col in bb_smiles_cols.items():
        if bb_col not in df_top.columns or smi_col not in df_top.columns:
            continue
        sub = df_top[[bb_col, smi_col]].dropna()
        mapping: Dict[str, str] = {}
        for bb, smi in sub.itertuples(index=False):
            bb_key = _strip_lib_suffix(str(bb))
            if bb_key in mapping:
                continue
            if isinstance(smi, str) and smi.strip():
                mapping[bb_key] = smi.strip()
        out[bb_col] = mapping
    return out


def build_interactive_html(df_top: pd.DataFrame, rec: pd.DataFrame, freq_df: pd.DataFrame,
                           out_html: str, score_col: str, top_n: int, rec_n: int,
                           bb_smiles_map: Optional[Dict[str, Dict[str, str]]] = None) -> None:
    if not _HAS_BOKEH:
        print("[WARN] bokeh is not available; skipping interactive HTML output.")
        return
    if not _HAS_RDKIT:
        print(f"[WARN] RDKit not available ({_RDKIT_IMPORT_ERR}); BB hover images will be blank.")

    cols = df_top.columns.tolist()
    y_col = _pick_y_column(cols, score_col)

    df_plot = df_top.copy()
    df_plot["rank"] = np.arange(1, len(df_plot) + 1, dtype=int)
    if "Consensus_hit" in df_plot.columns:
        df_plot["consensus_label"] = _to_int_series(df_plot["Consensus_hit"]).map({1: "pass"}).fillna("fail")
    else:
        df_plot["consensus_label"] = "na"

    for col in [score_col, y_col]:
        df_plot[col] = _safe_float_series(df_plot[col])

    if "mean_R1_norm" in df_plot.columns:
        size = np.log10(_safe_float_series(df_plot["mean_R1_norm"]).fillna(0.0) + 1.0) * 2.0 + 6.0
        df_plot["size"] = np.clip(size, 4.0, 12.0)
    else:
        df_plot["size"] = 8.0

    for col in ["ID_x", "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x", "LIB_ID_x"]:
        if col not in df_plot.columns:
            df_plot[col] = "NA"

    cds = ColumnDataSource(df_plot)
    hover = HoverTool(tooltips=[
        ("rank", "@rank"),
        ("ID", "@ID_x"),
        ("LIB", "@LIB_ID_x"),
        ("BB1", "@BB1_x"),
        ("BB2", "@BB2_x"),
        ("BB3", "@BB3_x"),
        ("BB4", "@BB4_x"),
        ("CP", "@CP_x"),
        (score_col, f"@{{{score_col}}}"),
        (y_col, f"@{{{y_col}}}"),
    ])

    p_scatter = figure(
        title=f"Top {top_n} hits: {score_col} vs {y_col}",
        x_axis_label=score_col,
        y_axis_label=y_col,
        height=420,
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )
    p_scatter.add_tools(hover)
    if "consensus_label" in df_plot.columns:
        colors = {"pass": "#2c7fb8", "fail": "#bdbdbd", "na": "#bdbdbd"}
        df_plot["color"] = df_plot["consensus_label"].map(colors).fillna("#bdbdbd")
        cds = ColumnDataSource(df_plot)
        p_scatter.scatter(x=score_col, y=y_col, size="size", marker="circle", color="color", alpha=0.75, source=cds)
    else:
        p_scatter.scatter(x=score_col, y=y_col, size="size", marker="circle", color="#2c7fb8", alpha=0.75, source=cds)

    # BB frequency bar charts
    bb_plots = []
    img_cache: Dict[str, str] = {}
    if not freq_df.empty and "bb_col" in freq_df.columns:
        for bb_col in freq_df["bb_col"].dropna().unique().tolist():
            sub = freq_df[freq_df["bb_col"] == bb_col].copy()
            if sub.empty:
                continue
            smi_map = (bb_smiles_map or {}).get(bb_col, {})
            sub["bb"] = sub["bb"].astype(str)
            sub["lib_ids"] = sub.get("lib_ids", "").astype(str)
            sub["lib_counts"] = sub.get("lib_counts", "").astype(str)
            sub = sub.sort_values("count", ascending=False)
            sub["smiles"] = sub["bb"].map(lambda bb: smi_map.get(bb, ""))
            def _img_for_smiles(smi: str) -> str:
                if not smi:
                    return _TRANSPARENT_PNG
                if smi in img_cache:
                    return img_cache[smi]
                img_cache[smi] = smiles_to_base64(smi)
                return img_cache[smi]
            sub["img"] = sub["smiles"].map(_img_for_smiles)
            x_range = sub["bb"].tolist()
            p = figure(
                title=f"{bb_col} frequency (top {len(x_range)})",
                x_range=x_range,
                height=260,
                tools="pan,reset,save",
            )
            src = ColumnDataSource(sub)
            p.vbar(x="bb", top="count", width=0.9, source=src, color="#74a9cf")
            p.xaxis.major_label_orientation = 1.0
            p.add_tools(HoverTool(tooltips="""
<div>
  <div><b>@bb</b> (@count, @pct{0.0}%)</div>
  <div>LIB_IDs: @lib_ids</div>
  <div>LIB counts: @lib_counts</div>
  <div style="max-width: 340px; word-break: break-all;">@smiles</div>
  <div><img src="@img" width="140"></div>
</div>
"""))
            bb_plots.append(p)

    # Recommended hits table
    table_cols = []
    table_fields = []
    for name in ["ID_x", "HitScore_GLM", "HitScore_RS", "BB1_x", "BB2_x", "BB3_x", "BB4_x",
                 "CP_x", "LFC_R1_vs_DEL2", "LFC_R2_vs_DEL2", "q_BEAD", "q_BEAD_R2", "q_BoostPaired"]:
        if name in rec.columns:
            table_fields.append(name)
    rec_table = rec[table_fields].copy() if table_fields else rec.copy()
    for col in rec_table.columns:
        if rec_table[col].dtype.kind in "if":
            rec_table[col] = rec_table[col].round(4)
            table_cols.append(TableColumn(field=col, title=col, formatter=NumberFormatter(format="0.0000")))
        else:
            table_cols.append(TableColumn(field=col, title=col))

    table = None
    if not rec_table.empty and table_cols:
        table = DataTable(
            source=ColumnDataSource(rec_table),
            columns=table_cols,
            height=280,
            width=1200,
            index_position=None,
        )

    info = Div(text=(
        f"<h2>Top {top_n} hit summary</h2>"
        f"<ul>"
        f"<li>Score column: <b>{score_col}</b></li>"
        f"<li>Recommended hits: <b>{len(rec)}</b> (target {rec_n})</li>"
        f"<li>Color: blue=Consensus_hit pass, gray=others</li>"
        f"</ul>"
    ))

    items = [info, p_scatter]
    if table is not None:
        items.append(Div(text="<h3>Recommended hits</h3>"))
        items.append(table)
    if bb_plots:
        items.append(Div(text="<h3>BB frequency in top hits</h3>"))
        items.extend(bb_plots)

    layout = column(*items, sizing_mode="stretch_width")
    output_file(out_html, title=f"Top {top_n} hits summary")
    save(layout)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True,
                    help="Run root (DELeGANce_out/<run>), 03_normalized/<preset>, or 05_hybrid_annot.tsv")
    ap.add_argument("--preset", default=None,
                    help="Optional subdir under 03_normalized (e.g., glm_full_dev_cuda_fp32)")
    ap.add_argument("--top-n", type=int, default=100, help="Top N hits to analyze")
    ap.add_argument("--recommend-n", type=int, default=10, help="Number of recommended hits to return")
    ap.add_argument("--bb-top-k", type=int, default=5, help="Top K BB frequencies per BB column")
    ap.add_argument("--score-col", default=None, help="Override score column (default: HitScore_GLM)")
    ap.add_argument("--out-prefix", default=None, help="Output prefix (default: derived from input dir)")
    ap.add_argument("--html-out", default=None, help="Write interactive HTML to this path")
    ap.add_argument("--no-html", action="store_true", help="Skip interactive HTML output")
    args = ap.parse_args()

    top_n = max(1, int(args.top_n))
    rec_n = max(1, int(args.recommend_n))
    bb_top_k = max(1, int(args.bb_top_k))

    hybrid = resolve_hybrid_path(args.output_dir, args.preset)
    out_dir = os.path.dirname(hybrid)
    prefix = args.out_prefix or os.path.join(out_dir, f"top{top_n}")

    # Determine columns
    header = pd.read_csv(hybrid, sep="\t", nrows=0).columns.tolist()
    score_col = args.score_col
    if not score_col:
        score_col = "HitScore_GLM" if "HitScore_GLM" in header else "HitScore_RS"
    if score_col not in header:
        raise SystemExit(f"[ERROR] score column not found: {score_col}")

    want = [
        "LIB_ID_x", "ID_x", "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x", "cycles",
        "bb1_smiles", "bb2_smiles", "bb3_smiles", "bb4_smiles",
        "HitScore_GLM", "HitScore_RS", "HitScore_pct", "SynthonScore",
        "GLM_hit", "RS_pass", "Consensus_hit", "NEG_hard_fail", "NEG_center_fail",
        "pass_filters", "fail_reasons",
        "mean_R1_norm", "mean_R2_norm", "DEL2_norm",
        "LFC_R1_vs_DEL2", "LFC_R2_vs_DEL2",
        "LFC_NEG_R1_vs_DEL2", "LFC_NEG_R2_vs_DEL2",
        "q_DEL2", "q_BEAD", "q_BEAD_R2", "q_BoostPaired",
    ]
    usecols = [c for c in want if c in header] + [score_col]
    usecols = list(dict.fromkeys(usecols))

    df = pd.read_csv(hybrid, sep="\t", usecols=usecols)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(float("-inf"))

    df_top = df.sort_values(score_col, ascending=False).head(top_n).copy()

    # BB frequencies (BB ID based; include LIB_ID list in metadata)
    bb_cols = [c for c in ["BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x"] if c in df_top.columns]
    freq_rows = []
    for col in bb_cols:
        tmp = pd.DataFrame({"bb": df_top[col].fillna("NA").astype(str)})
        tmp["bb"] = tmp["bb"].map(_strip_lib_suffix)
        if "LIB_ID_x" in df_top.columns:
            tmp["lib_id"] = df_top["LIB_ID_x"].fillna("NA").astype(str)
        else:
            tmp["lib_id"] = "NA"

        grp = tmp.groupby("bb", as_index=False).size().rename(columns={"size": "count"})
        grp = grp.sort_values(["count", "bb"], ascending=[False, True]).head(bb_top_k)

        lib_counts = tmp.groupby(["bb", "lib_id"]).size().reset_index(name="lib_count")
        lib_ids_map = (
            lib_counts.groupby("bb")["lib_id"]
            .apply(lambda s: ",".join(sorted(set(map(str, s)))))
            .to_dict()
        )
        lib_counts_map: Dict[str, str] = {}
        for bb_val, sub in lib_counts.groupby("bb", sort=False):
            sub = sub.sort_values("lib_count", ascending=False)
            pairs = [f"{r.lib_id}:{int(r.lib_count)}" for r in sub.itertuples(index=False)]
            lib_counts_map[str(bb_val)] = ", ".join(pairs)

        for _, row in grp.iterrows():
            bb = str(row["bb"])
            cnt = int(row["count"])
            freq_rows.append({
                "bb_col": col,
                "bb": bb,
                "count": cnt,
                "pct": (float(cnt) / float(len(df_top))) * 100.0,
                "lib_ids": lib_ids_map.get(bb, ""),
                "lib_counts": lib_counts_map.get(bb, ""),
            })
    freq_df = pd.DataFrame(freq_rows)
    if not freq_df.empty:
        freq_df = freq_df.sort_values(["bb_col", "count", "bb"], ascending=[True, False, True])
    freq_path = f"{prefix}_bb_frequency.tsv"
    freq_df.to_csv(freq_path, sep="\t", index=False)

    # Recommended hits (filter by available pass flags)
    rec = df_top.copy()
    if "Consensus_hit" in rec.columns:
        rec = rec[_to_int_series(rec["Consensus_hit"]) == 1]
    if "GLM_hit" in rec.columns:
        rec = rec[_to_int_series(rec["GLM_hit"]) == 1]
    if "RS_pass" in rec.columns:
        rec = rec[_to_int_series(rec["RS_pass"]) == 1]
    if "NEG_hard_fail" in rec.columns:
        rec = rec[_to_int_series(rec["NEG_hard_fail"]) == 0]

    rec = rec.sort_values(score_col, ascending=False).head(rec_n).copy()

    # Save tables
    top_path = f"{prefix}_hits.tsv"
    rec_path = f"{prefix}_recommended.tsv"
    df_top.to_csv(top_path, sep="\t", index=False)
    rec.to_csv(rec_path, sep="\t", index=False)

    html_path = args.html_out or f"{prefix}_interactive.html"
    if not args.no_html:
        bb_smiles_map = build_bb_smiles_map(df_top)
        build_interactive_html(df_top, rec, freq_df, html_path, score_col, top_n, rec_n, bb_smiles_map)

    # Console summary
    print(f"[INFO] input={hybrid}")
    print(f"[INFO] score_col={score_col}")
    print(f"[INFO] top_n={len(df_top)}  recommended={len(rec)}")
    print(f"[INFO] outputs: {top_path}, {rec_path}, {freq_path}" + ("" if args.no_html else f", {html_path}"))
    if bb_cols and not freq_df.empty:
        print("[INFO] most frequent BB per column:")
        for col in bb_cols:
            sub = freq_df[freq_df["bb_col"] == col].sort_values("count", ascending=False)
            if not sub.empty:
                row = sub.iloc[0]
                bb = row.get("bb", "NA")
                cnt = int(row.get("count", 0))
                pct = float(row.get("pct", 0.0))
                libs = row.get("lib_ids", "")
                lib_note = f" [LIB_IDs: {libs}]" if libs else ""
                print(f"  - {col}: {bb} ({cnt}, {pct:.1f}%)" + lib_note)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
