#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Selectivity report across active/inactive/both runs with diversity clustering.
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
    from bokeh.layouts import column, row
    from bokeh.models import (ColumnDataSource, CustomJS, DataTable, Div, HoverTool,
                              HTMLTemplateFormatter, NumberFormatter, Panel, TableColumn,
                              Tabs, TabPanel, Button, TextInput, MultiSelect)
    from bokeh.plotting import figure
except Exception:
    _HAS_BOKEH = False

_HAS_TABPANEL = False
try:
    from bokeh.models import TabPanel as _BokehTabPanel  # type: ignore
    _HAS_TABPANEL = True
except Exception:
    _HAS_TABPANEL = False

_HAS_RDKIT = True
_RDKIT_IMPORT_ERR = ""
try:
    from rdkit import Chem, DataStructs, RDLogger  # type: ignore
    from rdkit.Chem import AllChem, Draw  # type: ignore
    from rdkit.ML.Cluster import Butina  # type: ignore
    RDLogger.DisableLog("rdApp.warning")
except Exception as e:
    _HAS_RDKIT = False
    _RDKIT_IMPORT_ERR = str(e)

_TRANSPARENT_PNG = (
    "data:image/png;base64,"
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR4nGNgYAAAAAMA"
    "ASsJTYQAAAAASUVORK5CYII="
)


def _style_block() -> str:
    return """
<style>
.summary-cards { display:flex; flex-wrap:wrap; gap:12px; margin:6px 0 10px 0; }
.summary-card { background:#f7f7fb; border:1px solid #dedee6; border-radius:8px; padding:8px 12px; min-width:160px; }
.summary-label { font-size:11px; color:#666; letter-spacing:0.2px; text-transform:uppercase; }
.summary-value { font-size:18px; font-weight:600; }
.summary-table { border-collapse:collapse; margin:8px 0 12px 0; width:100%; max-width:520px; }
.summary-table th, .summary-table td { border:1px solid #e0e0e0; padding:6px 8px; font-size:12px; }
.summary-table th { background:#f5f5f7; text-align:left; }
.table-controls { display:flex; flex-wrap:wrap; gap:8px; align-items:center; margin:6px 0 8px 0; }
.table-count { font-size:12px; color:#555; }
.struct-thumb { position:relative; display:inline-block; }
.struct-thumb .struct-popup { display:none; position:absolute; z-index:10; left:70px; top:-10px; background:#fff; padding:6px; border:1px solid #ddd; box-shadow:0 2px 6px rgba(0,0,0,0.2); }
.struct-thumb:hover .struct-popup { display:block; }
.group-badge { display:inline-block; padding:2px 8px; border-radius:10px; font-size:11px; font-weight:600; color:#fff; }
</style>
"""


def _summary_cards(cards: List[Tuple[str, str]]) -> str:
    parts = ["<div class='summary-cards'>"]
    for label, value in cards:
        parts.append(
            "<div class='summary-card'>"
            f"<div class='summary-label'>{label}</div>"
            f"<div class='summary-value'>{value}</div>"
            "</div>"
        )
    parts.append("</div>")
    return "\n".join(parts)


def _group_color(code: str) -> str:
    return {
        "Active-specific": "#1b9e77",
        "Inactive-specific": "#d95f02",
        "Both-specific": "#7570b3",
        "Other": "#9e9e9e",
    }.get(code, "#9e9e9e")


def _add_group_badge(df: pd.DataFrame) -> pd.DataFrame:
    if "group" not in df.columns:
        return df
    out = df.copy()
    code_series = out["group_code"] if "group_code" in out.columns else out["group"]
    badges = []
    for label, code in zip(out["group"].astype(str).tolist(), code_series.astype(str).tolist()):
        color = _group_color(code)
        badges.append(f"<span class='group-badge' style='background:{color}'>{label}</span>")
    out["group_badge"] = badges
    return out


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
            if paths:
                return max(paths, key=lambda p: os.path.getmtime(p))
    raise FileNotFoundError(f"05_hybrid_annot.tsv not found under: {base}")


def _strip_lib_suffix(val: str) -> str:
    if val is None:
        return "NA"
    s = str(val).strip()
    if s in ("", "NA", "nan", "None"):
        return "NA"
    return re.sub(r"_LIB[\w\.-]+$", "", s)


def _to_int_series(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").fillna(0).astype(int)


def _make_compound_key(df: pd.DataFrame) -> pd.Series:
    bb_cols = [c for c in ["BB1_x", "BB2_x", "BB3_x", "BB4_x"] if c in df.columns]
    if not bb_cols:
        return df.get("ID_x", pd.Series(["NA"] * len(df))).astype(str)
    tmp = df[bb_cols].fillna("NA").astype(str).copy()
    for c in bb_cols:
        tmp[c] = tmp[c].map(_strip_lib_suffix)
    return tmp[bb_cols].agg("|".join, axis=1)


def _strip_bb_suffix_cols(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = out[col].fillna("NA").astype(str).map(_strip_lib_suffix)
    return out


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


def _infer_score_col(paths: List[str], preset: Optional[str], override: Optional[str]) -> str:
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


def _sample_cols_from_header(header: List[str]) -> List[str]:
    cpm_cols = [c for c in header if c.endswith("_CPM") and c[:-4] in header]
    raw_cols = [c[:-4] for c in cpm_cols]
    return raw_cols + cpm_cols


def _prefix_sample_cols(cols: List[str], prefix: str) -> Tuple[List[str], Dict[str, str]]:
    renamed = {c: f"{prefix}{c}" for c in cols}
    return [renamed[c] for c in cols], renamed


def _coalesce_unprefixed_samples(df: pd.DataFrame, base_cols: List[str],
                                 prefixes: List[str]) -> pd.DataFrame:
    out = df.copy()
    for col in base_cols:
        if col not in out.columns:
            continue
        series = pd.to_numeric(out[col], errors="coerce")
        for pre in prefixes:
            pref_col = f"{pre}{col}"
            if pref_col in out.columns:
                series = series.where(series.notna(), pd.to_numeric(out[pref_col], errors="coerce"))
        out[col] = series
    return out


def _sort_sample_cols_by_category(cols: List[str], prefixes: List[str]) -> List[str]:
    def _natural_key(s: str) -> List[object]:
        return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\\d+)", s)]

    def _sample_id(col: str, prefix: str) -> str:
        base = col[len(prefix):] if prefix and col.startswith(prefix) else col
        return base[:-4] if base.endswith("_CPM") else base

    def _order_group(group_cols: List[str], prefix: str) -> List[str]:
        sample_ids: List[str] = []
        for col in group_cols:
            sid = _sample_id(col, prefix)
            if sid not in sample_ids:
                sample_ids.append(sid)
        ordered_ids = sorted(sample_ids, key=_natural_key)
        raw_cols: List[str] = []
        cpm_cols: List[str] = []
        for sid in ordered_ids:
            raw = f"{prefix}{sid}" if prefix else sid
            cpm = f"{raw}_CPM"
            if raw in group_cols and raw not in raw_cols:
                raw_cols.append(raw)
            if cpm in group_cols and cpm not in cpm_cols:
                cpm_cols.append(cpm)
        extras = [c for c in group_cols if c not in raw_cols and c not in cpm_cols]
        extras = sorted(extras, key=_natural_key)
        return raw_cols + cpm_cols + extras

    seen = set()
    ordered: List[str] = []
    base_group = [c for c in cols if not any(c.startswith(p) for p in prefixes)]
    for col in _order_group(base_group, ""):
        if col in cols and col not in seen:
            ordered.append(col)
            seen.add(col)
    for prefix in prefixes:
        pref_group = [c for c in cols if c.startswith(prefix)]
        for col in _order_group(pref_group, prefix):
            if col in cols and col not in seen:
                ordered.append(col)
                seen.add(col)
    for col in cols:
        if col not in seen:
            ordered.append(col)
            seen.add(col)
    return ordered


def _dedupe_identical_sample_cols(df: pd.DataFrame, cols: List[str], prefixes: List[str]) -> List[str]:
    if df is None or df.empty or not cols or not prefixes:
        return cols
    prefix_order = {p: i for i, p in enumerate(prefixes)}
    def _col_prefix(col: str) -> str:
        for p in prefixes:
            if col.startswith(p):
                return p
        return ""
    base_groups: Dict[str, List[str]] = {}
    for col in cols:
        pref = _col_prefix(col)
        if not pref:
            continue
        base = col[len(pref):]
        base_groups.setdefault(base, []).append(col)
    drop: set = set()
    sentinel = "__NA__"
    for base, group_cols in base_groups.items():
        if len(group_cols) < 2:
            continue
        group_cols.sort(key=lambda c: prefix_order.get(_col_prefix(c), 0))
        ref = group_cols[0]
        ref_series = df[ref].fillna(sentinel)
        for col in group_cols[1:]:
            if ref_series.equals(df[col].fillna(sentinel)):
                drop.add(col)
    if not drop:
        return cols
    return [c for c in cols if c not in drop]


def _neg_thresholds_for_run(df_best: pd.DataFrame, sample_cols: List[str], prefix: str,
                            neg_samples: List[str], neg_pct: float) -> Dict[str, float]:
    thresholds: Dict[str, float] = {}
    if df_best is None or df_best.empty:
        return thresholds
    if not neg_samples:
        return thresholds
    for s in neg_samples:
        cpm = f"{s}_CPM"
        raw = s
        if cpm in sample_cols and cpm in df_best.columns:
            series = pd.to_numeric(df_best[cpm], errors="coerce")
            col_key = f"{prefix}{cpm}"
        elif raw in sample_cols and raw in df_best.columns:
            series = pd.to_numeric(df_best[raw], errors="coerce")
            col_key = f"{prefix}{raw}"
        else:
            continue
        if series.notna().any():
            thresholds[col_key] = float(np.nanpercentile(series, neg_pct))
    return thresholds


def _filter_enrich_cols(prefixed_cols: List[str], prefix: str,
                        exclude_bases: set) -> List[str]:
    out = []
    for col in prefixed_cols:
        base = col[len(prefix):] if prefix and col.startswith(prefix) else col
        base_no_cpm = base[:-4] if base.endswith("_CPM") else base
        if base_no_cpm.upper() in exclude_bases:
            continue
        out.append(col)
    return out


def _row_agg(df: pd.DataFrame, cols: List[str], agg: str) -> pd.Series:
    if not cols:
        return pd.Series([np.nan] * len(df), index=df.index)
    use = [c for c in cols if c in df.columns]
    if not use:
        return pd.Series([np.nan] * len(df), index=df.index)
    arr = df[use].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if agg == "mean":
        vals = np.nanmean(arr, axis=1)
    elif agg == "max":
        vals = np.nanmax(arr, axis=1)
    else:
        vals = np.nanmedian(arr, axis=1)
    return pd.Series(vals, index=df.index)


def _pick_enrich_cols(prefixed_cols: List[str]) -> List[str]:
    cpm = [c for c in prefixed_cols if c.endswith("_CPM")]
    if cpm:
        return cpm
    return [c for c in prefixed_cols if not c.endswith("_CPM")]


def _rank_maps(df: pd.DataFrame, score_col: str) -> Dict[str, Dict[str, float]]:
    df = df.copy()
    df["compound_key"] = _make_compound_key(df)
    score_map = df.groupby("compound_key")[score_col].max().to_dict()
    items = [(k, v) for k, v in score_map.items() if np.isfinite(v)]
    items.sort(key=lambda x: (-x[1], x[0]))
    rank_map = {k: i + 1 for i, (k, _) in enumerate(items)}
    denom = max(1, len(items) - 1)
    rank_pct_map = {k: 100.0 * (1.0 - (r - 1) / denom) for k, r in rank_map.items()}
    if "HitScore_pct" in df.columns:
        hitpct_map = df.groupby("compound_key")["HitScore_pct"].max().to_dict()
    else:
        hitpct_map = rank_pct_map
    scores = np.array([v for _, v in items], dtype=float)
    if scores.size == 0:
        mu = 0.0
        sd = 1.0
    else:
        mu = float(np.mean(scores))
        sd = float(np.std(scores))
        if sd == 0.0:
            sd = 1.0
    score_z_map = {k: ((score_map[k] - mu) / sd) if np.isfinite(score_map[k]) else np.nan for k in score_map}
    return {
        "score_map": score_map,
        "rank_map": rank_map,
        "rank_pct_map": rank_pct_map,
        "hitpct_map": hitpct_map,
        "score_z_map": score_z_map,
    }


def smiles_to_base64(smiles: str, size=(110, 110)) -> str:
    if not isinstance(smiles, str) or not smiles.strip():
        return _TRANSPARENT_PNG
    if not _HAS_RDKIT:
        return _TRANSPARENT_PNG
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return _TRANSPARENT_PNG
        img = Draw.MolToImage(mol, size=size)
        from io import BytesIO
        import base64
        buf = BytesIO()
        img.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{b64}"
    except Exception:
        return _TRANSPARENT_PNG


def _fp_from_smiles(smiles: str, radius: int, nbits: int):
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def _bbavg_similarity(fps_a: List[Optional[object]], fps_b: List[Optional[object]]) -> float:
    sims = []
    for fa, fb in zip(fps_a, fps_b):
        if fa is None or fb is None:
            continue
        sims.append(DataStructs.TanimotoSimilarity(fa, fb))
    if not sims:
        return 0.0
    return float(sum(sims) / len(sims))


def _tie_break_arrays(df: pd.DataFrame, primary_col: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    scores = pd.to_numeric(df.get(primary_col, np.nan), errors="coerce").fillna(float("-inf")).to_numpy()
    if "rank_pct" in df.columns:
        rank_pct = pd.to_numeric(df["rank_pct"], errors="coerce").fillna(0.0).to_numpy()
    elif "HitScore_pct" in df.columns:
        rank_pct = pd.to_numeric(df["HitScore_pct"], errors="coerce").fillna(0.0).to_numpy()
    else:
        rank_pct = np.zeros(len(df), dtype=float)

    cpm_cols = [c for c in df.columns if c.endswith("_CPM")]
    raw_cols = [c[:-4] for c in cpm_cols if c[:-4] in df.columns]
    if cpm_cols:
        cpm_mean = (
            df[cpm_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).mean(axis=1).to_numpy()
        )
    else:
        cpm_mean = np.zeros(len(df), dtype=float)
    if raw_cols:
        raw_sum = (
            df[raw_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).sum(axis=1).to_numpy()
        )
    else:
        raw_sum = np.zeros(len(df), dtype=float)
    return scores, rank_pct, cpm_mean, raw_sum


def _cluster_candidates(df: pd.DataFrame, rep_score_col: str, sim_cutoff: float,
                        radius: int, nbits: int, mode: str) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    if not _HAS_RDKIT:
        print(f"[WARN] RDKit not available ({_RDKIT_IMPORT_ERR}); clustering skipped.")
        out = df.copy()
        out["cluster_id"] = "NA"
        out["cluster_size"] = 0
        out["cluster_rep"] = 0
        return out, None

    smi_cols = [c for c in ["bb1_smiles", "bb2_smiles", "bb3_smiles", "bb4_smiles"] if c in df.columns]
    if not smi_cols:
        print("[WARN] No bb*_smiles columns available; clustering skipped.")
        out = df.copy()
        out["cluster_id"] = "NA"
        out["cluster_size"] = 0
        out["cluster_rep"] = 0
        return out, None

    fps = []
    valid_pos = []
    for i, row in df.reset_index(drop=True).iterrows():
        if mode == "compound_or":
            fp = None
            for col in smi_cols:
                smi = row.get(col, "")
                fpi = _fp_from_smiles(smi, radius, nbits)
                if fpi is None:
                    continue
                fp = fpi if fp is None else (fp | fpi)
            if fp is None:
                fps.append(None)
            else:
                fps.append(fp)
                valid_pos.append(i)
        else:
            row_fps = []
            for col in smi_cols:
                smi = row.get(col, "")
                row_fps.append(_fp_from_smiles(smi, radius, nbits))
            if any(fp is not None for fp in row_fps):
                fps.append(row_fps)
                valid_pos.append(i)
            else:
                fps.append(None)

    if not valid_pos:
        print("[WARN] No valid fingerprints; clustering skipped.")
        out = df.copy()
        out["cluster_id"] = "NA"
        out["cluster_size"] = 0
        out["cluster_rep"] = 0
        return out, None

    valid_fps = [fps[i] for i in valid_pos]
    dists = []
    if mode == "compound_or":
        for i in range(1, len(valid_fps)):
            sims = DataStructs.BulkTanimotoSimilarity(valid_fps[i], valid_fps[:i])
            dists.extend([1.0 - s for s in sims])
    else:
        for i in range(1, len(valid_fps)):
            fi = valid_fps[i]
            for j in range(i):
                sim = _bbavg_similarity(fi, valid_fps[j])
                dists.append(1.0 - sim)
    cutoff = max(0.0, min(1.0, 1.0 - float(sim_cutoff)))
    clusters = Butina.ClusterData(dists, len(valid_fps), cutoff, isDistData=True)

    out = df.reset_index(drop=True).copy()
    out["cluster_id"] = "NA"
    out["cluster_size"] = 0
    out["cluster_rep"] = 0
    out["cluster_medoid"] = 0

    scores, rank_pct, cpm_mean, raw_sum = _tie_break_arrays(out, rep_score_col)
    cluster_rows = []
    for cid, members in enumerate(clusters, start=1):
        members_pos = [valid_pos[i] for i in members]
        if not members_pos:
            continue
        size = len(members_pos)
        best_pos = max(
            members_pos,
            key=lambda idx: (scores[idx], rank_pct[idx], cpm_mean[idx], raw_sum[idx], -idx),
        )
        if len(members) == 1:
            medoid_valid = members[0]
            medoid_pos = members_pos[0]
            medoid_avg = 1.0
        else:
            medoid_scores: Dict[int, float] = {}
            for i in members:
                if mode == "compound_or":
                    others = [valid_fps[j] for j in members if j != i]
                    sims = DataStructs.BulkTanimotoSimilarity(valid_fps[i], others) if others else []
                    medoid_scores[i] = float(np.mean(sims)) if sims else 1.0
                else:
                    sims = []
                    for j in members:
                        if j == i:
                            continue
                        sims.append(_bbavg_similarity(valid_fps[i], valid_fps[j]))
                    medoid_scores[i] = float(np.mean(sims)) if sims else 1.0
            medoid_valid = max(
                members,
                key=lambda i: (
                    medoid_scores.get(i, 0.0),
                    scores[valid_pos[i]],
                    rank_pct[valid_pos[i]],
                    cpm_mean[valid_pos[i]],
                    raw_sum[valid_pos[i]],
                    -valid_pos[i],
                ),
            )
            medoid_pos = valid_pos[medoid_valid]
            medoid_avg = medoid_scores.get(medoid_valid, 0.0)
        for idx in members_pos:
            out.at[idx, "cluster_id"] = cid
            out.at[idx, "cluster_size"] = size
            out.at[idx, "cluster_rep"] = 1 if idx == best_pos else 0
            out.at[idx, "cluster_medoid"] = 1 if idx == medoid_pos else 0
        rep = out.loc[best_pos]
        medoid = out.loc[medoid_pos]
        cluster_rows.append({
            "cluster_id": cid,
            "cluster_size": size,
            "rep_index": int(best_pos),
            "rep_ID_x": rep.get("ID_x", "NA"),
            "rep_score": float(rep.get(rep_score_col, np.nan)) if pd.notna(rep.get(rep_score_col, np.nan)) else np.nan,
            "rep_compound_key": rep.get("compound_key", "NA"),
            "medoid_index": int(medoid_pos),
            "medoid_ID_x": medoid.get("ID_x", "NA"),
            "medoid_score": float(medoid.get(rep_score_col, np.nan)) if pd.notna(medoid.get(rep_score_col, np.nan)) else np.nan,
            "medoid_compound_key": medoid.get("compound_key", "NA"),
            "medoid_avg_sim": float(medoid_avg),
        })

    cluster_df = pd.DataFrame(cluster_rows)
    return out, cluster_df


def _pick_final_from_group(df_div: pd.DataFrame, df_all: pd.DataFrame, n: int) -> pd.DataFrame:
    if n <= 0:
        return df_all.head(0).copy() if df_all is not None else pd.DataFrame()
    df_div = df_div if df_div is not None else pd.DataFrame()
    df_all = df_all if df_all is not None else pd.DataFrame()
    if df_div.empty:
        return df_all.head(n).copy()
    picked = df_div.head(n).copy()
    if len(picked) >= n or df_all.empty:
        return picked
    key = "compound_key" if "compound_key" in df_all.columns else None
    if key and key in picked.columns:
        remainder = df_all[~df_all[key].isin(picked[key])].copy()
    else:
        remainder = df_all.copy()
    remainder = remainder.head(n - len(picked))
    return pd.concat([picked, remainder], ignore_index=True)


def _make_panel(child, title: str):
    if _HAS_TABPANEL:
        return _BokehTabPanel(child=child, title=title)
    return Panel(child=child, title=title)


def _col_width(col: str, series: pd.Series) -> int:
    if col.startswith("bb") and col.endswith("_img"):
        return 90
    if col in ("ID_x", "LIB_ID_x", "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x"):
        base = 150
    elif col in ("compound_key", "BB_SMILES_CONCAT"):
        base = 220
    else:
        base = 100
    try:
        vals = series.dropna().astype(str)
        if len(vals) > 500:
            vals = vals.head(500)
        max_len = max([len(str(col))] + [len(v) for v in vals])
    except Exception:
        max_len = len(str(col))
    max_len = min(max_len, 80)
    width = int(min(max(base, max_len * 8 + 28), 600))
    return width


def _img_formatter(width: int, zoom: int) -> HTMLTemplateFormatter:
    template = (
        "<div class='struct-thumb'>"
        f"<img src=\"<%= value %>\" width=\"{width}\">"
        f"<div class='struct-popup'><img src=\"<%= value %>\" width=\"{zoom}\"></div>"
        "</div>"
    )
    return HTMLTemplateFormatter(template=template)


def _make_table_controls(source: ColumnDataSource, full: ColumnDataSource,
                         filename: str, group_col: Optional[str],
                         page_state: ColumnDataSource, page_size: int) -> row:
    total = 0
    if full.data:
        first_key = next(iter(full.data.keys()))
        total = len(full.data[first_key])
    page_size = max(1, int(page_size))
    shown = min(total, page_size) if total else 0
    pages = (total + page_size - 1) // page_size if total else 0
    page_label = f"1/{pages}" if pages else "0/0"
    search = TextInput(title="Search", value="", placeholder="compound_key / ID / BB / CP")
    if shown:
        count_text = f"<span class='table-count'>Showing 1-{shown} of {total} (page {page_label})</span>"
    else:
        count_text = "<span class='table-count'>Showing 0 of 0 (page 0/0)</span>"
    count_div = Div(text=count_text)
    prev_btn = Button(label="Prev", button_type="default", width=60)
    next_btn = Button(label="Next", button_type="default", width=60)
    reset_btn = Button(label="Reset", button_type="default", width=70)
    download_btn = Button(label="Download CSV", button_type="primary")

    group_select = None
    if group_col and group_col in full.data:
        options = sorted({str(v) for v in full.data[group_col] if v is not None})
        group_select = MultiSelect(title="Filter group", value=options, options=options, size=6)

    filter_js = """
const data = full.data;
const cols = Object.keys(data);
const out = {};
for (const c of cols) { out[c] = []; }
const term = (search.value || '').toLowerCase();
const hasGroup = group !== null;
const groups = hasGroup ? new Set(group.value || []) : null;
const n = cols.length ? data[cols[0]].length : 0;
for (let i = 0; i < n; i++) {
  let keep = true;
  if (hasGroup && groups.size > 0) {
    const g = data[group_col][i];
    if (!groups.has(String(g))) {
      keep = false;
    }
  }
  if (keep && term) {
    let hit = false;
    for (const c of cols) {
      const v = data[c][i];
      if (v === null || v === undefined) continue;
      if (String(v).toLowerCase().includes(term)) { hit = true; break; }
    }
    keep = hit;
  }
  if (keep) {
    for (const c of cols) { out[c].push(data[c][i]); }
  }
}
const total = cols.length ? out[cols[0]].length : 0;
const pageSize = Math.max(1, page_size);
const pages = total > 0 ? Math.ceil(total / pageSize) : 0;
let page = 0;
page_state.data.page[0] = page;
page_state.change.emit();
const start = page * pageSize;
const end = Math.min(start + pageSize, total);
const view = {};
for (const c of cols) { view[c] = out[c].slice(start, end); }
source.data = view;
source.change.emit();
const shown = total === 0 ? '0' : `${start + 1}-${end}`;
const pageLabel = total === 0 ? '0/0' : `${page + 1}/${pages}`;
count_div.text = `<span class='table-count'>Showing ${shown} of ${total} (page ${pageLabel})</span>`;
"""
    args = dict(
        source=source, full=full, search=search, count_div=count_div,
        group=group_select, group_col=group_col or "", page_state=page_state, page_size=page_size
    )
    callback = CustomJS(args=args, code=filter_js)
    search.js_on_change("value", callback)
    if group_select is not None:
        group_select.js_on_change("value", callback)

    reset_js = """
const data = full.data;
search.value = '';
if (group !== null) { group.value = group.options; }
const cols = Object.keys(data);
const out = {};
for (const c of cols) { out[c] = data[c].slice(); }
const total = cols.length ? out[cols[0]].length : 0;
const pageSize = Math.max(1, page_size);
const pages = total > 0 ? Math.ceil(total / pageSize) : 0;
let page = 0;
page_state.data.page[0] = page;
page_state.change.emit();
const start = page * pageSize;
const end = Math.min(start + pageSize, total);
const view = {};
for (const c of cols) { view[c] = out[c].slice(start, end); }
source.data = view;
source.change.emit();
const shown = total === 0 ? '0' : `${start + 1}-${end}`;
const pageLabel = total === 0 ? '0/0' : `${page + 1}/${pages}`;
count_div.text = `<span class='table-count'>Showing ${shown} of ${total} (page ${pageLabel})</span>`;
"""
    reset_btn.js_on_click(CustomJS(args=args, code=reset_js))

    prev_js = """
const data = full.data;
const cols = Object.keys(data);
const out = {};
for (const c of cols) { out[c] = []; }
const term = (search.value || '').toLowerCase();
const hasGroup = group !== null;
const groups = hasGroup ? new Set(group.value || []) : null;
const n = cols.length ? data[cols[0]].length : 0;
for (let i = 0; i < n; i++) {
  let keep = true;
  if (hasGroup && groups.size > 0) {
    const g = data[group_col][i];
    if (!groups.has(String(g))) {
      keep = false;
    }
  }
  if (keep && term) {
    let hit = false;
    for (const c of cols) {
      const v = data[c][i];
      if (v === null || v === undefined) continue;
      if (String(v).toLowerCase().includes(term)) { hit = true; break; }
    }
    keep = hit;
  }
  if (keep) {
    for (const c of cols) { out[c].push(data[c][i]); }
  }
}
const total = cols.length ? out[cols[0]].length : 0;
const pageSize = Math.max(1, page_size);
const pages = total > 0 ? Math.ceil(total / pageSize) : 0;
let page = (page_state.data.page[0] || 0) - 1;
if (pages > 0) {
  page = Math.max(0, Math.min(page, pages - 1));
} else {
  page = 0;
}
page_state.data.page[0] = page;
page_state.change.emit();
const start = page * pageSize;
const end = Math.min(start + pageSize, total);
const view = {};
for (const c of cols) { view[c] = out[c].slice(start, end); }
source.data = view;
source.change.emit();
const shown = total === 0 ? '0' : `${start + 1}-${end}`;
const pageLabel = total === 0 ? '0/0' : `${page + 1}/${pages}`;
count_div.text = `<span class='table-count'>Showing ${shown} of ${total} (page ${pageLabel})</span>`;
"""
    next_js = """
const data = full.data;
const cols = Object.keys(data);
const out = {};
for (const c of cols) { out[c] = []; }
const term = (search.value || '').toLowerCase();
const hasGroup = group !== null;
const groups = hasGroup ? new Set(group.value || []) : null;
const n = cols.length ? data[cols[0]].length : 0;
for (let i = 0; i < n; i++) {
  let keep = true;
  if (hasGroup && groups.size > 0) {
    const g = data[group_col][i];
    if (!groups.has(String(g))) {
      keep = false;
    }
  }
  if (keep && term) {
    let hit = false;
    for (const c of cols) {
      const v = data[c][i];
      if (v === null || v === undefined) continue;
      if (String(v).toLowerCase().includes(term)) { hit = true; break; }
    }
    keep = hit;
  }
  if (keep) {
    for (const c of cols) { out[c].push(data[c][i]); }
  }
}
const total = cols.length ? out[cols[0]].length : 0;
const pageSize = Math.max(1, page_size);
const pages = total > 0 ? Math.ceil(total / pageSize) : 0;
let page = (page_state.data.page[0] || 0) + 1;
if (pages > 0) {
  page = Math.max(0, Math.min(page, pages - 1));
} else {
  page = 0;
}
page_state.data.page[0] = page;
page_state.change.emit();
const start = page * pageSize;
const end = Math.min(start + pageSize, total);
const view = {};
for (const c of cols) { view[c] = out[c].slice(start, end); }
source.data = view;
source.change.emit();
const shown = total === 0 ? '0' : `${start + 1}-${end}`;
const pageLabel = total === 0 ? '0/0' : `${page + 1}/${pages}`;
count_div.text = `<span class='table-count'>Showing ${shown} of ${total} (page ${pageLabel})</span>`;
"""
    prev_btn.js_on_click(CustomJS(args=args, code=prev_js))
    next_btn.js_on_click(CustomJS(args=args, code=next_js))

    download_js = """
function csvEscape(value) {
  if (value === null || value === undefined) return '';
  const s = String(value);
  if (s.includes('\"') || s.includes(',') || s.includes('\\n')) {
    return '\"' + s.replace(/\"/g, '\"\"') + '\"';
  }
  return s;
}
const data = full.data;
const cols = Object.keys(data);
const out = {};
for (const c of cols) { out[c] = []; }
const term = (search.value || '').toLowerCase();
const hasGroup = group !== null;
const groups = hasGroup ? new Set(group.value || []) : null;
const n = cols.length ? data[cols[0]].length : 0;
for (let i = 0; i < n; i++) {
  let keep = true;
  if (hasGroup && groups.size > 0) {
    const g = data[group_col][i];
    if (!groups.has(String(g))) {
      keep = false;
    }
  }
  if (keep && term) {
    let hit = false;
    for (const c of cols) {
      const v = data[c][i];
      if (v === null || v === undefined) continue;
      if (String(v).toLowerCase().includes(term)) { hit = true; break; }
    }
    keep = hit;
  }
  if (keep) {
    for (const c of cols) { out[c].push(data[c][i]); }
  }
}
const nout = cols.length ? out[cols[0]].length : 0;
let csv = cols.join(',') + '\\n';
for (let i = 0; i < nout; i++) {
  const row = cols.map(c => csvEscape(out[c][i]));
  csv += row.join(',') + '\\n';
}
const blob = new Blob([csv], {type: 'text/csv;charset=utf-8;'});
const url = URL.createObjectURL(blob);
const a = document.createElement('a');
a.href = url;
a.download = filename;
a.style.display = 'none';
document.body.appendChild(a);
a.click();
document.body.removeChild(a);
URL.revokeObjectURL(url);
"""
    download_btn.js_on_click(CustomJS(
        args=dict(source=source, full=full, search=search, group=group_select,
                  group_col=group_col or "", filename=filename),
        code=download_js
    ))

    controls = [search]
    if group_select is not None:
        controls.append(group_select)
    controls.extend([prev_btn, next_btn, reset_btn, download_btn, count_div])
    return row(*controls, css_classes=["table-controls"])

def _make_table(df: pd.DataFrame, with_images: bool, max_rows: int,
                with_structure: bool, sample_cols: List[str], page_size: int,
                download_name: str) -> Tuple[Optional[object], Optional[Div]]:
    if df is None or df.empty:
        return None, None
    page_size = max(1, int(page_size))
    if max_rows and len(df) > max_rows:
        df = df.head(max_rows).copy()
    if with_images or with_structure:
        df = _add_images(df)
    df = _add_group_badge(df)
    table_cols = []
    table_fields = [
        "group_badge", "group_rank", "specificity_score", "group_rank_score",
        "active_enrich", "inactive_enrich", "both_enrich", "selectivity_score",
        "active_rank_pct", "inactive_rank_pct", "both_rank_pct",
        "rank", "cluster_id", "cluster_size", "cluster_rep", "cluster_medoid", "ID_x", "LIB_ID_x",
        "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x",
    ]
    if with_images:
        table_fields = [
            "group_badge", "group_rank", "specificity_score", "group_rank_score",
            "active_enrich", "inactive_enrich", "both_enrich", "selectivity_score",
            "active_rank_pct", "inactive_rank_pct", "both_rank_pct",
            "rank", "cluster_id", "cluster_size", "cluster_rep", "cluster_medoid", "ID_x", "LIB_ID_x",
            "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x", "compound_img",
            "bb1_img", "bb2_img", "bb3_img", "bb4_img",
        ]
    if sample_cols:
        table_fields.extend(sample_cols)
    table_fields = [c for c in table_fields if c in df.columns]
    table_df = df[table_fields].copy()
    if with_structure and "structure_html" in df.columns:
        table_df["structure_html"] = df["structure_html"].values
    if "group" in df.columns and "group" not in table_df.columns:
        table_df["group"] = df["group"].values
    if "group_code" in df.columns and "group_code" not in table_df.columns:
        table_df["group_code"] = df["group_code"].values

    img_formatter = _img_formatter(70, 220)
    compound_formatter = _img_formatter(55, 220)
    badge_formatter = HTMLTemplateFormatter(template="<%= value %>")
    nan_fmt = "-"

    for col in table_df.columns:
        if col == "group_badge":
            table_cols.append(TableColumn(field=col, title="group", formatter=badge_formatter, width=_col_width(col, table_df[col])))
            continue
        if col.startswith("bb") and col.endswith("_img"):
            table_cols.append(TableColumn(field=col, title=col, formatter=img_formatter, width=_col_width(col, table_df[col])))
            continue
        if col == "compound_img":
            table_cols.append(TableColumn(field=col, title="compound_img", formatter=compound_formatter, width=_col_width(col, table_df[col])))
            continue
        if col == "structure_html":
            continue
        if table_df[col].dtype.kind in "if":
            if col in ("rank", "cluster_size", "cluster_rep", "cluster_medoid", "group_rank"):
                table_df[col] = pd.to_numeric(table_df[col], errors="coerce").fillna(0).astype(int)
                table_cols.append(TableColumn(field=col, title=col, formatter=NumberFormatter(format="0", nan_format=nan_fmt), width=_col_width(col, table_df[col])))
            elif col in sample_cols and not col.endswith("_CPM"):
                table_df[col] = pd.to_numeric(table_df[col], errors="coerce")
                table_cols.append(TableColumn(field=col, title=col, formatter=NumberFormatter(format="0", nan_format=nan_fmt), width=_col_width(col, table_df[col])))
            else:
                table_df[col] = pd.to_numeric(table_df[col], errors="coerce").round(4)
                table_cols.append(TableColumn(field=col, title=col, formatter=NumberFormatter(format="0.0000", nan_format=nan_fmt), width=_col_width(col, table_df[col])))
        else:
            series = table_df[col].astype(object)
            table_df[col] = series.where(series.notna(), nan_fmt)
            table_cols.append(TableColumn(field=col, title=col, width=_col_width(col, table_df[col])))

    if not table_cols:
        return None, None
    page_df = table_df.head(page_size).copy()
    source = ColumnDataSource(page_df)
    full_source = ColumnDataSource(table_df.copy())
    page_state = ColumnDataSource({"page": [0]})
    row_height = 80 if with_images else 30
    shown_rows = min(len(table_df), page_size) if page_size else len(table_df)
    table_height = 40 + row_height * max(1, shown_rows)
    total_width = sum([c.width or 0 for c in table_cols]) + 40
    table_width = min(5000, max(1200, total_width))
    freeze_col = None
    for col in ["CP_x", "BB4_x", "BB3_x", "BB2_x", "BB1_x", "LIB_ID_x", "ID_x"]:
        if col in table_df.columns:
            freeze_col = col
            break
    frozen_columns = 0
    if freeze_col:
        frozen_columns = int(table_df.columns.get_loc(freeze_col)) + 1
    table = DataTable(
        source=source,
        columns=table_cols,
        height=table_height,
        row_height=row_height,
        width=table_width,
        index_position=None,
        frozen_columns=frozen_columns,
    )
    group_col = "group" if "group" in table_df.columns else None
    controls = _make_table_controls(source, full_source, download_name, group_col, page_state, page_size)
    table_layout = column(controls, table, sizing_mode="stretch_width")
    struct_div = None
    if with_structure and "structure_html" in table_df.columns:
        struct_div = Div(text="<i>Click a row to view structure</i>")
        callback = CustomJS(args=dict(source=source, div=struct_div), code="""
const inds = source.selected.indices;
if (inds.length === 0) {
  div.text = '<i>Click a row to view structure</i>';
  return;
}
const i = inds[0];
const html = source.data['structure_html'][i];
div.text = html || '<i>No structure</i>';
""")
        source.selected.js_on_change("indices", callback)
    return table_layout, struct_div


def _add_images(df: pd.DataFrame) -> pd.DataFrame:
    if not _HAS_RDKIT:
        return df
    out = df.copy()
    for i, col in enumerate(["bb1_smiles", "bb2_smiles", "bb3_smiles", "bb4_smiles"], start=1):
        if col in out.columns:
            out[f"bb{i}_img"] = out[col].map(lambda s: smiles_to_base64(s, size=(90, 90)))
        else:
            out[f"bb{i}_img"] = _TRANSPARENT_PNG
    if "BB_SMILES_CONCAT" in out.columns:
        out["compound_img"] = out["BB_SMILES_CONCAT"].map(lambda s: smiles_to_base64(s, size=(220, 220)))
    else:
        def _fallback_concat(row):
            parts = []
            for col in ["bb1_smiles", "bb2_smiles", "bb3_smiles", "bb4_smiles"]:
                v = row.get(col, "")
                if isinstance(v, str) and v.strip():
                    parts.append(v.strip())
            if parts:
                return smiles_to_base64(".".join(parts), size=(220, 220))
            return _TRANSPARENT_PNG
        out["compound_img"] = out.apply(_fallback_concat, axis=1)

    def _structure_html(row):
        main = row.get("compound_img", _TRANSPARENT_PNG)
        if not main or main == _TRANSPARENT_PNG:
            main_html = "<div style='color:#888;font-size:12px'>No combined structure</div>"
        else:
            main_html = f"<img src=\"{main}\" width=\"220\">"
        parts = []
        for i in range(1, 5):
            img = row.get(f"bb{i}_img", _TRANSPARENT_PNG)
            parts.append(
                "<div style='display:inline-block;margin-right:6px;text-align:center'>"
                f"<img src=\"{img}\" width=\"70\"><div style='font-size:10px'>BB{i}</div></div>"
            )
        return (
            "<div style='display:flex;gap:12px;align-items:center'>"
            f"{main_html}<div>{''.join(parts)}</div></div>"
        )

    out["structure_html"] = out.apply(_structure_html, axis=1)
    return out


def _auto_range_pct(values: pd.Series, pad: float = 2.0, min_span: float = 5.0) -> Tuple[float, float]:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float, copy=False)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return (0.0, 100.0)
    vmin = float(np.min(arr))
    vmax = float(np.max(arr))
    if vmin == vmax:
        vmin -= min_span / 2.0
        vmax += min_span / 2.0
    pad_val = max(pad, (vmax - vmin) * 0.05)
    vmin -= pad_val
    vmax += pad_val
    vmin = max(0.0, vmin)
    vmax = min(100.0, vmax)
    if vmax - vmin < min_span:
        mid = (vmin + vmax) / 2.0
        vmin = max(0.0, mid - min_span / 2.0)
        vmax = min(100.0, mid + min_span / 2.0)
    return (vmin, vmax)


def _parse_range(arg: Optional[str], values: pd.Series) -> Tuple[float, float]:
    if arg is None:
        return _auto_range_pct(values)
    s = str(arg).strip().lower()
    if s in ("", "auto"):
        return _auto_range_pct(values)
    parts = [p for p in re.split(r"[,:\\s]+", s) if p]
    if len(parts) != 2:
        raise ValueError(f"Invalid range: {arg}")
    a, b = float(parts[0]), float(parts[1])
    if a > b:
        a, b = b, a
    return (a, b)


def build_html(df_all: pd.DataFrame, group_tables: Dict[str, Dict[str, pd.DataFrame]],
               out_html: str, title: str, max_table: int,
               plot_x_range: Optional[str], plot_y_range: Optional[str],
               sample_cols: List[str], group_display: Dict[str, str],
               active_label: str, inactive_label: str,
               table_page_size: int,
               final_df: Optional[pd.DataFrame] = None,
               final_title: Optional[str] = None) -> None:
    if not _HAS_BOKEH:
        print("[WARN] bokeh is not available; skipping interactive HTML output.")
        return

    items = [Div(text=_style_block())]
    summary = group_tables.get("summary", {})
    group_defs = summary.get("group_defs", {})
    a_def = group_defs.get("Active-specific", "")
    i_def = group_defs.get("Inactive-specific", "")
    b_def = group_defs.get("Both-specific", "")
    disp_active = group_display.get("Active-specific", "Active-specific")
    disp_inactive = group_display.get("Inactive-specific", "Inactive-specific")
    disp_both = group_display.get("Both-specific", "Both-specific")
    cards = [
        ("Candidates", str(summary.get("n_candidates", 0))),
        (disp_active, str(summary.get("n_active", 0))),
        (disp_inactive, str(summary.get("n_inactive", 0))),
        (disp_both, str(summary.get("n_both", 0))),
        ("Unclassified", str(summary.get("n_other", 0))),
    ]
    summary_lines = [
        f"<h2>{title}</h2>",
        _summary_cards(cards),
        "<table class='summary-table'>"
        "<tr><th>Group</th><th>All hits</th><th>Diverse hits</th></tr>"
        f"<tr><td>{disp_active}</td><td>{summary.get('n_active', 0)}</td><td>{summary.get('n_active_diverse', 0)}</td></tr>"
        f"<tr><td>{disp_inactive}</td><td>{summary.get('n_inactive', 0)}</td><td>{summary.get('n_inactive_diverse', 0)}</td></tr>"
        f"<tr><td>{disp_both}</td><td>{summary.get('n_both', 0)}</td><td>{summary.get('n_both_diverse', 0)}</td></tr>"
        "</table>",
        "<ul>",
        f"<li>Candidates: <b>{summary.get('n_candidates', 0)}</b></li>",
        f"<li>{disp_active}: <b>{summary.get('n_active', 0)}</b>{a_def}</li>",
        f"<li>{disp_inactive}: <b>{summary.get('n_inactive', 0)}</b>{i_def}</li>",
        f"<li>{disp_both}: <b>{summary.get('n_both', 0)}</b>{b_def}</li>",
        f"<li>Unclassified: <b>{summary.get('n_other', 0)}</b></li>",
        f"<li>Cluster reps ({disp_active}/{disp_inactive}/{disp_both}): <b>{summary.get('n_diverse', 0)}</b></li>",
        "<li>Cutoffs use rank_pct percentiles per run (e.g., 99 = top 1% within that run; higher is better).</li>",
        "</ul>",
    ]
    items.append(Div(text="\n".join(summary_lines)))

    df_plot = df_all.copy()
    df_plot["group_code"] = df_plot.get("group_code", df_plot.get("group", "Other"))
    df_plot["active_rank_pct"] = pd.to_numeric(df_plot.get("active_rank_pct", 0), errors="coerce").fillna(0.0)
    df_plot["inactive_rank_pct"] = pd.to_numeric(df_plot.get("inactive_rank_pct", 0), errors="coerce").fillna(0.0)
    df_plot["both_rank_pct"] = pd.to_numeric(df_plot.get("both_rank_pct", 0), errors="coerce").fillna(0.0)
    df_plot["selectivity_score"] = pd.to_numeric(df_plot.get("selectivity_score", 0), errors="coerce").fillna(0.0)

    colors = {
        "Active-specific": "#1b9e77",
        "Inactive-specific": "#d95f02",
        "Both-specific": "#7570b3",
        "Other": "#bdbdbd",
    }
    df_plot["color"] = df_plot["group_code"].map(colors).fillna("#bdbdbd")
    cds = ColumnDataSource(df_plot)

    hover = HoverTool(tooltips=[
        ("group", "@group"),
        ("group_rank", "@group_rank"),
        ("specificity", "@specificity_score{0.000}"),
        ("selectivity", "@selectivity_score{0.000}"),
        (f"{active_label}_pct", "@active_rank_pct{0.0}"),
        (f"{inactive_label}_pct", "@inactive_rank_pct{0.0}"),
        ("both_pct", "@both_rank_pct{0.0}"),
        ("cluster", "@cluster_id"),
        ("ID", "@ID_x"),
        ("BB1", "@BB1_x"),
        ("BB2", "@BB2_x"),
        ("BB3", "@BB3_x"),
        ("BB4", "@BB4_x"),
    ])

    x_min, x_max = _parse_range(plot_x_range, df_plot["active_rank_pct"])
    y_min, y_max = _parse_range(plot_y_range, df_plot["inactive_rank_pct"])
    p = figure(
        title=f"{active_label} vs {inactive_label} rank_pct (lower {inactive_label} is better)",
        x_axis_label=f"{active_label}_rank_pct",
        y_axis_label=f"{inactive_label}_rank_pct",
        height=380,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        x_range=(x_min, x_max),
        y_range=(y_min, y_max),
    )
    p.add_tools(hover)
    p.scatter(x="active_rank_pct", y="inactive_rank_pct", size=7, color="color", alpha=0.7, source=cds)
    items.append(p)

    summary_panel = _make_panel(column(*items, sizing_mode="stretch_width"), "Specific/Common Summary")

    tabs = [summary_panel]
    for group in ["Active-specific", "Inactive-specific", "Both-specific"]:
        group_info = group_tables.get(group, {})
        df_all_tier = group_info.get("all", pd.DataFrame())
        df_div = group_info.get("diverse", pd.DataFrame())
        display_name = group_display.get(group, group)
        panel_items = [Div(text=f"<h3>{display_name}</h3>")]
        div_table, div_struct = _make_table(
            df_div, with_images=True, max_rows=max_table, with_structure=True,
            sample_cols=sample_cols, page_size=table_page_size,
            download_name=f"{group.lower().replace(' ', '_')}_diverse.csv"
        )
        if div_table is not None:
            panel_items.append(Div(text=f"<p>Diverse hits (cluster reps): {len(df_div)}</p>"))
            if div_struct is not None:
                panel_items.append(div_struct)
            panel_items.append(div_table)
        all_table, all_struct = _make_table(
            df_all_tier, with_images=True, max_rows=max_table, with_structure=True,
            sample_cols=sample_cols, page_size=table_page_size,
            download_name=f"{group.lower().replace(' ', '_')}_all.csv"
        )
        if all_table is not None:
            panel_items.append(Div(text=f"<p>All hits (top {min(len(df_all_tier), max_table)} shown)</p>"))
            if all_struct is not None:
                panel_items.append(all_struct)
            panel_items.append(all_table)
        tabs.append(_make_panel(column(*panel_items, sizing_mode="stretch_width"), display_name))

    if final_df is not None and not final_df.empty:
        display_final = final_title or "Final hits"
        final_sample_cols = list(sample_cols)
        panel_items = [Div(text=f"<h3>{display_final}</h3>")]
        final_table, final_struct = _make_table(
            final_df, with_images=True, max_rows=max_table, with_structure=True,
            sample_cols=final_sample_cols, page_size=table_page_size,
            download_name="final_hits.csv"
        )
        if final_table is not None:
            panel_items.append(Div(text=f"<p>Final hits: {len(final_df)}</p>"))
            if final_struct is not None:
                panel_items.append(final_struct)
            panel_items.append(final_table)
        tabs.append(_make_panel(column(*panel_items, sizing_mode="stretch_width"), display_final))

    output_file(out_html, title=title)
    save(Tabs(tabs=tabs))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--active-run", required=True, help="Active run root or 05_hybrid_annot.tsv path")
    ap.add_argument("--inactive-run", required=True, help="Inactive run root or 05_hybrid_annot.tsv path")
    ap.add_argument("--both-run", default=None, help="Both run root or 05_hybrid_annot.tsv path")
    ap.add_argument("--labels", nargs="*", default=None,
                    help="Optional labels for active/inactive/(both) runs (same order)")
    ap.add_argument("--preset", default=None, help="Optional preset under 03_normalized")
    ap.add_argument("--score-col", default=None, help="Override score column")
    ap.add_argument("--top-n", type=int, default=1000, help="Top N from active run to consider")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    ap.add_argument("--out-prefix", default=None, help="Output prefix (default: tier_report)")
    ap.add_argument("--no-html", action="store_true", help="Skip HTML output")
    ap.add_argument("--max-table", type=int, default=200, help="Max rows per HTML table")
    ap.add_argument("--table-page-size", type=int, default=50, help="Rows per HTML table page")
    ap.add_argument("--plot-x-range", default="auto",
                    help="Summary plot X range 'min,max' or 'auto' (default: auto)")
    ap.add_argument("--plot-y-range", default="auto",
                    help="Summary plot Y range 'min,max' or 'auto' (default: auto)")
    ap.add_argument("--sample-cols-mode", choices=["coalesced", "prefixed", "both"], default="prefixed",
                    help="Sample columns in HTML tables: coalesced, prefixed, or both (default: prefixed)")
    ap.add_argument("--rank-by", choices=["enrichment", "specificity"], default="enrichment",
                    help="Group ranking basis (default: enrichment)")
    ap.add_argument("--enrich-agg", choices=["median", "mean", "max"], default="median",
                    help="Aggregate per-group enrichment across samples (default: median)")
    ap.add_argument("--common-enrich", choices=["min", "mean", "gmean"], default="min",
                    help="Common-group enrichment from active/inactive (default: min)")

    ap.add_argument("--active-spec-min", type=float, default=99.0,
                    help="Active-specific: min active rank_pct")
    ap.add_argument("--active-spec-max-inactive", type=float, default=50.0,
                    help="Active-specific: max inactive rank_pct")
    ap.add_argument("--active-spec-min-both", type=float, default=90.0,
                    help="Active-specific: min both rank_pct (ignored when using active/inactive-only)")
    ap.add_argument("--inactive-spec-min", type=float, default=99.0,
                    help="Inactive-specific: min inactive rank_pct")
    ap.add_argument("--inactive-spec-max-active", type=float, default=50.0,
                    help="Inactive-specific: max active rank_pct")
    ap.add_argument("--inactive-spec-min-both", type=float, default=90.0,
                    help="Inactive-specific: min both rank_pct (ignored when using active/inactive-only)")
    ap.add_argument("--both-spec-min", type=float, default=99.0,
                    help="Both-specific: min both rank_pct")
    ap.add_argument("--both-spec-min-active", type=float, default=90.0,
                    help="Both-specific: min active rank_pct")
    ap.add_argument("--both-spec-min-inactive", type=float, default=90.0,
                    help="Both-specific: min inactive rank_pct")
    ap.add_argument("--both-weight", type=float, default=0.3,
                    help="Weight for both_rank_pct in specificity scores (ignored when using active/inactive-only)")

    ap.add_argument("--neg-samples", nargs="+", default=None,
                    help="Negative control sample names (base names, e.g. K_R2C5 K_R3C5)")
    ap.add_argument("--neg-max-pct", type=float, default=None,
                    help="Exclude candidates if any NEG sample CPM is >= this percentile (e.g. 99)")
    ap.add_argument("--final-active-n", type=int, default=0,
                    help="Final hits: number of active-specific compounds")
    ap.add_argument("--final-inactive-n", type=int, default=0,
                    help="Final hits: number of inactive-specific compounds")
    ap.add_argument("--final-common-n", type=int, default=0,
                    help="Final hits: number of common (both-specific) compounds")
    ap.add_argument("--final-hits", default=None,
                    help="Optional final hits TSV to embed as a tab in the HTML report")

    ap.add_argument("--cluster", type=int, choices=[0, 1], default=1)
    ap.add_argument("--cluster-mode", choices=["bbavg", "compound_or"], default="bbavg",
                    help="Clustering similarity: bbavg (per-BB avg) or compound_or (OR fingerprint)")
    ap.add_argument("--cluster-sim", type=float, default=0.8)
    ap.add_argument("--cluster-radius", type=int, default=2)
    ap.add_argument("--cluster-nbits", type=int, default=2048)
    ap.add_argument("--cluster-rep", choices=["score", "medoid"], default="score",
                    help="Representative selection for diverse lists (default: score)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    if args.out_prefix:
        if os.path.isabs(args.out_prefix) or os.path.dirname(args.out_prefix):
            prefix = args.out_prefix
        else:
            prefix = os.path.join(args.out_dir, args.out_prefix)
    else:
        prefix = os.path.join(args.out_dir, "tier_report")

    run_paths = [args.active_run, args.inactive_run]
    if args.both_run:
        run_paths.append(args.both_run)
    labels = args.labels or []
    if labels and len(labels) != len(run_paths):
        raise SystemExit("[ERROR] --labels must match number of runs (active/inactive/(both))")
    if not labels:
        labels = [os.path.basename(os.path.normpath(p)) or f"run{i+1}" for i, p in enumerate(run_paths)]
    labels = [_sanitize_label(l) for l in labels]
    active_label = labels[0] if len(labels) > 0 else "sampleA"
    inactive_label = labels[1] if len(labels) > 1 else "sampleB"
    both_label = labels[2] if len(labels) > 2 else "both"
    group_display = {
        "Active-specific": f"{active_label}-specific",
        "Inactive-specific": f"{inactive_label}-specific",
        "Both-specific": f"Common ({active_label}&{inactive_label})",
        "Other": "Other",
    }
    rep_col = "cluster_medoid" if args.cluster_rep == "medoid" else "cluster_rep"
    score_col = _infer_score_col(run_paths, args.preset, args.score_col)

    active_path = resolve_hybrid_path(args.active_run, args.preset)
    inactive_path = resolve_hybrid_path(args.inactive_run, args.preset)
    both_path = resolve_hybrid_path(args.both_run, args.preset) if args.both_run else None

    def _load(path: str, want: List[str]) -> pd.DataFrame:
        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        usecols = [c for c in want if c in header]
        df = pd.read_csv(path, sep="\t", usecols=usecols)
        return df

    header_active = pd.read_csv(active_path, sep="\t", nrows=0).columns.tolist()
    header_inactive = pd.read_csv(inactive_path, sep="\t", nrows=0).columns.tolist()
    header_both = pd.read_csv(both_path, sep="\t", nrows=0).columns.tolist() if both_path else []

    active_sample_cols = _sample_cols_from_header(header_active)
    inactive_sample_cols = _sample_cols_from_header(header_inactive)
    both_sample_cols = _sample_cols_from_header(header_both) if header_both else []

    active_prefixed, active_rename = _prefix_sample_cols(active_sample_cols, f"{active_label}_")
    inactive_prefixed, inactive_rename = _prefix_sample_cols(inactive_sample_cols, f"{inactive_label}_")
    both_prefixed, both_rename = _prefix_sample_cols(both_sample_cols, f"{both_label}_")

    prefixed_sample_cols = []
    for col in active_prefixed + inactive_prefixed + both_prefixed:
        if col not in prefixed_sample_cols:
            prefixed_sample_cols.append(col)

    want_base = [
        "LIB_ID_x", "ID_x", "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x", "cycles",
        "bb1_smiles", "bb2_smiles", "bb3_smiles", "bb4_smiles", "BB_SMILES_CONCAT",
        "HitScore_GLM", "HitScore_RS", "HitScore_pct", "SynthonScore",
        "GLM_hit", "RS_pass", "Consensus_hit", "NEG_hard_fail", "NEG_center_fail",
        "pass_filters", "fail_reasons",
        "mean_R1_norm", "mean_R2_norm", "DEL2_norm",
        "LFC_R1_vs_DEL2", "LFC_R2_vs_DEL2",
        "LFC_NEG_R1_vs_DEL2", "LFC_NEG_R2_vs_DEL2",
        "q_DEL2", "q_BEAD", "q_BEAD_R2", "q_BoostPaired",
    ]
    want_active = list(want_base)
    for col in active_sample_cols:
        if col not in want_active:
            want_active.append(col)
    if score_col not in want_active:
        want_active.append(score_col)

    want_other = list(want_base)
    for col in inactive_sample_cols:
        if col not in want_other:
            want_other.append(col)
    for col in both_sample_cols:
        if col not in want_other:
            want_other.append(col)
    if score_col not in want_other:
        want_other.append(score_col)

    df_active = _load(active_path, want_active)
    df_inactive = _load(inactive_path, want_other)
    df_both = _load(both_path, want_other) if both_path else None

    for df in [df_active, df_inactive] + ([df_both] if df_both is not None else []):
        if df is None:
            continue
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(float("-inf"))
        if "HitScore_pct" in df.columns:
            df["HitScore_pct"] = pd.to_numeric(df["HitScore_pct"], errors="coerce")
        df["compound_key"] = _make_compound_key(df)

    maps_active = _rank_maps(df_active, score_col)
    maps_inactive = _rank_maps(df_inactive, score_col)
    maps_both = _rank_maps(df_both, score_col) if df_both is not None else None

    df_active_best = df_active.sort_values(score_col, ascending=False).drop_duplicates("compound_key")
    df_inactive_best = df_inactive.sort_values(score_col, ascending=False).drop_duplicates("compound_key")
    df_both_best = df_both.sort_values(score_col, ascending=False).drop_duplicates("compound_key") if df_both is not None else None

    def _sample_block(df_best: pd.DataFrame, cols: List[str], rename: Dict[str, str]) -> Optional[pd.DataFrame]:
        if df_best is None or df_best.empty or not cols:
            return None
        use = [c for c in cols if c in df_best.columns]
        if not use:
            return None
        block = df_best[["compound_key"] + use].copy()
        block = block.rename(columns=rename)
        return block

    active_samples = _sample_block(df_active_best, active_sample_cols, active_rename)
    inactive_samples = _sample_block(df_inactive_best, inactive_sample_cols, inactive_rename)
    both_samples = _sample_block(df_both_best, both_sample_cols, both_rename) if df_both_best is not None else None

    top_n = max(1, int(args.top_n))
    top_frames = []
    for df_best in (df_active_best, df_inactive_best, df_both_best):
        if df_best is None or df_best.empty:
            continue
        df_top = df_best.head(top_n).copy()
        df_top = _apply_recommend_filter(df_top)
        if not df_top.empty:
            top_frames.append(df_top)
    if not top_frames:
        raise SystemExit("[ERROR] No candidates after filtering; check top-n or filters.")
    df_base = pd.concat(top_frames, ignore_index=True).drop_duplicates("compound_key").copy()
    for block in (active_samples, inactive_samples, both_samples):
        if block is not None:
            df_base = df_base.merge(block, on="compound_key", how="left")
    for col in prefixed_sample_cols:
        if col in df_base.columns:
            df_base[col] = pd.to_numeric(df_base[col], errors="coerce")
    base_sample_cols = []
    for col in active_sample_cols + inactive_sample_cols + both_sample_cols:
        if col not in base_sample_cols:
            base_sample_cols.append(col)
    prefixes = [f"{active_label}_", f"{inactive_label}_"]
    if both_label and both_label not in (active_label, inactive_label):
        prefixes.append(f"{both_label}_")
    df_base = _coalesce_unprefixed_samples(df_base, base_sample_cols, prefixes)
    display_sample_cols = []
    if args.sample_cols_mode in ("coalesced", "both"):
        for col in base_sample_cols:
            if col in df_base.columns and col not in display_sample_cols:
                display_sample_cols.append(col)
    if args.sample_cols_mode in ("prefixed", "both"):
        for col in prefixed_sample_cols:
            if col in df_base.columns and col not in display_sample_cols:
                display_sample_cols.append(col)
    if display_sample_cols:
        display_sample_cols = _sort_sample_cols_by_category(display_sample_cols, prefixes)
        display_sample_cols = _dedupe_identical_sample_cols(df_base, display_sample_cols, prefixes)

    def _split_tokens(items):
        out = []
        for item in items or []:
            for tok in str(item).replace(",", " ").split():
                tok = tok.strip()
                if tok:
                    out.append(tok)
        return out

    neg_samples = _split_tokens(args.neg_samples)

    def _neg_filter(df: pd.DataFrame) -> pd.DataFrame:
        if args.neg_max_pct is None:
            return df
        if not neg_samples:
            return df
        thresholds: Dict[str, float] = {}
        thresholds.update(_neg_thresholds_for_run(
            df_active_best, active_sample_cols, f"{active_label}_", neg_samples, float(args.neg_max_pct)
        ))
        thresholds.update(_neg_thresholds_for_run(
            df_inactive_best, inactive_sample_cols, f"{inactive_label}_", neg_samples, float(args.neg_max_pct)
        ))
        if df_both_best is not None and both_sample_cols:
            thresholds.update(_neg_thresholds_for_run(
                df_both_best, both_sample_cols, f"{both_label}_", neg_samples, float(args.neg_max_pct)
            ))
        if not thresholds:
            print("[WARN] NEG filter requested but no NEG sample columns found; skipping.")
            return df
        mask = pd.Series(True, index=df.index)
        for col, thr in thresholds.items():
            if col not in df.columns:
                continue
            mask &= pd.to_numeric(df[col], errors="coerce").fillna(0.0) < thr
        kept = int(mask.sum())
        print(f"[INFO] NEG filter: pct={args.neg_max_pct}, cols={len(thresholds)}, kept={kept}/{len(df)}")
        return df[mask].copy()

    df_base = _neg_filter(df_base)

    df_base = _strip_bb_suffix_cols(df_base, ["BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x"])
    df_base = df_base.sort_values(score_col, ascending=False)
    df_base["rank"] = np.arange(1, len(df_base) + 1, dtype=int)

    df_base["active_score"] = df_base["compound_key"].map(maps_active["score_map"])
    df_base["active_rank"] = df_base["compound_key"].map(maps_active["rank_map"])
    df_base["active_rank_pct"] = df_base["compound_key"].map(maps_active["rank_pct_map"]).fillna(0.0)
    df_base["active_score_z"] = df_base["compound_key"].map(maps_active["score_z_map"]).fillna(0.0)
    df_base["active_hit_pct"] = df_base["compound_key"].map(maps_active["hitpct_map"]).fillna(0.0)

    df_base["inactive_score"] = df_base["compound_key"].map(maps_inactive["score_map"])
    df_base["inactive_rank_pct"] = df_base["compound_key"].map(maps_inactive["rank_pct_map"]).fillna(0.0)

    if maps_both:
        df_base["both_score"] = df_base["compound_key"].map(maps_both["score_map"])
        df_base["both_rank_pct"] = df_base["compound_key"].map(maps_both["rank_pct_map"]).fillna(0.0)
    else:
        df_base["both_score"] = np.nan
        df_base["both_rank_pct"] = np.nan

    df_base["selectivity_score"] = df_base["active_rank_pct"] - df_base["inactive_rank_pct"]
    df_base["inactive_selectivity_score"] = df_base["inactive_rank_pct"] - df_base["active_rank_pct"]
    df_base["both_specific_score"] = (df_base["active_rank_pct"] + df_base["inactive_rank_pct"]) / 2.0

    exclude_bases = {s.upper() for s in neg_samples} | {"DEL234"}
    active_enrich_cols = _pick_enrich_cols(
        _filter_enrich_cols([c for c in active_prefixed if c in df_base.columns], f"{active_label}_", exclude_bases)
    )
    inactive_enrich_cols = _pick_enrich_cols(
        _filter_enrich_cols([c for c in inactive_prefixed if c in df_base.columns], f"{inactive_label}_", exclude_bases)
    )
    df_base["active_enrich"] = _row_agg(df_base, active_enrich_cols, args.enrich_agg)
    df_base["inactive_enrich"] = _row_agg(df_base, inactive_enrich_cols, args.enrich_agg)
    if args.common_enrich == "mean":
        df_base["both_enrich"] = (df_base["active_enrich"] + df_base["inactive_enrich"]) / 2.0
    elif args.common_enrich == "gmean":
        df_base["both_enrich"] = np.sqrt(df_base["active_enrich"] * df_base["inactive_enrich"])
    else:
        df_base["both_enrich"] = np.fmin(df_base["active_enrich"], df_base["inactive_enrich"])

    def _assign_group(row: pd.Series) -> str:
        a_pass = (
            row["active_rank_pct"] >= float(args.active_spec_min)
            and row["inactive_rank_pct"] <= float(args.active_spec_max_inactive)
        )
        if a_pass:
            return "Active-specific"
        i_pass = (
            row["inactive_rank_pct"] >= float(args.inactive_spec_min)
            and row["active_rank_pct"] <= float(args.inactive_spec_max_active)
        )
        if i_pass:
            return "Inactive-specific"
        min_active = max(float(args.both_spec_min), float(args.both_spec_min_active))
        min_inactive = max(float(args.both_spec_min), float(args.both_spec_min_inactive))
        b_pass = (
            row["active_rank_pct"] >= min_active
            and row["inactive_rank_pct"] >= min_inactive
        )
        if b_pass:
            return "Both-specific"
        return "Other"

    df_base["group_code"] = df_base.apply(_assign_group, axis=1)
    df_base["specificity_score"] = np.select(
        [
            df_base["group_code"] == "Active-specific",
            df_base["group_code"] == "Inactive-specific",
            df_base["group_code"] == "Both-specific",
        ],
        [
            df_base["selectivity_score"],
            df_base["inactive_selectivity_score"],
            df_base["both_specific_score"],
        ],
        default=df_base["selectivity_score"],
    )

    if int(args.cluster) == 1:
        df_base, cluster_df = _cluster_candidates(
            df_base,
            rep_score_col="specificity_score",
            sim_cutoff=float(args.cluster_sim),
            radius=int(args.cluster_radius),
            nbits=int(args.cluster_nbits),
            mode=str(args.cluster_mode),
        )
    else:
        cluster_df = None

    group_order = ["Active-specific", "Inactive-specific", "Both-specific", "Other"]
    df_base["group_code"] = pd.Categorical(df_base["group_code"], categories=group_order, ordered=True)
    group_active_enrich = df_base["active_enrich"].fillna(float("-inf"))
    group_inactive_enrich = df_base["inactive_enrich"].fillna(float("-inf"))
    group_both_enrich = df_base["both_enrich"].fillna(float("-inf"))
    df_base["group_rank_score"] = np.select(
        [
            df_base["group_code"] == "Active-specific",
            df_base["group_code"] == "Inactive-specific",
            df_base["group_code"] == "Both-specific",
        ],
        [
            group_active_enrich if args.rank_by == "enrichment" else df_base["specificity_score"],
            group_inactive_enrich if args.rank_by == "enrichment" else df_base["specificity_score"],
            group_both_enrich if args.rank_by == "enrichment" else df_base["specificity_score"],
        ],
        default=df_base["specificity_score"] if args.rank_by == "specificity" else group_active_enrich,
    )
    df_base["group_tiebreak_score"] = np.select(
        [
            df_base["group_code"] == "Active-specific",
            df_base["group_code"] == "Inactive-specific",
            df_base["group_code"] == "Both-specific",
        ],
        [
            df_base["active_score"],
            df_base["inactive_score"],
            df_base["both_score"].fillna(df_base["both_specific_score"]),
        ],
        default=df_base["active_score"],
    )
    if args.rank_by == "enrichment":
        sort_cols = ["group_code", "group_rank_score", "specificity_score", "group_tiebreak_score"]
        sort_asc = [True, False, False, False]
    else:
        sort_cols = ["group_code", "group_rank_score", "group_tiebreak_score"]
        sort_asc = [True, False, False]
    df_base = df_base.sort_values(
        sort_cols,
        ascending=sort_asc,
    )
    df_base = df_base.drop(columns=["group_tiebreak_score"])
    df_base["group_rank"] = df_base.groupby("group_code", observed=False).cumcount() + 1
    group_code_str = df_base["group_code"].astype(str)
    df_base["group"] = group_code_str.map(group_display).fillna(group_code_str)

    out_all = f"{prefix}_all_candidates.tsv"
    df_base.to_csv(out_all, sep="\t", index=False)

    groups = {}
    diverse_total = 0
    group_map = {
        "Active-specific": "active_specific",
        "Inactive-specific": "inactive_specific",
        "Both-specific": "both_specific",
    }
    for group, tag in group_map.items():
        df_group = df_base[df_base["group_code"] == group].copy()
        df_div = df_group.copy()
        if rep_col in df_div.columns:
            df_div = df_div[_to_int_series(df_div[rep_col]) == 1].copy()
        df_group_path = f"{prefix}_{tag}.tsv"
        df_div_path = f"{prefix}_{tag}_diverse.tsv"
        df_group.to_csv(df_group_path, sep="\t", index=False)
        df_div.to_csv(df_div_path, sep="\t", index=False)
        diverse_total += len(df_div)
        groups[group] = {"all": df_group, "diverse": df_div, "path": df_group_path, "diverse_path": df_div_path}

    df_other = df_base[df_base["group_code"] == "Other"].copy()
    other_path = None
    if not df_other.empty:
        other_path = f"{prefix}_other.tsv"
        df_other.to_csv(other_path, sep="\t", index=False)

    if cluster_df is not None:
        cluster_path = f"{prefix}_clusters.tsv"
        cluster_df.to_csv(cluster_path, sep="\t", index=False)

    min_active = max(float(args.both_spec_min), float(args.both_spec_min_active))
    min_inactive = max(float(args.both_spec_min), float(args.both_spec_min_inactive))
    group_defs = {
        "Active-specific": f" ({active_label}≥{args.active_spec_min}, {inactive_label}≤{args.active_spec_max_inactive})",
        "Inactive-specific": f" ({inactive_label}≥{args.inactive_spec_min}, {active_label}≤{args.inactive_spec_max_active})",
        "Both-specific": f" ({active_label}≥{min_active}, {inactive_label}≥{min_inactive})",
    }
    summary = {
        "n_candidates": len(df_base),
        "n_active": len(groups.get("Active-specific", {}).get("all", [])),
        "n_inactive": len(groups.get("Inactive-specific", {}).get("all", [])),
        "n_both": len(groups.get("Both-specific", {}).get("all", [])),
        "n_other": len(df_other),
        "n_diverse": diverse_total,
        "n_active_diverse": len(groups.get("Active-specific", {}).get("diverse", [])),
        "n_inactive_diverse": len(groups.get("Inactive-specific", {}).get("diverse", [])),
        "n_both_diverse": len(groups.get("Both-specific", {}).get("diverse", [])),
        "group_defs": group_defs,
    }
    group_tables = {"summary": summary, **groups}

    final_df = None
    final_title = None
    final_counts = {
        "Active-specific": max(0, int(args.final_active_n)),
        "Inactive-specific": max(0, int(args.final_inactive_n)),
        "Both-specific": max(0, int(args.final_common_n)),
    }
    final_total = sum(final_counts.values())
    if final_total > 0:
        final_parts = []
        for group in ["Active-specific", "Inactive-specific", "Both-specific"]:
            n = final_counts.get(group, 0)
            if n <= 0:
                continue
            df_group = groups.get(group, {}).get("all", pd.DataFrame())
            df_div = groups.get(group, {}).get("diverse", pd.DataFrame())
            df_pick = _pick_final_from_group(df_div, df_group, n)
            if df_pick.empty:
                continue
            df_pick = df_pick.copy()
            df_pick["final_group_code"] = group
            df_pick["final_group"] = df_pick["final_group_code"].map(group_display).fillna(df_pick["final_group_code"])
            df_pick["final_group_rank"] = np.arange(1, len(df_pick) + 1, dtype=int)
            final_parts.append(df_pick)
        if final_parts:
            final_df = pd.concat(final_parts, ignore_index=True)
            final_title = (
                f"Final hits ({active_label}={final_counts['Active-specific']}, "
                f"{inactive_label}={final_counts['Inactive-specific']}, "
                f"Common={final_counts['Both-specific']})"
            )
            final_prefix = os.path.join(args.out_dir, "final_hits")
            final_df.to_csv(f"{final_prefix}.tsv", sep="\t", index=False)
            try:
                final_df.to_excel(f"{final_prefix}.xlsx", index=False)
            except Exception as exc:
                print(f"[WARN] Failed to write final hits Excel ({exc})")
            print(f"[INFO] final_hits: {final_prefix}.tsv")
        else:
            print("[WARN] final_hits requested but no rows matched group cutoffs.")
    else:
        final_hits_path = args.final_hits
        if final_hits_path is None:
            for candidate in [
                os.path.join(args.out_dir, "final_hits.tsv"),
            ]:
                if os.path.exists(candidate):
                    final_hits_path = candidate
                    break
        if final_hits_path:
            if os.path.exists(final_hits_path):
                try:
                    final_df = pd.read_csv(final_hits_path, sep="\t")
                    final_title = "Final hits"
                except Exception as exc:
                    print(f"[WARN] Failed to read final hits: {final_hits_path} ({exc})")
            else:
                print(f"[WARN] --final-hits not found: {final_hits_path}")

    html_path = f"{prefix}_interactive.html"
    if not args.no_html:
        build_html(df_base, group_tables, html_path, title="Selectivity & Diversity Report",
                   max_table=int(args.max_table),
                   plot_x_range=args.plot_x_range,
                   plot_y_range=args.plot_y_range,
                   sample_cols=display_sample_cols,
                   group_display=group_display,
                   active_label=active_label,
                   inactive_label=inactive_label,
                   table_page_size=int(args.table_page_size),
                   final_df=final_df,
                   final_title=final_title)

    print(f"[INFO] score_col={score_col}")
    print(
        "[INFO] candidates="
        f"{len(df_base)} active={summary['n_active']} inactive={summary['n_inactive']} "
        f"both={summary['n_both']} other={summary['n_other']}"
    )
    base_outputs = [
        out_all,
        f"{prefix}_active_specific.tsv",
        f"{prefix}_inactive_specific.tsv",
        f"{prefix}_both_specific.tsv",
    ]
    if other_path:
        base_outputs.append(other_path)
    print(f"[INFO] outputs: {', '.join(base_outputs)}")
    print(
        "[INFO] outputs: "
        f"{prefix}_active_specific_diverse.tsv, "
        f"{prefix}_inactive_specific_diverse.tsv, "
        f"{prefix}_both_specific_diverse.tsv"
    )
    if not args.no_html:
        print(f"[INFO] html: {html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
