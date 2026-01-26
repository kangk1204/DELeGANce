#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified report for up to 3 runs:
- Per-run top hit summaries
- Cross-run comparison + Venn cells
- Optional specificity report (active/inactive/both roles)
"""

import argparse
import base64
import os
import re
from io import BytesIO
from itertools import combinations
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


def _sample_cols_from_header(header: List[str]) -> List[str]:
    cpm_cols = [c for c in header if c.endswith("_CPM") and c[:-4] in header]
    raw_cols = [c[:-4] for c in cpm_cols]
    return raw_cols + cpm_cols


def _prefix_sample_cols(cols: List[str], prefix: str) -> Tuple[List[str], Dict[str, str]]:
    renamed = {c: f"{prefix}{c}" for c in cols}
    return [renamed[c] for c in cols], renamed


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


def _strip_bb_suffix_cols(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = out[col].fillna("NA").astype(str).map(_strip_lib_suffix)
    return out


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


def _compound_key_from_row(row: pd.Series) -> str:
    vals = []
    for col in ["BB1_x", "BB2_x", "BB3_x", "BB4_x"]:
        if col in row:
            vals.append(_strip_lib_suffix(row.get(col, "NA")))
    if not vals:
        return "NA"
    return "|".join(vals)


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


def _cluster_top_hits(df_top: pd.DataFrame, score_col: str, sim_cutoff: float,
                      radius: int, nbits: int, mode: str) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    if not _HAS_RDKIT:
        print(f"[WARN] RDKit not available ({_RDKIT_IMPORT_ERR}); clustering skipped.")
        out = df_top.copy()
        out["cluster_id"] = "NA"
        out["cluster_size"] = 0
        out["cluster_rep"] = 0
        return out, None

    smi_cols = [c for c in ["bb1_smiles", "bb2_smiles", "bb3_smiles", "bb4_smiles"] if c in df_top.columns]
    if not smi_cols:
        print("[WARN] No bb*_smiles columns available; clustering skipped.")
        out = df_top.copy()
        out["cluster_id"] = "NA"
        out["cluster_size"] = 0
        out["cluster_rep"] = 0
        return out, None

    fps = []
    valid_pos = []
    for i, row in df_top.reset_index(drop=True).iterrows():
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
        out = df_top.copy()
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

    out = df_top.reset_index(drop=True).copy()
    out["cluster_id"] = "NA"
    out["cluster_size"] = 0
    out["cluster_rep"] = 0
    out["cluster_medoid"] = 0

    scores, rank_pct, cpm_mean, raw_sum = _tie_break_arrays(out, score_col)
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
            "rep_score": float(rep.get(score_col, np.nan)) if pd.notna(rep.get(score_col, np.nan)) else np.nan,
            "rep_BB1_x": rep.get("BB1_x", "NA"),
            "rep_BB2_x": rep.get("BB2_x", "NA"),
            "rep_BB3_x": rep.get("BB3_x", "NA"),
            "rep_BB4_x": rep.get("BB4_x", "NA"),
            "rep_CP_x": rep.get("CP_x", "NA"),
            "rep_compound_key": _compound_key_from_row(rep),
            "medoid_index": int(medoid_pos),
            "medoid_ID_x": medoid.get("ID_x", "NA"),
            "medoid_score": float(medoid.get(score_col, np.nan)) if pd.notna(medoid.get(score_col, np.nan)) else np.nan,
            "medoid_compound_key": _compound_key_from_row(medoid),
            "medoid_avg_sim": float(medoid_avg),
        })

    cluster_df = pd.DataFrame(cluster_rows)
    return out, cluster_df


def _build_bb_frequency_df(df_top: pd.DataFrame, bb_top_k: int) -> pd.DataFrame:
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
    return freq_df


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


def _table_width(columns: List[TableColumn], min_width: int = 900, max_width: int = 4200) -> int:
    total = sum([c.width or 0 for c in columns]) + 40
    if total <= 0:
        total = min_width
    return int(min(max(total, min_width), max_width))


def _make_table_controls(source: ColumnDataSource, full: ColumnDataSource,
                         filename: str, group_col: Optional[str]) -> row:
    total = 0
    if full.data:
        first_key = next(iter(full.data.keys()))
        total = len(full.data[first_key])
    search = TextInput(title="Search", value="", placeholder="compound_key / ID / BB / CP")
    count_div = Div(text=f"<span class='table-count'>Showing {total} of {total}</span>")
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
source.data = out;
source.change.emit();
count_div.text = `<span class='table-count'>Showing ${out[cols[0]].length} of ${n}</span>`;
"""
    args = dict(source=source, full=full, search=search, count_div=count_div, group=group_select, group_col=group_col or "")
    callback = CustomJS(args=args, code=filter_js)
    search.js_on_change("value", callback)
    if group_select is not None:
        group_select.js_on_change("value", callback)

    reset_js = """
const data = full.data;
const out = {};
for (const c in data) { out[c] = data[c].slice(); }
source.data = out;
source.change.emit();
search.value = '';
if (group !== null) { group.value = group.options; }
const n = out[Object.keys(out)[0]].length;
count_div.text = `<span class='table-count'>Showing ${n} of ${n}</span>`;
"""
    reset_btn.js_on_click(CustomJS(args=args, code=reset_js))

    download_js = """
function csvEscape(value) {
  if (value === null || value === undefined) return '';
  const s = String(value);
  if (s.includes('\"') || s.includes(',') || s.includes('\\n')) {
    return '\"' + s.replace(/\"/g, '\"\"') + '\"';
  }
  return s;
}
const data = source.data;
const cols = Object.keys(data);
const n = cols.length ? data[cols[0]].length : 0;
let csv = cols.join(',') + '\\n';
for (let i = 0; i < n; i++) {
  const row = cols.map(c => csvEscape(data[c][i]));
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
    download_btn.js_on_click(CustomJS(args=dict(source=source, filename=filename), code=download_js))

    controls = [search]
    if group_select is not None:
        controls.append(group_select)
    controls.extend([reset_btn, download_btn, count_div])
    return row(*controls, css_classes=["table-controls"])


def _frozen_columns_for(cols: List[str]) -> int:
    for key in ["compound_key", "ID_x", "LIB_ID_x", "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x"]:
        if key in cols:
            return int(cols.index(key)) + 1
    return 0


def _table_layout(df: pd.DataFrame, cols: List[str], filename: str,
                  height: int = 320, group_col: Optional[str] = None) -> object:
    if df is None or df.empty:
        return Div(text="<i>No rows to display</i>")
    table_df = df[cols].copy()
    source = ColumnDataSource(table_df)
    full_source = ColumnDataSource(table_df.copy())
    table_cols = _table_columns(table_df, cols)
    frozen_columns = _frozen_columns_for(cols)
    table = DataTable(
        source=source,
        columns=table_cols,
        height=height,
        width=_table_width(table_cols),
        index_position=None,
        frozen_columns=frozen_columns,
    )
    controls = _make_table_controls(source, full_source, filename, group_col)
    return column(controls, table, sizing_mode="stretch_width")

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


def _make_table(df: pd.DataFrame, with_images: bool, max_rows: int,
                with_structure: bool, sample_cols: List[str],
                download_name: str) -> Tuple[Optional[object], Optional[Div]]:
    if df is None or df.empty:
        return None, None
    if max_rows and len(df) > max_rows:
        df = df.head(max_rows).copy()
    if with_images or with_structure:
        df = _add_images(df)
    df = _add_group_badge(df)
    table_cols = []
    table_fields = [
        "group_badge", "group_rank", "specificity_score", "selectivity_score",
        "active_rank_pct", "inactive_rank_pct", "both_rank_pct",
        "rank", "cluster_id", "cluster_size", "cluster_rep", "cluster_medoid", "ID_x", "LIB_ID_x",
        "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x",
    ]
    if with_images:
        table_fields = [
            "group_badge", "group_rank", "specificity_score", "selectivity_score",
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
                table_cols.append(TableColumn(field=col, title=col, formatter=NumberFormatter(format="0"), width=_col_width(col, table_df[col])))
            elif col in sample_cols and not col.endswith("_CPM"):
                table_df[col] = pd.to_numeric(table_df[col], errors="coerce").fillna(0).astype(int)
                table_cols.append(TableColumn(field=col, title=col, formatter=NumberFormatter(format="0"), width=_col_width(col, table_df[col])))
            else:
                table_df[col] = table_df[col].round(4)
                table_cols.append(TableColumn(field=col, title=col, formatter=NumberFormatter(format="0.0000"), width=_col_width(col, table_df[col])))
        else:
            table_cols.append(TableColumn(field=col, title=col, width=_col_width(col, table_df[col])))

    if not table_cols:
        return None, None
    source = ColumnDataSource(table_df)
    full_source = ColumnDataSource(table_df.copy())
    row_height = 80 if with_images else 30
    shown_rows = len(table_df)
    if max_rows:
        shown_rows = min(shown_rows, max_rows)
    table_height = min(800, 40 + row_height * max(1, shown_rows))
    table_width = _table_width(table_cols, min_width=1200, max_width=5000)
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
    controls = _make_table_controls(source, full_source, download_name, group_col)
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


def _table_columns(df: pd.DataFrame, cols: List[str]) -> List[TableColumn]:
    out = []
    img_small = _img_formatter(70, 220)
    img_compound = _img_formatter(55, 220)
    badge_fmt = HTMLTemplateFormatter(template="<%= value %>")
    for c in cols:
        if c not in df.columns:
            continue
        if c == "group_badge":
            out.append(TableColumn(field=c, title="group", formatter=badge_fmt, width=_col_width(c, df[c])))
            continue
        if c.endswith("_img"):
            fmt = img_compound if c == "compound_img" else img_small
            out.append(TableColumn(field=c, title=c, formatter=fmt, width=_col_width(c, df[c])))
            continue
        width = _col_width(c, df[c])
        if df[c].dtype.kind in "if":
            fmt = "0.0000"
            if c.endswith("_rank_pct") or c == "rank_pct" or c.endswith("_HitScore_pct") or c == "HitScore_pct":
                fmt = "0.0"
            if c.endswith("_rank") or c == "rank":
                fmt = "0"
            out.append(TableColumn(field=c, title=c, formatter=NumberFormatter(format=fmt), width=width))
        else:
            out.append(TableColumn(field=c, title=c, width=width))
    return out


def _make_panel(child, title: str):
    if _HAS_TABPANEL:
        return _BokehTabPanel(child=child, title=title)
    return Panel(child=child, title=title)


def _pattern_order(labels: List[str]) -> List[str]:
    out = []
    for r in range(1, len(labels) + 1):
        for combo in combinations(labels, r):
            out.append("+".join(combo))
    return out


def _role_label_map(roles: List[str], labels: List[str]) -> Dict[str, str]:
    role_to_label: Dict[str, str] = {}
    for role, label in zip(roles, labels):
        if role and role not in role_to_label:
            role_to_label[role] = label
    return role_to_label


def build_specificity_panels(df_all: pd.DataFrame, group_tables: Dict[str, Dict[str, pd.DataFrame]],
                             title: str, max_table: int, plot_x_range: Optional[str],
                             plot_y_range: Optional[str], sample_cols: List[str],
                             group_display: Dict[str, str], active_label: str,
                             inactive_label: str) -> List[Panel]:
    if not _HAS_BOKEH:
        return []

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
            df_div, with_images=True, max_rows=max_table, with_structure=True, sample_cols=sample_cols,
            download_name=f"{group.lower().replace(' ', '_')}_diverse.csv"
        )
        if div_table is not None:
            panel_items.append(Div(text=f"<p>Diverse hits (cluster reps): {len(df_div)}</p>"))
            if div_struct is not None:
                panel_items.append(div_struct)
            panel_items.append(div_table)
        all_table, all_struct = _make_table(
            df_all_tier, with_images=True, max_rows=max_table, with_structure=True, sample_cols=sample_cols,
            download_name=f"{group.lower().replace(' ', '_')}_all.csv"
        )
        if all_table is not None:
            panel_items.append(Div(text=f"<p>All hits (top {min(len(df_all_tier), max_table)} shown)</p>"))
            if all_struct is not None:
                panel_items.append(all_struct)
            panel_items.append(all_table)
        tabs.append(_make_panel(column(*panel_items, sizing_mode="stretch_width"), display_name))

    return tabs


def compute_specificity_report(active_info: Dict, inactive_info: Dict, both_info: Optional[Dict],
                               score_col: str, args, out_dir: str,
                               role_labels: Optional[Dict[str, str]] = None) -> Tuple[List[Panel], Dict[str, str]]:
    df_active = active_info["df"].copy()
    df_inactive = inactive_info["df"].copy()
    df_both = both_info["df"].copy() if both_info else None

    role_labels = role_labels or {}
    active_label = _sanitize_label(role_labels.get("active", "sampleA"))
    inactive_label = _sanitize_label(role_labels.get("inactive", "sampleB"))
    both_label = _sanitize_label(role_labels.get("both", "both"))
    group_display = {
        "Active-specific": f"{active_label}-specific",
        "Inactive-specific": f"{inactive_label}-specific",
        "Both-specific": f"Common ({active_label}&{inactive_label})",
        "Other": "Other",
    }
    rep_col = "cluster_medoid" if getattr(args, "cluster_rep", "score") == "medoid" else "cluster_rep"

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

    top_n = int(args.spec_top_n) if args.spec_top_n else int(args.top_n)

    def _sample_block(df_best: pd.DataFrame, cols: List[str], rename: Dict[str, str]) -> Optional[pd.DataFrame]:
        if df_best is None or df_best.empty or not cols:
            return None
        use = [c for c in cols if c in df_best.columns]
        if not use:
            return None
        block = df_best[["compound_key"] + use].copy()
        block = block.rename(columns=rename)
        return block

    frames = []
    for df_best in [active_info["df_best"], inactive_info["df_best"], both_info["df_best"] if both_info else None]:
        if df_best is None or df_best.empty:
            continue
        df_top = df_best.head(top_n).copy()
        df_top = _apply_recommend_filter(df_top)
        if not df_top.empty:
            frames.append(df_top)
    if not frames:
        return [], {}

    df_base = pd.concat(frames, ignore_index=True).drop_duplicates("compound_key").copy()

    sample_cols = []
    active_prefixed, active_rename = _prefix_sample_cols(active_info["sample_cols"], f"{active_label}_")
    inactive_prefixed, inactive_rename = _prefix_sample_cols(inactive_info["sample_cols"], f"{inactive_label}_")
    both_prefixed, both_rename = _prefix_sample_cols(both_info["sample_cols"], f"{both_label}_") if both_info else ([], {})

    blocks = [
        _sample_block(active_info["df_best"], active_info["sample_cols"], active_rename),
        _sample_block(inactive_info["df_best"], inactive_info["sample_cols"], inactive_rename),
        _sample_block(both_info["df_best"], both_info["sample_cols"], both_rename) if both_info else None,
    ]
    for block in blocks:
        if block is not None:
            df_base = df_base.merge(block, on="compound_key", how="left")
    for col in active_prefixed + inactive_prefixed + both_prefixed:
        if col in df_base.columns:
            df_base[col] = pd.to_numeric(df_base[col], errors="coerce").fillna(0)
            sample_cols.append(col)
    prefixes = [f"{active_label}_", f"{inactive_label}_"]
    if both_label and both_label not in (active_label, inactive_label):
        prefixes.append(f"{both_label}_")
    if sample_cols:
        sample_cols = _sort_sample_cols_by_category(sample_cols, prefixes)
        sample_cols = _dedupe_identical_sample_cols(df_base, sample_cols, prefixes)

    df_base = _strip_bb_suffix_cols(df_base, ["BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x"])
    df_base = df_base.sort_values(score_col, ascending=False)
    df_base["rank"] = np.arange(1, len(df_base) + 1, dtype=int)

    df_base["active_score"] = df_base[score_col]
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
    df_base = df_base.sort_values(["group_code", "specificity_score", "active_score"], ascending=[True, False, False])
    df_base["group_rank"] = df_base.groupby("group_code", observed=False).cumcount() + 1
    group_code_str = df_base["group_code"].astype(str)
    df_base["group"] = group_code_str.map(group_display).fillna(group_code_str)

    prefix = os.path.join(out_dir, args.spec_prefix)
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
        "group_defs": group_defs,
    }
    group_tables = {"summary": summary, **groups}

    spec_tabs = build_specificity_panels(
        df_base,
        group_tables,
        title="Selectivity & Diversity Report",
        max_table=int(args.max_table),
        plot_x_range=args.plot_x_range,
        plot_y_range=args.plot_y_range,
        sample_cols=sample_cols,
        group_display=group_display,
        active_label=active_label,
        inactive_label=inactive_label,
    )

    html_path = f"{prefix}_interactive.html"
    if spec_tabs and not getattr(args, "no_html", False):
        output_file(html_path, title="Selectivity & Diversity Report")
        save(Tabs(tabs=spec_tabs))

    output_meta = {
        "prefix": prefix,
        "all": out_all,
        "other": other_path or "",
        "html": html_path,
    }
    return spec_tabs, output_meta


def build_interactive_html(summary_html: str, pattern_df: pd.DataFrame,
                           overlap_df: pd.DataFrame, run_tabs: List[Panel],
                           venn_tabs: List[Panel], spec_tabs: List[Panel],
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
        table = _table_layout(overlap_df, ["compound_key"] + cols, filename="overlap.csv", height=360)
        items.append(Div(text="<h3>Overlapping compounds (>=2 runs)</h3>"))
        items.append(table)

    comp_panel = _make_panel(column(*items, sizing_mode="stretch_width"), "Comparison")
    tabs = Tabs(tabs=[comp_panel] + spec_tabs + venn_tabs + run_tabs)
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
    ap.add_argument("--fill-scores", type=int, choices=[0,1], default=1,
                    help="1이면 top-N에 없더라도 전체 테이블에서 score를 찾아 채웁니다.")
    ap.add_argument("--score-col", default=None, help="Override score column")
    ap.add_argument("--out-dir", default=".", help="Output directory (default: current dir)")
    ap.add_argument("--out-prefix", default=None, help="Output prefix (default: compare_topN)")
    ap.add_argument("--min-overlap", type=int, default=2, help="Min runs to consider overlap")
    ap.add_argument("--roles", nargs="*", default=None,
                    help="Optional roles per run (active/inactive/both/other; same order as --runs)")
    ap.add_argument("--include-summary", type=int, choices=[0, 1], default=1,
                    help="Include per-run summary panels and outputs (default: 1)")
    ap.add_argument("--include-specificity", type=int, choices=[0, 1], default=1,
                    help="Include specificity report when roles are provided (default: 1)")
    ap.add_argument("--spec-top-n", type=int, default=None,
                    help="Top N per run for specificity candidate pool (default: --top-n)")
    ap.add_argument("--spec-prefix", default="tier_report",
                    help="Prefix for specificity outputs under --out-dir")
    ap.add_argument("--bb-top-k", type=int, default=5, help="Top K BB frequencies per BB column")
    ap.add_argument("--cluster", type=int, choices=[0, 1], default=1,
                    help="Enable clustering from BB SMILES (requires RDKit).")
    ap.add_argument("--cluster-mode", choices=["bbavg", "compound_or"], default="bbavg",
                    help="Clustering similarity: bbavg (per-BB avg) or compound_or (OR fingerprint)")
    ap.add_argument("--cluster-sim", type=float, default=0.8,
                    help="Tanimoto similarity cutoff (default: 0.8)")
    ap.add_argument("--cluster-radius", type=int, default=2,
                    help="Morgan fingerprint radius (default: 2)")
    ap.add_argument("--cluster-nbits", type=int, default=2048,
                    help="Morgan fingerprint nBits (default: 2048)")
    ap.add_argument("--cluster-rep", choices=["score", "medoid"], default="score",
                    help="Representative selection for diverse lists (default: score)")
    ap.add_argument("--max-table", type=int, default=200, help="Max rows per HTML table (specificity)")
    ap.add_argument("--plot-x-range", default="auto",
                    help="Specificity plot X range 'min,max' or 'auto' (default: auto)")
    ap.add_argument("--plot-y-range", default="auto",
                    help="Specificity plot Y range 'min,max' or 'auto' (default: auto)")
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
    ap.add_argument("--no-html", action="store_true", help="Skip interactive HTML output")
    args = ap.parse_args()

    run_paths = args.runs
    if len(run_paths) > 3:
        raise SystemExit("[ERROR] --runs supports up to 3 runs.")
    labels = args.labels or []
    if labels and len(labels) != len(run_paths):
        raise SystemExit("[ERROR] --labels must match number of --runs")
    if not labels:
        labels = [os.path.basename(os.path.normpath(p)) or f"run{i+1}" for i, p in enumerate(run_paths)]
    labels = [_sanitize_label(l) for l in labels]

    roles = args.roles or []
    if roles and len(roles) != len(run_paths):
        raise SystemExit("[ERROR] --roles must match number of --runs")
    roles = [r.strip().lower() for r in roles] if roles else []
    for r in roles:
        if r not in ("active", "inactive", "both", "other"):
            raise SystemExit("[ERROR] --roles must be one of: active, inactive, both, other")
    role_labels = _role_label_map(roles, labels)

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
    rep_col = "cluster_medoid" if args.cluster_rep == "medoid" else "cluster_rep"

    out_dir = args.out_dir or "."
    os.makedirs(out_dir, exist_ok=True)
    name_top = max(top_n_list) if top_n_list else top_n
    prefix = args.out_prefix or os.path.join(out_dir, f"compare_top{name_top}")
    html_path = f"{prefix}_interactive.html"

    run_info = []
    for idx, (path, label) in enumerate(zip(run_paths, labels)):
        hybrid = resolve_hybrid_path(path, args.preset)
        header = pd.read_csv(hybrid, sep="\t", nrows=0).columns.tolist()
        sample_cols = _sample_cols_from_header(header)
        want = [
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
        usecols = [c for c in want if c in header] + [score_col]
        for col in sample_cols:
            if col not in usecols and col in header:
                usecols.append(col)
        usecols = list(dict.fromkeys(usecols))
        df = pd.read_csv(hybrid, sep="\t", usecols=usecols)
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(float("-inf"))
        df["compound_key"] = _make_compound_key(df)
        score_map = df.groupby("compound_key")[score_col].max().to_dict()
        items = [(k, v) for k, v in score_map.items() if v == v]
        items.sort(key=lambda x: (-x[1], x[0]))
        rank_map_all = {k: i + 1 for i, (k, _) in enumerate(items)}
        denom = max(1, len(items) - 1)
        rank_pct_map_all = {k: 100.0 * (1.0 - (r - 1) / denom) for k, r in rank_map_all.items()}
        hitpct_source = "rank_pct"
        hitpct_map_all = rank_pct_map_all.copy()
        if "HitScore_pct" in df.columns:
            df["HitScore_pct"] = pd.to_numeric(df["HitScore_pct"], errors="coerce")
            hitpct_map_all = df.groupby("compound_key")["HitScore_pct"].max().to_dict()
            hitpct_source = "HitScore_pct"
        scores = np.array([v for v in score_map.values() if np.isfinite(v)], dtype=float)
        if scores.size == 0:
            mu = 0.0
            sd = 1.0
        else:
            mu = float(np.mean(scores))
            sd = float(np.std(scores))
            if sd == 0.0:
                sd = 1.0
        score_z_map_all = {k: ((v - mu) / sd) if np.isfinite(v) else np.nan for k, v in score_map.items()}
        run_top_n = top_n_list[idx] if top_n_list else top_n
        run_rec_n = rec_n_list[idx] if rec_n_list else rec_n
        df_best = df.sort_values(score_col, ascending=False).drop_duplicates("compound_key")
        df_top = df_best.head(run_top_n).copy()
        df_top["compound_key"] = _make_compound_key(df_top)
        if rank_map_all:
            df_top["rank"] = df_top["compound_key"].map(rank_map_all)
            df_top["rank_pct"] = df_top["compound_key"].map(rank_pct_map_all)
            df_top["score_z"] = df_top["compound_key"].map(score_z_map_all)
            df_top["HitScore_pct"] = df_top["compound_key"].map(hitpct_map_all)

        cluster_df = None
        if int(args.include_summary) == 1 and int(args.cluster) == 1:
            df_top, cluster_df = _cluster_top_hits(
                df_top,
                score_col=score_col,
                sim_cutoff=float(args.cluster_sim),
                radius=int(args.cluster_radius),
                nbits=int(args.cluster_nbits),
                mode=str(args.cluster_mode),
            )
        df_rec = _apply_recommend_filter(df_top)
        df_rec = df_rec.sort_values(score_col, ascending=False)
        df_rec = df_rec.drop_duplicates("compound_key").head(run_rec_n).copy()
        if rank_map_all:
            df_rec["rank"] = df_rec["compound_key"].map(rank_map_all)
            df_rec["rank_pct"] = df_rec["compound_key"].map(rank_pct_map_all)
            df_rec["score_z"] = df_rec["compound_key"].map(score_z_map_all)
            df_rec["HitScore_pct"] = df_rec["compound_key"].map(hitpct_map_all)
        top_score_map = df_top.groupby("compound_key")[score_col].max().to_dict()
        diverse = df_rec.copy()
        if rep_col in diverse.columns:
            diverse = diverse[_to_int_series(diverse[rep_col]) == 1].copy()
        diverse = diverse.sort_values(score_col, ascending=False).copy()

        freq_df = _build_bb_frequency_df(df_top, int(args.bb_top_k)) if int(args.include_summary) == 1 else pd.DataFrame()

        run_info.append({
            "label": label,
            "hybrid": hybrid,
            "df": df,
            "df_best": df_best,
            "top": df_top,
            "rec": df_rec,
            "diverse": diverse,
            "freq": freq_df,
            "cluster_df": cluster_df,
            "sample_cols": sample_cols,
            "role": roles[idx] if roles else "",
            "top_n": run_top_n,
            "rec_n": run_rec_n,
            "score_map": score_map,
            "top_score_map": top_score_map,
            "rank_map_all": rank_map_all,
            "rank_pct_map_all": rank_pct_map_all,
            "score_z_map_all": score_z_map_all,
            "hitpct_map_all": hitpct_map_all,
            "hitpct_source": hitpct_source,
        })

    for info in run_info:
        label = info["label"]
        top_path = f"{prefix}_{label}_top{info['top_n']}.tsv"
        rec_path = f"{prefix}_{label}_recommended{info['rec_n']}.tsv"
        div_path = f"{prefix}_{label}_diverse.tsv"
        freq_path = f"{prefix}_{label}_bb_frequency.tsv"
        cluster_path = f"{prefix}_{label}_clusters.tsv"
        info["top"].to_csv(top_path, sep="\t", index=False)
        info["rec"].to_csv(rec_path, sep="\t", index=False)
        if not info["diverse"].empty:
            info["diverse"].to_csv(div_path, sep="\t", index=False)
        if isinstance(info.get("freq"), pd.DataFrame):
            info["freq"].to_csv(freq_path, sep="\t", index=False)
        if info.get("cluster_df") is not None and not info["cluster_df"].empty:
            info["cluster_df"].to_csv(cluster_path, sep="\t", index=False)
        info["top_path"] = top_path
        info["rec_path"] = rec_path
        info["diverse_path"] = div_path
        info["freq_path"] = freq_path
        info["cluster_path"] = cluster_path

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
        for info in run_info:
            lbl = info["label"]
            score = info.get("top_score_map", {}).get(key, np.nan)
            if (score != score) and int(args.fill_scores) == 1:
                score = info.get("score_map", {}).get(key, np.nan)
            row[f"{lbl}_score"] = float(score) if score == score else np.nan
            score_z = info.get("score_z_map_all", {}).get(key, np.nan)
            row[f"{lbl}_score_z"] = float(score_z) if score == score and score_z == score_z else np.nan
            hitpct = info.get("hitpct_map_all", {}).get(key, np.nan)
            row[f"{lbl}_HitScore_pct"] = float(hitpct) if score == score and hitpct == hitpct else np.nan
            rank = info.get("rank_map_all", {}) or {}
            r = rank.get(key, np.nan)
            row[f"{lbl}_rank"] = int(r) if r == r else np.nan
            rank_pct = info.get("rank_pct_map_all", {}) or {}
            rp = rank_pct.get(key, np.nan)
            row[f"{lbl}_rank_pct"] = float(rp) if rp == rp else np.nan
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

    venn_tabs: List[Panel] = []
    spec_tabs: List[Panel] = []
    spec_outputs: Dict[str, str] = {}

    if int(args.include_specificity) == 1 and roles:
        role_to_info = {}
        for info in run_info:
            role = (info.get("role") or "").strip().lower()
            if not role or role == "other":
                continue
            if role not in role_to_info:
                role_to_info[role] = info
        if "active" in role_to_info and "inactive" in role_to_info:
            both_info = role_to_info.get("both")
            spec_tabs, spec_outputs = compute_specificity_report(
                role_to_info["active"],
                role_to_info["inactive"],
                both_info,
                score_col,
                args,
                out_dir,
                role_labels=role_labels,
            )

    if not args.no_html:
        summary_lines = [_style_block(), "<h2>Top-hit comparison</h2>"]
        cards = [
            ("Runs", str(len(run_paths))),
            ("Top-N", str(top_n)),
            ("Recommend-N", str(rec_n)),
            ("Union", str(len(union_df))),
            ("Overlap", str(len(overlap_df))),
        ]
        summary_lines.append(_summary_cards(cards))
        summary_lines.append("<ul>")
        summary_lines.append(f"<li>Score column: <b>{score_col}</b></li>")
        summary_lines.append(f"<li>Top-N default: <b>{top_n}</b></li>")
        summary_lines.append(f"<li>Min overlap: <b>{min_overlap}</b> runs</li>")
        summary_lines.append("<li>HitScore_pct: prefer column if present, else fallback to rank_pct</li>")
        for info in run_info:
            summary_lines.append(
                f"<li>{info['label']}: top={len(info['top'])} (target {info['top_n']}), "
                f"recommended={len(info['rec'])} (target {info['rec_n']})</li>"
            )
        summary_lines.append(f"<li>Union size: <b>{len(union_df)}</b></li>")
        summary_lines.append(f"<li>Overlap (>= {min_overlap}): <b>{len(overlap_df)}</b></li>")
        if spec_outputs:
            summary_lines.append(f"<li>Specificity outputs: <b>{spec_outputs.get('prefix','')}</b>*</li>")
        summary_lines.append("<li>Venn cell tabs show each intersection pattern</li>")
        summary_lines.append("</ul>")

        run_tabs = []
        for info in run_info:
            top_df = info["top"].copy()
            rec_df = info["rec"].copy()
            div_df = info.get("diverse", pd.DataFrame()).copy()
            freq_df = info.get("freq", pd.DataFrame())
            if _HAS_RDKIT:
                top_df = _add_images(top_df)
                rec_df = _add_images(rec_df)
                if not div_df.empty:
                    div_df = _add_images(div_df)
            for col in ["BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x"]:
                if col in top_df.columns:
                    top_df[col] = top_df[col].fillna("NA")
                if col in rec_df.columns:
                    rec_df[col] = rec_df[col].fillna("NA")
                if col in div_df.columns:
                    div_df[col] = div_df[col].fillna("NA")
            cols_top = ["compound_key", "rank", "rank_pct", "HitScore_pct", score_col, "score_z",
                        "BB1_x", "BB2_x", "BB3_x", "BB4_x", "CP_x"]
            cols_rec = cols_top
            cols_top = [c for c in cols_top if c in top_df.columns]
            cols_rec = [c for c in cols_rec if c in rec_df.columns]
            cols_div = [c for c in cols_rec if c in div_df.columns]
            if "compound_img" in top_df.columns:
                cols_top = ["compound_img"] + cols_top
            if "compound_img" in rec_df.columns:
                cols_rec = ["compound_img"] + cols_rec
            if not div_df.empty and "compound_img" in div_df.columns:
                cols_div = ["compound_img"] + cols_div

            top_table = _table_layout(
                top_df,
                cols_top,
                filename=f"{info['label']}_top.csv",
                height=320,
            )
            rec_table = _table_layout(
                rec_df,
                cols_rec,
                filename=f"{info['label']}_recommended.csv",
                height=260,
            )
            items = [
                Div(text=f"<h3>{info['label']}</h3>"),
                Div(text=f"<p>Top-N table (n={len(top_df)})</p>"),
                top_table,
                Div(text=f"<p>Recommended (n={len(rec_df)})</p>"),
                rec_table,
            ]

            if not div_df.empty:
                div_table = _table_layout(
                    div_df,
                    cols_div,
                    filename=f"{info['label']}_diverse.csv",
                    height=260,
                )
                items.append(Div(text=f"<p>Diverse (cluster reps) (n={len(div_df)})</p>"))
                items.append(div_table)

            if int(args.include_summary) == 1 and not top_df.empty and _HAS_BOKEH:
                df_plot = top_df.copy()
                if "rank" not in df_plot.columns:
                    df_plot["rank"] = np.arange(1, len(df_plot) + 1, dtype=int)
                else:
                    df_plot["rank"] = pd.to_numeric(df_plot["rank"], errors="coerce").fillna(0).astype(int)
                if "Consensus_hit" in df_plot.columns:
                    df_plot["consensus_label"] = _to_int_series(df_plot["Consensus_hit"]).map({1: "pass"}).fillna("fail")
                else:
                    df_plot["consensus_label"] = "na"
                y_col = _pick_y_column(df_plot.columns.tolist(), score_col)
                for col in [score_col, y_col]:
                    df_plot[col] = _safe_float_series(df_plot[col])
                if "mean_R1_norm" in df_plot.columns:
                    size = np.log10(_safe_float_series(df_plot["mean_R1_norm"]).fillna(0.0) + 1.0) * 2.0 + 6.0
                    df_plot["size"] = np.clip(size, 4.0, 12.0)
                else:
                    df_plot["size"] = 8.0
                colors = {"pass": "#2c7fb8", "fail": "#bdbdbd", "na": "#bdbdbd"}
                df_plot["color"] = df_plot["consensus_label"].map(colors).fillna("#bdbdbd")
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
                    title=f"Top {len(df_plot)} hits: {score_col} vs {y_col}",
                    x_axis_label=score_col,
                    y_axis_label=y_col,
                    height=380,
                    tools="pan,wheel_zoom,box_zoom,reset,save",
                )
                p_scatter.add_tools(hover)
                p_scatter.scatter(x=score_col, y=y_col, size="size", color="color", alpha=0.75, source=cds)
                items.append(Div(text="<h4>Score scatter</h4>"))
                items.append(p_scatter)

                if isinstance(freq_df, pd.DataFrame) and not freq_df.empty:
                    bb_smiles_map = build_bb_smiles_map(top_df)
                    bb_plots = []
                    img_cache: Dict[str, str] = {}
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
                    if bb_plots:
                        items.append(Div(text="<h4>BB frequency in top hits</h4>"))
                        items.extend(bb_plots)

            panel = _make_panel(column(*items, sizing_mode="stretch_width"), info["label"])
            run_tabs.append(panel)

        # Build venn cell tabs
        if not union_df.empty:
            score_cols = [f"{lbl}_score" for lbl in labels if f"{lbl}_score" in union_df.columns]
            rank_cols = [f"{lbl}_rank" for lbl in labels if f"{lbl}_rank" in union_df.columns]
            rank_pct_cols = [f"{lbl}_rank_pct" for lbl in labels if f"{lbl}_rank_pct" in union_df.columns]
            hitpct_cols = [f"{lbl}_HitScore_pct" for lbl in labels if f"{lbl}_HitScore_pct" in union_df.columns]
            base_cols = ["compound_key", "pattern"] + rank_cols + rank_pct_cols + hitpct_cols + score_cols
            base_cols = [c for c in base_cols if c in union_df.columns]

            pattern_lookup = {p: i for i, p in enumerate(pattern_df["pattern"].tolist())} if not pattern_df.empty else {}
            for pattern in _pattern_order(labels):
                if pattern not in pattern_lookup:
                    continue
                sub = union_df[union_df["pattern"] == pattern].copy()
                if sub.empty:
                    continue
                if score_cols:
                    sub_scores = sub[score_cols].apply(pd.to_numeric, errors="coerce")
                    sub["best_score"] = sub_scores.max(axis=1)
                    sub = sub.sort_values("best_score", ascending=False)
                    cols = base_cols + ["best_score"] if "best_score" not in base_cols else base_cols
                else:
                    cols = base_cols
                cols = [c for c in cols if c in sub.columns]
                table = _table_layout(sub, cols, filename=f"venn_{pattern}.csv", height=360)
                panel = _make_panel(
                    column(
                        Div(text=f"<h3>Venn cell: {pattern} (n={len(sub)})</h3>"),
                        table,
                        sizing_mode="stretch_width",
                    ),
                    f"Venn {pattern}",
                )
                venn_tabs.append(panel)

        spec_tabs_main = []
        build_interactive_html("\n".join(summary_lines), pattern_df, overlap_df, run_tabs, venn_tabs, spec_tabs_main, html_path)

    print(f"[INFO] score_col={score_col}")
    print(f"[INFO] outputs: {union_path}, {overlap_path}, {pattern_path}")
    for info in run_info:
        extra = f", {info['diverse_path']}, {info['freq_path']}"
        if info.get("cluster_df") is not None and not info["cluster_df"].empty:
            extra += f", {info['cluster_path']}"
        print(f"[INFO] {info['label']}: {info['top_path']}, {info['rec_path']}{extra}")
    if spec_outputs:
        print(f"[INFO] specificity outputs: {spec_outputs.get('prefix','')}_*.tsv")
        if spec_outputs.get("html"):
            print(f"[INFO] specificity html: {spec_outputs.get('html')}")
    if not args.no_html:
        print(f"[INFO] html: {html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
