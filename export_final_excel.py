#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


def pick_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def strip_lib_suffix(bb: str) -> str:
    if bb is None:
        return "NA"
    s = str(bb)
    if s in ("", "NA", "nan", "None"):
        return "NA"
    # Anchored token: strip only a trailing "_LIB<lib>" namespace suffix
    return re.sub(r"_LIB[^_]+$", "", s)


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


def load_params_from_annot(annot_path: Path) -> Dict:
    hp = annot_path.parent / "hit_params.json"
    if hp.exists():
        return json.loads(hp.read_text())
    return {}


def sample_columns_from_params(params: Dict) -> List[str]:
    cols = []
    norm = params.get("normalized_columns", {})
    del2 = norm.get("del2")
    if del2:
        cols.append(del2)
    for key in ("r1", "r2", "neg"):
        cols.extend(norm.get(key, []) or [])
    return cols


def sample_columns_from_header(columns: List[str]) -> List[str]:
    """Fallback when hit_params.json is absent (03_call_hits.py run standalone):
    a sample column is any <name> that also has a <name>_CPM counterpart."""
    return [c[:-4] for c in columns if c.endswith("_CPM") and c[:-4] in columns]


def build_core(df: pd.DataFrame, sample_cols: List[str]) -> pd.DataFrame:
    lib_col = pick_col(df, ["LIB_ID", "LIB_ID_x", "LIB_ID_y"])
    bb1_col = pick_col(df, ["BB1", "BB1_x", "BB1_y"])
    bb2_col = pick_col(df, ["BB2", "BB2_x", "BB2_y"])
    bb3_col = pick_col(df, ["BB3", "BB3_x", "BB3_y"])
    bb4_col = pick_col(df, ["BB4", "BB4_x", "BB4_y"])

    sm1 = pick_col(df, ["bb1_smiles", "BB1_smiles", "SMILES1"])
    sm2 = pick_col(df, ["bb2_smiles", "BB2_smiles", "SMILES2"])
    sm3 = pick_col(df, ["bb3_smiles", "BB3_smiles", "SMILES3"])
    sm4 = pick_col(df, ["bb4_smiles", "BB4_smiles", "SMILES4"])

    # Keep df's index so scalar/Series assignments stay aligned even when LIB_ID is missing
    out = pd.DataFrame(index=df.index)
    out["LibID"] = df[lib_col] if lib_col else ""
    out["cycles"] = df["cycles"] if "cycles" in df.columns else pd.NA

    def _bb_series(col: Optional[str]) -> pd.Series:
        # Missing BB column -> Series of "NA" (a bare str has no .map)
        return df[col] if col else pd.Series(["NA"] * len(df), index=df.index)

    bb1_raw = _bb_series(bb1_col)
    bb2_raw = _bb_series(bb2_col)
    bb3_raw = _bb_series(bb3_col)
    bb4_raw = _bb_series(bb4_col)

    out["BB1"] = bb1_raw.map(strip_lib_suffix)
    out["BB2"] = bb2_raw.map(strip_lib_suffix)
    out["BB3"] = bb3_raw.map(strip_lib_suffix)
    out["BB4"] = bb4_raw.map(strip_lib_suffix)

    out["ID"] = (
        out["LibID"].astype(str) + "_" +
        out["BB1"].astype(str) + "_" +
        out["BB2"].astype(str) + "_" +
        out["BB3"].astype(str) + "_" +
        out["BB4"].astype(str)
    )
    # Backward-compatible alias
    out["ID_display"] = out["ID"]

    out["SMILES1"] = df[sm1] if sm1 else ""
    out["SMILES2"] = df[sm2] if sm2 else ""
    out["SMILES3"] = df[sm3] if sm3 else ""
    out["SMILES4"] = df[sm4] if sm4 else ""

    def add_if_exists(cols: List[str]):
        for c in cols:
            if c in df.columns and c not in out.columns:
                out[c] = df[c]

    add_if_exists([
        "HitScore_GLM", "HitScore_RS", "HitScore_pct",
        "GLM_hit", "RS_pass", "Consensus_hit", "NEG_hard_fail",
        "LFC_R1_vs_DEL2_used", "LFC_R2_vs_DEL2_used",
        "LFC_NEG_R1_vs_DEL2_used", "LFC_NEG_R2_vs_DEL2_used",
        "LFC_NEG_vs_DEL2_used", "LFC_NEG_centered",
        "LFC_NEG_centered_R1", "LFC_NEG_centered_R2",
        "NEG_center_shift_R1", "NEG_center_shift_R2",
        "NEG_control_used",
        "log2Boost_R2vsR1", "mean_log2Boost_R2vsR1_paired",
        "avg_R1", "avg_R2", "mean_R1_norm", "mean_R2_norm", "DEL2_norm",
        "E_component", "W_count", "SynthonScore", "SynthonBonus",
        "NEG_penalty_used", "Penalty_NEG", "Penalty",
        "fail_reasons", "pass_filters",
    ])

    for s in sample_cols:
        if s in df.columns:
            out[s] = df[s]
        cpm = f"{s}_CPM"
        if cpm in df.columns:
            out[cpm] = df[cpm]

    preferred = [
        "LibID", "ID", "cycles", "BB1", "BB2", "BB3", "BB4", "ID_display",
    ]
    front = [c for c in preferred if c in out.columns]
    rest = [c for c in out.columns if c not in front]
    return out[front + rest]


def build_summary(df: pd.DataFrame, run_name: str) -> pd.DataFrame:
    def count_true(col: str) -> int:
        if col not in df.columns:
            return 0
        s = df[col]
        if s.dtype == bool:
            return int(s.sum())
        return int((s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "y", "t"]).sum()))

    summary = {
        "run": run_name,
        "rows": len(df),
        "GLM_hit": count_true("GLM_hit"),
        "RS_pass": count_true("RS_pass"),
        "Consensus_hit": count_true("Consensus_hit"),
        "NEG_hard_fail": count_true("NEG_hard_fail"),
    }
    return pd.DataFrame([summary])


def build_guide() -> pd.DataFrame:
    items = [
        ("ID", "LibID + BB1~BB4 (LIB 접미사 제거)"),
        ("LibID", "라이브러리 ID"),
        ("cycles", "사이클 수(3 또는 4)"),
        ("BB1~BB4", "BB 코드(표시용, LIB 접미사 제거)"),
        ("SMILES1~SMILES4", "BB 구조 (SMILES)"),
        ("HitScore_GLM", "GLM 기반 히트 스코어"),
        ("HitScore_RS", "ReadScaler 기반 히트 스코어"),
        ("HitScore_pct", "GLM 스코어 백분위"),
        ("GLM_hit", "GLM 기준 히트 여부"),
        ("RS_pass", "ReadScaler 필터 통과 여부"),
        ("Consensus_hit", "GLM+RS 합의 히트"),
        ("NEG_hard_fail", "NEG 하드 게이트 실패"),
        ("LFC_R1_vs_DEL2_used", "R1 vs DEL2 LFC"),
        ("LFC_R2_vs_DEL2_used", "R2 vs DEL2 LFC"),
        ("LFC_NEG_R1_vs_DEL2_used", "NEG_R1 vs DEL2 LFC"),
        ("LFC_NEG_R2_vs_DEL2_used", "NEG_R2 vs DEL2 LFC"),
        ("LFC_NEG_centered", "NEG vs DEL2 (센터링)"),
        ("LFC_NEG_centered_R1", "NEG_R1 vs DEL2 (센터링, QC용)"),
        ("LFC_NEG_centered_R2", "NEG_R2 vs DEL2 (센터링, QC용)"),
        ("NEG_center_shift_R1", "NEG_R1 센터링 시프트"),
        ("NEG_center_shift_R2", "NEG_R2 센터링 시프트"),
        ("log2Boost_R2vsR1", "R2/R1 log2 boost"),
        ("mean_log2Boost_R2vsR1_paired", "paired log2 boost 평균"),
        ("avg_R1/avg_R2", "R1/R2 평균 카운트"),
        ("mean_R1_norm/mean_R2_norm", "R1/R2 정규화 평균"),
        ("DEL2_norm", "DEL2 정규화 카운트"),
        ("E_component", "GLM E 항"),
        ("W_count", "가중치 항"),
        ("SynthonScore/SynthonBonus", "Synthon 기반 보정"),
        ("NEG_penalty_used", "NEG 패널티 적용 방식"),
        ("Penalty_NEG/Penalty", "NEG/총 패널티"),
        ("<sample>", "샘플 raw count"),
        ("<sample>_CPM", "샘플 CPM"),
    ]
    return pd.DataFrame(items, columns=["Field", "Meaning"])

def bool_to_english(v):
    if v is None:
        return v
    try:
        if isinstance(v, float) and pd.isna(v):
            return v
    except Exception:
        pass
    if isinstance(v, bool):
        return "TRUE" if v else "FALSE"
    s = str(v).strip().lower()
    if s in ("true", "1", "yes", "y"):
        return "TRUE"
    if s in ("false", "0", "no", "n"):
        return "FALSE"
    return v

def normalize_bool_cols(df: pd.DataFrame) -> pd.DataFrame:
    bool_cols = ["GLM_hit", "RS_pass", "Consensus_hit", "NEG_hard_fail", "pass_filters"]
    for c in bool_cols:
        if c in df.columns:
            df[c] = df[c].map(bool_to_english)
    return df


_SHEET_BAD_CHARS = re.compile(r"[\[\]:*?/\\]")


def sheet_name(run_name: str, suffix: str, used: set) -> str:
    """Excel sheet names are limited to 31 chars and must be unique.
    Keep the suffix intact, truncate the run part, and add ~N when a collision remains
    (plain [:31] truncation made All_Core/TopN/Consensus/Params identical for long run names)."""
    run = _SHEET_BAD_CHARS.sub("_", str(run_name))
    base = run[: max(1, 31 - len(suffix))] + suffix
    name, k = base, 1
    while name in used:
        tag = f"~{k}"
        name = base[: 31 - len(tag)] + tag
        k += 1
    used.add(name)
    return name


def truthy_mask(series: pd.Series, length: int) -> pd.Series:
    if series is None:
        return pd.Series([False] * length)
    if getattr(series, "dtype", None) == bool:
        return series.fillna(False)
    s = series.astype(str).str.lower()
    return s.isin(["true", "1", "yes", "y", "t"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run_root", action="append", default=[])
    p.add_argument("--annot_tsv", action="append", default=[])
    p.add_argument("--prefer_dir", default="",
                   help="Preferred subdir under run_root for 05_hybrid_annot.tsv")
    p.add_argument("--out", required=True)
    p.add_argument("--top_n", type=int, default=1000)
    args = p.parse_args()

    out_path = Path(args.out)
    if not args.run_root and not args.annot_tsv:
        raise SystemExit("[ERROR] --run_root or --annot_tsv is required.")

    engine = None
    try:
        import openpyxl  # noqa: F401
        engine = "openpyxl"
    except Exception:
        try:
            import xlsxwriter  # noqa: F401
            engine = "xlsxwriter"
        except Exception:
            raise SystemExit("openpyxl 또는 xlsxwriter가 필요합니다.") from None

    out_path.parent.mkdir(parents=True, exist_ok=True)
    used_sheets: set = set()
    with pd.ExcelWriter(out_path, engine=engine) as writer:
        guide = build_guide()
        guide.to_excel(writer, sheet_name="Guide", index=False)

        summary_rows = []

        annot_jobs: List[Tuple[Path, str]] = []
        for run_root in args.run_root:
            rr = Path(run_root)
            annot_path = _resolve_annot(rr, args.prefer_dir)
            annot_jobs.append((annot_path, rr.name))
        for annot_path in args.annot_tsv:
            apath = Path(annot_path)
            if not apath.exists():
                raise SystemExit(f"[ERROR] annot_tsv not found: {apath}")
            annot_jobs.append((apath, _guess_run_name(apath)))

        for annot_path, run_name in annot_jobs:
            params = load_params_from_annot(annot_path)
            sample_cols = sample_columns_from_params(params)

            df = pd.read_csv(annot_path, sep="\t", low_memory=False)
            if not sample_cols:
                sample_cols = sample_columns_from_header(list(df.columns))
                print(f"[WARN] hit_params.json not found next to {annot_path}; "
                      f"sample columns inferred from header ({len(sample_cols)} found)")
            summary_rows.append(build_summary(df, run_name))

            core = build_core(df, sample_cols)

            if "HitScore_GLM" in core.columns:
                # stable sort + ID tie-breaker so Rank_GLM is reproducible for equal scores
                core = core.sort_values(["HitScore_GLM", "ID"], ascending=[False, True], kind="mergesort")

            core.insert(0, "Rank_GLM", range(1, len(core) + 1))

            top_n = min(args.top_n, len(core))
            top = core.head(top_n).copy()

            consensus_mask = truthy_mask(core.get("Consensus_hit"), len(core))
            consensus = core[consensus_mask].copy()

            core = normalize_bool_cols(core)
            top = normalize_bool_cols(top)
            consensus = normalize_bool_cols(consensus)

            core_sheet = sheet_name(run_name, "_All_Core", used_sheets)
            top_sheet = sheet_name(run_name, f"_Top{top_n}", used_sheets)
            cons_sheet = sheet_name(run_name, "_Consensus", used_sheets)
            param_sheet = sheet_name(run_name, "_Params", used_sheets)

            core.to_excel(writer, sheet_name=core_sheet, index=False)
            top.to_excel(writer, sheet_name=top_sheet, index=False)
            consensus.to_excel(writer, sheet_name=cons_sheet, index=False)

            if params:
                flat = []
                for k, v in params.get("args", {}).items():
                    flat.append((k, json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v))
                p_df = pd.DataFrame(flat, columns=["param", "value"])
                p_df.to_excel(writer, sheet_name=param_sheet, index=False)

        summary = pd.concat(summary_rows, ignore_index=True)
        summary.to_excel(writer, sheet_name="Summary", index=False)

    print(f"[OK] Excel written: {out_path}")


if __name__ == "__main__":
    main()
