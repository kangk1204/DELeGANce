#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Optional

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
    return re.sub(r"_LIB[\w\.-]+$", "", s)


def load_params(run_root: Path) -> Dict:
    hp = run_root / "03_normalized" / "glm_full_dev_cpu_fp64" / "hit_params.json"
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


def build_core(df: pd.DataFrame, sample_cols: List[str]) -> pd.DataFrame:
    lib_col = pick_col(df, ["LIB_ID", "LIB_ID_x", "LIB_ID_y"])
    id_col = pick_col(df, ["ID", "ID_x", "ID_y"])
    bb1_col = pick_col(df, ["BB1", "BB1_x", "BB1_y"])
    bb2_col = pick_col(df, ["BB2", "BB2_x", "BB2_y"])
    bb3_col = pick_col(df, ["BB3", "BB3_x", "BB3_y"])
    bb4_col = pick_col(df, ["BB4", "BB4_x", "BB4_y"])

    sm1 = pick_col(df, ["bb1_smiles", "BB1_smiles", "SMILES1"])
    sm2 = pick_col(df, ["bb2_smiles", "BB2_smiles", "SMILES2"])
    sm3 = pick_col(df, ["bb3_smiles", "BB3_smiles", "SMILES3"])
    sm4 = pick_col(df, ["bb4_smiles", "BB4_smiles", "SMILES4"])

    out = pd.DataFrame()
    out["LibID"] = df[lib_col] if lib_col else ""
    out["cycles"] = df["cycles"] if "cycles" in df.columns else pd.NA

    bb1_raw = df[bb1_col] if bb1_col else "NA"
    bb2_raw = df[bb2_col] if bb2_col else "NA"
    bb3_raw = df[bb3_col] if bb3_col else "NA"
    bb4_raw = df[bb4_col] if bb4_col else "NA"

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
        return int((s.astype(str).str.lower().isin(["true", "1", "yes", "y"]).sum()))

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


def truthy_mask(series: pd.Series, length: int) -> pd.Series:
    if series is None:
        return pd.Series([False] * length)
    if getattr(series, "dtype", None) == bool:
        return series.fillna(False)
    s = series.astype(str).str.lower()
    return s.isin(["true", "1", "yes", "y", "t"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run_root", action="append", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--top_n", type=int, default=1000)
    args = p.parse_args()

    out_path = Path(args.out)

    engine = None
    try:
        import openpyxl  # noqa: F401
        engine = "openpyxl"
    except Exception:
        try:
            import xlsxwriter  # noqa: F401
            engine = "xlsxwriter"
        except Exception:
            raise SystemExit("openpyxl 또는 xlsxwriter가 필요합니다.")

    with pd.ExcelWriter(out_path, engine=engine) as writer:
        guide = build_guide()
        guide.to_excel(writer, sheet_name="Guide", index=False)

        summary_rows = []

        for run_root in args.run_root:
            run_root = Path(run_root)
            run_name = run_root.name
            hybrid = run_root / "03_normalized" / "glm_full_dev_cpu_fp64" / "05_hybrid_annot.tsv"
            if not hybrid.exists():
                raise SystemExit(f"missing: {hybrid}")

            params = load_params(run_root)
            sample_cols = sample_columns_from_params(params)

            df = pd.read_csv(hybrid, sep="\t", low_memory=False)
            summary_rows.append(build_summary(df, run_name))

            core = build_core(df, sample_cols)

            if "HitScore_GLM" in core.columns:
                core = core.sort_values("HitScore_GLM", ascending=False)

            core.insert(0, "Rank_GLM", range(1, len(core) + 1))

            top_n = min(args.top_n, len(core))
            top = core.head(top_n).copy()

            consensus_mask = truthy_mask(core.get("Consensus_hit"), len(core))
            consensus = core[consensus_mask].copy()

            core = normalize_bool_cols(core)
            top = normalize_bool_cols(top)
            consensus = normalize_bool_cols(consensus)

            core_sheet = f"{run_name}_All_Core"
            top_sheet = f"{run_name}_Top{top_n}"
            cons_sheet = f"{run_name}_Consensus"
            param_sheet = f"{run_name}_Params"

            core.to_excel(writer, sheet_name=core_sheet[:31], index=False)
            top.to_excel(writer, sheet_name=top_sheet[:31], index=False)
            consensus.to_excel(writer, sheet_name=cons_sheet[:31], index=False)

            if params:
                flat = []
                for k, v in params.get("args", {}).items():
                    flat.append((k, json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v))
                p_df = pd.DataFrame(flat, columns=["param", "value"])
                p_df.to_excel(writer, sheet_name=param_sheet[:31], index=False)

        summary = pd.concat(summary_rows, ignore_index=True)
        summary.to_excel(writer, sheet_name="Summary", index=False)

    print(f"[OK] Excel written: {out_path}")


if __name__ == "__main__":
    main()
