#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd


# Anchored: strip only a trailing "_LIB<lib>" namespace token from a single BB value
LIB_SUFFIX_RE = re.compile(r"_LIB[^_]+$")
# Token-anchored variant for full tag IDs (cycles_BB1_BB2_BB3[_BB4]); the old unanchored
# pattern (\w includes "_") consumed everything after the first "_LIB" and truncated the ID.
LIB_ANY_RE = re.compile(r"_LIB[^_]+(?=_|$)")


def strip_lib_suffix(value) -> str:
    if value is None:
        return "NA"
    s = str(value)
    if s in ("", "NA", "nan", "None"):
        return "NA"
    return LIB_SUFFIX_RE.sub("", s)


def strip_lib_anywhere(value) -> str:
    if value is None:
        return ""
    return LIB_ANY_RE.sub("", str(value))


def pick_col(cols: List[str], candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in cols:
            return c
    return None

def pick_score_col(cols: List[str]) -> Optional[str]:
    for c in ("HitScore_GLM", "HitScore", "HitScore_RS"):
        if c in cols:
            return c
    return None


def parse_id_fields(id_str: str) -> Tuple[str, str, str, str]:
    s = str(id_str or "").strip()
    raw = [t for t in re.split(r"[\|_,:;/\s]+", s) if t != ""]
    if raw and re.fullmatch(r"\d+", raw[0]):
        raw = raw[1:]
    parts: List[str] = []
    i = 0
    while i < len(raw) and len(parts) < 4:
        t = raw[i]
        # Same guard as 03_call_hits.parse_id_fields: a token that itself starts with "LIB"
        # is a BB id, not a namespace fragment to be glued onto the previous token.
        if (
            i + 1 < len(raw)
            and raw[i + 1].startswith("LIB")
            and (t not in ("NA", "") and not t.startswith("LIB"))
        ):
            parts.append(f"{t}_{raw[i + 1]}")
            i += 2
        else:
            parts.append(t)
            i += 1
    while len(parts) < 4:
        parts.append("NA")
    return parts[0], parts[1], parts[2], parts[3]


def compute_display_fields(df: pd.DataFrame) -> pd.DataFrame:
    cols = list(df.columns)
    lib_col = pick_col(cols, ["LIB_ID", "LIB_ID_x", "LIB_ID_y", "LibID", "lib_id", "lib_id_x", "lib_id_y"])
    id_col = pick_col(cols, ["ID", "ID_x", "ID_y", "id", "id_x", "id_y"])

    bb1_col = pick_col(cols, ["BB1", "BB1_x", "BB1_y", "bb1_id"])
    bb2_col = pick_col(cols, ["BB2", "BB2_x", "BB2_y", "bb2_id"])
    bb3_col = pick_col(cols, ["BB3", "BB3_x", "BB3_y", "bb3_id"])
    bb4_col = pick_col(cols, ["BB4", "BB4_x", "BB4_y", "bb4_id"])

    if bb1_col and bb2_col and bb3_col and bb4_col:
        bb1 = df[bb1_col].map(strip_lib_suffix)
        bb2 = df[bb2_col].map(strip_lib_suffix)
        bb3 = df[bb3_col].map(strip_lib_suffix)
        bb4 = df[bb4_col].map(strip_lib_suffix)
    else:
        b1, b2, b3, b4 = [], [], [], []
        for s in df[id_col].astype(str):
            p1, p2, p3, p4 = parse_id_fields(s)
            b1.append(strip_lib_suffix(p1))
            b2.append(strip_lib_suffix(p2))
            b3.append(strip_lib_suffix(p3))
            b4.append(strip_lib_suffix(p4))
        bb1, bb2, bb3, bb4 = map(pd.Series, (b1, b2, b3, b4))

    lib = df[lib_col].astype(str) if lib_col else pd.Series([""] * len(df))
    id_display = lib + "_" + bb1.astype(str) + "_" + bb2.astype(str) + "_" + bb3.astype(str) + "_" + bb4.astype(str)
    valid = ~lib.str.lower().isin(["", "na", "nan", "none"])
    if id_col:
        fallback = df[id_col].astype(str).map(strip_lib_anywhere)
        id_display = id_display.where(valid, fallback)

    out = pd.DataFrame({
        "LibID": lib,
        "ID": id_display,
        "BB1": bb1,
        "BB2": bb2,
        "BB3": bb3,
        "BB4": bb4,
    })
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="Compare top hits and output display-normalized report")
    p.add_argument("--current", required=True)
    p.add_argument("--previous", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--top_n", type=int, default=200)
    p.add_argument("--top_summary_out", default="", help="(optional) write top-N shift summary TSV")
    p.add_argument("--top_summary_out_current", default="", help="(optional) write top-N by current score TSV")
    p.add_argument("--top_summary_n", type=int, default=50)
    p.add_argument("--top_summary_xlsx", default="", help="(optional) write Excel with full + summaries")
    args = p.parse_args()

    cur_df = pd.read_csv(args.current, sep="\t", low_memory=False)
    prev_df = pd.read_csv(args.previous, sep="\t", low_memory=False)

    score_cur_col = pick_score_col(cur_df.columns)
    score_prev_col = pick_score_col(prev_df.columns)
    if not score_cur_col or not score_prev_col:
        raise SystemExit("[ERROR] Missing HitScore columns in inputs.")

    if score_cur_col == score_prev_col:
        score = score_cur_col
        cur_df[score] = pd.to_numeric(cur_df[score], errors="coerce")
        prev_df[score] = pd.to_numeric(prev_df[score], errors="coerce")
    else:
        score = "_SCORE"
        cur_df[score] = pd.to_numeric(cur_df[score_cur_col], errors="coerce")
        prev_df[score] = pd.to_numeric(prev_df[score_prev_col], errors="coerce")

    cur_id = pick_col(cur_df.columns, ["ID", "ID_x", "ID_y", "id", "id_x", "id_y"])
    cur_lib = pick_col(cur_df.columns, ["LIB_ID", "LIB_ID_x", "LIB_ID_y"])
    prev_id = pick_col(prev_df.columns, ["ID", "ID_x", "ID_y", "id", "id_x", "id_y"])
    prev_lib = pick_col(prev_df.columns, ["LIB_ID", "LIB_ID_x", "LIB_ID_y"])

    if not cur_id or not cur_lib or not prev_id or not prev_lib:
        raise SystemExit("[ERROR] Missing ID/LIB columns in inputs.")

    cur_df["_KEY"] = cur_df[cur_lib].astype(str) + "|" + cur_df[cur_id].astype(str)
    prev_df["_KEY"] = prev_df[prev_lib].astype(str) + "|" + prev_df[prev_id].astype(str)

    # Stable sort with _KEY tie-breaker so top-N membership/ranks are reproducible for equal scores
    cur_top = cur_df.sort_values([score, "_KEY"], ascending=[False, True], kind="mergesort").head(args.top_n)
    prev_top = prev_df.sort_values([score, "_KEY"], ascending=[False, True], kind="mergesort").head(args.top_n)

    cur_set = set(cur_top["_KEY"])
    prev_set = set(prev_top["_KEY"])
    union_keys = list(cur_set | prev_set)

    cols = [
        score, "LFC_R1_vs_DEL2", "LFC_R2_vs_DEL2", "LFC_NEG_centered",
        "LFC_NEG_vs_DEL2_used", "LFC_NEG_vs_DEL2",
        "E_component", "Penalty", "NEG_hard_fail", "GLM_hit", "RS_pass", "Consensus_hit",
    ]
    cols = [c for c in cols if c in cur_df.columns and c in prev_df.columns]

    # Restrict to the union of the two top-N sets first, then OUTER-merge so a compound that is
    # absent from one input file is kept (the previous inner merge silently dropped it).
    cur_sub = cur_df.loc[cur_df["_KEY"].isin(union_keys), ["_KEY", *cols]].drop_duplicates("_KEY").copy()
    prev_sub = prev_df.loc[prev_df["_KEY"].isin(union_keys), ["_KEY", *cols]].drop_duplicates("_KEY").copy()

    m = cur_sub.merge(prev_sub, on="_KEY", how="outer", suffixes=("_cur", "_prev"))

    # rank info
    cur_rank = {k: i + 1 for i, k in enumerate(cur_top["_KEY"].tolist())}
    prev_rank = {k: i + 1 for i, k in enumerate(prev_top["_KEY"].tolist())}
    m["rank_cur"] = m["_KEY"].map(cur_rank)
    m["rank_prev"] = m["_KEY"].map(prev_rank)
    m["status"] = "both"
    m.loc[m["rank_cur"].isna(), "status"] = "only_prev"
    m.loc[m["rank_prev"].isna(), "status"] = "only_cur"
    # Distinguish "not in the other run's top-N" from "absent from the other run's file"
    in_cur_file = m["_KEY"].isin(set(cur_df["_KEY"]))
    in_prev_file = m["_KEY"].isin(set(prev_df["_KEY"]))
    m["in_cur_file"] = in_cur_file
    m["in_prev_file"] = in_prev_file

    # diffs
    for c in cols:
        if c in ("NEG_hard_fail", "GLM_hit", "RS_pass", "Consensus_hit"):
            continue
        m[c + "_diff"] = pd.to_numeric(m[c + "_cur"], errors="coerce") - pd.to_numeric(m[c + "_prev"], errors="coerce")

    # display fields (no LIB suffix)
    disp_cur = compute_display_fields(cur_df)
    disp_prev = compute_display_fields(prev_df)
    disp_cur["_KEY"] = cur_df["_KEY"]
    disp_prev["_KEY"] = prev_df["_KEY"]
    disp_cur = disp_cur[disp_cur["_KEY"].isin(union_keys)].drop_duplicates("_KEY")
    disp_prev = disp_prev[disp_prev["_KEY"].isin(union_keys)].drop_duplicates("_KEY")
    disp = disp_cur.merge(disp_prev, on="_KEY", how="outer", suffixes=("_cur", "_prev"))
    # prefer current display values; fall back to previous when the compound is only in prev
    disp_out = pd.DataFrame({
        "_KEY": disp["_KEY"],
        "LibID": disp["LibID_cur"].combine_first(disp["LibID_prev"]),
        "ID": disp["ID_cur"].combine_first(disp["ID_prev"]),
        "BB1": disp["BB1_cur"].combine_first(disp["BB1_prev"]),
        "BB2": disp["BB2_cur"].combine_first(disp["BB2_prev"]),
        "BB3": disp["BB3_cur"].combine_first(disp["BB3_prev"]),
        "BB4": disp["BB4_cur"].combine_first(disp["BB4_prev"]),
    })

    out = disp_out.merge(m, on="_KEY", how="left")
    # Add NEG center shift (estimated) for each run when columns exist
    if ("LFC_NEG_vs_DEL2_used_cur" in out.columns) and ("LFC_NEG_centered_cur" in out.columns):
        out["NEG_center_shift_cur"] = pd.to_numeric(out["LFC_NEG_vs_DEL2_used_cur"], errors="coerce") - pd.to_numeric(out["LFC_NEG_centered_cur"], errors="coerce")
    elif ("LFC_NEG_vs_DEL2_cur" in out.columns) and ("LFC_NEG_centered_cur" in out.columns):
        out["NEG_center_shift_cur"] = pd.to_numeric(out["LFC_NEG_vs_DEL2_cur"], errors="coerce") - pd.to_numeric(out["LFC_NEG_centered_cur"], errors="coerce")
    if ("LFC_NEG_vs_DEL2_used_prev" in out.columns) and ("LFC_NEG_centered_prev" in out.columns):
        out["NEG_center_shift_prev"] = pd.to_numeric(out["LFC_NEG_vs_DEL2_used_prev"], errors="coerce") - pd.to_numeric(out["LFC_NEG_centered_prev"], errors="coerce")
    elif ("LFC_NEG_vs_DEL2_prev" in out.columns) and ("LFC_NEG_centered_prev" in out.columns):
        out["NEG_center_shift_prev"] = pd.to_numeric(out["LFC_NEG_vs_DEL2_prev"], errors="coerce") - pd.to_numeric(out["LFC_NEG_centered_prev"], errors="coerce")
    if ("NEG_center_shift_cur" in out.columns) and ("NEG_center_shift_prev" in out.columns):
        out["NEG_center_shift_diff"] = pd.to_numeric(out["NEG_center_shift_cur"], errors="coerce") - pd.to_numeric(out["NEG_center_shift_prev"], errors="coerce")
    out = out.drop(columns=["_KEY"])

    if score == "_SCORE":
        out["Score_col_cur"] = score_cur_col
        out["Score_col_prev"] = score_prev_col
        out = out.rename(columns={
            "_SCORE_cur": "Score_cur",
            "_SCORE_prev": "Score_prev",
            "_SCORE_diff": "Score_diff",
        })
        score_cur_out = "Score_cur"
        score_diff_out = "Score_diff"
    else:
        score_cur_out = f"{score}_cur"
        score_diff_out = f"{score}_diff"

    out = out.sort_values(["status", "rank_cur", "rank_prev"]).reset_index(drop=True)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False, na_rep="NA")
    print(f"[OK] wrote {out_path}")

    top_df_diff = None
    top_df_cur = None

    if args.top_summary_out:
        if score_diff_out in out.columns:
            top_df = out.reindex(out[score_diff_out].abs().sort_values(ascending=False).index)
        else:
            top_df = out.copy()
        top_n = max(1, int(args.top_summary_n))
        top_df = top_df.head(top_n)
        top_df_diff = top_df
        top_path = Path(args.top_summary_out)
        top_path.parent.mkdir(parents=True, exist_ok=True)
        top_df.to_csv(top_path, sep="\t", index=False, na_rep="NA")
        print(f"[OK] wrote {top_path}")

    if args.top_summary_out_current:
        if score_cur_out in out.columns:
            top_df = out.reindex(pd.to_numeric(out[score_cur_out], errors="coerce").sort_values(ascending=False).index)
        else:
            top_df = out.copy()
        top_n = max(1, int(args.top_summary_n))
        top_df = top_df.head(top_n)
        top_df_cur = top_df
        top_path = Path(args.top_summary_out_current)
        top_path.parent.mkdir(parents=True, exist_ok=True)
        top_df.to_csv(top_path, sep="\t", index=False, na_rep="NA")
        print(f"[OK] wrote {top_path}")

    if args.top_summary_xlsx:
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
        xlsx_path = Path(args.top_summary_xlsx)
        xlsx_path.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(xlsx_path, engine=engine) as writer:
            out.to_excel(writer, sheet_name="Full", index=False)
            if top_df_diff is not None:
                top_df_diff.to_excel(writer, sheet_name="Top_By_Diff", index=False)
            if top_df_cur is not None:
                top_df_cur.to_excel(writer, sheet_name="Top_By_Current", index=False)
        print(f"[OK] wrote {xlsx_path}")


if __name__ == "__main__":
    main()
