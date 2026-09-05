#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import List

import pandas as pd

# Anchored: strip only a trailing "_LIB<lib>" namespace token (same rule as 06/07 _strip_lib_suffix)
LIB_SUFFIX_RE = re.compile(r"_LIB[^_]+$")


# Missing-BB placeholder inside compound_key. NOTE: the reviewed 06/07 sources emit "NA" here;
# the approved fix plan specifies "" — keep this constant in sync with 06/07 when their patch lands.
KEY_NA = "NA"   # must match 06_compare_top_hits.py / 07_tiered_report.py (_make_compound_key uses fillna("NA"))


def strip_lib_suffix(value) -> str:
    s = "" if value is None else str(value).strip()
    if s in ("", "NA", "nan", "None", "<NA>"):
        return KEY_NA
    return LIB_SUFFIX_RE.sub("", s)


def resolve_hybrid_path(base: str, preset: str | None) -> Path:
    p = Path(base)
    if p.is_file():
        return p
    if p.is_dir():
        direct = p / "05_hybrid_annot.tsv"
        if direct.exists():
            return direct
        norm = p / "03_normalized"
        if norm.exists():
            if preset:
                cand = norm / preset / "05_hybrid_annot.tsv"
                if cand.exists():
                    return cand
            # rglob order is filesystem-dependent: require a unique candidate instead of taking the first
            cands = sorted(norm.rglob("05_hybrid_annot.tsv"))
            if len(cands) == 1:
                return cands[0]
            if len(cands) > 1:
                listing = "\n".join(str(c) for c in cands)
                raise FileNotFoundError(
                    f"Multiple 05_hybrid_annot.tsv under {norm}; specify --preset:\n{listing}"
                )
    raise FileNotFoundError(f"05_hybrid_annot.tsv not found under {base}")


def sample_cols_from_header(header: List[str]) -> List[str]:
    cpm_cols = [c for c in header if c.endswith("_CPM") and c[:-4] in header]
    raw_cols = [c[:-4] for c in cpm_cols]
    return raw_cols + cpm_cols


def add_compound_key(df: pd.DataFrame) -> pd.DataFrame:
    """compound_key = BB1|BB2|BB3|BB4 with the LIB namespace suffix removed and NaN -> "" —
    the same definition as 06_compare_top_hits.py / 07_tiered_report.py and the README, so
    merged_all.tsv can be joined with those outputs. (Previously CP_x was included, the LIB
    suffix was kept and NaN became "nan", which made the key incompatible.)"""
    out = df.copy()
    bb_cols = ["BB1_x", "BB2_x", "BB3_x", "BB4_x"]
    for c in bb_cols:
        if c not in out.columns:
            out[c] = pd.NA
    key_parts = out[bb_cols].astype(object).apply(lambda col: col.map(strip_lib_suffix))
    out["compound_key"] = key_parts.agg("|".join, axis=1)
    return out


def reorder_cols(df: pd.DataFrame) -> pd.DataFrame:
    preferred = ["compound_key", "CP_x", "BB1_x", "BB2_x", "BB3_x", "BB4_x", "LIB_ID_x", "ID_x"]
    front = [c for c in preferred if c in df.columns]
    rest = [c for c in df.columns if c not in front]
    return df[front + rest]


def load_full(run_root: str, preset: str | None) -> tuple[pd.DataFrame, List[str]]:
    path = resolve_hybrid_path(run_root, preset)
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = add_compound_key(df)
    df = reorder_cols(df)
    sample_cols = sample_cols_from_header(df.columns.tolist())
    return df, sample_cols


def coalesce_columns(df: pd.DataFrame, base: str, suffix_a: str, suffix_b: str) -> pd.DataFrame:
    a = f"{base}{suffix_a}"
    b = f"{base}{suffix_b}"
    if a in df.columns and b in df.columns:
        series = df[a].combine_first(df[b])
        df[base] = series
        df = df.drop(columns=[a, b])
    elif a in df.columns:
        df = df.rename(columns={a: base})
    elif b in df.columns:
        df = df.rename(columns={b: base})
    return df


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--active-run", required=True)
    ap.add_argument("--inactive-run", required=True)
    ap.add_argument("--preset", default=None)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--active-label", default="active")
    ap.add_argument("--inactive-label", default="inactive")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_active, active_samples = load_full(args.active_run, args.preset)
    df_inactive, inactive_samples = load_full(args.inactive_run, args.preset)

    # Save full per-run TSVs
    active_path = out_dir / f"{args.active_label}_all.tsv"
    inactive_path = out_dir / f"{args.inactive_label}_all.tsv"
    df_active.to_csv(active_path, sep="\t", index=False)
    df_inactive.to_csv(inactive_path, sep="\t", index=False)

    # Merge by compound_key
    df_a = df_active.drop_duplicates("compound_key")
    df_b = df_inactive.drop_duplicates("compound_key")
    for label, full, dedup in ((args.active_label, df_active, df_a), (args.inactive_label, df_inactive, df_b)):
        n_dup = len(full) - len(dedup)
        if n_dup:
            print(f"[WARN] {label}: {n_dup} duplicate compound_key rows dropped before merge (first kept)")
    merged = df_a.merge(df_b, on="compound_key", how="outer", suffixes=(f"_{args.active_label}", f"_{args.inactive_label}"))

    # Coalesce identifier and shared sample columns
    id_cols = ["CP_x", "BB1_x", "BB2_x", "BB3_x", "BB4_x", "LIB_ID_x", "ID_x", "cycles"]
    shared_samples = sorted(set(active_samples) & set(inactive_samples))
    for base in id_cols + shared_samples:
        merged = coalesce_columns(merged, base, f"_{args.active_label}", f"_{args.inactive_label}")

    merged = reorder_cols(merged)
    merged_path = out_dir / "merged_all.tsv"
    merged.to_csv(merged_path, sep="\t", index=False)

    print("[INFO] wrote:")
    print(f"  {active_path}")
    print(f"  {inactive_path}")
    print(f"  {merged_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
