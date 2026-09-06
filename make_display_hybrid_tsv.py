#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd


# Anchored: strip only a trailing "_LIB<lib>" namespace token from a single BB value
LIB_SUFFIX_RE = re.compile(r"_LIB[\w.-]+$")
# Fallback for full tag IDs (cycles_BB1_BB2_BB3[_BB4]) when no lib_id is known. It is ambiguous for
# lib_ids containing "_" (01_preprocess allows [A-Za-z0-9_.-]), so strip_lib_anywhere() removes exact
# "_LIB<lib_id>" tokens whenever the run's lib_ids are known and uses this regex only otherwise.
LIB_ANY_RE = re.compile(r"_LIB[^_]+(?=_|$)")
NA_LIB_VALUES = {"", "na", "nan", "none", "<na>"}


def strip_lib_suffix(value) -> str:
    if value is None:
        return "NA"
    s = str(value)
    if s in ("", "NA", "nan", "None"):
        return "NA"
    return LIB_SUFFIX_RE.sub("", s)


def known_libs_from(series) -> Tuple[str, ...]:
    """Distinct real lib_ids of a LIB_ID column (NA-like values dropped), longest first."""
    if series is None:
        return ()
    vals = {str(v).strip() for v in series.dropna().unique()}
    vals = {v for v in vals if v.lower() not in NA_LIB_VALUES}
    return tuple(sorted(vals, key=len, reverse=True))


def strip_lib_anywhere(value, known_libs: Tuple[str, ...] = ()) -> str:
    """Remove "_LIB<lib>" tokens from a full tag ID: exact tokens for the known lib_ids (followed by "_" or
    end of string); the generic LIB_ANY_RE only when no lib_id is known."""
    if value is None:
        return ""
    s = str(value)
    if known_libs:
        for lib in known_libs:
            s = re.sub(r"_LIB" + re.escape(lib) + r"(?=_|$)", "", s)
        return s
    return LIB_ANY_RE.sub("", s)


def pick_col(cols: List[str], candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in cols:
            return c
    return None


def pick_bb_col(cols: List[str], base: str) -> Optional[str]:
    for c in (base, f"{base}_x", f"{base}_y", f"{base.lower()}_id"):
        if c in cols:
            return c
    return None


def normalize_display(df: pd.DataFrame) -> pd.DataFrame:
    cols = list(df.columns)
    lib_col = pick_col(cols, ["LIB_ID", "LIB_ID_x", "LIB_ID_y", "LibID", "lib_id", "lib_id_x", "lib_id_y"])
    id_cols = [c for c in cols if c in ("ID", "id", "ID_x", "ID_y", "id_x", "id_y")]

    bb_cols = []
    for c in cols:
        if re.fullmatch(r"BB[1-4](_[xy])?", c) or re.fullmatch(r"bb[1-4]_id", c):
            bb_cols.append(c)

    # Strip LIB suffix from BB columns (in place; values below reuse the stripped columns)
    for c in bb_cols:
        df[c] = df[c].map(strip_lib_suffix)

    b1_col = pick_bb_col(cols, "BB1")
    b2_col = pick_bb_col(cols, "BB2")
    b3_col = pick_bb_col(cols, "BB3")
    b4_col = pick_bb_col(cols, "BB4")

    def _bb(col: Optional[str]) -> pd.Series:
        # Missing BB column -> Series of "NA" (a bare str has no .astype)
        return df[col] if col else pd.Series(["NA"] * len(df), index=df.index)

    b1, b2, b3, b4 = _bb(b1_col), _bb(b2_col), _bb(b3_col), _bb(b4_col)

    if lib_col:
        libs = df[lib_col].fillna("").astype(str).str.strip()
        id_display = libs + "_" + b1.astype(str) + "_" + b2.astype(str) + "_" + b3.astype(str) + "_" + b4.astype(str)
        valid = ~libs.str.lower().isin(sorted(NA_LIB_VALUES))
        if id_cols:
            # rows without a lib_id: strip only the exact "_LIB<lib_id>" tokens of the libs seen in this table
            known = known_libs_from(df[lib_col])
            fallback = df[id_cols[0]].astype(str).map(lambda v: strip_lib_anywhere(v, known))
            id_display = id_display.where(valid, fallback)
        else:
            id_display = id_display.where(valid, "")
    else:
        # no LIB_ID column at all: no lib_id is known, generic regex fallback
        id_display = df[id_cols[0]].astype(str).map(strip_lib_anywhere) if id_cols else ""

    for c in id_cols:
        df[c] = id_display
    for c in ("DisplayID", "CombinedID"):
        if c in df.columns:
            df[c] = id_display

    return df


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
        "Specify --in_tsv or --prefer_dir.\n" + msg
    )


def process_one(in_tsv: Path, out_tsv: Path) -> None:
    df = pd.read_csv(in_tsv, sep="\t", low_memory=False)
    df = normalize_display(df)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False, na_rep="NA")
    print(f"[OK] wrote {out_tsv}")


def main() -> None:
    p = argparse.ArgumentParser(description="Create display-only hybrid TSV with LIB suffix removed")
    p.add_argument("--run_root", action="append", help="DELeGANce_out/<RUN> root (repeatable)")
    p.add_argument("--prefer_dir", default="",
                   help="Preferred subdir under run_root for 05_hybrid_annot.tsv")
    p.add_argument("--in_tsv", help="Direct input TSV path (optional)")
    p.add_argument("--out_tsv", help="Direct output TSV path (optional)")
    args = p.parse_args()

    if args.in_tsv:
        in_tsv = Path(args.in_tsv)
        if args.out_tsv:
            out_tsv = Path(args.out_tsv)
        else:
            out_tsv = in_tsv.with_name(in_tsv.stem + "_display.tsv")
        process_one(in_tsv, out_tsv)
        return

    if not args.run_root:
        raise SystemExit("[ERROR] --run_root or --in_tsv is required.")

    for run_root_str in args.run_root:
        run_root = Path(run_root_str)
        in_tsv = _resolve_annot(run_root, args.prefer_dir)
        out_tsv = in_tsv.with_name("05_hybrid_annot_display.tsv")
        process_one(in_tsv, out_tsv)


if __name__ == "__main__":
    main()
