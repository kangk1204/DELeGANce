#!/usr/bin/env python3
import argparse
import csv
import random
import re
from pathlib import Path
from typing import Dict, List, Tuple


LIB_RE = re.compile(r"(?:^|[_])LIB[A-Za-z0-9]+", re.IGNORECASE)


def has_lib_suffix(text: str) -> bool:
    return bool(LIB_RE.search(text or ""))


def reservoir_update(sample: List[str], rng: random.Random, i: int, value: str, n: int) -> None:
    if len(sample) < n:
        sample.append(value)
        return
    j = rng.randint(1, i)
    if j <= n:
        sample[j - 1] = value


def analyze_file(path: Path, sample_size: int, rng: random.Random) -> Tuple[Dict[str, int], Dict[str, List[str]]]:
    counts = {
        "total": 0,
        "id_has_lib": 0,
        "id_missing_lib": 0,
        "bb_has_lib": 0,
        "bb_missing_lib": 0,
        "id_missing_bb_has_lib": 0,
        "id_missing_bb_missing_lib": 0,
        "id_has_lib_bb_missing_lib": 0,
        "lib_id_present": 0,
        "id_missing_lib_with_lib_id": 0,
    }
    examples = {
        "id_missing_bb_has_lib": [],
        "id_missing_bb_missing_lib": [],
        "id_has_lib_bb_missing_lib": [],
    }
    idx = {k: 0 for k in examples.keys()}

    with path.open("r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            counts["total"] += 1
            rid = row.get("id", "")
            lib_id = (row.get("lib_id", "") or "").strip()
            if lib_id:
                counts["lib_id_present"] += 1

            id_has = has_lib_suffix(rid)
            bb_has = any(
                has_lib_suffix(row.get(k, ""))
                for k in ("C1_bb_id", "C2_bb_id", "C3_bb_id", "C4_bb_id")
            )

            if id_has:
                counts["id_has_lib"] += 1
            else:
                counts["id_missing_lib"] += 1
                if lib_id:
                    counts["id_missing_lib_with_lib_id"] += 1

            if bb_has:
                counts["bb_has_lib"] += 1
            else:
                counts["bb_missing_lib"] += 1

            if not id_has and bb_has:
                counts["id_missing_bb_has_lib"] += 1
                idx["id_missing_bb_has_lib"] += 1
                reservoir_update(examples["id_missing_bb_has_lib"], rng, idx["id_missing_bb_has_lib"], rid, sample_size)
            elif not id_has and not bb_has:
                counts["id_missing_bb_missing_lib"] += 1
                idx["id_missing_bb_missing_lib"] += 1
                reservoir_update(examples["id_missing_bb_missing_lib"], rng, idx["id_missing_bb_missing_lib"], rid, sample_size)
            elif id_has and not bb_has:
                counts["id_has_lib_bb_missing_lib"] += 1
                idx["id_has_lib_bb_missing_lib"] += 1
                reservoir_update(examples["id_has_lib_bb_missing_lib"], rng, idx["id_has_lib_bb_missing_lib"], rid, sample_size)

    return counts, examples


def pct(num: int, denom: int) -> float:
    if denom <= 0:
        return 0.0
    return 100.0 * float(num) / float(denom)


def main() -> None:
    p = argparse.ArgumentParser(description="Analyze lib suffix missing patterns in decoded_reads_*")
    p.add_argument("--run_root", action="append", required=True)
    p.add_argument("--out_tsv", default="DELeGANce_out/lib_suffix_missing_report.tsv")
    p.add_argument("--out_examples", default="DELeGANce_out/lib_suffix_missing_examples.txt")
    p.add_argument("--sample_size", type=int, default=8)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    out_rows = []
    example_lines = []

    rng = random.Random(args.seed)

    for run_root in [Path(r) for r in args.run_root]:
        dec_dir = run_root / "02_decoded"
        files = sorted(dec_dir.glob("decoded_reads_*.tsv"))
        if not files:
            raise SystemExit(f"[ERROR] No decoded_reads_*.tsv under {dec_dir}")

        for f in files:
            sample = f.stem.replace("decoded_reads_", "")
            counts, examples = analyze_file(f, args.sample_size, rng)

            total = counts["total"]
            row = {
                "run": run_root.name,
                "sample": sample,
                "total_reads": total,
                "id_has_lib": counts["id_has_lib"],
                "id_missing_lib": counts["id_missing_lib"],
                "id_missing_lib_pct": pct(counts["id_missing_lib"], total),
                "bb_has_lib": counts["bb_has_lib"],
                "bb_missing_lib": counts["bb_missing_lib"],
                "id_missing_bb_has_lib": counts["id_missing_bb_has_lib"],
                "id_missing_bb_missing_lib": counts["id_missing_bb_missing_lib"],
                "id_has_lib_bb_missing_lib": counts["id_has_lib_bb_missing_lib"],
                "lib_id_present": counts["lib_id_present"],
                "id_missing_lib_with_lib_id": counts["id_missing_lib_with_lib_id"],
            }
            out_rows.append(row)

            example_lines.append(f"=== {run_root.name} | {sample} ===")
            for key, vals in examples.items():
                example_lines.append(f"- {key}: {len(vals)} example(s)")
                for v in vals:
                    example_lines.append(f"  {v}")
            example_lines.append("")

    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(out_path, sep="\t", index=False)

    ex_path = Path(args.out_examples)
    ex_path.parent.mkdir(parents=True, exist_ok=True)
    ex_path.write_text("\n".join(example_lines), encoding="utf-8")

    print(f"[OK] wrote {out_path}")
    print(f"[OK] wrote {ex_path}")


if __name__ == "__main__":
    import pandas as pd  # local import to keep top minimal
    main()
