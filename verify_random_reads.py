#!/usr/bin/env python3
import argparse
import csv
import random
import re
from collections import Counter
from typing import Dict, List, Tuple

import pandas as pd


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def parse_id_fields(id_str: str) -> Tuple[int, str, str, str, str]:
    s = str(id_str or "").strip()
    if s == "":
        return 3, "NA", "NA", "NA", "NA"
    raw = [t for t in re.split(r"[\|_,:;/\s]+", s) if t != ""]

    cyc = None
    if raw and re.fullmatch(r"\d+", raw[0]):
        try:
            cyc = int(raw[0])
            raw = raw[1:]
        except Exception:
            cyc = None

    parts: List[str] = []
    i = 0
    while i < len(raw) and len(parts) < 4:
        t = raw[i]
        if t == "":
            i += 1
            continue
        # Same guard as 03_call_hits.parse_id_fields: a token starting with "LIB" is a BB id,
        # not a namespace fragment to glue onto the previous token.
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
    bb1, bb2, bb3, bb4 = parts[:4]
    if cyc is None:
        cyc = 4 if str(bb4) != "NA" else 3
    return int(cyc), str(bb1), str(bb2), str(bb3), str(bb4)


def extract_lib_from_bb(bb: str) -> str:
    if not bb or bb == "NA":
        return ""
    m = re.search(r"_LIB([A-Za-z0-9\.-]+)$", bb)
    return m.group(1) if m else ""


def safe_int(v):
    try:
        return int(v)
    except Exception:
        return None


def seq_match(read_seq: str, direction: str, pos, length, expected: str) -> bool:
    if expected in (None, "", "NA"):
        return True
    ipos = safe_int(pos)
    ilen = safe_int(length)
    if ipos is None or ilen is None:
        return True
    oriented = revcomp(read_seq) if str(direction).lower() == "revcomp" else read_seq
    sub = oriented[ipos:ipos + ilen]
    return sub == expected


def reservoir_update(sample: List[Dict[str, str]], rng: random.Random, i: int, row: Dict[str, str], n: int) -> None:
    if len(sample) < n:
        sample.append(row)
        return
    j = rng.randint(1, i)
    if j <= n:
        sample[j - 1] = row


def analyze_samples(sample_rows: List[Dict[str, str]]) -> Tuple[Counter, List[Tuple], List[str]]:
    stats = Counter()
    mismatch_examples: List[Tuple] = []
    ids: List[str] = []

    for row in sample_rows:
        rid = row.get("id", "")
        lib_id = row.get("lib_id", "")
        cycles = safe_int(row.get("cycles", "")) or 0
        bb1 = row.get("C1_bb_id", "NA")
        bb2 = row.get("C2_bb_id", "NA")
        bb3 = row.get("C3_bb_id", "NA")
        bb4 = row.get("C4_bb_id", "NA")

        ids.append(rid)

        cyc_p, p1, p2, p3, p4 = parse_id_fields(rid)
        if cycles and cyc_p != cycles:
            stats["cycles_mismatch"] += 1
            if len(mismatch_examples) < 5:
                mismatch_examples.append(("cycles", rid, cycles, cyc_p))

        if p1 != bb1 or p2 != bb2 or p3 != bb3 or (p4 != bb4 and p4 != "NA"):
            stats["bb_id_mismatch"] += 1
            if len(mismatch_examples) < 5:
                mismatch_examples.append(("bb_id", rid, (bb1, bb2, bb3, bb4), (p1, p2, p3, p4)))

        libs = {extract_lib_from_bb(x) for x in (p1, p2, p3, p4) if extract_lib_from_bb(x)}
        if libs:
            if len(libs) != 1 or (lib_id and lib_id not in libs):
                stats["lib_id_mismatch"] += 1
                if len(mismatch_examples) < 5:
                    mismatch_examples.append(("lib_id", rid, lib_id, list(libs)))
        else:
            stats["lib_suffix_missing"] += 1

        read_seq = row.get("read_seq", "")
        direction = row.get("direction", "")
        for i in (1, 2, 3, 4):
            if not seq_match(
                read_seq,
                direction,
                row.get(f"C{i}_pos", None),
                row.get(f"C{i}_len", None),
                row.get(f"C{i}_seq", None),
            ):
                stats[f"seq_mismatch_C{i}"] += 1
                if len(mismatch_examples) < 5:
                    mismatch_examples.append((f"seq_C{i}", rid, row.get(f"C{i}_seq"), row.get(f"C{i}_pos")))

    return stats, mismatch_examples, ids


def parse_seeds(args) -> List[int]:
    if args.seeds:
        out = []
        for s in re.split(r"[\s,]+", args.seeds.strip()):
            if not s:
                continue
            out.append(int(s))
        return out
    return [args.seed]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run_root", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--n_reads", type=int, default=1000)
    p.add_argument("--n_ids", type=int, default=50)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--seeds", default="")
    args = p.parse_args()

    decoded_path = f"{args.run_root}/02_decoded/decoded_reads_{args.sample}.tsv"
    matrix_path = f"{args.run_root}/02_decoded/raw_counts_matrix.tsv"

    seeds = parse_seeds(args)
    rngs = {seed: random.Random(seed) for seed in seeds}
    samples = {seed: [] for seed in seeds}

    counts_all = Counter()  # keyed by (lib_id, id) — raw_counts_matrix rows are (lib_id, id) pairs
    total_reads = 0

    with open(decoded_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader, 1):
            rid = row.get("id", "")
            if rid:
                counts_all[(row.get("lib_id", "") or "", rid)] += 1
            total_reads += 1
            for seed in seeds:
                reservoir_update(samples[seed], rngs[seed], i, row, args.n_reads)

    # raw_counts_matrix.tsv is keyed by (lib_id, id); the same id can occur under several lib_ids,
    # so a plain id index would return a Series (int() TypeError) and mis-compare counts.
    df_counts = pd.read_csv(matrix_path, sep="\t", usecols=["lib_id", "id", args.sample],
                            dtype={"lib_id": str, "id": str}, keep_default_na=False)
    matrix_counts = {
        (str(lib), str(rid)): int(v)
        for lib, rid, v in zip(df_counts["lib_id"], df_counts["id"], df_counts[args.sample])
    }

    print("=== Random Read Verification ===")
    print(f"run_root: {args.run_root}")
    print(f"sample: {args.sample}")
    print(f"decoded_path: {decoded_path}")
    print(f"matrix_path: {matrix_path}")
    print(f"seeds: {seeds}")
    print(f"sampled_reads_per_seed: {args.n_reads}")
    print(f"count_check_ids_per_seed: {args.n_ids}")
    print(f"total_reads_scanned: {total_reads}")

    for seed in seeds:
        stats, mismatch_examples, ids = analyze_samples(samples[seed])
        # pair each sampled id with its lib_id (analyze_samples preserves row order)
        keys = [(row.get("lib_id", "") or "", rid) for row, rid in zip(samples[seed], ids)]
        uniq_ids = list(dict.fromkeys(keys))
        rngc = random.Random(seed + 100003)
        rngc.shuffle(uniq_ids)
        check_ids = uniq_ids[: min(args.n_ids, len(uniq_ids))]

        count_mismatch = 0
        count_examples = []
        for key in check_ids:
            expected = matrix_counts.get(key)
            got = int(counts_all.get(key, 0))
            if expected is None or expected != got:
                count_mismatch += 1
                if len(count_examples) < 5:
                    count_examples.append((key, expected, got))

        print(f"\n[seed={seed}] checks:")
        print(f"  cycles_mismatch: {stats.get('cycles_mismatch', 0)}")
        print(f"  bb_id_mismatch: {stats.get('bb_id_mismatch', 0)}")
        print(f"  lib_id_mismatch: {stats.get('lib_id_mismatch', 0)}")
        print(f"  lib_suffix_missing: {stats.get('lib_suffix_missing', 0)}")
        for i in (1, 2, 3, 4):
            print(f"  seq_mismatch_C{i}: {stats.get(f'seq_mismatch_C{i}', 0)}")
        print(f"  count_check_ids: {len(check_ids)}")
        print(f"  count_mismatch: {count_mismatch}")
        if mismatch_examples:
            print("\n  examples (decode/parse):")
            for ex in mismatch_examples:
                print("   ", ex)
        if count_examples:
            print("\n  examples (count mismatches):")
            for ex in count_examples:
                print("   ", ex)


if __name__ == "__main__":
    main()
