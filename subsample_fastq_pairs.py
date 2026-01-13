#!/usr/bin/env python3
import argparse
import gzip
import random
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

READ_RE = re.compile(
    r"^(?P<base>.+)(?P<read>_R?1|_R?2|_1|_2|\.R1|\.R2|\.1|\.2)(?:_[0-9]{3})?(?P<ext>\.f(?:ast)?q(?:\.gz)?)$",
    re.IGNORECASE,
)


def open_maybe_gz(path: Path, mode: str):
    if path.name.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8", errors="replace")


def read_record(handle) -> Optional[Tuple[str, str, str, str]]:
    h = handle.readline()
    if not h:
        return None
    s = handle.readline()
    p = handle.readline()
    q = handle.readline()
    if not q:
        raise ValueError("FASTQ record truncated")
    return h, s, p, q


def detect_pairs(input_dir: Path) -> Dict[str, Dict[str, Path]]:
    pairs: Dict[str, Dict[str, Path]] = {}
    for path in sorted(input_dir.iterdir()):
        if not path.is_file():
            continue
        m = READ_RE.match(path.name)
        if not m:
            continue
        base = m.group("base")
        read = m.group("read")
        side = "R1" if read.endswith("1") else "R2"
        pairs.setdefault(base, {})
        if side in pairs[base]:
            raise ValueError(f"Duplicate {side} for base '{base}': {pairs[base][side]} vs {path}")
        pairs[base][side] = path
    return {k: v for k, v in pairs.items() if "R1" in v and "R2" in v}


def load_pairs_from_tsv(path: Path, input_dir: Path) -> Dict[str, Dict[str, Path]]:
    pairs: Dict[str, Dict[str, Path]] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\t,]", line)
            parts = [p for p in parts if p != ""]
            if len(parts) < 2:
                continue
            if parts[0].lower() in ("sample", "name") and "r1" in parts[1].lower():
                continue
            if len(parts) == 2:
                sample = Path(parts[0]).stem
                r1, r2 = parts
            else:
                sample, r1, r2 = parts[0], parts[1], parts[2]
            r1_path = Path(r1)
            r2_path = Path(r2)
            if not r1_path.is_absolute():
                r1_path = input_dir / r1_path
            if not r2_path.is_absolute():
                r2_path = input_dir / r2_path
            pairs[sample] = {"R1": r1_path, "R2": r2_path}
    return pairs


def build_out_paths(base: str, r1_in: Path, r2_in: Path, out_dir: Path, suffix: str) -> Tuple[Path, Path]:
    m1 = READ_RE.match(r1_in.name)
    m2 = READ_RE.match(r2_in.name)
    if not m1 or not m2:
        raise ValueError(f"Output naming failed for base '{base}'")
    out_r1 = f"{m1.group('base')}{suffix}{m1.group('read')}{m1.group('ext')}"
    out_r2 = f"{m2.group('base')}{suffix}{m2.group('read')}{m2.group('ext')}"
    return out_dir / out_r1, out_dir / out_r2


def subsample_pair(r1_path: Path, r2_path: Path, out_r1: Path, out_r2: Path,
                   n_pairs: int, mode: str, seed: int) -> Tuple[int, int]:
    rng = random.Random(seed)
    total = 0
    kept = 0

    if mode == "head":
        with open_maybe_gz(r1_path, "rt") as r1, open_maybe_gz(r2_path, "rt") as r2, \
                open_maybe_gz(out_r1, "wt") as w1, open_maybe_gz(out_r2, "wt") as w2:
            while kept < n_pairs:
                rec1 = read_record(r1)
                rec2 = read_record(r2)
                if rec1 is None and rec2 is None:
                    break
                if rec1 is None or rec2 is None:
                    raise ValueError("R1/R2 length mismatch")
                w1.writelines(rec1)
                w2.writelines(rec2)
                kept += 1
                total += 1
            for _ in r1:
                break
        return total, kept

    reservoir: List[Tuple[int, Tuple[str, str, str, str], Tuple[str, str, str, str]]] = []
    with open_maybe_gz(r1_path, "rt") as r1, open_maybe_gz(r2_path, "rt") as r2:
        while True:
            rec1 = read_record(r1)
            rec2 = read_record(r2)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise ValueError("R1/R2 length mismatch")
            total += 1
            if len(reservoir) < n_pairs:
                reservoir.append((total, rec1, rec2))
            else:
                j = rng.randint(1, total)
                if j <= n_pairs:
                    reservoir[j - 1] = (total, rec1, rec2)

    reservoir.sort(key=lambda x: x[0])
    kept = len(reservoir)
    with open_maybe_gz(out_r1, "wt") as w1, open_maybe_gz(out_r2, "wt") as w2:
        for _, rec1, rec2 in reservoir:
            w1.writelines(rec1)
            w2.writelines(rec2)
    return total, kept


def main() -> int:
    p = argparse.ArgumentParser(description="Subsample paired-end FASTQ files from a directory")
    p.add_argument("--input-dir", required=True, help="Directory containing paired FASTQs")
    p.add_argument("--output-dir", required=True, help="Output directory for subsampled FASTQs")
    p.add_argument("--n-pairs", type=int, required=True, help="Number of read pairs to keep per pair file")
    p.add_argument("--mode", choices=["random", "head"], default="random", help="Sampling mode")
    p.add_argument("--seed", type=int, default=42, help="Random seed for sampling")
    p.add_argument("--pairs", default="", help="Optional TSV/CSV with sample,r1,r2 paths")
    p.add_argument("--suffix", default="_subsampled", help="Suffix inserted before read token")
    p.add_argument("--overwrite", action="store_true", help="Overwrite output files if present")
    args = p.parse_args()

    in_dir = Path(args.input_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.pairs:
        pairs = load_pairs_from_tsv(Path(args.pairs), in_dir)
    else:
        pairs = detect_pairs(in_dir)
    if not pairs:
        raise SystemExit("[ERROR] No paired FASTQ files found.")

    for base, files in pairs.items():
        r1_path = files["R1"]
        r2_path = files["R2"]
        if not r1_path.exists() or not r2_path.exists():
            raise SystemExit(f"[ERROR] Missing FASTQ: {r1_path} or {r2_path}")
        out_r1, out_r2 = build_out_paths(base, r1_path, r2_path, out_dir, args.suffix)
        if (out_r1.exists() or out_r2.exists()) and not args.overwrite:
            raise SystemExit(f"[ERROR] Output exists: {out_r1} or {out_r2} (use --overwrite)")
        total, kept = subsample_pair(r1_path, r2_path, out_r1, out_r2, args.n_pairs, args.mode, args.seed)
        print(f"[OK] {base}: total_pairs={total}, kept={kept} → {out_r1.name}, {out_r2.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
