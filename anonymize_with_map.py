#!/usr/bin/env python3
"""
Local-only helper: apply or reverse anonymized names using local_target_map.tsv.
This script is safe to keep local; the mapping file is gitignored.

Usage:
  python3 anonymize_with_map.py --map local_target_map.tsv --mode anonymize --root .
  python3 anonymize_with_map.py --map local_target_map.tsv --mode deanonymize --root .

Notes:
- Only renames paths that exist.
- Renames are performed in a safe order (longer names first).
"""
import argparse
from pathlib import Path


def read_map(path: Path):
    rows = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.strip() or line.strip().startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        category, anon_name, real_name = (p.strip() for p in parts[:3])
        if not anon_name or not real_name:
            continue
        rows.append((category, anon_name, real_name))
    return rows


def build_pairs(rows, mode: str):
    pairs = []
    for _, anon, real in rows:
        if mode == "anonymize":
            pairs.append((real, anon))
        else:
            pairs.append((anon, real))
    # rename longer names first to avoid partial replacement
    pairs.sort(key=lambda x: len(x[0]), reverse=True)
    return pairs


def rename_paths(root: Path, pairs):
    renamed = []
    for src, dst in pairs:
        for path in sorted(root.rglob("*")):
            if src not in path.name:
                continue
            new_name = path.name.replace(src, dst)
            if new_name == path.name:
                continue
            target = path.with_name(new_name)
            if target.exists():
                continue
            path.rename(target)
            renamed.append((path, target))
    return renamed


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--map", default="local_target_map.tsv")
    p.add_argument("--mode", choices=["anonymize", "deanonymize"], required=True)
    p.add_argument("--root", default=".")
    args = p.parse_args()

    map_path = Path(args.map)
    if not map_path.exists():
        raise SystemExit(f"[ERROR] mapping file not found: {map_path}")

    rows = read_map(map_path)
    pairs = build_pairs(rows, args.mode)
    if not pairs:
        print("[INFO] No mappings to apply.")
        return

    renamed = rename_paths(Path(args.root), pairs)
    print(f"[OK] renamed {len(renamed)} path(s).")


if __name__ == "__main__":
    main()
