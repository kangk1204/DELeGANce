#!/usr/bin/env python3
"""
Local-only helper: apply or reverse anonymized names using local_target_map.tsv.
This script is safe to keep local; the mapping file is gitignored.

Usage:
  python3 anonymize_with_map.py --map local_target_map.tsv --mode anonymize --root .
  python3 anonymize_with_map.py --map local_target_map.tsv --mode deanonymize --root .

Notes:
- Only renames paths that exist.
- Renames are performed in a safe order (longer names first, deepest paths first).
- .git/, .venv/, __pycache__/ and the mapping file itself are never renamed.
- Use --dry-run to preview the renames without touching the filesystem (paths are printed as they
  are at preview time, i.e. before any parent directory would have been renamed).
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


_SKIP_DIRS = {".git", ".venv", "venv", ".conda", "__pycache__", ".idea", ".vscode"}


def _skip(path: Path, root: Path, protected: set) -> bool:
    rel_parts = path.relative_to(root).parts
    if any(part in _SKIP_DIRS for part in rel_parts):
        return True
    return path.resolve() in protected


def rename_paths(root: Path, pairs, protected=None, dry_run: bool = False):
    renamed = []
    failed = []
    protected = protected or set()
    for src, dst in pairs:
        # Bottom-up (deepest first): renaming a parent directory before its children invalidated the
        # pre-collected child paths (FileNotFoundError, tree left half-renamed). rglob is re-run per pair.
        paths = sorted(root.rglob("*"), key=lambda p: (len(p.parts), str(p)), reverse=True)
        for path in paths:
            if src not in path.name or _skip(path, root, protected):
                continue
            new_name = path.name.replace(src, dst)
            if new_name == path.name:
                continue
            target = path.with_name(new_name)
            if target.exists():
                print(f"[SKIP] target exists: {target}")
                continue
            if dry_run:
                print(f"[DRY] {path} -> {target}")
                renamed.append((path, target))
                continue
            try:
                path.rename(target)
            except OSError as e:
                failed.append((path, target, str(e)))
                print(f"[FAIL] {path} -> {target}: {e}")
                continue
            renamed.append((path, target))
    return renamed, failed


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--map", default="local_target_map.tsv")
    p.add_argument("--mode", choices=["anonymize", "deanonymize"], required=True)
    p.add_argument("--root", default=".")
    p.add_argument("--dry-run", action="store_true", help="Print planned renames without applying them")
    args = p.parse_args()

    map_path = Path(args.map)
    if not map_path.exists():
        raise SystemExit(f"[ERROR] mapping file not found: {map_path}")

    rows = read_map(map_path)
    pairs = build_pairs(rows, args.mode)
    if not pairs:
        print("[INFO] No mappings to apply.")
        return

    # never rename the mapping file itself, nor its JSON twin (both are gitignored local-only files)
    protected = {map_path.resolve(), map_path.with_suffix(".json").resolve(), map_path.with_suffix(".tsv").resolve()}
    renamed, failed = rename_paths(Path(args.root), pairs, protected=protected, dry_run=args.dry_run)
    verb = "would rename" if args.dry_run else "renamed"
    if args.dry_run and renamed:
        print("[DRY] note: child paths are shown under their current (not yet renamed) parent directories")
    print(f"[OK] {verb} {len(renamed)} path(s); {len(failed)} failure(s).")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
