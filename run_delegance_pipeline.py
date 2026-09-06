#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DELeGANce — End‑to‑End Pipeline Runner

Runs the full process in order:
  1) Preprocess FASTQ (01_preprocess_reads.pl)
  2) Decode merged reads (02_decode_reads.pl)
  3) Call hits + reports (03_call_hits.py)

Defaults align with this repository:
  - RUN_ROOT:      DELeGANce_out/<run_name>
  - FASTQ_DIR:     00_original_files
  - BBINFO:        auto‑detect or explicitly pass via --bbinfo
  - All‑in‑one:    --run_root DELeGANce_out/<run_name>

Examples
  # Full run (requires explicit fastq_dir, bbinfo, output_dir, and column names)
  python3 run_delegance_pipeline.py \
      --fastq-dir 00_original_files \
      --bbinfo 00_BB_information.txt \
      --output-dir DELeGANce_out/example_run \
      --r1 R1C1,R1C2,R1C3 \
      --r2 R2C1,R2C2,R2C3 \
      --neg NEG_R1 \
      --del2 DEL2 \
      --threads 6 --mismatch hp_op_cp

  # Hit stage only on an existing run (still requires column names)
  python3 run_delegance_pipeline.py \
      --only hit \
      --output-dir DELeGANce_out/example_run \
      --r1 R1C1,R1C2,R1C3 --r2 R2C1,R2C2,R2C3 --neg NEG_R1 --del2 DEL2
"""

import os
import sys
import argparse
import shutil
import subprocess
import datetime as _dt
import json
import hashlib
import glob
import re
from html import escape as _html_escape

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
MERGED_FASTQ_RE = re.compile(r'(?:\.fpmerged\.fq(?:\.gz)?|_merged\.(?:fq|fastq)(?:\.gz)?)$', re.IGNORECASE)


def _timestamp() -> str:
    return _dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def ensure_dir(p: str):
    if p:
        os.makedirs(p, exist_ok=True)


def which(prog: str) -> str:
    return shutil.which(prog) or ""


def run_cmd(cmd, log_file: str, env=None) -> int:
    """Run a command, streaming stdout/stderr to both console and a log file."""
    print(f"[CMD] {' '.join(map(str, cmd))}")
    ensure_dir(os.path.dirname(log_file))
    with open(log_file, 'a', encoding='utf-8', errors='replace') as lf:
        lf.write(f"\n[{_timestamp()}] RUN: {' '.join(map(str, cmd))}\n")
        lf.flush()
        # errors='replace': a single non-UTF-8 byte from a child must not kill the orchestrator mid-stage
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                encoding='utf-8', errors='replace', env=env)
        for line in proc.stdout:
            sys.stdout.write(line)
            lf.write(line)
        proc.wait()
        rc = proc.returncode
        lf.write(f"[{_timestamp()}] EXIT CODE: {rc}\n")
        lf.flush()
    return rc


def _has_interactive_deps() -> bool:
    try:
        import bokeh  # noqa: F401
        return True
    except Exception:
        return False


def _file_or_gz_exists(path: str) -> bool:
    return os.path.isfile(path) or os.path.isfile(path + ".gz")


def _log_done_after_last_start(path: str, start_markers, done_marker: str) -> bool:
    """True only if the LAST run recorded in an append-mode stage log completed.

    The Perl stages append to their logs, so an 'All done.' line from an earlier successful run
    survives a later failed rerun. We therefore require the last `done_marker` line to appear
    AFTER the last run-start line (any of `start_markers`). Logs without any start line (very old
    runs) fall back to plain presence of the done marker.
    """
    last_start = -1
    last_done = -1
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            for i, line in enumerate(f):
                if any(m in line for m in start_markers):
                    last_start = i
                if done_marker in line:
                    last_done = i
    except Exception:
        return False
    if last_done < 0:
        return False
    if last_start < 0:
        return True
    return last_done > last_start


def _env_int(name: str, default: int) -> int:
    """Read an integer environment variable; fall back to default on missing/invalid values."""
    try:
        return int(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default


# Stage completion markers. 01/02 Perl scripts open their logs in append mode, so an
# 'All done.' line may be left over from an earlier successful run. Prefer explicit marker
# files (written by the stage on success); accept the log marker only as legacy fallback.
PREPROCESS_DONE_MARKER = os.path.join('01_fastp_out', '.preprocess_done')
DECODE_DONE_MARKER = os.path.join('02_decoded', '.decode_done')
# First log line each Perl stage writes at run start (present in legacy versions as well).
PREPROCESS_START_MARKERS = ("Starting preprocess.",)
DECODE_START_MARKERS = ("Merged dir = ", "Effective config:")


def _stage_done(marker_path: str, log_path: str, start_markers) -> bool:
    """True if the stage marker file exists, or (legacy runs) the log's LAST run ends with 'All done.'."""
    if os.path.isfile(marker_path):
        return True
    return _log_done_after_last_start(log_path, start_markers, "All done.")


def _clear_marker(marker_path: str):
    try:
        if os.path.isfile(marker_path):
            os.remove(marker_path)
    except Exception as e:
        print(f"[WARN] Failed to remove stale marker {marker_path}: {e}")


def _write_marker(marker_path: str):
    """Write the completion marker (ISO timestamp). Idempotent; the stage itself may also write it."""
    try:
        ensure_dir(os.path.dirname(marker_path))
        with open(marker_path, 'w', encoding='utf-8') as f:
            f.write(_dt.datetime.now().isoformat(timespec='seconds') + "\n")
    except Exception as e:
        print(f"[WARN] Failed to write marker {marker_path}: {e}")


def _list_merged_fastqs(dir_path: str):
    if not os.path.isdir(dir_path):
        return []
    try:
        items = os.listdir(dir_path)
    except Exception:
        return []
    return [s for s in items if MERGED_FASTQ_RE.search(s)]


def _has_fastp_outputs(run_root: str) -> bool:
    """Return True when merged FASTQs and fixed BB info already exist."""
    merged_dir = os.path.join(run_root, '01_fastp_out')
    fixed_bb  = os.path.join(run_root, 'BB_information_fixed.tsv')
    log_path = os.path.join(run_root, '01_preprocess_reads.log')
    marker = os.path.join(run_root, PREPROCESS_DONE_MARKER)
    merged = _list_merged_fastqs(merged_dir)
    return bool(merged) and _file_or_gz_exists(fixed_bb) and _stage_done(marker, log_path, PREPROCESS_START_MARKERS)


def _has_merged_fastqs(run_root: str) -> bool:
    """Return True when merged FASTQs exist (regardless of BB metadata)."""
    merged_dir = os.path.join(run_root, '01_fastp_out')
    return bool(_list_merged_fastqs(merged_dir))


def _has_decoded_outputs(run_root: str) -> bool:
    """Return True when raw_counts_matrix.tsv exists (required for hit calling)."""
    decoded_dir = os.path.join(run_root, '02_decoded')
    if not os.path.isdir(decoded_dir):
        return False
    raw = os.path.join(decoded_dir, 'raw_counts_matrix.tsv')
    log_path = os.path.join(decoded_dir, '02_decode_reads.log')
    marker = os.path.join(run_root, DECODE_DONE_MARKER)
    return _file_or_gz_exists(raw) and _stage_done(marker, log_path, DECODE_START_MARKERS)


def _pick_decoded_matrix(run_root: str) -> str:
    """Pick the best available decoded matrix for sizing; prefer raw counts, fallback to scaled."""
    decoded_dir = os.path.join(run_root, '02_decoded')
    for name in ('raw_counts_matrix.tsv', 'scaled_counts_matrix.tsv'):
        p = os.path.join(decoded_dir, name)
        if os.path.isfile(p):
            return p
        if os.path.isfile(p + '.gz'):
            return p + '.gz'
    return os.path.join(decoded_dir, 'raw_counts_matrix.tsv')


_ESTIMATE_ROWS_CACHE = {}


def _count_lines_gz(path: str) -> int:
    """Count lines of a gzipped file via `gzip -cd | wc -l` (no shell), falling back to Python."""
    import gzip
    if shutil.which('gzip') and shutil.which('wc'):
        gz = subprocess.Popen(['gzip', '-cd', path], stdout=subprocess.PIPE)
        try:
            out = subprocess.check_output(['wc', '-l'], stdin=gz.stdout, text=True)
        finally:
            gz.stdout.close()
            gz.wait()
        if gz.returncode == 0:
            return int(out.strip().split()[0])
    cnt = 0
    with gzip.open(path, 'rb') as f:
        for _ in f:
            cnt += 1
    return cnt


def _estimate_rows(path: str) -> int:
    """Heuristically estimate row count for (possibly gzipped) TSV, ignoring header.

    Result is memoised per path: a multi-GB matrix must not be re-scanned for every
    auto-optimisation question asked in one run. Returns -1 when the file does not exist, or
    (with a warning) when counting fails.
    """
    if path in _ESTIMATE_ROWS_CACHE:
        return _ESTIMATE_ROWS_CACHE[path]
    n = -1
    try:
        target = None
        if os.path.isfile(path):
            target = path
        elif not path.endswith('.gz') and os.path.isfile(path + '.gz'):
            target = path + '.gz'
        if target is not None:
            if target.endswith('.gz'):
                n = max(0, _count_lines_gz(target) - 1)
            elif shutil.which('wc'):
                out = subprocess.check_output(['wc', '-l', target], text=True).strip().split()[0]
                n = max(0, int(out) - 1)
            else:
                # Fallback: Python loop (may be slow on huge files)
                cnt = 0
                with open(target, 'rb') as f:
                    for _ in f:
                        cnt += 1
                n = max(0, cnt - 1)
    except Exception as e:
        print(f"[WARN] Could not estimate row count of {path}: {e} (auto-opt falls back to defaults)")
        n = -1
    _ESTIMATE_ROWS_CACHE[path] = n
    return n


def _find_existing_hybrid_dir(norm_root: str) -> str:
    """Return the most recent directory under 03_normalized that has 05_hybrid_annot.tsv, else ''. """
    best_dir = ''
    best_mt = -1.0
    # Base (legacy) file
    base = os.path.join(norm_root, '05_hybrid_annot.tsv')
    if os.path.isfile(base):
        best_dir = norm_root; best_mt = os.path.getmtime(base)
    try:
        for name in os.listdir(norm_root):
            p = os.path.join(norm_root, name)
            if not os.path.isdir(p):
                continue
            f = os.path.join(p, '05_hybrid_annot.tsv')
            if os.path.isfile(f):
                mt = os.path.getmtime(f)
                if mt > best_mt:
                    best_mt = mt; best_dir = p
    except Exception:
        pass
    return best_dir


def _jsonable(o):
    if isinstance(o, (str, int, float)) or o is None:
        return o
    if isinstance(o, (list, tuple)):
        return [_jsonable(x) for x in o]
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    return str(o)


# Only these CLI knobs influence the hit-caller output. Hashing all of vars(args) made the
# cache miss on --only/--force-hit/--threads/--fastq-dir changes that do not affect results.
HIT_CACHE_ARG_KEYS = (
    'neg_gate_mode', 'neg_centering', 'topk', 'hard_filter', 'preset',
    'glm_mode', 'glm_top_pct', 'glm_top_k', 'force_gpu_top', 'device', 'dtype',
    'prefilter_del2_q', 'prefilter_min_del2', 'prefilter_min_total',
    'streaming_agg', 'streaming_chunk_rows', 'validate_smiles',
)


def _build_hit_payload(args, r1_cols, r2_cols, neg_cols, del2_col, run_root_abs, hit_out_abs, hitter_path):
    payload = {
        "version": 2,
        "script": os.path.basename(hitter_path),
        "script_mtime": os.path.getmtime(hitter_path) if os.path.isfile(hitter_path) else None,
        "run_root": run_root_abs,
        "hit_out": hit_out_abs,
        "normalized_columns": {
            "r1": r1_cols,
            "r2": r2_cols,
            "neg": neg_cols,
            "del2": del2_col,
        },
        "args": {k: _jsonable(getattr(args, k, None)) for k in HIT_CACHE_ARG_KEYS},
    }
    return payload


def _hash_payload(payload: dict) -> str:
    s = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(s.encode('utf-8')).hexdigest()


def _load_cached_hit_params(hit_out_abs: str):
    path = os.path.join(hit_out_abs, 'hit_params.json')
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return None


def _write_hit_params(hit_out_abs: str, payload: dict):
    try:
        ensure_dir(hit_out_abs)
        path = os.path.join(hit_out_abs, 'hit_params.json')
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"[WARN] Failed to write hit params: {e}")


def _file_fingerprint(path: str):
    try:
        st = os.stat(path)
        return {"path": os.path.abspath(path), "size": st.st_size, "mtime": st.st_mtime}
    except Exception:
        return None


def _fastq_dir_fingerprint(fastq_dir: str):
    pats = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
    items = []
    for pat in pats:
        for p in glob.glob(os.path.join(fastq_dir, pat)):
            fp = _file_fingerprint(p)
            if fp:
                items.append(fp)
    items = sorted(items, key=lambda x: x["path"])
    return items


def _merged_dir_fingerprint(merged_dir: str):
    items = []
    for name in _list_merged_fastqs(merged_dir):
        fp = _file_fingerprint(os.path.join(merged_dir, name))
        if fp:
            items.append(fp)
    return sorted(items, key=lambda x: x["path"])


def _build_preproc_payload(args, preproc_path: str, fastq_dir: str, bbinfo: str, merged_dir: str):
    skip_fastp = bool(args.skip_fastp)
    return {
        "version": 1,
        "script": os.path.basename(preproc_path),
        "script_mtime": os.path.getmtime(preproc_path) if os.path.isfile(preproc_path) else None,
        "fastq_dir": os.path.abspath(fastq_dir) if fastq_dir else "",
        # With --skip-fastp the raw FASTQs are never read (01 skips run_fastp); the merged FASTQs
        # already present in 01_fastp_out are the real inputs, so fingerprint those instead.
        "fastq_fingerprint": _fastq_dir_fingerprint(fastq_dir) if (fastq_dir and not skip_fastp) else [],
        "merged_fingerprint": _merged_dir_fingerprint(merged_dir) if skip_fastp else [],
        "bbinfo": _file_fingerprint(bbinfo),
        "mismatch": args.mismatch,
        "skip_fastp": skip_fastp,
        # NOTE: fastp thread count does not change outputs; deliberately not part of the hash.
    }


def _build_decode_payload(decode_path: str, merged_dir: str, fixed_bb: str, mismatch: str):
    # v2: the input dependency is the fingerprint of the merged FASTQs the decoder reads (plus the
    # fixed BB table). v1 chained the preprocess hash, which depended on raw-FASTQ mtimes, so
    # deleting the raw FASTQ directory after a successful run invalidated the decode cache.
    return {
        "version": 2,
        "script": os.path.basename(decode_path),
        "script_mtime": os.path.getmtime(decode_path) if os.path.isfile(decode_path) else None,
        "merged_dir": os.path.abspath(merged_dir),
        "merged_fingerprint": _merged_dir_fingerprint(merged_dir),
        "fixed_bb": _file_fingerprint(fixed_bb),
        "mismatch": mismatch,
    }


def _load_json(path: str):
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return None


def _write_json(path: str, payload: dict):
    try:
        ensure_dir(os.path.dirname(path))
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print(f"[WARN] Failed to write {path}: {e}")


def build_parser():
    p = argparse.ArgumentParser(
        description="DELeGANce — End‑to‑End pipeline runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Core inputs (REQUIRED by default)
    p.add_argument('--fastq-dir', required=False,
                   help='Directory containing paired FASTQs (required if running preprocess without --skip-fastp)')
    p.add_argument('--bbinfo', required=False,
                   help='BB information file (required if running preprocess)')
    p.add_argument('--output-dir', required=True,
                   help='Output run root directory (absolute or relative)')
    p.add_argument('--threads', type=int, default=_env_int('FASTP_THREADS', 4),
                   help='Threads for fastp')
    p.add_argument('--mismatch', choices=['none', 'hp_op_cp'], default='hp_op_cp',
                   help='Preprocess mismatch mode for HP/OP/CP')
    p.add_argument('--skip-fastp', action='store_true', help='Skip fastp stage (use existing merged fastqs if any)')
    p.add_argument('--force-preprocess', action='store_true', help='Force rerun preprocess even if outputs/cache exist')

    # All-in-one output location
    p.add_argument('--hit-out', default=None,
                   help='All-in-one output dir (default: RUN_ROOT/03_normalized/<preset>_<glm>_<dev>_<prefilter> tokens). '
                        'Do not point this at RUN_ROOT/03_normalized itself; per-run outputs live in subdirectories.')

    # Phase toggles
    p.add_argument('--only', choices=['all', 'preprocess', 'decode', 'hit'], default='all',
                   help='Run only a specific phase (or all)')
    p.add_argument('--dry-run', action='store_true', help='Print commands without executing')
    # Toggle-able stop-on-error
    p.add_argument('--stop-on-error', dest='stop_on_error', action='store_true', default=True,
                   help='Stop pipeline on first error (default true)')
    p.add_argument('--no-stop-on-error', dest='stop_on_error', action='store_false',
                   help='Do not stop on the first error')

    # Matrix column names for hit-caller
    # Accept comma- or space-separated lists for R1/R2; NEG supports 1 or 2 names (R1, optional R2)
    p.add_argument('--r1', nargs='+', required=False, help='R1 columns (one or more)')
    p.add_argument('--r2', nargs='+', required=False, help='R2 columns (optional; omit for R1-only)')
    p.add_argument('--neg', nargs='+', required=False, help='NEG columns: R1 [R2 optional] (at most two)')
    p.add_argument('--del2', required=False, help='DEL2 column name')

    # Advanced: pass-through knobs for all-in-one (optional; else defaults apply)
    p.add_argument('--neg-gate-mode', choices=['none', 'soft', 'hard'], default=None)
    p.add_argument('--neg-centering', type=int, choices=[0, 1], default=None)
    p.add_argument('--topk', type=int, default=None)
    p.add_argument('--hard-filter', action='store_true', help='Apply aggressive NEG-heavy preset in hit-caller')
    p.add_argument('--preset', choices=['auto','balanced','medium','neg_heavy','very_hard','lenient'], default=None,
                   help='Hit-caller preset (overrides individual knobs)')
    p.add_argument('--auto-opt', type=int, choices=[0,1], default=1,
                   help='Auto-detect CUDA and data size to choose device/dtype and GLM mode')

    # GLM performance knobs (pass-through)
    p.add_argument('--glm-mode', choices=['full','top','skip'], default=None,
                   help='GLM strategy: full=all IDs, top=fit top by baseline, skip=baseline only')
    p.add_argument('--glm-top-pct', type=float, default=None,
                   help='When --glm-mode=top, percentage of IDs to GLM (e.g., 1=top 1%%)')
    p.add_argument('--glm-top-k', type=int, default=None,
                   help='When --glm-mode=top, cap absolute number of IDs to GLM')
    p.add_argument('--force-gpu-top', type=int, choices=[0,1], default=None,
                   help='Force GPU for glm_mode=top via DELEGANCE_FORCE_GPU_TOP (0/1)')

    # Device/dtype pass-through for GLM acceleration
    p.add_argument('--device', default=None, help='Device for GLM: cpu, cuda, cuda:0 ...')
    p.add_argument('--dtype', choices=['auto','float32','float64'], default=None,
                   help='GLM dtype: auto=float32 on CUDA, float64 on CPU')

    # Prefilter pass-through (03 honors these)
    p.add_argument('--prefilter-del2-q', type=float, default=None,
                   help='Quantile gate for DEL2 in prefilter (0..1). Higher = fewer IDs.')
    p.add_argument('--prefilter-min-del2', type=int, default=None,
                   help='Absolute minimum DEL2 counts to keep (after quantile).')
    p.add_argument('--prefilter-min-total', type=int, default=None,
                   help='Absolute minimum total counts across all samples to keep.')

    # Force rerun of HitCaller even if outputs exist
    p.add_argument('--force-hit', action='store_true', help='Force rerun hit-caller even if outputs exist')
    p.add_argument('--force-decode', action='store_true', help='Force rerun decode even if outputs/cache exist')
    # Streaming + SMILES validation passthrough to hit-caller
    p.add_argument('--streaming-agg', type=int, choices=[0,1], default=None,
                   help='Pass to 03_call_hits.py: enable streaming aggregation to reduce memory')
    p.add_argument('--streaming-chunk-rows', type=int, default=None,
                   help='Pass to 03_call_hits.py: chunk size for streaming aggregation')
    p.add_argument('--validate-smiles', type=int, choices=[0,1], default=None,
                   help='Pass to 03_call_hits.py: validate SMILES of consensus hits (requires RDKit)')

    return p


def main():
    ap = build_parser(); args = ap.parse_args()

    # Resolve paths (relative output paths are resolved from the current working directory)
    base_dir = os.getcwd()
    run_root_abs = os.path.abspath(os.path.join(base_dir, args.output_dir)) if not os.path.isabs(args.output_dir) else os.path.abspath(args.output_dir)
    run_root_rel = os.path.relpath(run_root_abs, base_dir)
    # Default normalized root; the final per-run output directory will be set after auto-opt
    norm_root_abs = os.path.join(run_root_abs, '03_normalized')
    if not args.dry_run:
        ensure_dir(run_root_abs)
        ensure_dir(norm_root_abs)
    hit_out = None  # decide after auto-optimization so that GLM mode suffix is reflected
    hit_out_abs = None

    # Script paths (support legacy filenames; prefer the new simplified names)
    # New names appear first so modern installs pick them up, while older repos still work.
    def _pick_script(candidates):
        for name in candidates:
            p = os.path.join(THIS_DIR, name)
            if os.path.isfile(p):
                return p
        return os.path.join(THIS_DIR, candidates[0])

    preproc = _pick_script(['01_preprocess_reads.pl', '01_preprocess_fastq.pl', '01_preprocess_fastq_250816p2.pl'])
    decode  = _pick_script(['02_decode_reads.pl', '02_decode_fastq.pl', '02_decode_fastq_250816p2.pl'])
    hitter  = _pick_script(['03_call_hits.py', '03_all_in_one_hitcaller.py', '03_all_in_one_hitcaller_250828.py'])
    for f in (preproc, decode, hitter):
        if not os.path.isfile(f):
            print(f"[ERROR] Script not found: {f}")
            return 2

    # Dependency checks (fastp only if we will run preprocess)
    if not which('perl'):
        print('[ERROR] perl not found in PATH'); return 2

    # Logs
    master_log = os.path.join(run_root_abs, '00_pipeline.log')
    print(f"[INFO] Pipeline log → {master_log}")

    # Normalize R1/R2/NEG lists (split on commas and flatten)
    def _split_list(lst):
        out = []
        for item in lst or []:
            out.extend([t for t in str(item).split(',') if t != ''])
        return [x.strip() for x in out if x.strip()]
    r1_cols = _split_list(args.r1)
    r2_cols = _split_list(args.r2)
    neg_cols = _split_list(args.neg)
    del2_col = (args.del2 or "").strip()

    # Only require hit-caller columns if we will run hit stage
    if args.only in ('all', 'hit'):
        if len(r1_cols) == 0:
            print('[ERROR] --r1 must specify at least one column (R1)')
            return 2
        if len(neg_cols) == 0:
            print('[ERROR] --neg must specify at least one column (NEG for R1; optional second for R2)')
            return 2
        if len(neg_cols) > 2:
            print(f"[ERROR] --neg accepts at most two columns (NEG for R1, optional NEG for R2); got {len(neg_cols)}: {' '.join(neg_cols)}")
            return 2
        if del2_col == '':
            print('[ERROR] --del2 must not be empty')
            return 2

    # Preprocess cache/hash setup (computed BEFORE input validation so that an existing,
    # cache-matching preprocess output can satisfy a run that omits --fastq-dir/--bbinfo).
    preproc_out_dir = os.path.join(run_root_abs, '01_fastp_out')
    preproc_params_path = os.path.join(preproc_out_dir, 'preprocess_params.json')
    preproc_marker = os.path.join(run_root_abs, PREPROCESS_DONE_MARKER)
    merged_dir = preproc_out_dir
    cached_preproc = _load_json(preproc_params_path)
    # Effective inputs: explicit CLI values win; otherwise fall back to what the cached run used.
    fastq_dir_eff = args.fastq_dir or ""
    bbinfo_eff = args.bbinfo or ""
    if isinstance(cached_preproc, dict):
        if not fastq_dir_eff:
            fastq_dir_eff = str(cached_preproc.get("fastq_dir") or "")
        if not bbinfo_eff:
            bb_fp = cached_preproc.get("bbinfo")
            if isinstance(bb_fp, dict):
                bbinfo_eff = str(bb_fp.get("path") or "")
    preprocess_done = _has_fastp_outputs(run_root_abs)
    preproc_payload = _build_preproc_payload(args, preproc, fastq_dir_eff, bbinfo_eff, merged_dir)
    preproc_hash = _hash_payload(preproc_payload)
    preproc_payload["hash"] = preproc_hash
    preproc_cache_ok = preprocess_done and isinstance(cached_preproc, dict) and str(cached_preproc.get("hash", "")) == preproc_hash
    will_run_preprocess = args.only in ('all', 'preprocess') and (args.force_preprocess or not preproc_cache_ok)
    preproc_skip_note = f"outputs + cache match (hash={preproc_hash[:12]})"
    if preprocess_done and not preproc_cache_ok and args.only in ('all', 'preprocess'):
        cached_hash = cached_preproc.get("hash") if isinstance(cached_preproc, dict) else None
        if cached_hash:
            print(f"[INFO] Preprocess outputs found but parameter hash differs (cached={cached_hash[:12]}, current={preproc_hash[:12]}). Will rerun.")
        else:
            print("[INFO] Preprocess outputs found but cache missing/invalid. Will rerun.")
    if will_run_preprocess and (not args.skip_fastp) and not which('fastp'):
        if not _has_merged_fastqs(run_root_abs):
            print('[ERROR] fastp not found in PATH and no merged FASTQs exist under RUN_ROOT/01_fastp_out.')
            print('        Install fastp, or provide pre-merged FASTQs in 01_fastp_out and pass --skip-fastp.')
            return 2
        print('[WARN] fastp not found; switching to --skip-fastp for preprocess')
        args.skip_fastp = True
        preproc_payload["skip_fastp_effective"] = True  # informational only (hash already computed)

    # Validate required inputs only if we will actually run preprocess. With --skip-fastp the
    # preprocess stage never reads the raw FASTQ directory (01 skips run_fastp), so it is optional.
    if will_run_preprocess:
        if not args.skip_fastp and (not fastq_dir_eff or not os.path.isdir(fastq_dir_eff)):
            if preprocess_done and not args.force_preprocess:
                # Outputs + completion marker exist and only the raw FASTQs are gone (e.g. deleted to
                # save space): keep the existing preprocess outputs instead of failing the whole run.
                print('[WARN] Raw FASTQ directory is missing but complete preprocess outputs exist under RUN_ROOT; '
                      'keeping them (preprocess cache not refreshed). Use --force-preprocess with --fastq-dir to redo.')
                will_run_preprocess = False
                preproc_skip_note = "existing outputs kept (raw FASTQ dir missing; cache not refreshed)"
            else:
                print('[ERROR] --fastq-dir is required for preprocess and must exist (or provide existing, cache-matching outputs under RUN_ROOT/01_fastp_out)')
                return 2
        if will_run_preprocess and (not bbinfo_eff or not os.path.isfile(bbinfo_eff)):
            print('[ERROR] --bbinfo is required for preprocess and must point to a file (or provide existing, cache-matching outputs under RUN_ROOT)')
            return 2

    # Decode paths (hash is computed AFTER the preprocess stage, see below)
    decode_out_dir = os.path.join(run_root_abs, '02_decoded')
    decode_params_path = os.path.join(decode_out_dir, 'decode_params.json')
    decode_marker = os.path.join(run_root_abs, DECODE_DONE_MARKER)
    fixed_bb   = os.path.join(run_root_abs, 'BB_information_fixed.tsv')

    # If we are skipping fastp, merged FASTQs must already exist for preprocess/decode (the hit
    # stage only needs the decoded matrix). BB_information_fixed.tsv is produced by the preprocess
    # stage itself (also with --skip-fastp), so it is only required when preprocess will NOT run.
    if args.skip_fastp and args.only in ('all', 'decode'):
        if not _has_merged_fastqs(run_root_abs):
            print('[ERROR] --skip-fastp requested but no merged FASTQs found under RUN_ROOT/01_fastp_out.')
            print('        Either run fastp (omit --skip-fastp) or provide merged FASTQs in 01_fastp_out first.')
            return 2
        if not will_run_preprocess and not _file_or_gz_exists(fixed_bb):
            print('[ERROR] --skip-fastp requested but BB_information_fixed.tsv is missing under RUN_ROOT.')
            print('        Provide RUN_ROOT/BB_information_fixed.tsv or rerun without --skip-fastp.')
            return 2

    # Hit-only requires decoded outputs to already exist.
    if args.only == 'hit' and not _has_decoded_outputs(run_root_abs):
        print('[ERROR] Hit-only requested but decoded outputs are missing (RUN_ROOT/02_decoded/raw_counts_matrix.tsv).')
        print('        Run --only decode (or full pipeline) first, or point --output-dir to an existing run.')
        return 2

    # 1) Preprocess
    rc = 0
    if args.only in ('all', 'preprocess'):
        if not will_run_preprocess:
            print(f"[INFO] Preprocess {preproc_skip_note}; skipping preprocess.")
        else:
            cmd = ['perl', preproc]
            cmd += ['-b', bbinfo_eff]
            # pass absolute run_root to preprocess; it detects absolute and uses it directly
            if fastq_dir_eff:
                cmd += ['-f', fastq_dir_eff]
            cmd += ['-o', run_root_abs, '-t', str(int(args.threads)), '--mismatch', args.mismatch]
            if args.skip_fastp:
                cmd += ['--skip-fastp']
            print(f"[INFO] Preprocess → RUN_ROOT={run_root_abs}")
            if args.dry_run:
                print('DRY-RUN:', ' '.join(cmd))
            else:
                # A failed rerun must not leave the previous run's marker/cache behind.
                _clear_marker(preproc_marker)
                _clear_marker(preproc_params_path)
                rc = run_cmd(cmd, master_log)
                if rc != 0 and args.stop_on_error:
                    print('[ERROR] Preprocess failed; aborting.')
                    return rc
                if rc == 0:
                    _write_marker(preproc_marker)
                    _write_json(preproc_params_path, preproc_payload)
                    print(f"[OK] Stored preprocess parameter hash: {preproc_hash[:12]}")

    # Decode cache/hash setup — computed here, after preprocess has (re)written
    # BB_information_fixed.tsv, so the stored fingerprint matches the file the decode actually used.
    decode_done = _has_decoded_outputs(run_root_abs)
    decode_payload = _build_decode_payload(decode, merged_dir, fixed_bb, args.mismatch)
    decode_hash = _hash_payload(decode_payload)
    decode_payload["hash"] = decode_hash
    cached_decode = _load_json(decode_params_path)
    decode_cache_ok = decode_done and isinstance(cached_decode, dict) and str(cached_decode.get("hash", "")) == decode_hash
    if decode_done and not decode_cache_ok and args.only in ('all', 'decode'):
        cached_hash = cached_decode.get("hash") if isinstance(cached_decode, dict) else None
        if cached_hash:
            print(f"[INFO] Decode outputs found but parameter hash differs (cached={cached_hash[:12]}, current={decode_hash[:12]}). Will rerun.")
        else:
            print("[INFO] Decode outputs found but cache missing/invalid. Will rerun.")
    if not decode_cache_ok:
        decode_done = False

    # 2) Decode
    if rc == 0 and args.only in ('all', 'decode'):
        if decode_done and not args.force_decode:
            print(f"[INFO] Decode outputs + cache match (hash={decode_hash[:12]}); skipping decode.")
        else:
            out_dir    = decode_out_dir
            cmd = ['perl', decode, '--merged-dir', merged_dir, '--fixed-bb-file', fixed_bb, '--out-dir', out_dir, '--mismatch', args.mismatch]
            print(f"[INFO] Decode   → MERGED={merged_dir}  FIXED_BB={fixed_bb}  OUT={out_dir}")
            if args.dry_run:
                print('DRY-RUN:', ' '.join(cmd))
            else:
                _clear_marker(decode_marker)
                _clear_marker(decode_params_path)
                rc = run_cmd(cmd, master_log)
                if rc != 0 and args.stop_on_error:
                    print('[ERROR] Decode failed; aborting.')
                    return rc
                if rc == 0:
                    _write_marker(decode_marker)
                    _write_json(decode_params_path, decode_payload)
                    print(f"[OK] Stored decode parameter hash: {decode_hash[:12]}")

    # Auto-optimization: choose device/dtype/GLM mode heuristically if not specified.
    # Only when the hit stage will run: estimating the matrix size scans the (possibly multi-GB) file.
    if int(getattr(args, 'auto_opt', 1)) == 1 and args.only in ('all', 'hit'):
        # Device/dtype
        # Honour DELEGANCE_DISABLE_TORCH (03_call_hits.py does): otherwise the output dir would be
        # labelled dev_cuda_fp32 while the hit-caller actually runs its CPU/float64 fallback.
        torch_disabled = os.environ.get("DELEGANCE_DISABLE_TORCH", "").strip().lower() in ("1", "true", "yes")
        if (args.device is None) and not torch_disabled:
            try:
                import torch  # type: ignore
                if hasattr(torch, 'cuda') and torch.cuda.is_available():
                    args.device = 'cuda'
                    print('[AUTO] CUDA detected → using --device cuda')
            except Exception:
                pass
        if (args.dtype is None) and (args.device is not None and str(args.device).startswith('cuda')):
            args.dtype = 'auto'  # float32 on CUDA
        # GLM mode by rough data size
        if args.glm_mode is None:
            mat = _pick_decoded_matrix(run_root_abs)
            nrows = _estimate_rows(mat)
            if nrows >= 5_000_000:
                args.glm_mode = 'top'
                if args.glm_top_pct is None: args.glm_top_pct = 0.5
                if args.glm_top_k is None: args.glm_top_k = 500_000
                print(f"[AUTO] Large matrix (~{nrows:,} rows) → GLM top mode: pct={args.glm_top_pct}%, cap={args.glm_top_k}")
            elif nrows >= 1_000_000:
                args.glm_mode = 'top'
                if args.glm_top_pct is None: args.glm_top_pct = 2.0
                if args.glm_top_k is None: args.glm_top_k = 200_000
                print(f"[AUTO] Medium matrix (~{nrows:,} rows) → GLM top mode: pct={args.glm_top_pct}%, cap={args.glm_top_k}")
            elif nrows > 0:
                print(f"[AUTO] Matrix size ~{nrows:,} rows → GLM full mode")
        # Prefilter heuristics (only if not explicitly provided)
        if args.prefilter_del2_q is None or args.prefilter_min_total is None or args.prefilter_min_del2 is None:
            mat = _pick_decoded_matrix(run_root_abs)
            nrows = _estimate_rows(mat)
            q = args.prefilter_del2_q
            min_tot = args.prefilter_min_total
            min_del2 = args.prefilter_min_del2
            if nrows >= 50_000_000:
                q = q if q is not None else 0.90
                min_tot = min_tot if min_tot is not None else 5
                min_del2 = min_del2 if min_del2 is not None else 2
            elif nrows >= 20_000_000:
                q = q if q is not None else 0.85
                min_tot = min_tot if min_tot is not None else 4
                min_del2 = min_del2 if min_del2 is not None else 2
            elif nrows >= 5_000_000:
                q = q if q is not None else 0.75
                min_tot = min_tot if min_tot is not None else 3
                min_del2 = min_del2 if min_del2 is not None else 1
            elif nrows >= 1_000_000:
                q = q if q is not None else 0.50
                min_tot = min_tot if min_tot is not None else 2
                min_del2 = min_del2 if min_del2 is not None else 1
            else:
                q = q if q is not None else 0.0
                min_tot = min_tot if min_tot is not None else 0
                min_del2 = min_del2 if min_del2 is not None else 0
            args.prefilter_del2_q = q
            args.prefilter_min_total = min_tot
            args.prefilter_min_del2 = min_del2
            print(f"[AUTO] Prefilter tuned: del2_q={q}, min_total={min_tot}, min_del2={min_del2}")
        # 03_call_hits.py forces glm_mode=top back to CPU unless the GPU is explicitly forced
        # (--force_gpu_top 1 or DELEGANCE_FORCE_GPU_TOP). Mirror that here so the output directory
        # token (dev_*) and hit_params.json describe the device actually used.
        env_force_gpu = os.environ.get("DELEGANCE_FORCE_GPU_TOP", "").strip().lower() in ("1", "true", "yes", "y")
        if (args.glm_mode == 'top' and str(args.device or '').startswith('cuda')
                and int(args.force_gpu_top or 0) != 1 and not env_force_gpu):
            print('[AUTO] glm_mode=top on CUDA is forced to CPU by 03_call_hits.py; using --device cpu '
                  '(pass --force-gpu-top 1 to keep CUDA)')
            args.device = 'cpu'
            args.dtype = None

    # Decide final output directory name (after auto-opt so we know glm_mode)
    if args.hit_out:
        hit_out = args.hit_out
    else:
        preset_name = args.preset if args.preset else None
        glm_mode = (args.glm_mode or 'full')
        # 1) GLM mode token with topK if applicable
        mode_part = f"glm_{glm_mode}"
        if glm_mode == 'top':
            k_label = None
            if args.glm_top_k is not None:
                try:
                    k_label = f"topK{int(args.glm_top_k)}"
                except Exception:
                    k_label = None
            if k_label is None:
                # Estimate K from pct and matrix size if possible
                try:
                    mat = _pick_decoded_matrix(run_root_abs)
                    nrows = _estimate_rows(mat)
                    pct = float(args.glm_top_pct) if args.glm_top_pct is not None else 0.0
                    if nrows > 0 and pct > 0.0:
                        k_est = max(1, int(round(nrows * (pct / 100.0))))
                        k_label = f"topK{k_est}"
                    else:
                        k_label = f"topPct{pct}"
                except Exception:
                    pct = float(args.glm_top_pct) if args.glm_top_pct is not None else 0.0
                    k_label = f"topPct{pct}"
            mode_part = f"glm_top_{k_label}"

        # 2) R1-only token
        r1_only = (len(r2_cols) == 0)
        rtoken = 'r1only' if r1_only else None

        # 3) Device/dtype token (effective dtype: auto-> float32 on CUDA, else float64)
        dev = (args.device or 'cpu')
        dev_token = dev.replace(':', '')
        eff_dtype = (args.dtype or 'auto')
        if eff_dtype == 'auto':
            eff_dtype = 'float32' if (dev and str(dev).startswith('cuda')) else 'float64'
        dtype_token = 'fp32' if eff_dtype == 'float32' else 'fp64'
        dd_part = f"dev_{dev_token}_{dtype_token}"

        # 4) Prefilter token
        pf_parts = []
        if args.prefilter_del2_q is not None:
            try:
                q = float(args.prefilter_del2_q)
                if 0.0 < q < 1.0:
                    pf_parts.append(f"q{int(round(q*100))}")
                elif q > 0:
                    pf_parts.append(f"q{q}")
            except Exception:
                pass
        if args.prefilter_min_total is not None and int(args.prefilter_min_total) > 0:
            pf_parts.append(f"mt{int(args.prefilter_min_total)}")
        if args.prefilter_min_del2 is not None and int(args.prefilter_min_del2) > 0:
            pf_parts.append(f"md{int(args.prefilter_min_del2)}")
        pf_part = ("pf_" + "_".join(pf_parts)) if pf_parts else None

        tokens = []
        if preset_name:
            tokens.append(preset_name)
        tokens.append(mode_part)
        if rtoken:
            tokens.append(rtoken)
        tokens.append(dd_part)
        if pf_part:
            tokens.append(pf_part)
        dir_name = "_".join(tokens)
        hit_out = os.path.join(norm_root_abs, dir_name)
    hit_out_abs = hit_out if os.path.isabs(hit_out) else os.path.abspath(os.path.join(base_dir, hit_out))
    # NOTE: hit_out_abs is created only when the hit stage actually runs (see below), so
    # --only preprocess/decode and --dry-run do not leave empty preset directories behind.
    if os.path.realpath(hit_out_abs) == os.path.realpath(norm_root_abs):
        print("[WARN] --hit-out points at RUN_ROOT/03_normalized itself; per-run outputs should live in a subdirectory. Legacy migration is disabled for this run.")

    # Build current hit parameters payload + hash (used for cache validation).
    # Input dependency = fingerprint of the decoded matrix actually consumed by the hit-caller
    # (size/mtime), not the decode parameter hash (which in turn depended on fastq mtimes/threads).
    hit_payload = _build_hit_payload(args, r1_cols, r2_cols, neg_cols, del2_col, run_root_abs, hit_out_abs, hitter)
    hit_payload["decoded_matrix"] = _file_fingerprint(_pick_decoded_matrix(run_root_abs))
    current_hash = _hash_payload(hit_payload)
    hit_payload["hash"] = current_hash

    cached_payload = _load_cached_hit_params(hit_out_abs)
    has_hybrid = os.path.isfile(os.path.join(hit_out_abs, '05_hybrid_annot.tsv'))
    cached_hit_ok = has_hybrid and isinstance(cached_payload, dict) and str(cached_payload.get("hash", "")) == current_hash
    if has_hybrid and not cached_hit_ok:
        cached_hash = (cached_payload.get("hash") if isinstance(cached_payload, dict) else None)
        if cached_hash:
            print(f"[INFO] Cached HitCaller outputs found but parameter hash differs (cached={cached_hash[:12]}, current={current_hash[:12]}). Recomputing.")
        else:
            print("[INFO] Cached HitCaller outputs found but hit_params.json missing or invalid. Recomputing.")
    used_hit_dir = None

    # 3) All‑in‑one hit caller
    if rc == 0 and args.only in ('all', 'hit'):
        if cached_hit_ok and not getattr(args, 'force_hit', False):
            used_hit_dir = hit_out_abs
            print(f"[INFO] Cached HitCaller outputs match parameters (hash={current_hash[:12]}). Skipping recompute.")
        else:
            cmd = [sys.executable, hitter, '--run_root', run_root_abs, '--outdir', hit_out]
            # Column names
            cmd += ['--del2_col', del2_col]
            cmd += ['--r1_cols', *r1_cols]
            if len(r2_cols) > 0:
                cmd += ['--r2_cols', *r2_cols]
            if len(neg_cols) >= 1:
                cmd += ['--neg_r1_col', neg_cols[0]]
            if len(neg_cols) >= 2:
                cmd += ['--neg_r2_col', neg_cols[1]]
            # Optional passthroughs
            if args.neg_gate_mode: cmd += ['--neg_gate_mode', args.neg_gate_mode]
            if args.neg_centering is not None: cmd += ['--neg_centering', str(int(args.neg_centering))]
            if args.topk is not None: cmd += ['--topk', str(int(args.topk))]
            if args.hard_filter: cmd += ['--hard_filter', '1']
            if args.preset: cmd += ['--preset', args.preset]
            # GLM performance pass-throughs
            if args.glm_mode: cmd += ['--glm_mode', args.glm_mode]
            if args.glm_top_pct is not None: cmd += ['--glm_top_pct', str(float(args.glm_top_pct))]
            if args.glm_top_k is not None: cmd += ['--glm_top_k', str(int(args.glm_top_k))]
            if args.force_gpu_top is not None: cmd += ['--force_gpu_top', str(int(args.force_gpu_top))]
            if args.device: cmd += ['--device', args.device]
            if args.dtype: cmd += ['--dtype', args.dtype]
            # Prefilter pass-throughs
            if args.prefilter_del2_q is not None: cmd += ['--prefilter_del2_q', str(float(args.prefilter_del2_q))]
            if args.prefilter_min_total is not None: cmd += ['--prefilter_min_total', str(int(args.prefilter_min_total))]
            if args.prefilter_min_del2 is not None: cmd += ['--prefilter_min_del2', str(int(args.prefilter_min_del2))]
            # Streaming/SMILES validation pass-throughs
            if args.streaming_agg is not None: cmd += ['--streaming_agg', str(int(args.streaming_agg))]
            if args.streaming_chunk_rows is not None: cmd += ['--streaming_chunk_rows', str(int(args.streaming_chunk_rows))]
            if args.validate_smiles is not None: cmd += ['--validate_smiles', str(int(args.validate_smiles))]
            print(f"[INFO] HitCaller → RUN_ROOT={run_root_abs}, OUT={hit_out}")
            if args.dry_run:
                print('DRY-RUN:', ' '.join(cmd))
            else:
                ensure_dir(hit_out_abs)
                rc = run_cmd(cmd, master_log)
                if rc != 0 and args.stop_on_error:
                    print('[ERROR] HitCaller failed; aborting.')
                    return rc
                if rc == 0:
                    used_hit_dir = hit_out_abs
                    _write_hit_params(hit_out_abs, hit_payload)
                    print(f"[OK] Stored hit parameter hash: {current_hash[:12]}")

    # 4) Interactive hits (Bokeh) + index.html aggregator
    if rc == 0 and args.only in ('all', 'hit'):
        try:
            # Migrate any legacy outputs directly under 03_normalized into a tagged folder once
            try:
                legacy_markers = [
                    os.path.join(norm_root_abs, '05_hybrid_annot.tsv'),
                    os.path.join(norm_root_abs, 'report.html'),
                ]
                # Never move files in --dry-run, and never when --hit-out IS 03_normalized
                # (the "legacy" files would be the outputs we just produced).
                hit_out_is_norm_root = (os.path.realpath(hit_out_abs) == os.path.realpath(norm_root_abs))
                if any(os.path.isfile(p) for p in legacy_markers) and not args.dry_run and not hit_out_is_norm_root:
                    import time, shutil as _sh
                    stamp = time.strftime('%Y%m%d_%H%M%S')
                    legacy_dir = os.path.join(norm_root_abs, f'legacy_{stamp}')
                    ensure_dir(legacy_dir)
                    # Move known result files and plots
                    for fname in (
                        '01_raw_by_id.tsv','02_cpm_by_id.tsv','03_glm_results.csv','04_read_scaler.tsv',
                        '05_hybrid_annot.tsv','06_topk_glm.tsv','07_topk_rs.tsv','08_topk_consensus.tsv','report.html'):
                        src = os.path.join(norm_root_abs, fname)
                        if os.path.isfile(src):
                            _sh.move(src, os.path.join(legacy_dir, fname))
                    plots_src = os.path.join(norm_root_abs, 'plots')
                    plots_dst = os.path.join(legacy_dir, 'plots')
                    if os.path.isdir(plots_src):
                        if os.path.isdir(plots_dst):
                            _sh.rmtree(plots_dst)
                        _sh.move(plots_src, plots_dst)
                    print(f"[INFO] Migrated legacy outputs to {legacy_dir}")
            except Exception as _e:
                print(f"[WARN] Legacy migration skipped: {_e}")

            # Choose source dir: prefer current run output; fallback to latest existing under norm_root
            cand1 = hit_out_abs
            cand2 = _find_existing_hybrid_dir(norm_root_abs)
            def _has_hybrid(d):
                try:
                    return os.path.isfile(os.path.join(d, '05_hybrid_annot.tsv'))
                except Exception:
                    return False
            if used_hit_dir and _has_hybrid(used_hit_dir):
                src_dir = used_hit_dir
            else:
                src_dir = cand1 if _has_hybrid(cand1) else (cand2 if _has_hybrid(cand2) else cand1)
            interactive_html = os.path.join(src_dir, 'interactive_hits.html')
            hybrid_tsv = os.path.join(src_dir, '05_hybrid_annot.tsv')
            bbinfo_fixed = os.path.join(run_root_abs, 'BB_information_fixed.tsv')
            if os.path.isfile(hybrid_tsv):
                if not _has_interactive_deps():
                    print("[WARN] Interactive HTML skipped: missing bokeh (add narwhals if bokeh requests it).")
                else:
                    cmd = [
                        sys.executable, _pick_script(['04_build_interactive_report.py','04_hit_finder.py','04_hit_finder_250816p2.py']),
                        '--master_tsv', hybrid_tsv,
                        '--bbinfo', bbinfo_fixed,
                        '--out', interactive_html,
                        '--top_hitscore', '10000',
                        '--plot_height', '260',
                    ]
                    print(f"[INFO] Interactive (Bokeh) → {interactive_html}")
                    if args.dry_run:
                        print('DRY-RUN:', ' '.join(cmd))
                    else:
                        rc2 = run_cmd(cmd, master_log)
                        if rc2 != 0:
                            print('[WARN] Interactive HTML generation failed (continuing).')
            else:
                print(f"[WARN] Skipping interactive HTML; not found: {hybrid_tsv}")
        except Exception as e:
            print(f"[WARN] Interactive generation error: {e}")

        # Write aggregator index.html at run_root
        try:
            index_html = os.path.join(run_root_abs, 'index.html')
            rel_norm_root = os.path.relpath(norm_root_abs, run_root_abs)
            # Collect preset subdirs (and include root if it has results)
            # include base 03_normalized if it has outputs
            def _has_outputs(d):
                return os.path.isfile(os.path.join(d, '05_hybrid_annot.tsv')) or \
                       os.path.isfile(os.path.join(d, 'report.html'))
            presets = []
            try:
                for name in sorted(os.listdir(norm_root_abs)):
                    p = os.path.join(norm_root_abs, name)
                    # Only directories that actually hold results (skips empty/aborted preset dirs)
                    if os.path.isdir(p) and _has_outputs(p):
                        presets.append((name, p))
            except Exception:
                pass
            include_base = _has_outputs(norm_root_abs)
            # Build HTML
            def _h(s: str) -> str:
                try:
                    return _html_escape(str(s), quote=True)
                except Exception:
                    return str(s)
            # Known per-run result files: (filename, label). Links are emitted only for files that exist.
            RESULT_LINKS = (
                ('report.html', 'All-in-one report'),
                ('beginner_qc_report.html', 'Beginner QC report'),
                ('05_hybrid_annot.tsv', '05_hybrid_annot.tsv'),
                ('08_topk_consensus.tsv', '08_topk_consensus.tsv'),
                ('06_topk_glm.tsv', '06_topk_glm.tsv'),
                ('07_topk_rs.tsv', '07_topk_rs.tsv'),
            )
            def _result_items(rel_path: str, abs_dir: str):
                items = []
                for fname, label in RESULT_LINKS:
                    if os.path.isfile(os.path.join(abs_dir, fname)):
                        items.append(f"<li><a href='{_h(rel_path)}/{fname}'>{_h(label)}</a></li>")
                return items
            def _interactive_items(rel_path: str, abs_dir: str):
                items = []
                iframe = ""
                ih = os.path.join(abs_dir, 'interactive_hits.html')
                if os.path.isfile(ih):
                    items.append(f"<li><a href='{_h(rel_path)}/interactive_hits.html'>Interactive (embedded below)</a></li>")
                    popouts = []
                    for fname, label in (
                        ('interactive_p1.html', 'BB1×BB2'),
                        ('interactive_p2.html', 'BB1×BB3'),
                        ('interactive_p3.html', 'BB2×BB3'),
                        ('interactive_p4.html', 'BB1×BB4'),
                        ('interactive_table.html', 'Top table'),
                    ):
                        if os.path.isfile(os.path.join(abs_dir, fname)):
                            popouts.append(f"<a href='{_h(rel_path)}/{fname}'>{label}</a>")
                    if popouts:
                        items.append(f"<li>Pop-outs: {' · '.join(popouts)}</li>")
                    iframe = f"<iframe src='{_h(rel_path)}/interactive_hits.html'></iframe>"
                return items, iframe

            html = [
                "<!DOCTYPE html>",
                "<html><head><meta charset='utf-8'><title>DELeGANce Results Index</title>",
                "<style>",
                "body{font-family:system-ui,-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;line-height:1.5;color:#222}",
                "section{margin:16px 0;padding:16px;border:1px solid #e5e7eb;border-radius:8px}",
                "ul{margin:8px 0 0 18px} a{text-decoration:none;color:#1f2937} a:hover{text-decoration:underline}",
                "iframe{width:100%;height:640px;border:1px solid #e5e7eb;border-radius:6px}",
                ".preset{background:#fafafa}",
                ".tags{font-weight:normal;color:#555;font-size:0.9em;margin-left:6px}",
                "</style></head><body>",
                "<h1>DELeGANce — Run Results</h1>",
                f"<p><b>Run root:</b> {_h(run_root_rel)}</p>",
                f"<p>Normalized outputs are under <code>{_h(rel_norm_root)}</code></p>",
            ]
            # Base section (no preset)
            if include_base:
                rel_base = os.path.relpath(norm_root_abs, run_root_abs)
                inter_items, inter_iframe = _interactive_items(rel_base, norm_root_abs)
                html += [
                    "<section class='preset'><h2>Preset: (default)</h2>",
                    "<ul>",
                    *_result_items(rel_base, norm_root_abs),
                    *inter_items,
                    "</ul>",
                    inter_iframe,
                    "</section>",
                ]
            # Each preset subdir
            for name, absdir in presets:
                rel = os.path.relpath(absdir, run_root_abs)
                inter_items, inter_iframe = _interactive_items(rel, absdir)
                # Derive preset and tags from folder name convention. Directories created without
                # --preset start directly with a tag token (glm_/r1only/dev_/pf_/legacy_).
                first_tok = name.split('_', 1)[0]
                if first_tok in ('glm', 'r1only', 'dev', 'pf'):
                    preset_name = '(none)'
                    tag_str = name
                elif '_' in name:
                    preset_name = first_tok
                    tag_str = name[len(preset_name)+1:]
                else:
                    preset_name = name
                    tag_str = ''
                tag_html = f" <span class='tags'>({_h(tag_str)})</span>" if tag_str else ''
                html += [
                    f"<section class='preset'><h2>Preset: {_h(preset_name)}{tag_html}</h2>",
                    "<ul>",
                    *_result_items(rel, absdir),
                    *inter_items,
                    "</ul>",
                    inter_iframe,
                    "</section>",
                ]
            html += ["</body></html>"]
            if args.dry_run:
                print(f"DRY-RUN: would write {index_html}")
            else:
                with open(index_html, 'w', encoding='utf-8') as f:
                    f.write("\n".join(html))
                print(f"[OK] Wrote index.html → {index_html}")
        except Exception as e:
            print(f"[WARN] Failed to write index.html: {e}")

    if rc != 0:
        # Reached only with --no-stop-on-error: a stage failed and later stages were skipped.
        print(f'\n[FAIL] Pipeline finished with errors (last stage rc={rc}); later stages were skipped.')
        print(f" - RUN_ROOT:  {run_root_abs}")
        print(f" - Master log: {master_log}")
        return 1
    print('\n[OK] Pipeline completed')
    print(f" - RUN_ROOT:  {run_root_abs}")
    print(f" - HIT_OUT:   {hit_out_abs}")
    print(f" - Master log: {master_log}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
