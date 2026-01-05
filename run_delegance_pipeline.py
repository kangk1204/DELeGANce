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
      --bbinfo 00_DELeGANce_KRASMAT2A_BB_information.txt \
      --output-dir DELeGANce_out/KRAS_example_run \
      --r1 D_R1C1,D_R1C2,D_R1C3 \
      --r2 D_R2C1,D_R2C2,D_R2C3 \
      --neg D_R1C4 \
      --del2 DEL2 \
      --threads 6 --mismatch hp_op_cp

  # Hit stage only on an existing run (still requires column names)
  python3 run_delegance_pipeline.py \
      --only hit \
      --output-dir DELeGANce_out/KRAS_example_run \
      --r1 D_R1C1,D_R1C2,D_R1C3 --r2 D_R2C1,D_R2C2,D_R2C3 --neg D_R1C4 --del2 DEL2
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
from html import escape as _html_escape

THIS_DIR = os.path.abspath(os.path.dirname(__file__))


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
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
        for line in proc.stdout:
            sys.stdout.write(line)
            lf.write(line)
        proc.wait()
        rc = proc.returncode
        lf.write(f"[{_timestamp()}] EXIT CODE: {rc}\n")
        lf.flush()
    return rc


def _file_or_gz_exists(path: str) -> bool:
    return os.path.isfile(path) or os.path.isfile(path + ".gz")


def _has_fastp_outputs(run_root: str) -> bool:
    """Return True when merged FASTQs and fixed BB info already exist."""
    merged_dir = os.path.join(run_root, '01_fastp_out')
    fixed_bb  = os.path.join(run_root, 'BB_information_fixed.tsv')
    if not os.path.isdir(merged_dir):
        return False
    try:
        items = os.listdir(merged_dir)
    except Exception:
        items = []
    has_fastq = any(
        s.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')) for s in items
    )
    return has_fastq and _file_or_gz_exists(fixed_bb)


def _has_merged_fastqs(run_root: str) -> bool:
    """Return True when merged FASTQs exist (regardless of BB metadata)."""
    merged_dir = os.path.join(run_root, '01_fastp_out')
    if not os.path.isdir(merged_dir):
        return False
    try:
        items = os.listdir(merged_dir)
    except Exception:
        items = []
    return any(s.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')) for s in items)


def _has_decoded_outputs(run_root: str) -> bool:
    """Return True when raw_counts_matrix.tsv exists (required for hit calling)."""
    decoded_dir = os.path.join(run_root, '02_decoded')
    if not os.path.isdir(decoded_dir):
        return False
    raw = os.path.join(decoded_dir, 'raw_counts_matrix.tsv')
    return _file_or_gz_exists(raw)


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


def _estimate_rows(path: str) -> int:
    """Heuristically estimate row count for (possibly gzipped) TSV, ignoring header."""
    try:
        import subprocess, gzip, shlex
        if os.path.isfile(path):
            if shutil.which('wc'):
                out = subprocess.check_output(['wc', '-l', path], text=True).strip().split()[0]
                return max(0, int(out) - 1)
            # Fallback: Python loop (may be slow on huge files)
            cnt = 0
            with open(path, 'rb') as f:
                for _ in f:
                    cnt += 1
            return max(0, cnt - 1)
        gz = path if path.endswith('.gz') and os.path.isfile(path) else path + '.gz'
        if os.path.isfile(gz):
            if shutil.which('gzip') and shutil.which('wc'):
                cmd = f"gzip -cd {shlex.quote(gz)} | wc -l"
                out = subprocess.check_output(['bash', '-lc', cmd], text=True).strip()
                return max(0, int(out) - 1)
            cnt = 0
            with gzip.open(gz, 'rt', encoding='utf-8', errors='ignore') as f:
                for _ in f:
                    cnt += 1
            return max(0, cnt - 1)
    except Exception:
        return -1
    return -1


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


def _build_hit_payload(args, r1_cols, r2_cols, neg_cols, del2_col, run_root_abs, hit_out_abs, hitter_path):
    payload = {
        "version": 1,
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
        "args": {k: _jsonable(v) for k, v in sorted(vars(args).items())},
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


def _build_preproc_payload(args, preproc_path: str, fastq_dir: str, bbinfo: str):
    return {
        "version": 1,
        "script": os.path.basename(preproc_path),
        "script_mtime": os.path.getmtime(preproc_path) if os.path.isfile(preproc_path) else None,
        "fastq_dir": os.path.abspath(fastq_dir) if fastq_dir else "",
        "fastq_fingerprint": _fastq_dir_fingerprint(fastq_dir) if fastq_dir else [],
        "bbinfo": _file_fingerprint(bbinfo),
        "mismatch": args.mismatch,
        "skip_fastp": bool(args.skip_fastp),
        "threads": int(args.threads),
    }


def _build_decode_payload(decode_path: str, merged_dir: str, fixed_bb: str, preproc_hash: str):
    return {
        "version": 1,
        "script": os.path.basename(decode_path),
        "script_mtime": os.path.getmtime(decode_path) if os.path.isfile(decode_path) else None,
        "merged_dir": os.path.abspath(merged_dir),
        "fixed_bb": _file_fingerprint(fixed_bb),
        "preprocess_hash": preproc_hash,
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


# (deduplicated helper definitions)


def build_parser():
    p = argparse.ArgumentParser(
        description="DELeGANce — End‑to‑End pipeline runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Core inputs (REQUIRED by default)
    p.add_argument('--fastq-dir', required=False,
                   help='Directory containing paired FASTQs (required if running preprocess)')
    p.add_argument('--bbinfo', required=False,
                   help='BB information file (required if running preprocess)')
    p.add_argument('--output-dir', required=True,
                   help='Output run root directory (absolute or relative)')
    p.add_argument('--threads', type=int, default=int(os.environ.get('FASTP_THREADS', 4)),
                   help='Threads for fastp')
    p.add_argument('--mismatch', choices=['none', 'hp_op_cp'], default='hp_op_cp',
                   help='Preprocess mismatch mode for HP/OP/CP')
    p.add_argument('--skip-fastp', action='store_true', help='Skip fastp stage (use existing merged fastqs if any)')
    p.add_argument('--force-preprocess', action='store_true', help='Force rerun preprocess even if outputs/cache exist')

    # All-in-one output location
    p.add_argument('--hit-out', default=None,
                   help='All-in-one output dir (default: RUN_ROOT/03_normalized)')

    # Phase toggles
    p.add_argument('--only', choices=['all', 'preprocess', 'decode', 'hit'], default='all',
                   help='Run only a specific phase (or all)')
    p.add_argument('--dry-run', action='store_true', help='Print commands without executing')
    # Toggle-able stop-on-error
    p.add_argument('--stop-on-error', dest='stop_on_error', action='store_true', default=True,
                   help='Stop pipeline on first error (default true)')
    p.add_argument('--no-stop-on-error', dest='stop_on_error', action='store_false',
                   help='Do not stop on the first error')

    # REQUIRED: Matrix column names for hit-caller
    # Accept comma- or space-separated lists for R1/R2; NEG supports 1 or 2 names (R1, optional R2)
    p.add_argument('--r1', nargs='+', required=True, help='R1 columns (one or more)')
    p.add_argument('--r2', nargs='+', required=False, help='R2 columns (optional; omit for R1-only)')
    p.add_argument('--neg', nargs='+', required=True, help='NEG columns: R1 [R2 optional]')
    p.add_argument('--del2', required=True, help='DEL2 column name')

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

    # Resolve paths
    run_root_abs = os.path.abspath(os.path.join(THIS_DIR, args.output_dir)) if not os.path.isabs(args.output_dir) else args.output_dir
    run_root_rel = os.path.relpath(run_root_abs, THIS_DIR)
    ensure_dir(run_root_abs)
    # Default normalized root; the final per-run output directory will be set after auto-opt
    norm_root_abs = os.path.join(run_root_abs, '03_normalized')
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

    # Init flags before dependency checks
    preprocess_done = False
    decode_done = False

    # Dependency checks (fastp only if we will run preprocess)
    if not which('perl'):
        print('[ERROR] perl not found in PATH'); return 2
    if not which('python3'):
        print('[ERROR] python3 not found in PATH'); return 2
    if (args.only in ('all', 'preprocess')) and (not preprocess_done) and (not args.skip_fastp) and not which('fastp'):
        print('[WARN] fastp not found; switching to --skip-fastp for preprocess')
        args.skip_fastp = True

    # Validate required inputs only if we will run preprocess
    if args.only in ('all', 'preprocess') and not preprocess_done:
        if not args.fastq_dir or not os.path.isdir(args.fastq_dir):
            print('[ERROR] --fastq-dir is required for preprocess and must exist (or provide existing outputs under RUN_ROOT/01_fastp_out)')
            return 2
        if not args.bbinfo or not os.path.isfile(args.bbinfo):
            print('[ERROR] --bbinfo is required for preprocess and must point to a file (or provide existing RUN_ROOT/BB_information_fixed.tsv)')
            return 2

    # Normalize R1/R2/NEG lists (split on commas and flatten)
    def _split_list(lst):
        out = []
        for item in lst or []:
            out.extend([t for t in str(item).split(',') if t != ''])
        return [x.strip() for x in out if x.strip()]
    r1_cols = _split_list(args.r1)
    r2_cols = _split_list(args.r2)
    neg_cols = _split_list(args.neg)
    if len(r1_cols) == 0:
        print('[ERROR] --r1 must specify at least one column (R1)')
        return 2
    if len(neg_cols) == 0:
        print('[ERROR] --neg must specify at least one column (NEG for R1; optional second for R2)')
        return 2
    del2_col = args.del2.strip()
    if del2_col == '':
        print('[ERROR] --del2 must not be empty')
        return 2

    # Logs
    master_log = os.path.join(run_root_abs, '00_pipeline.log')
    print(f"[INFO] Pipeline log → {master_log}")

    # Preprocess cache/hash setup
    preproc_out_dir = os.path.join(run_root_abs, '01_fastp_out')
    preproc_params_path = os.path.join(preproc_out_dir, 'preprocess_params.json')
    preprocess_done = _has_fastp_outputs(run_root_abs)
    preproc_payload = _build_preproc_payload(args, preproc, args.fastq_dir or "", args.bbinfo or "")
    preproc_hash = _hash_payload(preproc_payload)
    preproc_payload["hash"] = preproc_hash
    cached_preproc = _load_json(preproc_params_path)
    preproc_cache_ok = preprocess_done and isinstance(cached_preproc, dict) and str(cached_preproc.get("hash", "")) == preproc_hash
    if preprocess_done and not preproc_cache_ok:
        cached_hash = cached_preproc.get("hash") if isinstance(cached_preproc, dict) else None
        if cached_hash:
            print(f"[INFO] Preprocess outputs found but parameter hash differs (cached={cached_hash[:12]}, current={preproc_hash[:12]}). Will rerun.")
        else:
            print("[INFO] Preprocess outputs found but cache missing/invalid. Will rerun.")
        preprocess_done = False

    # Decode cache/hash setup
    decode_out_dir = os.path.join(run_root_abs, '02_decoded')
    decode_params_path = os.path.join(decode_out_dir, 'decode_params.json')
    decode_done = _has_decoded_outputs(run_root_abs)
    merged_dir = os.path.join(run_root_abs, '01_fastp_out')
    fixed_bb   = os.path.join(run_root_abs, 'BB_information_fixed.tsv')
    preproc_hash_for_decode = preproc_hash if preproc_cache_ok else None
    decode_payload = _build_decode_payload(decode, merged_dir, fixed_bb, preproc_hash_for_decode)
    decode_hash = _hash_payload(decode_payload)
    decode_payload["hash"] = decode_hash
    cached_decode = _load_json(decode_params_path)
    decode_cache_ok = decode_done and isinstance(cached_decode, dict) and str(cached_decode.get("hash", "")) == decode_hash
    if decode_done and not decode_cache_ok:
        cached_hash = cached_decode.get("hash") if isinstance(cached_decode, dict) else None
        if cached_hash:
            print(f"[INFO] Decode outputs found but parameter hash differs (cached={cached_hash[:12]}, current={decode_hash[:12]}). Will rerun.")
        else:
            print("[INFO] Decode outputs found but cache missing/invalid. Will rerun.")
        decode_done = False

    # If we are skipping fastp, require merged FASTQs to already exist when decode/hit will run.
    if args.skip_fastp and args.only in ('all', 'decode', 'hit'):
        if not _has_merged_fastqs(run_root_abs):
            print('[ERROR] --skip-fastp requested but no merged FASTQs found under RUN_ROOT/01_fastp_out.')
            print('        Either run fastp (omit --skip-fastp) or provide merged FASTQs in 01_fastp_out first.')
            return 2

    # 1) Preprocess
    rc = 0
    if args.only in ('all', 'preprocess'):
        if preprocess_done and not args.force_preprocess:
            print(f"[INFO] Preprocess outputs + cache match (hash={preproc_hash[:12]}); skipping preprocess.")
        else:
            cmd = ['perl', preproc]
            cmd += ['-b', args.bbinfo]
            # pass absolute run_root to preprocess; it detects absolute and uses it directly
            cmd += ['-f', args.fastq_dir, '-o', run_root_abs, '-t', str(int(args.threads)), '--mismatch', args.mismatch]
            if args.skip_fastp:
                cmd += ['--skip-fastp']
            print(f"[INFO] Preprocess → RUN_ROOT={run_root_abs}")
            if args.dry_run:
                print('DRY-RUN:', ' '.join(cmd))
            else:
                rc = run_cmd(cmd, master_log)
                if rc != 0 and args.stop_on_error:
                    print('[ERROR] Preprocess failed; aborting.')
                    return rc
                _write_json(preproc_params_path, preproc_payload)
                print(f"[OK] Stored preprocess parameter hash: {preproc_hash[:12]}")

    # 2) Decode
    if rc == 0 and args.only in ('all', 'decode'):
        if decode_done and not args.force_decode:
            print(f"[INFO] Decode outputs + cache match (hash={decode_hash[:12]}); skipping decode.")
        else:
            out_dir    = decode_out_dir
            cmd = ['perl', decode, '--merged-dir', merged_dir, '--fixed-bb-file', fixed_bb, '--out-dir', out_dir]
            print(f"[INFO] Decode   → MERGED={merged_dir}  FIXED_BB={fixed_bb}  OUT={out_dir}")
            if args.dry_run:
                print('DRY-RUN:', ' '.join(cmd))
            else:
                rc = run_cmd(cmd, master_log)
                if rc != 0 and args.stop_on_error:
                    print('[ERROR] Decode failed; aborting.')
                    return rc
                _write_json(decode_params_path, decode_payload)
                print(f"[OK] Stored decode parameter hash: {decode_hash[:12]}")

    # Auto-optimization: choose device/dtype/GLM mode heuristically if not specified
    if int(getattr(args, 'auto_opt', 1)) == 1:
        # Device/dtype
        if (args.device is None):
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
                args.glm_mode = 'top'; args.glm_top_pct = args.glm_top_pct or 0.5; args.glm_top_k = args.glm_top_k or 500_000
                print(f"[AUTO] Large matrix (~{nrows:,} rows) → GLM top mode: pct={args.glm_top_pct}%, cap={args.glm_top_k}")
            elif nrows >= 1_000_000:
                args.glm_mode = 'top'; args.glm_top_pct = args.glm_top_pct or 2.0; args.glm_top_k = args.glm_top_k or 200_000
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
    hit_out_abs = hit_out if os.path.isabs(hit_out) else os.path.join(THIS_DIR, hit_out)
    ensure_dir(hit_out_abs)

    # Build current hit parameters payload + hash (used for cache validation)
    hit_payload = _build_hit_payload(args, r1_cols, r2_cols, neg_cols, del2_col, run_root_abs, hit_out_abs, hitter)
    hit_payload["decode_hash"] = decode_hash
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
            print(f"[INFO] Cached HitCaller outputs found but hit_params.json missing or invalid. Recomputing.")
    used_hit_dir = None

    # 3) All‑in‑one hit caller
    if rc == 0 and args.only in ('all', 'hit'):
        if cached_hit_ok and not getattr(args, 'force_hit', False):
            used_hit_dir = hit_out_abs
            print(f"[INFO] Cached HitCaller outputs match parameters (hash={current_hash[:12]}). Skipping recompute.")
        else:
            cmd = ['python3', hitter, '--run_root', run_root_abs, '--outdir', hit_out]
        # Column names
            cmd += ['--del2_col', del2_col]
            cmd += ['--r1_cols'] + r1_cols
            if len(r2_cols) > 0:
                cmd += ['--r2_cols'] + r2_cols
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
                rc = run_cmd(cmd, master_log)
                if rc != 0 and args.stop_on_error:
                    print('[ERROR] HitCaller failed; aborting.')
                    return rc
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
                if any(os.path.isfile(p) for p in legacy_markers):
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
                cmd = [
                    'python3', _pick_script(['04_build_interactive_report.py','04_hit_finder.py','04_hit_finder_250816p2.py']),
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
            presets = []
            try:
                for name in sorted(os.listdir(norm_root_abs)):
                    p = os.path.join(norm_root_abs, name)
                    if os.path.isdir(p):
                        presets.append((name, p))
            except Exception:
                pass
            # include base 03_normalized if it has outputs
            def _has_outputs(d):
                return os.path.isfile(os.path.join(d, '05_hybrid_annot.tsv')) or \
                       os.path.isfile(os.path.join(d, 'report.html'))
            include_base = _has_outputs(norm_root_abs)
            # Build HTML
            def _h(s: str) -> str:
                try:
                    return _html_escape(str(s), quote=True)
                except Exception:
                    return str(s)

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
                f"<h1>DELeGANce — Run Results</h1>",
                f"<p><b>Run root:</b> {_h(run_root_rel)}</p>",
                f"<p>Normalized outputs are under <code>{_h(rel_norm_root)}</code></p>",
            ]
            # Base section (no preset)
            if include_base:
                rel_base = os.path.relpath(norm_root_abs, run_root_abs)
                html += [
                    f"<section class='preset'><h2>Preset: (default)</h2>",
                    "<ul>",
                    f"<li><a href='{_h(rel_base)}/report.html'>All-in-one report</a></li>",
                    f"<li><a href='{_h(rel_base)}/05_hybrid_annot.tsv'>05_hybrid_annot.tsv</a></li>",
                    f"<li><a href='{_h(rel_base)}/08_topk_consensus.tsv'>08_topk_consensus.tsv</a></li>",
                    f"<li><a href='{_h(rel_base)}/interactive_hits.html'>Interactive (embedded below)</a></li>",
                    f"<li>Pop-outs: <a href='{_h(rel_base)}/interactive_p1.html'>BB1×BB2</a> · <a href='{_h(rel_base)}/interactive_p2.html'>BB1×BB3</a> · <a href='{_h(rel_base)}/interactive_p3.html'>BB2×BB3</a> · <a href='{_h(rel_base)}/interactive_p4.html'>BB1×BB4</a> · <a href='{_h(rel_base)}/interactive_table.html'>Top table</a></li>",
                    "</ul>",
                    f"<p>Interactive (if present): <a href='{_h(rel_base)}/interactive_hits.html'>interactive_hits.html</a></p>",
                    f"<iframe src='{_h(rel_base)}/interactive_hits.html'></iframe>",
                    "</section>",
                ]
            # Each preset subdir
            for name, absdir in presets:
                rel = os.path.relpath(absdir, run_root_abs)
                # Derive preset and tags from folder name convention
                if '_' in name:
                    preset_name = name.split('_', 1)[0]
                    tag_str = name[len(preset_name)+1:]
                else:
                    preset_name = name
                    tag_str = ''
                tag_html = f" <span class='tags'>({_h(tag_str)})</span>" if tag_str else ''
                html += [
                    f"<section class='preset'><h2>Preset: {_h(preset_name)}{tag_html}</h2>",
                    "<ul>",
                    f"<li><a href='{_h(rel)}/report.html'>All-in-one report</a></li>",
                    f"<li><a href='{_h(rel)}/05_hybrid_annot.tsv'>05_hybrid_annot.tsv</a></li>",
                    f"<li><a href='{_h(rel)}/08_topk_consensus.tsv'>08_topk_consensus.tsv</a></li>",
                    f"<li><a href='{_h(rel)}/06_topk_glm.tsv'>06_topk_glm.tsv</a></li>",
                    f"<li><a href='{_h(rel)}/07_topk_rs.tsv'>07_topk_rs.tsv</a></li>",
                    f"<li><a href='{_h(rel)}/interactive_hits.html'>Interactive (embedded below)</a></li>",
                    f"<li>Pop-outs: <a href='{_h(rel)}/interactive_p1.html'>BB1×BB2</a> · <a href='{_h(rel)}/interactive_p2.html'>BB1×BB3</a> · <a href='{_h(rel)}/interactive_p3.html'>BB2×BB3</a> · <a href='{_h(rel)}/interactive_p4.html'>BB1×BB4</a> · <a href='{_h(rel)}/interactive_table.html'>Top table</a></li>",
                    "</ul>",
                    f"<p>Interactive: <a href='{_h(rel)}/interactive_hits.html'>interactive_hits.html</a></p>",
                    f"<iframe src='{_h(rel)}/interactive_hits.html'></iframe>",
                    "</section>",
                ]
            html += ["</body></html>"]
            with open(index_html, 'w', encoding='utf-8') as f:
                f.write("\n".join(html))
            print(f"[OK] Wrote index.html → {index_html}")
        except Exception as e:
            print(f"[WARN] Failed to write index.html: {e}")

    print('\n[OK] Pipeline completed')
    print(f" - RUN_ROOT:  {run_root_abs}")
    print(f" - HIT_OUT:   {hit_out_abs}")
    print(f" - Master log: {master_log}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
