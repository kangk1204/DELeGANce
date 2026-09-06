#!/usr/bin/env python3
import os
import gzip
import json
import argparse
import random
from collections import defaultdict, Counter
from datetime import datetime


def _open_text(path):
    """Open a (possibly gzipped) text file; tolerate non-UTF-8 bytes instead of aborting."""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'r', encoding='utf-8', errors='replace')


def one_bp_substitution(s):
    out = set()
    nucs = 'AGCT'
    for i, orig in enumerate(s):
        for n in nucs:
            if n == orig:
                continue
            out.add(s[:i] + n + s[i+1:])
    return out


def one_bp_insertion(s):
    # NOTE: intentionally mirrors 02_decode_reads.pl::one_bp_insertion, including its
    # unreachable `i == 0` / `i == L` guards (no insertions at either end). Keep in sync with
    # the decoder; do not "fix" here alone or validator/decoder equivalence breaks.
    out = set()
    nucs = 'AGCT'
    L = len(s)
    first = s[0]
    last = s[-1]
    for i in range(1, L):
        for n in nucs:
            if (i == 0 or i == 1) and n == first:
                continue
            if ((i == L) or (i == L-1)) and n == last:
                continue
            out.add(s[:i] + n + s[i:])
    return out


def one_bp_deletion(s):
    out = set()
    L = len(s)
    for i in range(L):
        out.add(s[:i] + s[i+1:])
    return out


def merge_seq_qual(seq_map, seq, qual):
    if seq in seq_map:
        if seq_map[seq] == 'miss' and qual == 'perf':
            seq_map[seq] = 'perf'
    else:
        seq_map[seq] = qual


def build_variant_map(seq, mismatch='hp_op_cp'):
    """Perfect sequence plus (if mismatch != 'none') all 1-bp sub/ins/del variants, as the decoder indexes them."""
    seq = seq.upper()
    out = {}
    merge_seq_qual(out, seq, 'perf')
    if mismatch == 'none':
        return out
    for s in one_bp_substitution(seq):
        merge_seq_qual(out, s, 'miss')
    for s in one_bp_insertion(seq):
        merge_seq_qual(out, s, 'miss')
    for s in one_bp_deletion(seq):
        merge_seq_qual(out, s, 'miss')
    return out


def build_variant_set(seq):
    return set(build_variant_map(seq).keys())


def _build_by_len_map(seq_qual_map):
    by_len = defaultdict(dict)
    for seq, qual in seq_qual_map.items():
        by_len[len(seq)][seq] = qual
    return by_len


def _merge_len_maps(a, b):
    out = {L: dict(d) for L, d in a.items()}
    for L, seqmap in b.items():
        dest = out.setdefault(L, {})
        for seq, qual in seqmap.items():
            if seq in dest:
                if dest[seq] == 'miss' and qual == 'perf':
                    dest[seq] = 'perf'
            else:
                dest[seq] = qual
    return out


def load_bb_sets(bb_path, mismatch='hp_op_cp'):
    hp_sets = defaultdict(set)
    op_sets = defaultdict(set)
    cp_sets = defaultdict(set)
    codon_sets = defaultdict(set)
    codon_sets_by_cycle = defaultdict(lambda: defaultdict(set))
    lib_cycles = defaultdict(set)
    hp_variants = defaultdict(dict)
    op_variants = defaultdict(dict)
    cp_variants = defaultdict(dict)
    cp_owner = {}
    with _open_text(bb_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 7:
                continue
            typ, seq_raw, bb_id_fixed, cycle, tag_id, lib_id, smiles = parts[:7]
            if typ.lower() in ('type', 'bb_type', 'bb-type'):
                continue
            seq = ''.join(seq_raw.split()).upper()
            if not seq or any(c not in 'ACGT' for c in seq):
                continue
            lib_norm = lib_id if lib_id is not None else ''
            if lib_norm.upper() in ('HP', 'OP', 'CP', 'NA', 'NONE'):
                lib_norm = ''
            if typ.lower().startswith('codon'):
                if lib_norm:
                    codon_sets[lib_norm].add(seq)
                    # Same rule as 02_decode_reads.pl (cycle =~ /^[1-4]$/): cycle values outside 1..4 are
                    # ignored for the expected-cycle decision instead of forcing a 4-cycle library.
                    if str(cycle) in ('1', '2', '3', '4'):
                        c = int(cycle)
                        codon_sets_by_cycle[lib_norm][c].add(seq)
                        lib_cycles[lib_norm].add(c)
            elif typ == 'HP':
                variants = build_variant_map(seq, mismatch)
                hp_sets[lib_norm].update(variants.keys())
                for s, q in variants.items():
                    merge_seq_qual(hp_variants[lib_norm], s, q)
            elif typ == 'OP':
                variants = build_variant_map(seq, mismatch)
                op_sets[lib_norm].update(variants.keys())
                for s, q in variants.items():
                    merge_seq_qual(op_variants[lib_norm], s, q)
            elif typ == 'CP':
                if lib_norm:
                    variants = build_variant_map(seq, mismatch)
                    cp_sets[lib_norm].update(variants.keys())
                    for s, q in variants.items():
                        if s in cp_owner and cp_owner[s] != lib_norm:
                            raise ValueError(f"CP collision across libs (variant={s}): {cp_owner[s]} vs {lib_norm}")
                        cp_owner[s] = lib_norm
                        merge_seq_qual(cp_variants[lib_norm], s, q)

    lib_expected_cycles = {}
    for lib, cycles in lib_cycles.items():
        if cycles:
            lib_expected_cycles[lib] = 4 if max(cycles) >= 4 else 3
        else:
            lib_expected_cycles[lib] = 3
    for lib in cp_sets.keys():
        lib_expected_cycles.setdefault(lib, 3)

    hp_by_len_global = _build_by_len_map(hp_variants.get('', {}))
    op_by_len_global = _build_by_len_map(op_variants.get('', {}))

    hp_by_len_merged = {}
    op_by_len_merged = {}
    for lib in lib_expected_cycles.keys():
        hp_by_len_merged[lib] = _merge_len_maps(hp_by_len_global, _build_by_len_map(hp_variants.get(lib, {})))
        op_by_len_merged[lib] = _merge_len_maps(op_by_len_global, _build_by_len_map(op_variants.get(lib, {})))

    cp_by_len = defaultdict(dict)
    for lib, seqmap in cp_variants.items():
        for seq, qual in seqmap.items():
            L = len(seq)
            if seq in cp_by_len[L]:
                prev_lib, prev_qual = cp_by_len[L][seq]
                if prev_lib != lib:
                    raise ValueError(f"CP collision across libs (variant={seq}): {prev_lib} vs {lib}")
                if prev_qual == 'miss' and qual == 'perf':
                    cp_by_len[L][seq] = (lib, qual)
            else:
                cp_by_len[L][seq] = (lib, qual)

    codon_by_len_by_cycle = defaultdict(lambda: defaultdict(dict))
    for lib, cycles in codon_sets_by_cycle.items():
        for c, seqs in cycles.items():
            by_len = defaultdict(set)
            for seq in seqs:
                by_len[len(seq)].add(seq)
            codon_by_len_by_cycle[lib][c] = by_len

    return {
        'hp_sets': hp_sets,
        'op_sets': op_sets,
        'cp_sets': cp_sets,
        'codon_sets': codon_sets,
        'codon_sets_by_cycle': codon_sets_by_cycle,
        'lib_expected_cycles': lib_expected_cycles,
        'hp_by_len_merged': hp_by_len_merged,
        'op_by_len_merged': op_by_len_merged,
        'cp_by_len': cp_by_len,
        'codon_by_len_by_cycle': codon_by_len_by_cycle,
    }


def hp_allowed(hp_sets, lib_id):
    s = set(hp_sets.get('', set()))
    s.update(hp_sets.get(lib_id, set()))
    return s


def op_allowed(op_sets, lib_id):
    s = set(op_sets.get('', set()))
    s.update(op_sets.get(lib_id, set()))
    return s


def rc_seq(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def parse_active_centers(raw):
    if raw is None:
        return []
    s = str(raw).strip()
    if not s or s.upper() == 'NA':
        return []
    out = []
    for part in s.split(','):
        part = part.strip()
        if not part:
            continue
        if part.isdigit():
            out.append(int(part))
    return out


def classify_len_py(L, exp3, exp4, tol, active_mode, active_centers):
    in3 = abs(L - exp3) <= tol
    in4 = abs(L - exp4) <= tol
    if active_mode == '3':
        in4 = False
    elif active_mode == '4':
        in3 = False

    if active_mode != 'multi':
        if in3 and in4:
            return 'within_both'
        if in3:
            return 'within_3'
        if in4:
            return 'within_4'
        return 'outside'

    in_multi = any(abs(L - c) <= tol for c in active_centers)
    if in_multi:
        if in3 and in4:
            return 'within_both'
        if in3:
            return 'within_3'
        if in4:
            return 'within_4'
        return 'within_multi'
    return 'outside'


def len_filter_reject_py(cls, policy):
    if policy == 'reject_outside':
        return cls == 'outside'
    if policy == 'reject_inside':
        return cls.startswith('within_')
    return cls == 'outside'


def qual_rank(q):
    return 1 if q == 'perf' else 0


def scan_hits(seq, by_len_map):
    hits = []
    n = len(seq)
    for L, seqmap in by_len_map.items():
        if L > n:
            continue
        for i in range(0, n - L + 1):
            sub = seq[i:i+L]
            if isinstance(seqmap, set):
                if sub in seqmap:
                    hits.append({'pos': i, 'end': i + L - 1, 'seq': sub, 'qual': 'perf'})
            else:
                q = seqmap.get(sub)
                if q:
                    hits.append({'pos': i, 'end': i + L - 1, 'seq': sub, 'qual': q})
    return hits


def scan_cp_hits(seq, cp_by_len):
    hits = []
    n = len(seq)
    for L, seqmap in cp_by_len.items():
        if L > n:
            continue
        for i in range(0, n - L + 1):
            sub = seq[i:i+L]
            val = seqmap.get(sub)
            if val:
                lib_id, q = val
                hits.append({'pos': i, 'end': i + L - 1, 'seq': sub, 'qual': q, 'lib_id': lib_id})
    return hits


def sort_hits(hits, max_n=None):
    # Same decision rule as 02_decode_reads.pl pick_best_*: perfect matches first, then larger position first.
    # Ties (same qual and pos, different match length) follow the decoder's trie scan order (shorter match
    # first) so that --max-cp-cands/--max-anchor-cands truncate the same candidate set as the decoder.
    hits.sort(key=lambda h: (-qual_rank(h.get('qual')), -h.get('pos', 0), len(h.get('seq', ''))))
    if max_n is not None and len(hits) > max_n:
        return hits[:max_n]
    return hits


def gap_ok(next_start, prev_end, tol):
    gap = next_start - (prev_end + 1)
    return abs(gap) <= tol


def codon_chain_exists(by_cycle, cycle, expected, prev_end, cp_pos, tol, memo):
    key = (cycle, prev_end)
    if key in memo:
        return memo[key]
    if cycle > expected:
        ok = gap_ok(cp_pos, prev_end, tol)
        memo[key] = ok
        return ok
    for cand in by_cycle.get(cycle, []):
        if not gap_ok(cand['pos'], prev_end, tol):
            continue
        if codon_chain_exists(by_cycle, cycle + 1, expected, cand['end'], cp_pos, tol, memo):
            memo[key] = True
            return True
    memo[key] = False
    return False


def pick_codons_by_adjacency_py(seq, op_end, cp_pos, expected_cycles, adj_tol, codon_by_len_by_cycle):
    by_cycle = {}
    for c in range(1, expected_cycles + 1):
        by_len = codon_by_len_by_cycle.get(c)
        if not by_len:
            return False
        hits = scan_hits(seq, by_len)
        filt = []
        for h in hits:
            if h['pos'] < op_end - adj_tol:
                continue
            if h['end'] > cp_pos + adj_tol:
                continue
            filt.append(h)
        if not filt:
            return False
        by_cycle[c] = sorted(filt, key=lambda h: h['pos'])

    memo = {}
    return codon_chain_exists(by_cycle, 1, expected_cycles, op_end, cp_pos, adj_tol, memo)


def _cap(n):
    # 0 (decoder default for --max-cp-cands/--max-anchor-cands) means "no limit"
    return None if not n or n <= 0 else n


def decode_attempt_orient(seq, adj_tol, bb, max_cp_cands=0, max_anchor_cands=0):
    cp_hits = scan_cp_hits(seq, bb['cp_by_len'])
    if not cp_hits:
        return False
    cp_hits = sort_hits(cp_hits, max_n=_cap(max_cp_cands))

    for cp in cp_hits:
        lib_id = cp['lib_id']
        expected_cycles = bb['lib_expected_cycles'].get(lib_id)
        if expected_cycles is None:
            continue
        op_map = bb['op_by_len_merged'].get(lib_id)
        hp_map = bb['hp_by_len_merged'].get(lib_id)
        codon_by_cycle = bb['codon_by_len_by_cycle'].get(lib_id)
        if not op_map or not hp_map or not codon_by_cycle:
            continue

        op_hits_all = scan_hits(seq, op_map)
        op_hits = [h for h in op_hits_all if h['pos'] < cp['pos']]
        if not op_hits:
            continue
        op_hits = sort_hits(op_hits, max_n=_cap(max_anchor_cands))

        hp_hits_all = scan_hits(seq, hp_map)
        if not hp_hits_all:
            continue

        for op in op_hits:
            hp_hits = [h for h in hp_hits_all if h['pos'] < op['pos']]
            if not hp_hits:
                continue
            hp_hits = sort_hits(hp_hits, max_n=_cap(max_anchor_cands))
            for hp in hp_hits:
                if not gap_ok(op['pos'], hp['end'], adj_tol):
                    continue
                op_end = op['end']
                if pick_codons_by_adjacency_py(seq, op_end, cp['pos'], expected_cycles, adj_tol, codon_by_cycle):
                    return True
    return False


def decode_attempt(seq, adj_tol, bb, max_cp_cands=0, max_anchor_cands=0):
    if decode_attempt_orient(seq, adj_tol, bb, max_cp_cands, max_anchor_cands):
        return True
    rc = rc_seq(seq)
    return decode_attempt_orient(rc, adj_tol, bb, max_cp_cands, max_anchor_cands)


def reservoir_sample(path, sample_n, seed):
    rng = random.Random(seed)
    sample = []
    with _open_text(path) as f:
        header = f.readline()
        for i, line in enumerate(f):
            if len(sample) < sample_n:
                sample.append(line)
            else:
                j = rng.randint(0, i)
                if j < sample_n:
                    sample[j] = line
    return header, sample


def validate_undecoded_sample(path, sample_n, seed, adj_tol, len_conf, bb, out_path=None,
                              max_cp_cands=0, max_anchor_cands=0):
    header, sample = reservoir_sample(path, sample_n, seed)
    if not sample:
        return {'sampled_n': 0}
    cols = header.rstrip('\n').split('\t')
    idx = {h: i for i, h in enumerate(cols)}
    counts = Counter()
    rows = []

    for line in sample:
        row = line.rstrip('\n').split('\t')
        if len(row) != len(cols):
            counts['row_len_mismatch'] += 1
            continue
        counts['sampled_n'] += 1
        read_id = row[idx['read_id']] if 'read_id' in idx else 'NA'
        if 'read_seq' not in idx:
            counts['missing_read_seq'] += 1
            continue
        read_seq = row[idx['read_seq']]
        if 'read_len' in idx:
            try:
                read_len = int(row[idx['read_len']])
            except Exception:
                counts['read_len_parse'] += 1
                read_len = len(read_seq)
        else:
            read_len = len(read_seq)
        reason = row[idx['reason']] if 'reason' in idx else ''

        len_reject = False
        if len_conf:
            cls = classify_len_py(read_len, len_conf['exp3'], len_conf['exp4'], len_conf['tol'],
                                  len_conf['active_mode'], len_conf['active_centers'])
            len_reject = len_filter_reject_py(cls, len_conf['policy'])

        decode_possible = False
        if len_reject:
            counts['len_reject'] += 1
        else:
            decode_possible = decode_attempt(read_seq, adj_tol, bb, max_cp_cands, max_anchor_cands)
            if decode_possible:
                counts['unexpected_decoded'] += 1
            else:
                counts['expected_fail'] += 1

        # reason_len_mismatch = decoder said len_out_of_range XOR the validator's re-classification rejects.
        # Without length_window_stats (len_conf None) len_reject is always False and every decoder length
        # rejection would be counted as a mismatch, so the comparison is skipped and flagged instead.
        if len_conf:
            reason_has_len = 'len_out_of_range' in reason
            if reason_has_len != len_reject:
                counts['reason_len_mismatch'] += 1
        else:
            counts['len_conf_missing'] = 1

        if out_path:
            rows.append({
                'read_id': read_id,
                'read_len': str(read_len),
                'reason': reason,
                'len_reject': '1' if len_reject else '0',
                'decode_possible': '1' if decode_possible else '0',
            })

    if out_path and rows:
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write('\t'.join(['read_id', 'read_len', 'reason', 'len_reject', 'decode_possible']) + '\n')
            for r in rows:
                f.write('\t'.join([r['read_id'], r['read_len'], r['reason'], r['len_reject'], r['decode_possible']]) + '\n')

    return dict(counts)


def read_single_row_tsv(path):
    with _open_text(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        row = f.readline().rstrip('\n').split('\t')
    return {header[i]: row[i] for i in range(min(len(header), len(row)))}


def read_table_by_sample(path):
    if not os.path.exists(path):
        return {}
    with _open_text(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {h: i for i, h in enumerate(header)}
        s_i = idx.get('sample')
        if s_i is None:
            return {}
        out = {}
        for line in f:
            row = line.rstrip('\n').split('\t')
            if len(row) <= s_i:
                continue
            sample = row[s_i]
            out[sample] = {header[i]: row[i] if i < len(row) else '' for i in range(len(header))}
        return out


def _strip_tsv_suffix(fn, prefix):
    """Return the sample name for '<prefix><sample>.tsv[.gz]', or None if fn does not match."""
    if not fn.startswith(prefix):
        return None
    if fn.endswith('.tsv.gz'):
        return fn[len(prefix):-7]
    if fn.endswith('.tsv'):
        return fn[len(prefix):-4]
    return None


def _find_sample_file(sdir, prefix, sample):
    """Locate '<prefix><sample>.tsv' (or .tsv.gz) in sdir; None if absent."""
    for suffix in ('.tsv', '.tsv.gz'):
        cand = os.path.join(sdir, f"{prefix}{sample}{suffix}")
        if os.path.exists(cand):
            return cand
    return None


def list_samples_from_decoded(decoded_dir):
    samples = []
    for fn in os.listdir(decoded_dir):
        name = _strip_tsv_suffix(fn, 'decoded_reads_')
        if name is not None:
            samples.append(name)
    return sorted(samples)


def read_reason_counts(path):
    counts = Counter()
    with _open_text(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {h:i for i,h in enumerate(header)}
        r_i = idx.get('reason')
        if r_i is None:
            return counts
        for line in f:
            row = line.rstrip('\n').split('\t')
            if len(row) <= r_i:
                continue
            counts[row[r_i]] += 1
    return counts


def validate_decoded_file(path, adj_tol, hp_sets, op_sets, cp_sets, codon_sets,
                          lib_expected_cycles, len_conf, progress_every=1000000, log=None,
                          max_lines=0, sample_n=0, sample_seed=2026):
    counts = Counter()
    with _open_text(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {h:i for i,h in enumerate(header)}
        line_no = 0
        processed = 0  # rows actually examined (excludes the row that triggered the max_lines break)
        if sample_n and sample_n > 0:
            rng = random.Random(sample_seed)
            sample_lines = []
            for i, line in enumerate(f, 1):
                if len(sample_lines) < sample_n:
                    sample_lines.append(line)
                else:
                    j = rng.randint(0, i - 1)
                    if j < sample_n:
                        sample_lines[j] = line
            if log:
                log(f"  sampled {len(sample_lines):,} decoded reads (reservoir)")
            lines_iter = enumerate(sample_lines, 1)
        else:
            lines_iter = enumerate(f, 1)

        for line_no, line in lines_iter:
            if not sample_n and max_lines and line_no > max_lines:
                break
            processed += 1
            row = line.rstrip('\n').split('\t')
            if len(row) != len(header):
                counts['row_len_mismatch'] += 1
                continue
            counts['n'] += 1
            errs = []

            try:
                read_len = int(row[idx['read_len']])
                read_seq_raw = row[idx['read_seq']]
            except Exception:
                counts['parse_error'] += 1
                continue

            if read_len != len(read_seq_raw):
                errs.append('read_len_mismatch')

            direction = row[idx['direction']]
            read_seq = rc_seq(read_seq_raw) if direction == 'revcomp' else read_seq_raw

            lib_id = row[idx['lib_id']]
            if len_conf:
                cls = classify_len_py(read_len, len_conf['exp3'], len_conf['exp4'], len_conf['tol'],
                                      len_conf['active_mode'], len_conf['active_centers'])
                if len_filter_reject_py(cls, len_conf['policy']):
                    errs.append('length_filter_reject')

            try:
                cycles = int(row[idx['cycles']])
            except Exception:
                errs.append('cycles_parse')
                cycles = 0
            if cycles < 0 or cycles > 4:
                # decoded header only carries C1..C4; guard against corrupt rows (would KeyError below)
                errs.append('cycles_out_of_range')
                cycles = 0
            if cycles:
                expected = lib_expected_cycles.get(lib_id)
                if expected is None:
                    errs.append('lib_unknown')
                elif cycles != expected:
                    errs.append('cycles_mismatch')

            # ruff B023: closes over `row`/`read_seq` of the current iteration but is called
            # immediately below within the same iteration, so late binding cannot occur.
            def check_anchor(tag):
                pos = row[idx[f'{tag}_pos']]
                ln = row[idx[f'{tag}_len']]
                seq = row[idx[f'{tag}_seq']]
                if pos == 'NA' or ln == 'NA' or seq == 'NA':
                    return [f'missing_{tag.lower()}'], None, None, None
                try:
                    pos_i = int(pos); ln_i = int(ln)
                except Exception:
                    return [f'parse_{tag.lower()}'], None, None, None
                if pos_i < 0 or pos_i + ln_i > len(read_seq):
                    return [f'{tag}_pos_oob'], pos_i, ln_i, seq
                if read_seq[pos_i:pos_i+ln_i] != seq:
                    return [f'{tag}_seq_mismatch'], pos_i, ln_i, seq
                return [], pos_i, ln_i, seq

            e, hp_pos, hp_len, hp_seq = check_anchor('HP'); errs += e
            e, op_pos, op_len, op_seq = check_anchor('OP'); errs += e
            e, cp_pos, cp_len, cp_seq = check_anchor('CP'); errs += e

            cpos = []
            clen = []
            cseq = []
            cidx = []  # cycle number of each retained codon (a missing/corrupt cycle leaves a hole)
            for i in range(1, cycles + 1):
                p = row[idx[f'C{i}_pos']]
                l = row[idx[f'C{i}_len']]
                s = row[idx[f'C{i}_seq']]
                if p == 'NA' or l == 'NA' or s == 'NA':
                    errs.append(f'missing_c{i}')
                    continue
                try:
                    p_i = int(p); l_i = int(l)
                except Exception:
                    errs.append(f'parse_c{i}')
                    continue
                if p_i < 0 or p_i + l_i > len(read_seq):
                    errs.append(f'c{i}_pos_oob')
                    continue
                if read_seq[p_i:p_i+l_i] != s:
                    errs.append(f'c{i}_seq_mismatch')
                cpos.append(p_i); clen.append(l_i); cseq.append(s); cidx.append(i)

            if hp_pos is not None and op_pos is not None and hp_pos > op_pos:
                errs.append('order_violation')
            prev = op_pos
            for p in cpos:
                if prev is not None and p < prev:
                    errs.append('order_violation')
                prev = p
            if cpos and cp_pos is not None and cpos[-1] > cp_pos:
                errs.append('order_violation')

            # ruff B023: closes over `errs` of the current iteration; called immediately below
            # (same iteration), so this is safe. Named gap_check to avoid shadowing module-level gap_ok().
            def gap_check(start_pos, start_len, next_pos, label):
                if start_pos is None or start_len is None or next_pos is None:
                    return
                gap = next_pos - (start_pos + start_len)
                if abs(gap) > adj_tol:
                    errs.append(label)

            if hp_pos is not None and hp_len is not None and op_pos is not None:
                gap_check(hp_pos, hp_len, op_pos, 'gap_HP_OP')
            if op_pos is not None and op_len is not None and cpos:
                gap_check(op_pos, op_len, cpos[0], 'gap_OP_C1')
            for i in range(1, len(cpos)):
                # label with the actual cycle numbers: if C2 is missing, C1->C3 is reported as gap_C1_C3
                gap_check(cpos[i-1], clen[i-1], cpos[i], f'gap_C{cidx[i-1]}_C{cidx[i]}')
            if cpos and cp_pos is not None:
                gap_check(cpos[-1], clen[-1], cp_pos, 'gap_C_last_CP')

            for l in clen:
                if l != 9:
                    errs.append('codon_len_not9')
                    break
            if cp_len is not None:
                if not (27 <= cp_len <= 29):
                    errs.append('cp_len_out_of_27_29')

            if hp_seq is not None:
                if hp_seq not in hp_allowed(hp_sets, lib_id):
                    errs.append('HP_not_in_lib')
            if op_seq is not None:
                if op_seq not in op_allowed(op_sets, lib_id):
                    errs.append('OP_not_in_lib')
            if cp_seq is not None:
                if cp_seq not in cp_sets.get(lib_id, set()):
                    errs.append('CP_not_in_lib')
            if lib_id in codon_sets:
                for s in cseq:
                    if s not in codon_sets[lib_id]:
                        errs.append('codon_not_in_lib')
                        break

            if errs:
                counts['fail'] += 1
                for e in errs:
                    counts[e] += 1
            else:
                counts['pass'] += 1

            if not sample_n and progress_every and log and line_no % progress_every == 0:
                log(f"  validated {line_no:,} lines")

        if sample_n and sample_n > 0:
            counts['validated_lines'] = len(sample_lines)
        else:
            counts['validated_lines'] = processed

    return counts


def parse_length_conf(row):
    if not row:
        return None
    try:
        return {
            'policy': row.get('policy', 'reject_outside'),
            'tol': int(row.get('tol', 0)),
            'exp3': int(row.get('exp3', 0)),
            'exp4': int(row.get('exp4', 0)),
            'active_mode': row.get('active_mode', 'both'),
            'active_centers': parse_active_centers(row.get('active_centers', '')),
        }
    except Exception:
        return None


def collect_run(run_name, run_dir, adj_tol, bb, progress_every, log, undecoded_sample=0,
                sample_seed=2026, decoded_max_lines=0, decoded_sample=0,
                max_cp_cands=0, max_anchor_cands=0):
    out = {
        'name': run_name,
        'dir': run_dir,
        'adj_tol': adj_tol,
        'samples': {},
        'totals': {
            'total_reads': 0,
            'length_passed_reads': 0,
            'decoded_reads': 0,
        },
        'decode_fail_totals': Counter(),
        'qc_totals': Counter(),
        'undecoded_reason_totals': Counter(),
        'validation_totals': Counter(),
        'undecoded_sample_totals': Counter(),
        'length_conf_missing_samples': 0,
    }

    if not os.path.isdir(run_dir):
        raise SystemExit(f"[ERROR] run '{run_name}': directory not found: {run_dir}")
    sample_dirs = sorted([d for d in os.listdir(run_dir) if d.startswith('02_decoded_')])
    decoded_root = os.path.join(run_dir, '02_decoded')
    use_legacy = bool(sample_dirs)
    if not use_legacy and not os.path.isdir(decoded_root):
        # Neither the shared 02_decoded/ nor legacy 02_decoded_<sample>/ layout exists: the run would
        # otherwise be reported with zero samples and no message.
        log(f"[{run_name}] WARN: neither {decoded_root} nor 02_decoded_<sample>/ found under {run_dir}; no samples to validate")
    if not use_legacy and os.path.isdir(decoded_root):
        sample_stats_map = read_table_by_sample(os.path.join(decoded_root, 'sample_stats.tsv'))
        decoding_summary_map = read_table_by_sample(os.path.join(decoded_root, 'decoding_summary.tsv'))
        qc_map = read_table_by_sample(os.path.join(decoded_root, 'qc_checks.tsv'))
        length_stats_map = read_table_by_sample(os.path.join(decoded_root, 'length_window_stats.tsv'))
        # Union of decoded files and stats tables: a sample present in sample_stats but lacking a
        # decoded_reads file must not silently vanish from the summary.
        samples = sorted(set(list_samples_from_decoded(decoded_root))
                         | set(sample_stats_map) | set(decoding_summary_map) | set(qc_map))
    else:
        samples = [d.replace('02_decoded_', '') for d in sample_dirs]

    for sample in samples:
        sdir = os.path.join(run_dir, f'02_decoded_{sample}') if use_legacy else decoded_root
        s = {}

        log(f"[{run_name}] sample={sample} start")

        ss = None
        if use_legacy:
            ss_path = os.path.join(sdir, 'sample_stats.tsv')
            if os.path.exists(ss_path):
                ss = read_single_row_tsv(ss_path)
        else:
            ss = sample_stats_map.get(sample)
        if ss:
            s['sample_stats'] = ss
            out['totals']['total_reads'] += int(ss['total_reads'])
            out['totals']['length_passed_reads'] += int(ss['length_passed_reads'])
            out['totals']['decoded_reads'] += int(ss['decoded_reads'])

        ds = None
        if use_legacy:
            ds_path = os.path.join(sdir, 'decoding_summary.tsv')
            if os.path.exists(ds_path):
                ds = read_single_row_tsv(ds_path)
        else:
            ds = decoding_summary_map.get(sample)
        if ds:
            s['decoding_summary'] = ds
            for k in ('no_cp', 'no_op', 'no_hp', 'codon_fail', 'len_out_of_range', 'order_violation'):
                out['decode_fail_totals'][k] += int(ds.get(k, 0))

        qc = None
        if use_legacy:
            qc_path = os.path.join(sdir, 'qc_checks.tsv')
            if os.path.exists(qc_path):
                qc = read_single_row_tsv(qc_path)
        else:
            qc = qc_map.get(sample)
        if qc:
            s['qc_checks'] = qc
            for k in ('codon_len_not9', 'order_violation', 'cp_len_out_of_27_29'):
                out['qc_totals'][k] += int(qc.get(k, 0))

        und_path = _find_sample_file(sdir, 'undecoded_reads_', sample)
        if und_path is None and use_legacy:
            # Legacy per-sample directories may use a differently named file; in the shared
            # 02_decoded layout this fallback would pick ANOTHER sample's file, so legacy-only.
            for fn in os.listdir(sdir):
                if _strip_tsv_suffix(fn, 'undecoded_reads_') is not None:
                    und_path = os.path.join(sdir, fn)
                    break
        if und_path is None:
            s['undecoded_missing'] = True
        if und_path:
            rc = read_reason_counts(und_path)
            s['undecoded_reason_counts'] = dict(rc)
            out['undecoded_reason_totals'].update(rc)

        len_conf = None
        lc = None
        if use_legacy:
            len_path = os.path.join(sdir, 'length_window_stats.tsv')
            if os.path.exists(len_path):
                lc = read_single_row_tsv(len_path)
        else:
            lc = length_stats_map.get(sample)
        if lc:
            len_conf = parse_length_conf(lc)
            s['length_window_stats'] = lc
            s['length_conf'] = len_conf
        else:
            s['length_conf_missing'] = True
            out['length_conf_missing_samples'] += 1

        dec_path = _find_sample_file(sdir, 'decoded_reads_', sample)
        if dec_path is None and use_legacy:
            for fn in os.listdir(sdir):
                if _strip_tsv_suffix(fn, 'decoded_reads_') is not None:
                    dec_path = os.path.join(sdir, fn)
                    break
        if dec_path is None:
            s['decoded_missing'] = True
            log(f"[{run_name}] sample={sample} WARN: decoded_reads file not found; validation skipped")
        if dec_path:
            vc = validate_decoded_file(
                dec_path, adj_tol,
                bb['hp_sets'], bb['op_sets'], bb['cp_sets'], bb['codon_sets'],
                bb['lib_expected_cycles'], len_conf,
                progress_every, log, max_lines=decoded_max_lines,
                sample_n=decoded_sample, sample_seed=sample_seed
            )
            s['validation'] = dict(vc)
            out['validation_totals'].update(vc)

        if und_path and undecoded_sample > 0:
            # File name carries the sample: in the shared 02_decoded layout a sample-less name
            # was overwritten by every subsequent sample.
            out_path = os.path.join(sdir, f"undecoded_sample_recheck_{sample}_{undecoded_sample}.tsv")
            uv = validate_undecoded_sample(und_path, undecoded_sample, sample_seed, adj_tol, len_conf, bb, out_path,
                                           max_cp_cands=max_cp_cands, max_anchor_cands=max_anchor_cands)
            s['undecoded_sample_validation'] = uv
            s['undecoded_sample_recheck'] = out_path
            out['undecoded_sample_totals'].update(uv)

        out['samples'][sample] = s
        log(f"[{run_name}] sample={sample} done")

    return out


def write_summary_tsv(out_path, runs):
    cols = [
        'run','sample','total_reads','length_passed_reads','decoded_reads','decode_rate_pct',
        'no_cp','no_op','no_hp','codon_fail','len_out_of_range','order_violation',
        'qc_codon_len_not9','qc_order_violation','qc_cp_len_out_of_27_29',
        'validation_sampled','validation_pass','validation_fail','validation_cycles_mismatch','validation_length_filter_reject','validation_lib_unknown',
        'length_conf_missing',
        'undecoded_sampled','undecoded_unexpected_decoded','undecoded_len_reject','undecoded_reason_len_mismatch'
    ]
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(cols) + '\n')
        for r in runs:
            for sample, s in sorted(r['samples'].items()):
                ss = s.get('sample_stats', {})
                ds = s.get('decoding_summary', {})
                qc = s.get('qc_checks', {})
                vc = s.get('validation', {})
                uv = s.get('undecoded_sample_validation', {})
                row = {
                    'run': r['name'],
                    'sample': sample,
                    'total_reads': ss.get('total_reads', '0'),
                    'length_passed_reads': ss.get('length_passed_reads', '0'),
                    'decoded_reads': ss.get('decoded_reads', '0'),
                    'decode_rate_pct': ss.get('decode_rate_pct', '0'),
                    'no_cp': ds.get('no_cp', '0'),
                    'no_op': ds.get('no_op', '0'),
                    'no_hp': ds.get('no_hp', '0'),
                    'codon_fail': ds.get('codon_fail', '0'),
                    'len_out_of_range': ds.get('len_out_of_range', '0'),
                    'order_violation': ds.get('order_violation', '0'),
                    'qc_codon_len_not9': qc.get('codon_len_not9', '0'),
                    'qc_order_violation': qc.get('order_violation', '0'),
                    'qc_cp_len_out_of_27_29': qc.get('cp_len_out_of_27_29', '0'),
                    'validation_sampled': str(vc.get('validated_lines', 0)),
                    'validation_pass': str(vc.get('pass', 0)),
                    'validation_fail': str(vc.get('fail', 0)),
                    'validation_cycles_mismatch': str(vc.get('cycles_mismatch', 0)),
                    'validation_length_filter_reject': str(vc.get('length_filter_reject', 0)),
                    'validation_lib_unknown': str(vc.get('lib_unknown', 0)),
                    'length_conf_missing': '1' if s.get('length_conf_missing') else '0',
                    'undecoded_sampled': str(uv.get('sampled_n', 0)),
                    'undecoded_unexpected_decoded': str(uv.get('unexpected_decoded', 0)),
                    'undecoded_len_reject': str(uv.get('len_reject', 0)),
                    'undecoded_reason_len_mismatch': str(uv.get('reason_len_mismatch', 0)),
                }
                f.write('\t'.join(row[c] for c in cols) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bb', required=True)
    ap.add_argument('--run', action='append', nargs=3, metavar=('NAME','DIR','ADJ_TOL'), required=True)
    ap.add_argument('--out-json', required=True)
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--progress-every', type=int, default=1000000)
    ap.add_argument('--undecoded-sample', type=int, default=0)
    ap.add_argument('--decoded-max-lines', type=int, default=0,
                    help="Max decoded lines to validate per sample (0 = all; for quick sampling)")
    ap.add_argument('--decoded-sample', type=int, default=0,
                    help="Randomly sample N decoded reads per sample (0 = disabled; uses reservoir sampling)")
    ap.add_argument('--sample-seed', type=int, default=2026)
    # Mirror the decoder's effective settings so the undecoded re-check reproduces 02_decode_reads.pl.
    ap.add_argument('--mismatch', choices=['none', 'hp_op_cp'], default='hp_op_cp',
                    help="Anchor (HP/OP/CP) 1-bp mismatch indexing mode used by the decode run (default hp_op_cp)")
    ap.add_argument('--max-cp-cands', type=int, default=0,
                    help="Cap on CP candidates per read in the re-check (0 = unlimited, decoder default)")
    ap.add_argument('--max-anchor-cands', type=int, default=0,
                    help="Cap on OP/HP anchor candidates in the re-check (0 = unlimited, decoder default)")
    ap.add_argument('--log', required=True)
    args = ap.parse_args()

    def log(msg):
        ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(args.log, 'a', encoding='utf-8') as lf:
            lf.write(f"[{ts}] {msg}\n")
        print(f"[{ts}] {msg}")

    log(f'Loading BB sets (mismatch={args.mismatch})')
    bb = load_bb_sets(args.bb, mismatch=args.mismatch)

    runs = []
    for name, rdir, adj_tol_s in args.run:
        adj_tol = int(adj_tol_s)
        log(f"Start run: {name} dir={rdir} adj_tol={adj_tol}")
        runs.append(collect_run(
            name, rdir, adj_tol, bb, args.progress_every, log,
            undecoded_sample=args.undecoded_sample, sample_seed=args.sample_seed,
            decoded_max_lines=args.decoded_max_lines, decoded_sample=args.decoded_sample,
            max_cp_cands=args.max_cp_cands, max_anchor_cands=args.max_anchor_cands
        ))
        log(f"End run: {name}")

    # convert Counters to dict for json
    for r in runs:
        for k in ('decode_fail_totals','qc_totals','undecoded_reason_totals','validation_totals','undecoded_sample_totals'):
            r[k] = dict(r[k])

    with open(args.out_json, 'w', encoding='utf-8') as f:
        json.dump(runs, f, indent=2, sort_keys=True)
    write_summary_tsv(args.out_tsv, runs)
    log(f"Wrote {args.out_json} and {args.out_tsv}")

if __name__ == '__main__':
    main()
