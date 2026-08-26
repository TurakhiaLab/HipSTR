#!/usr/bin/env python3
"""
Tolerant VCF comparator for HipSTR-MT correctness testing under precision changes
(e.g. double -> float). Unlike a byte/md5 diff, this treats small numeric
drift in float-typed fields as acceptable while still hard-failing on any
change to genotype calls, allele/ref content, or non-numeric text.

It doesn't hardcode HipSTR-MT's schema. Every field value (INFO subfields and
per-sample FORMAT subfields) is tokenized into numeric and literal segments.
Literal segments must match exactly. Numeric segments are compared as floats:
- if both sides look like integers (no '.', no exponent) and disagree, that's
  an "int_mismatch" (a discrete count/index changed -- usually a real
  regression, not rounding).
- otherwise, relative drift is computed and compared against --float-tol.

GT is always forced into the int_mismatch/critical bucket regardless of
--lenient-ints, since a genotype call change is never acceptable.

Usage:
  python3 compare_vcf_tolerant.py original.vcf.gz candidate.vcf.gz
  python3 compare_vcf_tolerant.py orig.vcf.gz cand.vcf.gz --float-tol 1e-3 \
      --lenient-ints --max-report 100 --out report.tsv

Exit status: 0 if no critical findings (genotype mismatches, structural
mismatches, missing/extra records, or float drift beyond --float-tol),
nonzero otherwise. Suitable for CI gating.
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict

NUM_RE = re.compile(r'[+-]?\d+\.\d+(?:[eE][+-]?\d+)?|[+-]?\d+(?:[eE][+-]?\d+)?')


def smart_open(path):
    with open(path, 'rb') as fh:
        magic = fh.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rt')
    return open(path, 'rt')


def tokenize(s):
    tokens = []
    last = 0
    for m in NUM_RE.finditer(s):
        if m.start() > last:
            tokens.append(('lit', s[last:m.start()]))
        tokens.append(('num', m.group()))
        last = m.end()
    if last < len(s):
        tokens.append(('lit', s[last:]))
    return tokens


def is_int_token(tok):
    return '.' not in tok and 'e' not in tok.lower()


def compare_values(a, b, float_tol, force_int=False):
    """Returns (status, max_rel_drift). status in:
    'exact', 'float_drift_ok', 'float_fail', 'int_mismatch', 'structural'
    """
    if a == b:
        return 'exact', 0.0

    ta, tb = tokenize(a), tokenize(b)
    if len(ta) != len(tb) or [t[0] for t in ta] != [t[0] for t in tb]:
        return 'structural', None

    max_drift = 0.0
    any_int_mismatch = False
    for (ka, va), (kb, vb) in zip(ta, tb):
        if ka == 'lit':
            if va != vb:
                return 'structural', None
            continue
        if va == vb:
            continue
        fa, fb = float(va), float(vb)
        if force_int or (is_int_token(va) and is_int_token(vb)):
            any_int_mismatch = True
        else:
            denom = max(abs(fa), abs(fb), 1e-9)
            max_drift = max(max_drift, abs(fa - fb) / denom)

    if any_int_mismatch:
        return 'int_mismatch', max_drift
    if max_drift > float_tol:
        return 'float_fail', max_drift
    if max_drift > 0:
        return 'float_drift_ok', max_drift
    return 'exact', 0.0


def parse_vcf(path):
    """Returns (samples, records) where records is a dict keyed by
    (chrom, pos) -> dict of fields, preserving raw INFO/FORMAT structure."""
    samples = None
    records = {}
    dup_keys = 0
    with smart_open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                cols = line.split('\t')
                samples = cols[9:]
                continue
            cols = line.split('\t')
            chrom, pos, rid, ref, alt, qual, filt, info, fmt = cols[:9]
            sample_vals = cols[9:]
            key = (chrom, pos)
            if key in records:
                dup_keys += 1
            records[key] = {
                'id': rid, 'ref': ref, 'alt': alt, 'qual': qual, 'filter': filt,
                'info': info, 'format': fmt, 'samples': sample_vals,
            }
    if samples is None:
        raise ValueError(f"{path}: no #CHROM header line found")
    return samples, records, dup_keys


def info_to_dict(info_str):
    d = {}
    if info_str == '.':
        return d
    for kv in info_str.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
        else:
            d[kv] = None
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('original', help='Reference/baseline VCF (.vcf or .vcf.gz)')
    ap.add_argument('candidate', help='VCF to check against the baseline')
    ap.add_argument('--float-tol', type=float, default=1e-3,
                     help='Max allowed relative drift for float-like fields (default: 1e-3)')
    ap.add_argument('--lenient-ints', action='store_true',
                     help='Downgrade non-GT integer mismatches (e.g. DAB, AC) from critical to warning')
    ap.add_argument('--max-report', type=int, default=50,
                     help='Max per-record findings printed to stdout (default: 50)')
    ap.add_argument('--out', help='Write full TSV of all findings to this path')
    args = ap.parse_args()

    samples_a, records_a, dup_a = parse_vcf(args.original)
    samples_b, records_b, dup_b = parse_vcf(args.candidate)

    if dup_a or dup_b:
        print(f"WARNING: duplicate (CHROM,POS) keys found -- original={dup_a} candidate={dup_b}; "
              f"only the last record per key was kept per file.", file=sys.stderr)

    if samples_a != samples_b:
        print(f"CRITICAL: sample columns differ.\n  original:  {samples_a}\n  candidate: {samples_b}",
              file=sys.stderr)
        sys.exit(2)

    keys_a, keys_b = set(records_a), set(records_b)
    missing_in_b = keys_a - keys_b
    missing_in_a = keys_b - keys_a
    shared = keys_a & keys_b

    findings = []  # (chrom, pos, field, sample_or_INFO, status, drift, a_val, b_val)

    def record_finding(chrom, pos, field, scope, status, drift, a_val, b_val):
        findings.append((chrom, pos, field, scope, status, drift, a_val, b_val))

    field_drift_stats = defaultdict(lambda: [0.0, 0.0, 0])  # field -> [max, sum, count]

    for key in sorted(shared):
        chrom, pos = key
        ra, rb = records_a[key], records_b[key]

        for field in ('ref', 'alt', 'id', 'filter'):
            if ra[field] != rb[field]:
                record_finding(chrom, pos, field.upper(), '-', 'structural', None, ra[field], rb[field])

        status, drift = compare_values(ra['qual'], rb['qual'], args.float_tol)
        if status not in ('exact',):
            record_finding(chrom, pos, 'QUAL', '-', status, drift, ra['qual'], rb['qual'])

        info_a, info_b = info_to_dict(ra['info']), info_to_dict(rb['info'])
        for k in sorted(set(info_a) | set(info_b)):
            va, vb = info_a.get(k), info_b.get(k)
            if va is None or vb is None:
                record_finding(chrom, pos, f'INFO/{k}', '-', 'structural', None, va, vb)
                continue
            status, drift = compare_values(va, vb, args.float_tol)
            if status == 'float_drift_ok' or status == 'float_fail':
                s = field_drift_stats[f'INFO/{k}']
                s[0] = max(s[0], drift); s[1] += drift; s[2] += 1
            if status != 'exact':
                record_finding(chrom, pos, f'INFO/{k}', '-', status, drift, va, vb)

        fmt_a, fmt_b = ra['format'].split(':'), rb['format'].split(':')
        for si, sample in enumerate(samples_a):
            vals_a = ra['samples'][si].split(':')
            vals_b = rb['samples'][si].split(':')
            da = dict(zip(fmt_a, vals_a))
            db = dict(zip(fmt_b, vals_b))
            for k in sorted(set(da) | set(db)):
                va, vb = da.get(k), db.get(k)
                if va is None or vb is None:
                    record_finding(chrom, pos, f'FORMAT/{k}', sample, 'structural', None, va, vb)
                    continue
                force_int = (k == 'GT')
                status, drift = compare_values(va, vb, args.float_tol, force_int=force_int)
                if status in ('float_drift_ok', 'float_fail'):
                    s = field_drift_stats[f'FORMAT/{k}']
                    s[0] = max(s[0], drift); s[1] += drift; s[2] += 1
                if status != 'exact':
                    record_finding(chrom, pos, f'FORMAT/{k}', sample, status, drift, va, vb)

    def severity(status, field):
        if status == 'structural':
            return 'CRITICAL'
        if status == 'int_mismatch':
            if field == 'FORMAT/GT':
                return 'CRITICAL'
            return 'WARN' if args.lenient_ints else 'CRITICAL'
        if status == 'float_fail':
            return 'FAIL'
        if status == 'float_drift_ok':
            return 'ok'
        return 'ok'

    critical = [f for f in findings if severity(f[4], f[2]) == 'CRITICAL']
    fail = [f for f in findings if severity(f[4], f[2]) == 'FAIL']
    warn = [f for f in findings if severity(f[4], f[2]) == 'WARN']

    print(f"Compared {len(shared)} shared records "
          f"({len(missing_in_b)} missing in candidate, {len(missing_in_a)} extra in candidate)")
    print(f"Findings: {len(critical)} CRITICAL, {len(fail)} FAIL (drift > {args.float_tol}), "
          f"{len(warn)} WARN, {sum(f[2]=='FORMAT/GT' for f in critical)} of which are GT mismatches")

    shown = 0
    for f in critical + fail + warn:
        if shown >= args.max_report:
            print(f"... ({len(critical) + len(fail) + len(warn) - shown} more findings suppressed, "
                  f"raise --max-report or use --out)")
            break
        chrom, pos, field, scope, status, drift, a_val, b_val = f
        sev = severity(status, field)
        drift_str = f" drift={drift:.3g}" if drift else ""
        print(f"[{sev}] {chrom}:{pos} {field} ({scope}) {status}{drift_str}: "
              f"{a_val!r} -> {b_val!r}")
        shown += 1

    if missing_in_b:
        print(f"\n{len(missing_in_b)} records present in original but missing from candidate (first 10):")
        for chrom, pos in sorted(missing_in_b)[:10]:
            print(f"  {chrom}:{pos}")
    if missing_in_a:
        print(f"\n{len(missing_in_a)} records present in candidate but not original (first 10):")
        for chrom, pos in sorted(missing_in_a)[:10]:
            print(f"  {chrom}:{pos}")

    if field_drift_stats:
        print("\nFloat drift by field (max / mean relative drift, over comparisons where values differed):")
        for k in sorted(field_drift_stats):
            mx, total, cnt = field_drift_stats[k]
            print(f"  {k}: max={mx:.3g} mean={total/cnt:.3g} n={cnt}")

    if args.out:
        with open(args.out, 'w') as fh:
            fh.write("chrom\tpos\tfield\tscope\tseverity\tstatus\tdrift\toriginal\tcandidate\n")
            for chrom, pos, field, scope, status, drift, a_val, b_val in findings:
                sev = severity(status, field)
                drift_str = f"{drift:.6g}" if drift else ""
                fh.write(f"{chrom}\t{pos}\t{field}\t{scope}\t{sev}\t{status}\t{drift_str}\t{a_val}\t{b_val}\n")
        print(f"\nFull findings written to {args.out}")

    ok = not critical and not fail and not missing_in_a and not missing_in_b
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
