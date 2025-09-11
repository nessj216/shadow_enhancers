#!/usr/bin/env python3
"""
Make a presence/absence matrix from a directory of PAF files.

Rows:   query names (qname, PAF col 1)
Cols:   one column per .paf/.paf.gz file (by basename)
Cells:  1 if any alignment for that query in that file passes `keep()`, else 0.

Usage:
  python3 paf_matrix.py /path/to/paf_dir presence_matrix.tsv
"""

import sys, os, gzip, csv, glob

def open_maybe_gzip(path):
    if path in ("-", "", None):
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path, "r")

def keep(fields):
    # Require the 12 mandatory PAF fields
    if len(fields) < 12:
        return False
    try:
        qlen   = int(fields[1])   # PAF col 2 (0-based idx 1)
        # PAF col 11 (0-based idx 10) = alignment block length
        aligned_length = int(fields[10])
        # PAF col 12 (0-based idx 11) = mapping quality
        mapq   = int(fields[11])
    except ValueError:
        return False

    if mapq < 10:
        return False

    # Your original logic
    t_span = aligned_length - 2000
    base = qlen - 2000
    lo = base * 0.3
    hi = base * 2

    if lo > hi:
        lo, hi = hi, lo
    lo = max(0, lo)
    if hi < 0:
        return False

    return lo <= t_span <= hi

def scan_paf(path):
    """Return set of query names that have ≥1 passing hit in this PAF."""
    hits = set()
    with open_maybe_gzip(path) as fin:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if keep(f):
                qname = f[0]
                hits.add(qname)
    return hits

def paf_basename(p):
    b = os.path.basename(p)
    if b.endswith(".paf.gz"):
        return b[:-7]  # strip .paf.gz
    if b.endswith(".paf"):
        return b[:-4]  # strip .paf
    return b

def main():
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} <paf_dir> <out_matrix.tsv>", file=sys.stderr)
        sys.exit(2)

    paf_dir = sys.argv[1]
    out_tsv = sys.argv[2]

    if not os.path.isdir(paf_dir):
        print(f"[ERROR] Not a directory: {paf_dir}", file=sys.stderr)
        sys.exit(1)

    # Collect .paf and .paf.gz files (non-recursive; make it recursive if needed)
    paf_files = sorted(
        glob.glob(os.path.join(paf_dir, "*.paf"))
        + glob.glob(os.path.join(paf_dir, "*.paf.gz"))
    )

    if not paf_files:
        print(f"[ERROR] No .paf or .paf.gz files found in {paf_dir}", file=sys.stderr)
        sys.exit(1)

    # For each PAF, compute hits = set of qnames passing keep()
    file_to_hits = {}
    all_queries = set()
    col_names = []

    for p in paf_files:
        sample = paf_basename(p)
        col_names.append(sample)
        hits = scan_paf(p)
        file_to_hits[sample] = hits
        all_queries.update(hits)

    # Build and write matrix
    with open(out_tsv, "w", newline="") as fout:
        w = csv.writer(fout, delimiter="\t")
        w.writerow(["query"] + col_names)
        for q in sorted(all_queries):
            row = [q]
            for sample in col_names:
                row.append(1 if q in file_to_hits.get(sample, set()) else 0)
            w.writerow(row)

    print(f"Wrote matrix for {len(all_queries)} queries across {len(col_names)} files → {out_tsv}")

if __name__ == "__main__":
    main()
































'''#!/usr/bin/env python3
import sys, gzip

def open_maybe_gzip(path):
    if path in ("-", "", None):
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)

def keep(fields):
    # PAF required 12 fields
    if len(fields) < 12:
        return False
    try:
        qlen   = int(fields[1])          # column 2
        #tstart = int(fields[7])
        aligned_length= int(fields[10])# column 8
        #tend   = int(fields[8])          # column 9
        mapq   = int(fields[11])         # column 12
    except ValueError:
        return False

    if mapq < 10:
        return False

    # Target match length as span on target
    t_span = aligned_length-2000

    base = qlen - 2000
    lo = base * 0.3
    hi = base * 2

    # normalize & clamp
    if lo > hi:
        lo, hi = hi, lo
    lo = max(0, lo)
    if hi < 0:
        return False  # nothing can match if both bounds are negative

    return lo <= t_span <= hi

def main():
    if len(sys.argv) < 2:
        print(f"usage: {sys.argv[0]} <in.paf[.gz]> [out.paf]", file=sys.stderr)
        sys.exit(2)
    in_path = sys.argv[1]
    out = open(sys.argv[2], "w") if len(sys.argv) > 2 else sys.stdout

    with open_maybe_gzip(in_path) as fin:
        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if keep(f):
                out.write(line)

if __name__ == "__main__":
    main()'''
