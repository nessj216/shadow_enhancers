#!/usr/bin/env python3
"""
Fill / update alignment columns in alignment_with_kmer.csv by scanning enhancer_pairs/
and running EMBOSS needle on each enhancer FASTA pair.

Key features:
- Scans: <root>/enhancer_pairs/<shadowset>/<pair_dir>/enhancers/*.fa
- Updates existing CSV rows (matching by shadowset + undirected pair)
- Adds NEW rows for pairs found on disk but missing from the CSV
- Writes alignment fields; if can't align, keeps blanks and sets alignment_status:
    OK, NO_FASTA, ALIGN_FAIL

Requires:
- EMBOSS needle on PATH

Example:
  python3 fill_alignments_from_enhancer_pairs.py \
    --root /path/to/project_root \
    --csv alignment_with_kmer.csv \
    --out alignment_with_kmer_updated.csv \
    --gapopen 10 --gapextend 0.5 \
    --only-missing
"""

import argparse
import csv
import os
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

FA_EXTS = [".fa", ".fasta", ".fna", ".fas", ".fsa", ".fa.gz", ".fasta.gz", ".fna.gz"]

# -------------------------
# CSV helpers
# -------------------------
def sniff_delimiter(path: Path) -> str:
    s = str(path).lower()
    return "\t" if (s.endswith(".tsv") or s.endswith(".tab")) else ","

def read_table(path: Path) -> Tuple[List[Dict[str, str]], List[str], str]:
    delim = sniff_delimiter(path)
    with open(path, newline="") as f:
        rdr = csv.DictReader(f, delimiter=delim)
        rows = list(rdr)
        cols = rdr.fieldnames or []
    return rows, cols, delim

def write_table(path: Path, rows: List[Dict[str, object]], cols: List[str], delim: str = ",") -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter=delim)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})

def ensure_cols(cols: List[str], wanted: List[str]) -> List[str]:
    seen = set(cols)
    for c in wanted:
        if c not in seen:
            cols.append(c)
            seen.add(c)
    return cols

# -------------------------
# FASTA + alignment helpers
# -------------------------
def normalize_name(s: str) -> str:
    return str(s).strip().replace(":", "_")

def name_variants(s: str) -> List[str]:
    s = str(s).strip()
    return list({s, s.replace(":", "_")})

def find_fasta(enh_dir: Path, name: str) -> Optional[Path]:
    # exact name.ext first, then prefix match
    for stem in name_variants(name):
        for ext in FA_EXTS:
            p = enh_dir / f"{stem}{ext}"
            if p.exists():
                return p
        for ext in FA_EXTS:
            hits = list(enh_dir.glob(f"{stem}*{ext}"))
            if hits:
                return hits[0]
    return None

def run_needle(fa1: Path, fa2: Path, gapopen: str, gapextend: str) -> Tuple[Optional[float], Optional[float], str]:
    cmd = [
        "needle",
        "-asequence", str(fa1),
        "-bsequence", str(fa2),
        "-gapopen", str(gapopen),
        "-gapextend", str(gapextend),
        "-auto",
        "-stdout",
    ]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        return None, None, (p.stderr or p.stdout or "").strip()

    out = p.stdout or ""
    m_id = re.search(r"Identity:\s+\d+/\d+\s+\(([\d\.]+)%\)", out)
    pid = float(m_id.group(1)) if m_id else None
    m_sc = re.search(r"Score:\s+([-\d\.]+)", out)
    score = float(m_sc.group(1)) if m_sc else None
    return pid, score, ""

def is_missing_value(x) -> bool:
    s = str(x).strip()
    return (s == "" or s.lower() == "nan" or s.lower() == "none")

# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Path containing enhancer_pairs/")
    ap.add_argument("--csv", required=True, help="Existing alignment_with_kmer CSV/TSV to update")
    ap.add_argument("--out", default=None, help="Output file (default: overwrite --csv)")
    ap.add_argument("--gapopen", default="10")
    ap.add_argument("--gapextend", default="0.5")
    ap.add_argument("--only-missing", action="store_true",
                    help="Only run needle for rows with missing percent_identity/alignment (recommended)")
    ap.add_argument("--limit", type=int, default=None, help="Process only first N pair dirs (debug)")
    args = ap.parse_args()

    root = Path(args.root)
    ep_root = root / "enhancer_pairs"
    if not ep_root.is_dir():
        raise SystemExit(f"ERROR: not found: {ep_root}")

    csv_path = Path(args.csv)
    out_path = Path(args.out) if args.out else csv_path

    rows, cols, _in_delim = read_table(csv_path)

    # Ensure output columns exist
    needed = [
        "set", "cohort", "pair_dir", "nameA", "nameB",
        "enhancer_pair",
        "percent_identity", "alignment", "alignment_score",
        "aligner", "gapopen", "gapextend",
        "alignment_status",
    ]
    cols = ensure_cols(cols, needed)

    # Build index from existing CSV: (shadowset, undirected_pair) -> row idx
    idx: Dict[Tuple[str, str], int] = {}
    for i, r in enumerate(rows):
        shadowset = (r.get("set") or "").strip()
        a = normalize_name(r.get("nameA") or "")
        b = normalize_name(r.get("nameB") or "")
        if not shadowset or not a or not b:
            continue
        und = "__".join(sorted([a, b]))
        idx[(shadowset, und)] = i

    # Gather all pair dirs on disk
    pair_dirs: List[Tuple[str, Path]] = []
    for shadowset in sorted(os.listdir(ep_root)):
        sdir = ep_root / shadowset
        if not sdir.is_dir():
            continue
        for pair_name in sorted(os.listdir(sdir)):
            pdir = sdir / pair_name
            if pdir.is_dir() and "__" in pair_name:
                pair_dirs.append((shadowset, pdir))

    if args.limit is not None:
        pair_dirs = pair_dirs[: args.limit]

    processed = 0
    updated = 0
    added = 0
    no_fasta = 0
    align_fail = 0
    ok = 0

    for shadowset, pdir in pair_dirs:
        processed += 1
        pair_name = pdir.name
        try:
            a_disk, b_disk = pair_name.split("__", 1)
        except ValueError:
            continue

        und = "__".join(sorted([a_disk, b_disk]))
        key = (shadowset, und)

        row_i = idx.get(key)
        if row_i is None:
            # ADD NEW ROW (your requested behavior)
            r_new = {c: "" for c in cols}
            r_new["set"] = shadowset
            r_new["cohort"] = ""
            r_new["pair_dir"] = str(pdir)
            r_new["nameA"] = a_disk
            r_new["nameB"] = b_disk
            r_new["enhancer_pair"] = f"{a_disk}__{b_disk}"
            r_new["aligner"] = "needle"
            r_new["gapopen"] = str(args.gapopen)
            r_new["gapextend"] = str(args.gapextend)
            r_new["alignment_status"] = "MISSING"
            rows.append(r_new)
            row_i = len(rows) - 1
            idx[key] = row_i
            added += 1

        r = rows[row_i]

        # Respect only-missing: skip if we already have alignment
        if args.only_missing:
            if not is_missing_value(r.get("percent_identity")) or not is_missing_value(r.get("alignment")):
                # if either is present, treat as filled
                continue

        enh_dir = pdir / "enhancers"
        fa1 = find_fasta(enh_dir, a_disk)
        fa2 = find_fasta(enh_dir, b_disk)

        # Always ensure these metadata fields are set/updated
        r["pair_dir"] = str(pdir)
        r["nameA"] = a_disk
        r["nameB"] = b_disk
        r["enhancer_pair"] = f"{a_disk}__{b_disk}"
        r["aligner"] = "needle"
        r["gapopen"] = str(args.gapopen)
        r["gapextend"] = str(args.gapextend)

        if not fa1 or not fa2:
            r["percent_identity"] = ""
            r["alignment"] = ""
            r["alignment_score"] = ""
            r["alignment_status"] = "NO_FASTA"
            no_fasta += 1
            updated += 1
            continue

        pid, score, err = run_needle(fa1, fa2, args.gapopen, args.gapextend)
        if pid is None:
            r["percent_identity"] = ""
            r["alignment"] = ""
            r["alignment_score"] = ""
            r["alignment_status"] = "ALIGN_FAIL"
            align_fail += 1
        else:
            r["percent_identity"] = pid
            r["alignment"] = pid  # alias for convenience
            r["alignment_score"] = score if score is not None else ""
            r["alignment_status"] = "OK"
            ok += 1

        updated += 1

    # Write as CSV (comma) by default for compatibility
    write_table(out_path, rows, cols, delim=",")

    print(f"Scanned pair dirs: {processed}")
    print(f"Rows added: {added}")
    print(f"Rows updated (touched): {updated}")
    print(f"OK: {ok} | NO_FASTA: {no_fasta} | ALIGN_FAIL: {align_fail}")
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
