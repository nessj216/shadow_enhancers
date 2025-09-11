#!/usr/bin/env python3
import os, re, glob, argparse, csv
import pandas as pd
from typing import Dict, List, Tuple, Optional

def parse_paf_line(line: str):
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 12:
        return None
    d = {
        "qname": parts[0],
        "qlen": int(parts[1]),
        "qstart": int(parts[2]),
        "qend": int(parts[3]),
        "strand": parts[4],
        "tname": parts[5],
        "tlen": int(parts[6]),
        "tstart": int(parts[7]),
        "tend": int(parts[8]),
        "nmatch": int(parts[9]),
        "block_len": int(parts[10]),
        "mapq": int(parts[11]),
        "tags_raw": parts[12:],
    }
    tags = {}
    for tag in d["tags_raw"]:
        if ":" in tag:
            k,t,v = tag.split(":",2)
            tags[k]=(t,v)
    d["tags"] = tags
    return d

def is_primary(tags: dict) -> bool:
    return tags.get("tp", (None,None))[1] in (None, "P")

def parse_cg(cg: str):
    # returns list of (op, length) with numeric length as int
    return [(op, int(n)) for n, op in re.findall(r"(\d+)([MIDNSHP=XB])", cg)]

def flank_coverages_from_cg_and_qstart(cg: str, qstart: int, flank_len: int = 500) -> Tuple[int,int]:
    """
    Compute aligned coverage (M/=X/I) for left [0,flank_len) and right [flank_len,2*flank_len) on the FULL query.
    cg only covers the aligned portion; qstart shifts the aligned region within the full query length.
    """
    ops = parse_cg(cg)
    qpos = qstart  # position in full query coordinates
    left = 0
    right = 0
    for op, length in ops:
        if op in ("M","=","X","I"):  # consume query
            q_start = qpos
            q_end = qpos + length
            if op in ("M","=","X","I"):
                # overlap with left and right windows
                ol_left = max(0, min(q_end, flank_len) - max(q_start, 0))
                left += ol_left
                ol_right = max(0, min(q_end, 2*flank_len) - max(q_start, flank_len))
                right += ol_right
            qpos = q_end
        else:
            # ref-only ops (D, N, H, P, B etc.) do not consume query
            pass
    return left, right

def n_total_from_cg(cg: str) -> int:
    return sum(int(n) for n,op in re.findall(r"(\d+)([A-Z=])", cg) if op == "N")

def choose_best_alignments(paf_path: str) -> Dict[str, dict]:
    """
    Keep the best/primary alignment per query name.
    Preference: primary (tp:P) > higher MAPQ > longer block_len.
    """
    best = {}
    with open(paf_path, "r") as fh:
        for line in fh:
            if not line.strip(): continue
            r = parse_paf_line(line)
            if not r: continue
            qn = r["qname"]
            prev = best.get(qn)
            if prev is None:
                best[qn] = r
            else:
                p_pri = is_primary(prev["tags"])
                r_pri = is_primary(r["tags"])
                if r_pri and not p_pri:
                    best[qn] = r
                elif r_pri == p_pri:
                    if r["mapq"] > prev["mapq"] or (r["mapq"] == prev["mapq"] and r["block_len"] > prev["block_len"]):
                        best[qn] = r
    return best

def call_te_status(record: dict, flank_len: int = 500, min_mapq: int = 20, min_flank_cov: float = 0.85):
    """
    Returns one of {"1","0","N/A"} based on user rules:
      - Require MAPQ >= min_mapq and at least one flank coverage >= min_flank_cov
      - Then TE present if cg contains N (=> "1"), else "0"
      - If cg missing or filters fail => "N/A"
    """
    if record["mapq"] < min_mapq:
        return "N/A"
    cg = record["tags"].get("cg", (None,None))[1] if "cg" in record["tags"] else None
    if not cg:
        return "N/A"
    left_aln, right_aln = flank_coverages_from_cg_and_qstart(cg, record["qstart"], flank_len=flank_len)
    left_cov = left_aln / float(flank_len)
    right_cov = right_aln / float(flank_len)
    if not (left_cov >= min_flank_cov or right_cov >= min_flank_cov):
        return "N/A"
    hasN = n_total_from_cg(cg) > 0
    return "1" if hasN else "0"

def build_matrix(paf_dir: str, out_tsv: str, flank_len: int = 500, min_mapq: int = 20, min_flank_cov: float = 0.85, glob_pat: str = "*.paf"):
    paf_files = sorted(glob.glob(os.path.join(paf_dir, glob_pat)))
    if not paf_files:
        raise SystemExit(f"No PAF files found in {paf_dir} (pattern {glob_pat})")
    # Derive genome names from file basenames (strip common suffixes)
    def genome_name(path: str) -> str:
        base = os.path.basename(path)
        for suf in [".paf", ".txt", ".tsv"]:
            if base.endswith(suf):
                base = base[: -len(suf)]
                break
        return base

    genome_cols = [genome_name(p) for p in paf_files]

    # Collect best alignment per TE per genome
    all_te_ids = set()
    genome_to_calls: Dict[str, Dict[str, str]] = {}
    for path, gcol in zip(paf_files, genome_cols):
        best = choose_best_alignments(path)
        calls = {}
        for te_id, rec in best.items():
            status = call_te_status(rec, flank_len=flank_len, min_mapq=min_mapq, min_flank_cov=min_flank_cov)
            calls[te_id] = status
        genome_to_calls[gcol] = calls
        all_te_ids.update(best.keys())

    # Build matrix
    te_ids_sorted = sorted(all_te_ids, key=lambda x: (x.split("|")[0], x))
    data = []
    header = ["TE_ID"] + genome_cols
    for te in te_ids_sorted:
        row = [te]
        for g in genome_cols:
            val = genome_to_calls[g].get(te, "N/A")
            row.append(val)
        data.append(row)

    # Write TSV
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        w.writerows(data)

    # Also return a DataFrame for convenience
    df = pd.DataFrame(data, columns=header)
    return df

def main():
    import argparse
    ap = argparse.ArgumentParser(description="Build TE presence/absence matrix from PAFs (simple rules).")
    ap.add_argument("paf_dir", help="Directory containing PAF files")
    ap.add_argument("--out", default="te_matrix.tsv", help="Output TSV path (default: te_matrix.tsv)")
    ap.add_argument("--flank-len", type=int, default=500, help="Length of each flank in query (default 500)")
    ap.add_argument("--min-mapq", type=int, default=20, help="Minimum MAPQ (default 20)")
    ap.add_argument("--min-flank-cov", type=float, default=0.85, help="Minimum coverage of at least one flank (default 0.85)")
    ap.add_argument("--glob", default="*.paf", help="Glob pattern for PAF files (default *.paf)")
    args = ap.parse_args()

    df = build_matrix(args.paf_dir, args.out, flank_len=args.flank_len, min_mapq=args.min_mapq, min_flank_cov=args.min_flank_cov, glob_pat=args.glob)
    print(f"Wrote matrix to {args.out}")
    # quick QC summary
    # Count per column values
    summary_rows = []
    for col in df.columns[1:]:
        vc = df[col].value_counts()
        present = int(vc.get("1", 0))
        absent  = int(vc.get("0", 0))
        na      = int(vc.get("N/A", 0))
        summary_rows.append([col, present, absent, na])
    summ = pd.DataFrame(summary_rows, columns=["genome","present_1","absent_0","NA"])
    summ.to_csv(os.path.splitext(args.out)[0] + ".summary.tsv", sep="\t", index=False)
    print(summ.to_string(index=False))


