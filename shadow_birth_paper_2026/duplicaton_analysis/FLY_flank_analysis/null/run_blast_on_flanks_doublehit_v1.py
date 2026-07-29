#!/usr/bin/env python3
"""
run_blast_on_flanks_doublehit_v1.py

BLAST-based flank comparison runner for null or observed enhancer pairs.

Comparison planning follows the same adjacency-aware rules as
run_lastz_on_flanks_v3.py:
- FAR (gap >= --gap-threshold): run {LL, RR, LR, RL}
- CLOSE (gap < --gap-threshold):
    * ADJACENT neighbors: run {LL, RR, LR}
    * non-adjacent: run {LL, RR, LR, RL}

For each comparison, the script runs blastn query-vs-subject on the two flank
FASTA files and writes a per-comparison hit table.

After the per-comparison pass, the script builds pair-level summaries:
- LL+RR = concordant double hit
- LR+RL = cross double hit
- BOTH if both patterns are present

Among positive double-hit pairs, it also writes a null cutoff table so you can
choose a percentile-based threshold.
"""

import argparse
import math
import re
import shlex
import subprocess as sp
from multiprocessing.dummy import Pool
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_BLASTN_ARGS = "-task blastn -strand both -dust no -soft_masking false"
DEFAULT_CUTOFF_PERCENTILES = "50,75,90,95,99"
BLAST_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore qcovhsp"
)
BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovhsp",
]

PAIR_RE = re.compile(
    r"^(?P<chrom1>[^_]+)_(?P<s1>\d+)-(?P<e1>\d+)__"
    r"(?P<chrom2>[^_]+)_(?P<s2>\d+)-(?P<e2>\d+)$"
)


def parse_args():
    p = argparse.ArgumentParser(
        description="Run BLAST on flank pairs and summarize double-hit null scores."
    )
    p.add_argument("--flanks-root", required=True, help="Root containing set/pair/flank FASTAs")
    p.add_argument("--manifest", required=True, help="Path to flank manifest TSV")
    p.add_argument("--parallel", type=int, default=4)
    p.add_argument("--gap-threshold", type=int, default=5000)
    p.add_argument("--blastn", default="blastn")
    p.add_argument("--blastn-args", default=DEFAULT_BLASTN_ARGS)
    p.add_argument(
        "--hit-evalue-max",
        type=float,
        default=1.0,
        help="Maximum e-value retained as a hit and passed to blastn.",
    )
    p.add_argument(
        "--cutoff-percentiles",
        default=DEFAULT_CUTOFF_PERCENTILES,
        help="Comma-separated percentiles for the positive-double-hit null cutoff table.",
    )
    p.add_argument(
        "--zero-evalue-floor",
        type=float,
        default=1e-300,
        help="Replacement value when BLAST reports evalue 0.0 for -log10 transforms.",
    )
    p.add_argument(
        "--pair-limit",
        type=int,
        default=0,
        help="Optional cap on pair_dir groups, useful for smoke tests.",
    )
    p.add_argument("--resume", action="store_true")
    return p.parse_args()


def parse_pair_name(pair_name: str):
    m = PAIR_RE.match(pair_name)
    if not m:
        return None
    c1, s1, e1 = m.group("chrom1"), int(m.group("s1")), int(m.group("e1"))
    c2, s2, e2 = m.group("chrom2"), int(m.group("s2")), int(m.group("e2"))
    if c1 == c2 and s1 > s2:
        (c1, s1, e1), (c2, s2, e2) = (c2, s2, e2), (c1, s1, e1)
    return (c1, s1, e1), (c2, s2, e2)


def enh_id_tuple_to_str(tup):
    return f"{tup[0]}:{tup[1]}-{tup[2]}"


def gap_bp(a, b):
    if a[0] != b[0]:
        return 10**12
    return max(0, b[1] - a[2])


def parse_percentiles(text: str):
    values = []
    for raw in str(text).split(","):
        raw = raw.strip()
        if not raw:
            continue
        pct = float(raw)
        if pct < 0 or pct > 100:
            raise SystemExit(f"[ERROR] Percentile must be between 0 and 100: {raw}")
        values.append(pct)
    if not values:
        raise SystemExit("[ERROR] No valid --cutoff-percentiles provided.")
    return values


def read_blast_stats(out_path: Path, hit_evalue_max: float):
    empty = {
        "hits_count": 0,
        "best_evalue": np.nan,
        "best_bitscore": np.nan,
        "best_pident": np.nan,
        "best_length": 0,
        "best_qcovhsp": np.nan,
    }
    if not out_path.exists() or out_path.stat().st_size == 0:
        return empty

    df = pd.read_csv(out_path, sep="\t", names=BLAST_COLUMNS)
    if df.empty:
        return empty

    df = df[pd.to_numeric(df["evalue"], errors="coerce") <= float(hit_evalue_max)].copy()
    if df.empty:
        return empty

    numeric_cols = [
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qcovhsp",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.sort_values(
        by=["evalue", "bitscore", "length", "pident"],
        ascending=[True, False, False, False],
        kind="mergesort",
    )
    best = df.iloc[0]
    return {
        "hits_count": int(len(df)),
        "best_evalue": float(best["evalue"]),
        "best_bitscore": float(best["bitscore"]),
        "best_pident": float(best["pident"]),
        "best_length": int(best["length"]),
        "best_qcovhsp": float(best["qcovhsp"]),
    }


def run_blast(
    subject_fa: Path,
    query_fa: Path,
    out_path: Path,
    log_path: Path,
    blastn_cmd: str,
    blastn_args: str,
    hit_evalue_max: float,
):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        blastn_cmd,
        "-subject",
        str(subject_fa),
        "-query",
        str(query_fa),
        "-evalue",
        str(hit_evalue_max),
        "-outfmt",
        BLAST_OUTFMT,
    ] + shlex.split(blastn_args)

    with log_path.open("w") as logfh, out_path.open("w") as outfh:
        try:
            proc = sp.run(cmd, stdout=outfh, stderr=logfh, check=False)
            return proc.returncode
        except FileNotFoundError:
            logfh.write(f"[ERROR] blast executable not found: {blastn_cmd}\n")
            return 127


def build_set_enhancers_from_pairs(df: pd.DataFrame):
    set_enh = {}
    for pair_dir in df["pair_dir"].drop_duplicates():
        parts = str(pair_dir).split("/", 1)
        set_label = parts[0]
        pair_name = parts[1] if len(parts) > 1 else ""
        parsed = parse_pair_name(pair_name)
        if not parsed:
            continue
        a, b = parsed
        d = set_enh.setdefault(set_label, {})
        d.setdefault(a[0], set()).add((a[1], a[2]))
        d.setdefault(b[0], set()).add((b[1], b[2]))
    for set_label, chrom_map in set_enh.items():
        for chrom, coords in list(chrom_map.items()):
            chrom_map[chrom] = sorted(coords)
    return set_enh


def are_adjacent_neighbors(set_enhancers, set_label, a, b):
    if a[0] != b[0]:
        return False
    ordered = set_enhancers.get(set_label, {}).get(a[0], [])
    try:
        ia = ordered.index((a[1], a[2]))
        ib = ordered.index((b[1], b[2]))
    except ValueError:
        return False
    return abs(ia - ib) == 1


def plan_comparisons(up, down, dpaths, gap_thr, is_adjacent):
    have = lambda enh, side: dpaths.get((enh_id_tuple_to_str(enh), side))
    comps = []
    ll = (have(up, "LEFT"), have(down, "LEFT"))
    rr = (have(up, "RIGHT"), have(down, "RIGHT"))
    lr = (have(up, "LEFT"), have(down, "RIGHT"))
    rl = (have(up, "RIGHT"), have(down, "LEFT"))
    is_close = gap_bp(up, down) < gap_thr

    if not is_close:
        allow = [("LL", ll), ("RR", rr), ("LR", lr), ("RL", rl)]
    else:
        allow = [("LL", ll), ("RR", rr), ("LR", lr)]
        if not is_adjacent:
            allow.append(("RL", rl))

    for tag, pair in allow:
        if pair[0] and pair[1]:
            comps.append((tag, pair[0], pair[1]))
    return comps, is_close


def worker(task):
    (
        pair_dir,
        comparison,
        subject_fa,
        query_fa,
        out_path,
        log_path,
        blastn_cmd,
        blastn_args,
        hit_evalue_max,
        resume,
        is_close_gap,
        is_adjacent_neighbors,
    ) = task

    out_path = Path(out_path)
    log_path = Path(log_path)
    if resume and out_path.exists():
        stats = read_blast_stats(out_path, hit_evalue_max)
        status = "SKIP" if out_path.stat().st_size > 0 else "SKIP_EMPTY"
    else:
        rc = run_blast(
            Path(subject_fa),
            Path(query_fa),
            out_path,
            log_path,
            blastn_cmd,
            blastn_args,
            hit_evalue_max,
        )
        stats = read_blast_stats(out_path, hit_evalue_max)
        status = "OK" if rc == 0 else f"WARN:rc={rc}"

    return {
        "pair_dir": pair_dir,
        "comparison": comparison,
        "status": status,
        "blast_path": str(out_path),
        "log_path": str(log_path),
        "is_close_gap": bool(is_close_gap),
        "is_adjacent_neighbors": bool(is_adjacent_neighbors),
        **stats,
    }


def safe_neg_log10(value: float, floor: float):
    if pd.isna(value):
        return np.nan
    return -math.log10(max(float(value), float(floor)))


def summarize_pairs(comp_df: pd.DataFrame, zero_evalue_floor: float):
    if comp_df.empty:
        return pd.DataFrame()

    rows = []
    for pair_dir, group in comp_df.groupby("pair_dir", sort=True):
        by_comp = {row["comparison"]: row for _, row in group.iterrows()}
        ll_hit = int(by_comp.get("LL", {}).get("hits_count", 0) > 0) if "LL" in by_comp else 0
        rr_hit = int(by_comp.get("RR", {}).get("hits_count", 0) > 0) if "RR" in by_comp else 0
        lr_hit = int(by_comp.get("LR", {}).get("hits_count", 0) > 0) if "LR" in by_comp else 0
        rl_hit = int(by_comp.get("RL", {}).get("hits_count", 0) > 0) if "RL" in by_comp else 0

        ll_rr_double = bool(ll_hit and rr_hit)
        lr_rl_double = bool(lr_hit and rl_hit)
        any_double = bool(ll_rr_double or lr_rl_double)

        if ll_rr_double and lr_rl_double:
            double_hit_type = "BOTH"
            selected = ["LL", "RR", "LR", "RL"]
        elif ll_rr_double:
            double_hit_type = "LL_RR"
            selected = ["LL", "RR"]
        elif lr_rl_double:
            double_hit_type = "LR_RL"
            selected = ["LR", "RL"]
        else:
            double_hit_type = "NONE"
            selected = []

        positive_rows = group[group["comparison"].isin(selected)].copy()
        if positive_rows.empty:
            best_double_evalue = np.nan
            best_double_bitscore = np.nan
            best_double_comparison = ""
        else:
            positive_rows = positive_rows.sort_values(
                by=["best_evalue", "best_bitscore", "best_length"],
                ascending=[True, False, False],
                na_position="last",
                kind="mergesort",
            )
            best = positive_rows.iloc[0]
            best_double_evalue = float(best["best_evalue"])
            best_double_bitscore = float(
                positive_rows["best_bitscore"].max(skipna=True)
            )
            best_double_comparison = str(best["comparison"])

        rows.append(
            {
                "pair_dir": pair_dir,
                "is_close_gap": bool(group["is_close_gap"].iloc[0]),
                "is_adjacent_neighbors": bool(group["is_adjacent_neighbors"].iloc[0]),
                "LL_hit": ll_hit,
                "RR_hit": rr_hit,
                "LR_hit": lr_hit,
                "RL_hit": rl_hit,
                "LL_RR_double_hit": int(ll_rr_double),
                "LR_RL_double_hit": int(lr_rl_double),
                "any_double_hit": int(any_double),
                "double_hit_type": double_hit_type,
                "best_double_hit_comparison": best_double_comparison,
                "best_double_hit_evalue": best_double_evalue,
                "best_double_hit_bitscore": best_double_bitscore,
                "best_double_hit_neg_log10_evalue": safe_neg_log10(
                    best_double_evalue, zero_evalue_floor
                ),
            }
        )
    return pd.DataFrame(rows)


def build_cutoff_table(pair_df: pd.DataFrame, percentiles, zero_evalue_floor: float):
    positive = pair_df[pair_df["any_double_hit"] == 1].copy()
    if positive.empty:
        return pd.DataFrame(
            columns=[
                "metric",
                "percentile",
                "cutoff_value",
                "stronger_values",
                "n_positive_double_hits",
                "notes",
            ]
        )

    rows = []
    n_positive = int(len(positive))

    neglog_values = positive["best_double_hit_neg_log10_evalue"].dropna().to_numpy()
    bitscore_values = positive["best_double_hit_bitscore"].dropna().to_numpy()
    evalue_values = positive["best_double_hit_evalue"].dropna().to_numpy()

    for pct in percentiles:
        q = pct / 100.0

        if len(neglog_values) > 0:
            cutoff = float(np.quantile(neglog_values, q))
            rows.append(
                {
                    "metric": "neg_log10_evalue",
                    "percentile": pct,
                    "cutoff_value": cutoff,
                    "stronger_values": "higher",
                    "n_positive_double_hits": n_positive,
                    "notes": f"Equivalent evalue ~= {10 ** (-cutoff):.6g}",
                }
            )

        if len(bitscore_values) > 0:
            cutoff = float(np.quantile(bitscore_values, q))
            rows.append(
                {
                    "metric": "bitscore",
                    "percentile": pct,
                    "cutoff_value": cutoff,
                    "stronger_values": "higher",
                    "n_positive_double_hits": n_positive,
                    "notes": "",
                }
            )

        if len(evalue_values) > 0:
            cutoff = float(np.quantile(evalue_values, 1.0 - q))
            rows.append(
                {
                    "metric": "evalue",
                    "percentile": pct,
                    "cutoff_value": cutoff,
                    "stronger_values": "lower",
                    "n_positive_double_hits": n_positive,
                    "notes": (
                        "For evalue, this is the lower-tail cutoff matching the same "
                        "strength percentile as the score-based rows."
                    ),
                }
            )

    return pd.DataFrame(rows)


def main():
    args = parse_args()
    cutoff_percentiles = parse_percentiles(args.cutoff_percentiles)

    root = Path(args.flanks_root).resolve()
    manifest = Path(args.manifest).resolve()
    df = pd.read_csv(manifest, sep="\t")

    required = {"gene_label", "pair_dir", "fasta_path"}
    if not required.issubset(df.columns):
        missing = sorted(required - set(df.columns))
        raise SystemExit(f"[ERROR] manifest missing columns: {missing}")

    set_enhancers = build_set_enhancers_from_pairs(df)

    tasks = []
    pair_groups = list(df.groupby("pair_dir", sort=True))
    if args.pair_limit and args.pair_limit > 0:
        pair_groups = pair_groups[: int(args.pair_limit)]

    for pair_dir, group in pair_groups:
        parts = str(pair_dir).split("/", 1)
        set_label = parts[0]
        pair_name = parts[1] if len(parts) > 1 else ""
        parsed = parse_pair_name(pair_name)
        if not parsed:
            print(f"[WARN] Cannot parse pair_dir: {pair_dir}")
            continue

        up, down = parsed
        dpaths = {}
        for _, row in group.iterrows():
            enh_id = str(row.get("enhancer_id") or row.get("flank_enhancer_id") or "")
            side = str(row.get("flank_side") or row.get("side") or "").upper()
            fasta_path = root / str(row["fasta_path"])
            if enh_id and side in {"LEFT", "RIGHT"} and fasta_path.exists():
                dpaths[(enh_id, side)] = fasta_path

        if not dpaths:
            print(f"[INFO] No flanks found for {pair_dir}")
            continue

        is_adjacent = (
            are_adjacent_neighbors(set_enhancers, set_label, up, down)
            if gap_bp(up, down) < args.gap_threshold
            else False
        )
        comps, is_close = plan_comparisons(
            up, down, dpaths, args.gap_threshold, is_adjacent
        )
        if not comps:
            print(f"[INFO] No valid comparisons for {pair_dir}")
            continue

        for comparison, subject_fa, query_fa in comps:
            blast_dir = root / pair_dir / "blast" / comparison
            out_path = blast_dir / "blast_hits.tsv"
            log_path = blast_dir / "blast.log"
            tasks.append(
                (
                    pair_dir,
                    comparison,
                    str(subject_fa),
                    str(query_fa),
                    str(out_path),
                    str(log_path),
                    args.blastn,
                    args.blastn_args,
                    args.hit_evalue_max,
                    args.resume,
                    is_close,
                    is_adjacent,
                )
            )

    if not tasks:
        print("[INFO] Planned comparisons: 0")
        return

    if int(args.parallel) <= 1:
        results = [worker(task) for task in tasks]
    else:
        with Pool(processes=int(args.parallel)) as pool:
            results = list(pool.map(worker, tasks))

    comp_df = pd.DataFrame(results).sort_values(
        by=["pair_dir", "comparison"], kind="mergesort"
    )
    pair_df = summarize_pairs(comp_df, args.zero_evalue_floor)
    cutoff_df = build_cutoff_table(
        pair_df, cutoff_percentiles, args.zero_evalue_floor
    )

    comp_path = root / "blast_comparisons.tsv"
    pair_path = root / "blast_pairs_summary.tsv"
    cutoff_path = root / "blast_null_cutoffs.tsv"

    comp_df.to_csv(comp_path, sep="\t", index=False)
    pair_df.to_csv(pair_path, sep="\t", index=False)
    cutoff_df.to_csv(cutoff_path, sep="\t", index=False)

    n_positive = int((pair_df["any_double_hit"] == 1).sum()) if not pair_df.empty else 0
    print(f"[DONE] Wrote per-comparison BLAST summary: {comp_path}")
    print(f"[DONE] Wrote pair-level double-hit summary: {pair_path}")
    print(f"[DONE] Wrote null cutoff table: {cutoff_path}")
    print(f"[DONE] Positive double-hit pairs: {n_positive}")


if __name__ == "__main__":
    main()
