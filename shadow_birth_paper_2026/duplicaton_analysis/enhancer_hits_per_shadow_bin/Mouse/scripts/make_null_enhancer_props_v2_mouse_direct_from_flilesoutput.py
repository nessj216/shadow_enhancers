#!/usr/bin/env python3
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import math
import multiprocessing as mp
from pathlib import Path
import random
import re
import subprocess
import tempfile
import zipfile

import numpy as np
import pandas as pd


WORKER_STATE = {}


# ---------------------------
# Basic file utilities
# ---------------------------

def is_enhancer_file(p: Path) -> bool:
    n = p.name.lower()
    return (
        p.is_file()
        and n.endswith(
            (
                ".fa",
                ".fna",
                ".fasta",
                ".fas",
                ".fa.gz",
                ".fasta.gz",
                ".fna.gz",
                ".txt",
            )
        )
    )


def enhancer_stem(p: Path) -> str:
    n = re.sub(r"\.gz$", "", p.name, flags=re.I)
    return re.sub(r"\.(fa|fna|fasta|fas|txt)$", "", n, flags=re.I)


def parse_fasta_like_text(txt: str) -> dict:
    records = {}
    seq = []
    hdr = None
    for line in txt.splitlines():
        if line.startswith(">"):
            if hdr is not None:
                records[hdr] = "".join(seq)
            hdr = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())
    if hdr is not None:
        records[hdr] = "".join(seq)
    return records


def read_sequence(p: Path) -> str:
    txt = p.read_text()
    recs = parse_fasta_like_text(txt)
    if not recs:
        return ""
    return next(iter(recs.values()))


def read_first_header(p: Path) -> str:
    txt = p.read_text()
    recs = parse_fasta_like_text(txt)
    if not recs:
        return ""
    return next(iter(recs.keys()))


def enhancer_id_from_header(header: str) -> str:
    token = str(header).split()[0]
    if ":" not in token or "-" not in token:
        raise ValueError(f"Cannot parse enhancer id from FASTA header: {header}")
    return token


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return (a_start < b_end) and (b_start < a_end)


# ---------------------------
# BLAST wrapper
# ---------------------------

def blast_has_any_hit(seq_a: str, seq_b: str, blastn_path: str, params: dict, dry_run: bool = False) -> bool:
    if dry_run:
        h = hash(seq_a[:50] + "|" + seq_b[:50])
        return (h % 100) < 3

    with tempfile.TemporaryDirectory() as td:
        q = Path(td) / "q.fa"
        s = Path(td) / "s.fa"
        q.write_text(">q\n" + seq_a + "\n")
        s.write_text(">s\n" + seq_b + "\n")

        cmd = [
            blastn_path,
            "-query",
            str(q),
            "-subject",
            str(s),
            "-evalue",
            str(params.get("evalue", "0.01")),
            "-word_size",
            str(params.get("word_size", "7")),
            "-gapopen",
            str(params.get("gapopen", "5")),
            "-gapextend",
            str(params.get("gapextend", "2")),
            "-reward",
            str(params.get("reward", "2")),
            "-penalty",
            str(params.get("penalty", "-3")),
            "-dust",
            str(params.get("dust", "yes")),
            "-outfmt",
            "6",
        ]
        out = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            timeout=params.get("timeout", 60),
        )
        return bool(out.stdout.strip())


# ---------------------------
# Build catalog of all enhancers
# ---------------------------

def build_catalog(shadow_zip: Path, working_dir: Path) -> pd.DataFrame:
    with zipfile.ZipFile(shadow_zip) as z:
        z.extractall(working_dir)

    candidates = [p for p in working_dir.iterdir() if p.is_dir() and p.name != "__MACOSX"]
    root = max(candidates, key=lambda p: len([q for q in p.iterdir() if q.is_dir()])) if len(candidates) > 1 else candidates[0]

    rows = []
    for gd in [p for p in root.iterdir() if p.is_dir() and p.name != "__MACOSX"]:
        files = [p for p in gd.iterdir() if is_enhancer_file(p)]
        if not files:
            for sub in gd.iterdir():
                if sub.is_dir():
                    files += [p for p in sub.iterdir() if is_enhancer_file(p)]
        for f in files:
            hdr = read_first_header(f)
            rows.append(
                {
                    "gene": gd.name,
                    "path": str(f),
                    "filename_stem": enhancer_stem(f),
                    "enhancer_id": enhancer_id_from_header(hdr),
                }
            )
    return pd.DataFrame(rows)


def build_catalog_from_root(root: Path) -> pd.DataFrame:
    rows = []
    for gd in sorted([p for p in root.iterdir() if p.is_dir() and p.name != "__MACOSX"]):
        files = [p for p in gd.iterdir() if is_enhancer_file(p)]
        if not files:
            for sub in gd.iterdir():
                if sub.is_dir():
                    files += [p for p in sub.iterdir() if is_enhancer_file(p)]
        for f in sorted(files):
            hdr = read_first_header(f)
            rows.append(
                {
                    "gene": gd.name,
                    "path": str(f),
                    "filename_stem": enhancer_stem(f),
                    "enhancer_id": enhancer_id_from_header(hdr),
                }
            )
    if not rows:
        raise SystemExit(f"[ERROR] No enhancer FASTA/TXT files found under {root}")
    return pd.DataFrame(rows)


def derive_gene_sizes_from_catalog(catalog: pd.DataFrame) -> pd.DataFrame:
    return (
        catalog.groupby("gene", as_index=False)
        .agg(num_shadow_enhancers=("filename_stem", "nunique"))
        .rename(columns={"gene": "Gene Name"})
    )


def read_flank_manifest(manifest_path: Path, flanks_root: Path, allowed_keys: set | None = None):
    need = {"gene_label", "enhancer_id", "flank_side", "fasta_path"}
    coord_cols = {"chrom", "start", "end"}
    usecols = list(need | coord_cols)
    df = pd.read_csv(manifest_path, sep="\t", usecols=lambda c: c in usecols)
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"[ERROR] flank manifest missing columns: {sorted(missing)}")

    flank_paths = defaultdict(dict)
    flank_coords = defaultdict(dict)
    allowed = set(allowed_keys or [])
    for row in df.itertuples(index=False):
        gene = str(row.gene_label)
        enh_id = str(row.enhancer_id)
        key = (gene, enh_id)
        if allowed and key not in allowed:
            continue

        side = str(row.flank_side).upper()
        if side not in {"LEFT", "RIGHT"}:
            continue

        fpath = flanks_root / str(row.fasta_path)
        if not fpath.exists():
            continue

        flank_paths[key][side] = str(fpath)
        chrom = str(getattr(row, "chrom", ""))
        start = pd.to_numeric(getattr(row, "start", None), errors="coerce")
        end = pd.to_numeric(getattr(row, "end", None), errors="coerce")
        if chrom and pd.notna(start) and pd.notna(end):
            flank_coords[key][side] = (chrom, int(start), int(end))

    return dict(flank_paths), dict(flank_coords)


def any_flank_overlap_for_pair(key_a, key_b, flank_coords) -> bool:
    if key_a not in flank_coords or key_b not in flank_coords:
        return False
    if not {"LEFT", "RIGHT"}.issubset(flank_coords[key_a]) or not {"LEFT", "RIGHT"}.issubset(flank_coords[key_b]):
        return False

    a_left, a_right = flank_coords[key_a]["LEFT"], flank_coords[key_a]["RIGHT"]
    b_left, b_right = flank_coords[key_b]["LEFT"], flank_coords[key_b]["RIGHT"]

    coords_a = [a_left, a_right]
    coords_b = [b_left, b_right]
    for chrom_a, start_a, end_a in coords_a:
        for chrom_b, start_b, end_b in coords_b:
            if chrom_a == chrom_b and intervals_overlap(start_a, end_a, start_b, end_b):
                return True
    return False


# ---------------------------
# Sampling a decoy set
# ---------------------------

def build_catalog_by_gene(catalog: pd.DataFrame) -> dict[str, list[dict]]:
    out = defaultdict(list)
    for row in catalog.to_dict("records"):
        out[row["gene"]].append(row)
    return dict(out)


def build_enhancer_conflict_graph(catalog: pd.DataFrame) -> dict[str, set[str]]:
    graph = defaultdict(set)
    for _, sub in catalog.groupby("gene"):
        enh_ids = sorted(set(sub["enhancer_id"].astype(str)))
        for enh_id in enh_ids:
            graph[enh_id].add(enh_id)
        for i, enh_a in enumerate(enh_ids):
            for enh_b in enh_ids[i + 1:]:
                graph[enh_a].add(enh_b)
                graph[enh_b].add(enh_a)
    return {k: set(v) for k, v in graph.items()}


def sample_decoy_set(
    catalog_by_gene: dict,
    size_n: int,
    forbid_names: set | None = None,
    forbid_enhancer_ids: set | None = None,
    used_enhancer_ids: set | None = None,
    conflict_graph: dict | None = None,
    rng: random.Random | None = None,
) -> list:
    rng = rng or random
    genes = list(catalog_by_gene.keys())
    rng.shuffle(genes)

    picked = []
    seen_stems = set()
    seen_enhancer_ids = set()
    forbid = set(forbid_names or [])
    blocked_enhancer_ids = set(forbid_enhancer_ids or set()) | set(used_enhancer_ids or set())
    conflicts = conflict_graph or {}
    blocked = seen_stems | forbid

    for g in genes:
        if len(picked) >= size_n:
            break
        candidates = []
        for row in catalog_by_gene[g]:
            enh_id = str(row["enhancer_id"])
            if row["filename_stem"] in blocked:
                continue
            if enh_id in blocked_enhancer_ids or enh_id in seen_enhancer_ids:
                continue
            if any(prev_id in conflicts.get(enh_id, set()) for prev_id in seen_enhancer_ids):
                continue
            candidates.append(row)
        if not candidates:
            continue
        row = candidates[rng.randrange(len(candidates))]
        picked.append(row)
        seen_stems.add(row["filename_stem"])
        seen_enhancer_ids.add(str(row["enhancer_id"]))
        blocked.add(row["filename_stem"])

    return picked if len(picked) == size_n else []


# ---------------------------
# Observed enhancer-level proportions
# ---------------------------

def split_pair_ids(s: str):
    parts = str(s).split("_")
    return (parts[0], parts[1]) if len(parts) >= 2 else (str(s), None)


def observed_enhancer_proportions(observed_pairs_csv: Path, gene_sizes_csv: Path) -> dict:
    sizes = pd.read_csv(gene_sizes_csv)
    if "Gene Name" not in sizes.columns and "gene" in sizes.columns:
        sizes = sizes.rename(columns={"gene": "Gene Name"})

    obs = pd.read_csv(observed_pairs_csv)

    hits_by_gene: dict[str, int] = {}
    for g, gdf in obs.groupby("Gene Name"):
        touched = set()
        for val in gdf["Comparisons"].dropna().astype(str):
            a, b = split_pair_ids(val)
            touched.add(a)
            if b is not None:
                touched.add(b)
        hits_by_gene[g] = len(touched)

    per_gene_enh = pd.Series(hits_by_gene, name="unique_hit_enhancers").reset_index().rename(columns={"index": "Gene Name"})
    df = sizes.merge(per_gene_enh, on="Gene Name", how="left")
    df["unique_hit_enhancers"] = df["unique_hit_enhancers"].fillna(0).astype(int)

    def bucket(n: int) -> str:
        if n == 2:
            return "2"
        if n == 3:
            return "3"
        if n >= 4:
            return ">=4"
        return "<2"

    df["bucket"] = df["num_shadow_enhancers"].apply(bucket)

    out = {}
    for b in ["2", "3", ">=4"]:
        sub = df[df["bucket"] == b]
        total_enh = int(sub["num_shadow_enhancers"].sum())
        hit_enh = int(sub["unique_hit_enhancers"].sum())
        out[b] = (hit_enh / total_enh) if total_enh > 0 else float("nan")
    return out


# ---------------------------
# Null enhancer-level metrics
# ---------------------------

def pooled_null_hit_props(
    decoy_sets: list,
    body_params: dict,
    flank_params: dict,
    seq_cache: dict,
    body_hit_cache: dict,
    flank_hit_cache: dict,
    flank_paths: dict,
    flank_coords: dict,
    dry_run: bool = False,
) -> dict:
    touched_body = set()
    touched_double_flank = set()
    touched_duplication = set()
    total_enh = 0

    def cached_seq(path_str: str) -> str:
        seq = seq_cache.get(path_str)
        if seq is None:
            seq = read_sequence(Path(path_str))
            seq_cache[path_str] = seq
        return seq

    def cached_body_hit(path_a: str, path_b: str) -> bool:
        key = tuple(sorted((path_a, path_b)))
        if key not in body_hit_cache:
            seq_a = cached_seq(path_a)
            seq_b = cached_seq(path_b)
            if not seq_a or not seq_b:
                body_hit_cache[key] = False
            else:
                body_hit_cache[key] = blast_has_any_hit(seq_a, seq_b, body_params.get("blastn", "blastn"), body_params, dry_run=dry_run)
        return body_hit_cache[key]

    def cached_flank_hit(path_a: str, path_b: str) -> bool:
        key = tuple(sorted((path_a, path_b)))
        if key not in flank_hit_cache:
            seq_a = cached_seq(path_a)
            seq_b = cached_seq(path_b)
            if not seq_a or not seq_b:
                flank_hit_cache[key] = False
            else:
                flank_hit_cache[key] = blast_has_any_hit(seq_a, seq_b, flank_params.get("blastn", "blastn"), flank_params, dry_run=dry_run)
        return flank_hit_cache[key]

    for decoy in decoy_sets:
        n = len(decoy)
        total_enh += n
        for i in range(n):
            for j in range(i + 1, n):
                a, b = decoy[i], decoy[j]
                pa, pb = a["path"], b["path"]
                key_a = (a["gene"], a["enhancer_id"])
                key_b = (b["gene"], b["enhancer_id"])

                body_hit = cached_body_hit(pa, pb)
                if body_hit:
                    touched_body.add(pa)
                    touched_body.add(pb)

                if (
                    key_a not in flank_paths
                    or key_b not in flank_paths
                    or not {"LEFT", "RIGHT"}.issubset(flank_paths[key_a])
                    or not {"LEFT", "RIGHT"}.issubset(flank_paths[key_b])
                    or any_flank_overlap_for_pair(key_a, key_b, flank_coords)
                ):
                    continue

                ll = cached_flank_hit(flank_paths[key_a]["LEFT"], flank_paths[key_b]["LEFT"])
                rr = cached_flank_hit(flank_paths[key_a]["RIGHT"], flank_paths[key_b]["RIGHT"])
                lr = cached_flank_hit(flank_paths[key_a]["LEFT"], flank_paths[key_b]["RIGHT"])
                rl = cached_flank_hit(flank_paths[key_a]["RIGHT"], flank_paths[key_b]["LEFT"])
                double_flank = (ll and rr) or (lr and rl)

                if double_flank:
                    touched_double_flank.add(pa)
                    touched_double_flank.add(pb)
                    if body_hit:
                        touched_duplication.add(pa)
                        touched_duplication.add(pb)

    denom = total_enh if total_enh > 0 else math.nan
    return {
        "body_hit_prop": (len(touched_body) / denom) if total_enh > 0 else float("nan"),
        "double_flank_hit_prop": (len(touched_double_flank) / denom) if total_enh > 0 else float("nan"),
        "duplication_hit_prop": (len(touched_duplication) / denom) if total_enh > 0 else float("nan"),
        "body_hit_enhancers": len(touched_body),
        "double_flank_hit_enhancers": len(touched_double_flank),
        "duplication_hit_enhancers": len(touched_duplication),
        "total_enhancers": total_enh,
    }


# ---------------------------
# Parallel rep workers
# ---------------------------

def init_rep_worker(
    catalog_by_gene: dict,
    gene_to_n: dict,
    target: dict,
    real_stems: dict,
    real_enhancer_ids: dict,
    enhancer_conflict_graph: dict,
    body_params: dict,
    flank_params: dict,
    flank_paths: dict,
    flank_coords: dict,
    dry_run: bool,
    base_seed: int,
):
    global WORKER_STATE
    WORKER_STATE = {
        "catalog_by_gene": catalog_by_gene,
        "gene_to_n": gene_to_n,
        "target": target,
        "real_stems": real_stems,
        "real_enhancer_ids": real_enhancer_ids,
        "enhancer_conflict_graph": enhancer_conflict_graph,
        "body_params": body_params,
        "flank_params": flank_params,
        "flank_paths": flank_paths,
        "flank_coords": flank_coords,
        "dry_run": dry_run,
        "base_seed": base_seed,
        "seq_cache": {},
        "body_hit_cache": {},
        "flank_hit_cache": {},
    }


def run_rep_task(task: tuple[str, int]) -> dict:
    bucket_label, rep_index = task
    st = WORKER_STATE
    rng = random.Random(st["base_seed"] + rep_index + (101 * sum(ord(c) for c in bucket_label)))
    genes = st["target"].get(bucket_label, [])
    decoys_for_bucket = []
    decoy_rows = []
    used_enhancer_ids = set()

    for g in genes:
        n = st["gene_to_n"].get(g, 0)
        if n < 2:
            continue
        forbid = st["real_stems"].get(g, set())
        forbid_enhancer_ids = st["real_enhancer_ids"].get(g, set()) | used_enhancer_ids
        picked = []
        for _ in range(50):
            cand = sample_decoy_set(
                st["catalog_by_gene"],
                n,
                forbid_names=forbid,
                forbid_enhancer_ids=forbid_enhancer_ids,
                used_enhancer_ids=used_enhancer_ids,
                conflict_graph=st["enhancer_conflict_graph"],
                rng=rng,
            )
            if cand:
                picked = cand
                break
        if not picked:
            continue

        decoys_for_bucket.append(picked)
        used_enhancer_ids.update(str(row["enhancer_id"]) for row in picked)
        for idx, row in enumerate(picked):
            decoy_rows.append(
                {
                    "bucket": bucket_label,
                    "rep": rep_index,
                    "gene_name": g,
                    "decoy_index_in_set": idx,
                    "decoy_gene_source": row["gene"],
                    "decoy_enhancer_path": row["path"],
                    "decoy_filename_stem": row["filename_stem"],
                }
            )

    pooled = (
        pooled_null_hit_props(
            decoys_for_bucket,
            st["body_params"],
            st["flank_params"],
            st["seq_cache"],
            st["body_hit_cache"],
            st["flank_hit_cache"],
            st["flank_paths"],
            st["flank_coords"],
            dry_run=st["dry_run"],
        )
        if decoys_for_bucket
        else {
            "body_hit_prop": float("nan"),
            "double_flank_hit_prop": float("nan"),
            "duplication_hit_prop": float("nan"),
            "body_hit_enhancers": 0,
            "double_flank_hit_enhancers": 0,
            "duplication_hit_enhancers": 0,
            "total_enhancers": 0,
        }
    )

    return {
        "null_row": {
            "bucket": bucket_label,
            "rep": rep_index,
            "pooled_null_prop_enh": pooled["body_hit_prop"],
            "pooled_null_prop_double_flank": pooled["double_flank_hit_prop"],
            "pooled_null_prop_duplication": pooled["duplication_hit_prop"],
            "body_hit_enhancers": pooled["body_hit_enhancers"],
            "double_flank_hit_enhancers": pooled["double_flank_hit_enhancers"],
            "duplication_hit_enhancers": pooled["duplication_hit_enhancers"],
            "total_enhancers": pooled["total_enhancers"],
        },
        "decoy_rows": decoy_rows,
        "bucket": bucket_label,
        "rep": rep_index,
    }


# ---------------------------
# Plotting
# ---------------------------

def write_null_summary_plot(summary_df: pd.DataFrame, out_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        print("WARNING: matplotlib not available; skipping plot:", e)
        return

    bucket_order = [b for b in ["2", "3", ">=4"] if b in set(summary_df["bucket"])]
    labels = {"2": "2 shadows/set", "3": "3 shadows/set", ">=4": ">=4 shadows/set"}

    dfp = summary_df.copy().set_index("bucket").reindex(bucket_order).reset_index()
    means = dfp["null_mean_duplication_prop"].to_numpy(dtype=float)
    sds = dfp["null_sd_duplication_prop"].to_numpy(dtype=float)
    nreps = dfp["n_reps"].to_numpy(dtype=float)

    fig = plt.figure(figsize=(7.2, 4.2))
    ax = fig.add_subplot(111)
    x = np.arange(len(bucket_order))

    ax.bar(x, means, edgecolor="white", linewidth=1.2, alpha=0.9)
    ax.errorbar(x, means, yerr=sds, fmt="none", capsize=4, linewidth=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels([labels[b] for b in bucket_order])
    ax.set_ylabel("Null expected fraction duplicated")
    ax.set_xlabel("Shadow enhancers per set")

    for i, _ in enumerate(bucket_order):
        reps_i = int(nreps[i]) if i < len(nreps) and not np.isnan(nreps[i]) else 0
        ax.text(i, -0.12, f"reps={reps_i}\nmean={means[i]:.3f}\nsd={sds[i]:.3f}", ha="center", va="top", transform=ax.get_xaxis_transform(), fontsize=9)

    ymax = np.nanmax(means + sds) if np.isfinite(np.nanmax(means + sds)) else 0.0
    ax.set_ylim(0, float(ymax) + 0.05)
    ax.set_title("Null control: expected duplication-hit fraction\n(mean +/- sd across reps)")
    plt.tight_layout()
    fig.savefig(out_dir / "null_dup_fraction_by_shadow_category.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "null_dup_fraction_by_shadow_category.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Enhancer-level null distributions per shadow-set bucket (2,3,>=4).")
    ap.add_argument("--no_plot", action="store_true", help="Do not write summary plot PNG/PDF")
    ap.add_argument("--shadow_zip", help="Zip with gene subdirs of enhancer FASTA/TXT files")
    ap.add_argument("--enhancer_root", help="Existing root directory with per-gene shadow enhancer FASTA/TXT files")
    ap.add_argument("--observed_pairs_csv", help="CSV with columns: Gene Name, Comparisons")
    ap.add_argument("--gene_sizes_csv", help="CSV with Gene Name,num_shadow_enhancers. If omitted, derive set sizes from enhancer_root/shadow_zip.")
    ap.add_argument("--out_dir", help="Directory for outputs. Defaults to the gene_sizes_csv parent, enhancer_root parent, or shadow_zip parent.")
    ap.add_argument("--flank_manifest", help="Manifest TSV describing flank FASTAs for the same shadow enhancers")
    ap.add_argument("--flanks_root", help="Root directory that manifest fasta_path values are relative to")
    ap.add_argument("--buckets", default="2,3,>=4", help="Comma-separated bucket subset to run: any of 2,3,>=4")
    ap.add_argument("--reps", type=int, default=200, help="Null repetitions per bucket")
    ap.add_argument("--rep_workers", type=int, default=1, help="Number of worker processes used to parallelize reps")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--blastn", default="blastn")
    ap.add_argument("--timeout", type=int, default=90)
    ap.add_argument("--evalue", default="0.01", help="Shared fallback BLAST e-value cutoff if body/flank-specific cutoffs are not provided")
    ap.add_argument("--enhancer_hit_evalue", help="BLAST e-value cutoff for enhancer-body hits")
    ap.add_argument("--double_flank_hit_evalue", help="BLAST e-value cutoff for flank-vs-flank hits used to call double-flank support")
    ap.add_argument("--dry_run", action="store_true", help="Practice mode: simulate hits, no BLAST")
    ap.add_argument("--verbose", action="store_true", help="Print progress")
    ap.add_argument("--skip_observed", action="store_true", help="Skip observed calc, only output null distributions")
    args = ap.parse_args()

    if not args.shadow_zip and not args.enhancer_root:
        raise SystemExit("[ERROR] Pass either --shadow_zip or --enhancer_root")
    if args.shadow_zip and args.enhancer_root:
        raise SystemExit("[ERROR] Pass only one of --shadow_zip or --enhancer_root")
    if (not args.skip_observed) and not args.observed_pairs_csv:
        raise SystemExit("[ERROR] Pass --observed_pairs_csv unless using --skip_observed")
    if args.rep_workers < 1:
        raise SystemExit("[ERROR] --rep_workers must be >= 1")

    selected_buckets = [b.strip() for b in str(args.buckets).split(",") if b.strip()]
    allowed_buckets = {"2", "3", ">=4"}
    invalid_buckets = [b for b in selected_buckets if b not in allowed_buckets]
    if invalid_buckets:
        raise SystemExit(f"[ERROR] Invalid bucket(s): {invalid_buckets}. Allowed: 2, 3, >=4")
    if not selected_buckets:
        raise SystemExit("[ERROR] --buckets must include at least one of: 2,3,>=4")

    if args.out_dir:
        out_dir = Path(args.out_dir).resolve()
    elif args.gene_sizes_csv:
        out_dir = Path(args.gene_sizes_csv).parent
    elif args.enhancer_root:
        out_dir = Path(args.enhancer_root).resolve().parent
    else:
        out_dir = Path(args.shadow_zip).resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.enhancer_root:
        catalog = build_catalog_from_root(Path(args.enhancer_root).resolve())
    else:
        work = out_dir / "_shadowsets_extracted"
        work.mkdir(parents=True, exist_ok=True)
        catalog = build_catalog(Path(args.shadow_zip), work)

    derived_gene_sizes_path = out_dir / "_derived_gene_sizes.csv"
    if args.gene_sizes_csv:
        sizes_df = pd.read_csv(args.gene_sizes_csv)
        if "Gene Name" not in sizes_df.columns and "gene" in sizes_df.columns:
            sizes_df = sizes_df.rename(columns={"gene": "Gene Name"})
    else:
        sizes_df = derive_gene_sizes_from_catalog(catalog)
        sizes_df.to_csv(derived_gene_sizes_path, index=False)

    def bucket(n: int) -> str:
        if n == 2:
            return "2"
        if n == 3:
            return "3"
        if n >= 4:
            return ">=4"
        return "<2"

    sizes_df["bucket"] = sizes_df["num_shadow_enhancers"].apply(bucket)
    target = {b: list(sizes_df[sizes_df["bucket"] == b]["Gene Name"]) for b in selected_buckets}
    gene_to_n = dict(zip(sizes_df["Gene Name"], sizes_df["num_shadow_enhancers"]))

    real_stems = defaultdict(set)
    real_enhancer_ids = defaultdict(set)
    for bucket_label in selected_buckets:
        for g in target[bucket_label]:
            real_stems[g] |= set(catalog[catalog["gene"] == g]["filename_stem"].tolist())
            real_enhancer_ids[g] |= set(catalog[catalog["gene"] == g]["enhancer_id"].astype(str).tolist())

    enhancer_hit_evalue = args.enhancer_hit_evalue or args.evalue
    double_flank_hit_evalue = args.double_flank_hit_evalue or args.evalue
    body_params = {
        "evalue": enhancer_hit_evalue,
        "word_size": "7",
        "gapopen": "5",
        "gapextend": "2",
        "reward": "2",
        "penalty": "-3",
        "dust": "yes",
        "timeout": args.timeout,
        "blastn": args.blastn,
    }
    flank_params = {
        "evalue": double_flank_hit_evalue,
        "word_size": "7",
        "gapopen": "5",
        "gapextend": "2",
        "reward": "2",
        "penalty": "-3",
        "dust": "yes",
        "timeout": args.timeout,
        "blastn": args.blastn,
    }

    flank_paths = {}
    flank_coords = {}
    if args.flank_manifest or args.flanks_root:
        if not args.flank_manifest or not args.flanks_root:
            raise SystemExit("[ERROR] Pass both --flank_manifest and --flanks_root")
        flank_paths, flank_coords = read_flank_manifest(
            Path(args.flank_manifest).resolve(),
            Path(args.flanks_root).resolve(),
            allowed_keys=set(zip(catalog["gene"], catalog["enhancer_id"])),
        )

    if args.verbose:
        print("Buckets (genes):", {k: len(v) for k, v in target.items()})
        print("Selected buckets:", selected_buckets)
        print("Hit cutoffs:", {"enhancer_hit_evalue": enhancer_hit_evalue, "double_flank_hit_evalue": double_flank_hit_evalue})
        print("Rep workers:", args.rep_workers)

    catalog_by_gene = build_catalog_by_gene(catalog)
    enhancer_conflict_graph = build_enhancer_conflict_graph(catalog)
    null_rows = []
    decoy_rows = []
    tasks = [(bucket_label, r) for bucket_label in selected_buckets for r in range(args.reps) if target[bucket_label]]

    if args.rep_workers == 1:
        init_rep_worker(
            catalog_by_gene,
            gene_to_n,
            target,
            dict(real_stems),
            dict(real_enhancer_ids),
            enhancer_conflict_graph,
            body_params,
            flank_params,
            flank_paths,
            flank_coords,
            args.dry_run,
            args.seed,
        )
        for bucket_label, rep_index in tasks:
            if args.verbose:
                print(f"[{bucket_label}] rep {rep_index + 1}/{args.reps} ...", end="", flush=True)
            result = run_rep_task((bucket_label, rep_index))
            null_rows.append(result["null_row"])
            decoy_rows.extend(result["decoy_rows"])
            if args.verbose:
                print(" done")
    else:
        max_workers = min(args.rep_workers, len(tasks))
        ctx = mp.get_context("fork")
        with ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=ctx,
            initializer=init_rep_worker,
            initargs=(
                catalog_by_gene,
                gene_to_n,
                target,
                dict(real_stems),
                dict(real_enhancer_ids),
                enhancer_conflict_graph,
                body_params,
                flank_params,
                flank_paths,
                flank_coords,
                args.dry_run,
                args.seed,
            ),
        ) as ex:
            future_map = {ex.submit(run_rep_task, task): task for task in tasks}
            for fut in as_completed(future_map):
                result = fut.result()
                null_rows.append(result["null_row"])
                decoy_rows.extend(result["decoy_rows"])
                if args.verbose:
                    print(f"[{result['bucket']}] rep {result['rep'] + 1}/{args.reps} done")

    null_df = pd.DataFrame(null_rows).sort_values(["bucket", "rep"]).reset_index(drop=True)
    null_path = out_dir / "null_distributions_enhancer.csv"
    null_df.to_csv(null_path, index=False)

    decoy_df = pd.DataFrame(decoy_rows).sort_values(["bucket", "rep", "gene_name", "decoy_index_in_set"]).reset_index(drop=True)
    decoy_path = out_dir / "null_decoy_sets_enhancer.csv"
    decoy_df.to_csv(decoy_path, index=False)

    observed = {b: float("nan") for b in selected_buckets}
    if not args.skip_observed:
        observed = observed_enhancer_proportions(
            Path(args.observed_pairs_csv),
            Path(args.gene_sizes_csv) if args.gene_sizes_csv else derived_gene_sizes_path,
        )

    summary_rows = []
    for b in selected_buckets:
        sub_body = null_df[null_df["bucket"] == b]["pooled_null_prop_enh"].dropna()
        sub_flank = null_df[null_df["bucket"] == b]["pooled_null_prop_double_flank"].dropna()
        sub_dup = null_df[null_df["bucket"] == b]["pooled_null_prop_duplication"].dropna()
        mu = sub_body.mean() if len(sub_body) > 0 else float("nan")
        sd = sub_body.std(ddof=1) if len(sub_body) > 1 else float("nan")
        obs = observed.get(b, float("nan"))
        z = (obs - mu) / sd if (not math.isnan(obs)) and (not math.isnan(mu)) and sd and sd > 0 else float("nan")
        summary_rows.append(
            {
                "bucket": b,
                "observed_enhancer_prop": obs,
                "null_mean_enhancer_prop": mu,
                "null_sd_enhancer_prop": sd,
                "null_mean_double_flank_prop": sub_flank.mean() if len(sub_flank) > 0 else float("nan"),
                "null_sd_double_flank_prop": sub_flank.std(ddof=1) if len(sub_flank) > 1 else float("nan"),
                "null_mean_duplication_prop": sub_dup.mean() if len(sub_dup) > 0 else float("nan"),
                "null_sd_duplication_prop": sub_dup.std(ddof=1) if len(sub_dup) > 1 else float("nan"),
                "z_score": z,
                "n_reps": int(len(sub_body)),
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_path = out_dir / "null_zscores_summary_enhancer.csv"
    summary_df.to_csv(summary_path, index=False)

    if not args.no_plot:
        write_null_summary_plot(summary_df, out_dir)

    if args.verbose:
        print("\nWrote:")
        print(" ", null_path)
        print(" ", decoy_path)
        print(" ", summary_path)


if __name__ == "__main__":
    main()
