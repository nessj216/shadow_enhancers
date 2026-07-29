#!/usr/bin/env python3

import argparse
import itertools
import math
import random
import re
import shlex
import shutil
import subprocess as sp
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock
from typing import Dict, List, NamedTuple, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd


PAIR_RE = re.compile(
    r"^(?P<chrom1>[^_]+)_(?P<s1>\d+)-(?P<e1>\d+)__(?P<chrom2>[^_]+)_(?P<s2>\d+)-(?P<e2>\d+)$"
)
RELAXED_LASTZ = (
    "--strand=both --seed=12of19 --transition=2 --step=1 "
    "--ambiguous=iupac --hspthresh=1500 --gappedthresh=1500 "
    "--ydrop=5400 --chain --format=maf+"
)
DEFAULT_FLANK_BLAST_ARGS = "-task blastn -strand both -dust no -soft_masking false"
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


class EnhancerRecord(NamedTuple):
    member_id: str
    enhancer_id: str
    set_label: str
    seq_path: Path
    seq_header: str
    enhancer_file: str


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Shadow-duplication null workflows. "
            "Default mode samples null enhancer sets matched to observed set sizes, "
            "runs BLAST on enhancer sequences, runs LASTZ on BLAST-positive flank "
            "pairs, and writes one null distribution for each set-size bin (2, 3, >=4). "
            "The flank-null-blast mode instead reads a flank manifest directly and "
            "builds fly-style BLAST double-hit null cutoff tables. "
            "The validate-flank-manifest mode only checks manifest/path consistency."
        ),
        epilog=(
            "SCC shortcut example:\n"
            "  python shadow_duplication_null_pipeline.py "
            "--base-dir /project/wunderl/Jillian --reps 100 --compute-observed --overwrite"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--mode",
        choices=["shadow-null-pipeline", "flank-null-blast", "validate-flank-manifest"],
        default="shadow-null-pipeline",
        help=(
            "shadow-null-pipeline: existing enhancer-set null workflow. "
            "flank-null-blast: fly-style BLAST double-hit scan over a flank manifest. "
            "validate-flank-manifest: only validate manifest rows and resolved flank paths."
        ),
    )
    p.add_argument("--base-dir",
                   help=(
                       "Base project directory. If set, defaults become: "
                       "<base-dir>/FlilesOutput, <base-dir>/shadow_pairs_flanks, "
                       "and <base-dir>/null_shadow_duplication_run"
                   ))
    p.add_argument("--enhancer-root",
                   help="Root directory of shadow-set enhancer FASTAs/TXT files (FlilesOutput)")
    p.add_argument("--source-flanks-root",
                   help="Root holding observed shadow pair flank FASTAs")
    p.add_argument("--source-manifest",
                   help="Manifest TSV describing the flank FASTAs. If omitted, autodetect inside source-flanks-root.")
    p.add_argument("--flanks-root",
                   help="Flank FASTA root for --mode flank-null-blast. Alias of --source-flanks-root in that mode.")
    p.add_argument("--manifest",
                   help="Flank manifest TSV for --mode flank-null-blast. Alias of --source-manifest in that mode.")
    p.add_argument("--outdir",
                   help="Output directory for summaries and optional pair tables")
    p.add_argument("--reps", type=int, default=100,
                   help="Null replicates per run")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--threads", type=int, default=8,
                   help="Concurrent BLAST/LASTZ workers")
    p.add_argument("--max-sample-attempts-per-set", type=int, default=200,
                   help="Attempts to sample a valid null set for one observed set")

    p.add_argument("--blastn", default="blastn")
    p.add_argument("--blast-evalue", type=float, default=1e-5)
    p.add_argument("--blast-word-size", type=int, default=7)
    p.add_argument("--blast-timeout", type=int, default=120)
    p.add_argument("--flank-blast-args", default=DEFAULT_FLANK_BLAST_ARGS,
                   help="Extra blastn args for --mode flank-null-blast")
    p.add_argument("--flank-hit-evalue-max", type=float, default=1.0,
                   help="Maximum e-value retained as a flank hit in --mode flank-null-blast")
    p.add_argument("--cutoff-percentiles", default=DEFAULT_CUTOFF_PERCENTILES,
                   help="Percentiles for fly-style null cutoff tables; comma, colon, semicolon, or whitespace separated")
    p.add_argument("--zero-evalue-floor", type=float, default=1e-300,
                   help="Replacement value when BLAST reports evalue 0.0 for -log10 transforms")
    p.add_argument("--pair-limit", type=int, default=0,
                   help="Optional cap on pair_dir groups in --mode flank-null-blast")
    p.add_argument("--lastz", default="lastz")
    p.add_argument("--lastz-args", default=RELAXED_LASTZ)
    p.add_argument("--lastz-timeout", type=int, default=120)
    p.add_argument("--gap-threshold", type=int, default=5000)
    p.add_argument("--len-threshold", type=int, default=40,
                   help="Require best_block_len > this")
    p.add_argument("--identity-threshold", type=float, default=65.0,
                   help="Require best_block_identity_pct > this")

    p.add_argument("--compute-observed", action="store_true",
                   help="Also compute the observed per-bin hit fractions from the real sets")
    p.add_argument("--without-replacement-across-sets", action="store_true",
                   help="Within a null replicate, do not reuse an enhancer across different sampled sets")
    p.add_argument("--allow-same-source-set-within-null-set", action="store_true",
                   help="Allow multiple enhancers from the same source set inside one sampled null set")
    p.add_argument("--exclude-observed-pairs-globally", action="store_true",
                   help="Do not allow a sampled null pair if it is one of the observed real shadow pairs")
    p.add_argument("--skip-overlapping-flanks", action="store_true",
                   help="Skip sampled same-chromosome pairs whose flank intervals overlap")
    p.add_argument("--write-null-pair-results", action="store_true",
                   help="Write the full null pair-level result table (can be large)")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def is_sequence_file(path: Path) -> bool:
    name = path.name.lower()
    if path.name.startswith("._"):
        return False
    if path.name.startswith("."):
        return False
    return path.is_file() and name.endswith((
        ".fa", ".fasta", ".fna", ".fas", ".txt"
    ))


def read_first_fasta_record(path: Path) -> Tuple[str, str]:
    header = None
    seq_chunks: List[str] = []
    encodings = ("utf-8", "latin-1")
    last_error = None
    for encoding in encodings:
        try:
            with path.open(encoding=encoding, errors="strict") as fh:
                for raw in fh:
                    line = raw.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if header is not None:
                            break
                        header = line[1:].strip()
                        continue
                    if header is not None:
                        seq_chunks.append(line)
            break
        except UnicodeDecodeError as exc:
            header = None
            seq_chunks = []
            last_error = exc
            continue
    if header is None:
        if last_error is not None:
            raise ValueError(f"Could not decode FASTA file {path}: {last_error}") from last_error
        raise ValueError(f"{path} does not contain a FASTA header")
    return header, "".join(seq_chunks)


def parse_enhancer_id(enhancer_id: str) -> Tuple[str, int, int]:
    chrom, rest = str(enhancer_id).split(":")
    start, end = rest.split("-")
    return chrom, int(start), int(end)


def enhancer_id_from_header(header: str) -> str:
    token = header.split()[0]
    if ":" not in token or "-" not in token:
        raise ValueError(f"Cannot parse enhancer id from FASTA header: {header}")
    return token


def bucket_for_size(n: int) -> str:
    if n == 2:
        return "2"
    if n == 3:
        return "3"
    if n >= 4:
        return ">=4"
    return "<2"


def make_member_id(set_label: str, seq_path: Path) -> str:
    return f"{set_label}::{seq_path.name}"


def manifest_member_key(set_label: str, enhancer_id: str) -> Tuple[str, str]:
    return (set_label, enhancer_id)


def member_manifest_key(member_id: str, enhancer_records: Dict[str, EnhancerRecord]) -> Tuple[str, str]:
    rec = enhancer_records[member_id]
    return manifest_member_key(rec.set_label, rec.enhancer_id)


def member_sort_key(member_id: str, enhancer_records: Dict[str, EnhancerRecord]):
    rec = enhancer_records[member_id]
    chrom, start, end = parse_enhancer_id(rec.enhancer_id)
    return (chrom, start, end, rec.set_label, rec.enhancer_file)


def member_biological_id(member_id: str, enhancer_records: Dict[str, EnhancerRecord]) -> str:
    return enhancer_records[member_id].enhancer_file


def member_pair_ids(a: str, b: str) -> Tuple[str, str]:
    return tuple(sorted((a, b)))


def genomic_pair_for_members(a: str, b: str, enhancer_records: Dict[str, EnhancerRecord]) -> Tuple[str, str]:
    return canonical_pair(enhancer_records[a].enhancer_id, enhancer_records[b].enhancer_id)


def pair_name_for_members(a: str, b: str, enhancer_records: Dict[str, EnhancerRecord]) -> str:
    return canonical_pair_name(enhancer_records[a].enhancer_id, enhancer_records[b].enhancer_id)


def canonical_pair(a: str, b: str) -> Tuple[str, str]:
    a_t = parse_enhancer_id(a)
    b_t = parse_enhancer_id(b)
    if a_t[0] == b_t[0]:
        return (a, b) if a_t[1] <= b_t[1] else (b, a)
    return (a, b) if a_t <= b_t else (b, a)


def canonical_pair_name(a: str, b: str) -> str:
    left, right = canonical_pair(a, b)
    lc, ls, le = parse_enhancer_id(left)
    rc, rs, re_ = parse_enhancer_id(right)
    return f"{lc}_{ls}-{le}__{rc}_{rs}-{re_}"


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def fasta_len(path: Path) -> int:
    n = 0
    encodings = ("utf-8", "latin-1")
    last_error = None
    for encoding in encodings:
        try:
            with path.open(encoding=encoding, errors="strict") as fh:
                for line in fh:
                    if line.startswith(">"):
                        continue
                    n += len(line.strip())
            return n
        except UnicodeDecodeError as exc:
            n = 0
            last_error = exc
            continue
    if last_error is not None:
        raise ValueError(f"Could not decode FASTA file {path}: {last_error}") from last_error
    return n


def infer_flank_coords(chrom: str, enh_start: int, enh_end: int, side: str, flank_len: int):
    if side == "LEFT":
        return chrom, max(0, enh_start - flank_len), enh_start
    if side == "RIGHT":
        return chrom, enh_end, enh_end + flank_len
    raise ValueError(side)


def read_enhancer_catalog(enhancer_root: Path):
    enhancer_records: Dict[str, EnhancerRecord] = {}
    set_to_enhancers: Dict[str, List[str]] = defaultdict(list)
    rows = []

    for set_dir in sorted([p for p in enhancer_root.iterdir() if p.is_dir()]):
        for seq_path in sorted([p for p in set_dir.iterdir() if is_sequence_file(p)]):
            header, _ = read_first_fasta_record(seq_path)
            enh_id = enhancer_id_from_header(header)
            member_id = make_member_id(set_dir.name, seq_path)
            if member_id in enhancer_records:
                raise SystemExit(
                    f"[ERROR] Duplicate set/file identity in enhancer catalog: {member_id}\n"
                    f"  seen in {enhancer_records[member_id].seq_path}\n"
                    f"  and    in {seq_path}"
                )
            rec = EnhancerRecord(
                member_id=member_id,
                enhancer_id=enh_id,
                set_label=set_dir.name,
                seq_path=seq_path.resolve(),
                seq_header=header,
                enhancer_file=seq_path.name,
            )
            enhancer_records[member_id] = rec
            set_to_enhancers[set_dir.name].append(member_id)
            rows.append({
                "member_id": member_id,
                "set_label": set_dir.name,
                "enhancer_id": enh_id,
                "enhancer_file": seq_path.name,
                "seq_path": str(seq_path.resolve()),
                "seq_header": header,
            })

    if not enhancer_records:
        raise SystemExit(f"[ERROR] No enhancer sequence files found under {enhancer_root}")

    return enhancer_records, set_to_enhancers, pd.DataFrame(rows)


def read_source_manifest(source_flanks_root: Path, source_manifest: Path):
    df = pd.read_csv(source_manifest, sep="\t")
    need = {"gene_label", "pair_dir", "enhancer_id", "flank_side", "fasta_path"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"[ERROR] source manifest missing columns: {sorted(missing)}")

    enh_to_flanks: Dict[Tuple[str, str], Dict[str, str]] = {}
    enh_to_source_set: Dict[Tuple[str, str], str] = {}
    enh_to_coords: Dict[Tuple[str, str], Tuple[str, int, int]] = {}
    flank_coords: Dict[Tuple[str, str], Dict[str, Tuple[str, int, int]]] = defaultdict(dict)
    observed_pairs: Set[Tuple[str, str]] = set()

    for _, row in df.iterrows():
        set_label = str(row["gene_label"])
        enh_id = str(row["enhancer_id"])
        member_key = manifest_member_key(set_label, enh_id)
        side = str(row["flank_side"]).upper()
        rel = str(row["fasta_path"])
        fpath = source_flanks_root / rel

        if side not in {"LEFT", "RIGHT"}:
            continue
        if not fpath.exists():
            continue

        enh_to_flanks.setdefault(member_key, {})[side] = rel
        enh_to_source_set.setdefault(member_key, set_label)

        if member_key not in enh_to_coords:
            enh_to_coords[member_key] = parse_enhancer_id(enh_id)

        chrom = str(row.get("chrom", ""))
        start = pd.to_numeric(row.get("start", None), errors="coerce")
        end = pd.to_numeric(row.get("end", None), errors="coerce")
        if chrom and pd.notna(start) and pd.notna(end):
            flank_coords[member_key][side] = (chrom, int(start), int(end))

    valid = {
        member_key for member_key, sides in enh_to_flanks.items()
        if {"LEFT", "RIGHT"}.issubset(set(sides))
    }
    for member_key in valid:
        chrom, start, end = enh_to_coords[member_key]
        for side in ("LEFT", "RIGHT"):
            if side not in flank_coords[member_key]:
                flank_len = fasta_len(source_flanks_root / enh_to_flanks[member_key][side])
                flank_coords[member_key][side] = infer_flank_coords(chrom, start, end, side, flank_len)

    for pair_dir in df["pair_dir"].astype(str).unique():
        parts = pair_dir.split("/", 1)
        pair_name = parts[1] if len(parts) > 1 else parts[0]
        match = PAIR_RE.match(pair_name)
        if not match:
            continue
        a = f"{match.group('chrom1')}:{match.group('s1')}-{match.group('e1')}"
        b = f"{match.group('chrom2')}:{match.group('s2')}-{match.group('e2')}"
        observed_pairs.add(canonical_pair(a, b))

    return enh_to_flanks, enh_to_source_set, enh_to_coords, flank_coords, observed_pairs


def autodetect_manifest(source_flanks_root: Path) -> Path:
    required = {"gene_label", "pair_dir", "enhancer_id", "flank_side", "fasta_path"}

    def looks_like_manifest(path: Path) -> bool:
        try:
            with path.open() as fh:
                header = fh.readline().strip().split("\t")
            return required.issubset(set(header))
        except OSError:
            return False

    preferred_names = [
        "manifest.tsv",
        "source_manifest.tsv",
        "flanks_manifest.tsv",
        "shadow_pairs_flanks_manifest.tsv",
    ]
    for name in preferred_names:
        candidate = source_flanks_root / name
        if candidate.exists() and looks_like_manifest(candidate):
            return candidate.resolve()

    shallow_candidates = []
    shallow_patterns = [
        "*.tsv",
        "*manifest*.tsv",
        "*/*.tsv",
        "*/*manifest*.tsv",
    ]
    seen = set()
    for pattern in shallow_patterns:
        for path in sorted(source_flanks_root.glob(pattern)):
            if path.is_file() and path not in seen:
                shallow_candidates.append(path)
                seen.add(path)

    manifest_like = [p for p in shallow_candidates if looks_like_manifest(p)]
    if len(manifest_like) == 1:
        return manifest_like[0].resolve()
    if len(manifest_like) > 1:
        raise SystemExit(
            "[ERROR] Multiple candidate source manifests were found near the flanks root. "
            "Please pass --source-manifest explicitly:\n"
            + "\n".join(f"  {p}" for p in manifest_like)
        )

    raise SystemExit(
        "[ERROR] Could not quickly detect a source manifest near the flanks root. "
        "Pass --source-manifest explicitly, for example:\n"
        f"  --source-manifest {source_flanks_root}/manifest.tsv"
    )


def resolve_input_paths(args):
    base_dir = Path(args.base_dir).resolve() if args.base_dir else None

    enhancer_root = Path(args.enhancer_root).resolve() if args.enhancer_root else None
    source_flanks_root = Path(args.source_flanks_root).resolve() if args.source_flanks_root else None
    outdir = Path(args.outdir).resolve() if args.outdir else None

    if base_dir is not None:
        if enhancer_root is None:
            enhancer_root = base_dir / "FlilesOutput"
        if source_flanks_root is None:
            source_flanks_root = base_dir / "shadow_pairs_flanks"
        if outdir is None:
            outdir = base_dir / "null_shadow_duplication_run"

    if enhancer_root is None:
        raise SystemExit("[ERROR] Pass --enhancer-root or --base-dir")
    if source_flanks_root is None:
        raise SystemExit("[ERROR] Pass --source-flanks-root or --base-dir")
    if outdir is None:
        raise SystemExit("[ERROR] Pass --outdir or --base-dir")

    source_manifest = (
        Path(args.source_manifest).resolve()
        if args.source_manifest
        else autodetect_manifest(source_flanks_root)
    )

    return enhancer_root.resolve(), source_flanks_root.resolve(), source_manifest.resolve(), outdir.resolve()


def resolve_flank_blast_inputs(args):
    flanks_root = (
        Path(args.flanks_root).resolve()
        if args.flanks_root
        else (Path(args.source_flanks_root).resolve() if args.source_flanks_root else None)
    )
    manifest = (
        Path(args.manifest).resolve()
        if args.manifest
        else (
            Path(args.source_manifest).resolve()
            if args.source_manifest
            else (autodetect_manifest(flanks_root) if flanks_root else None)
        )
    )
    outdir = Path(args.outdir).resolve() if args.outdir else flanks_root

    if flanks_root is None:
        raise SystemExit("[ERROR] Pass --flanks-root or --source-flanks-root for --mode flank-null-blast")
    if manifest is None:
        raise SystemExit("[ERROR] Pass --manifest or --source-manifest for --mode flank-null-blast")
    if outdir is None:
        raise SystemExit("[ERROR] Could not determine output directory for --mode flank-null-blast")
    return flanks_root.resolve(), manifest.resolve(), outdir.resolve()


def parse_percentiles(text: str) -> List[float]:
    values = []
    for raw in re.split(r"[\s,;:]+", str(text).strip()):
        raw = raw.strip()
        if not raw:
            continue
        value = float(raw)
        if value < 0 or value > 100:
            raise SystemExit(f"[ERROR] Percentile must be between 0 and 100: {raw}")
        values.append(value)
    if not values:
        raise SystemExit("[ERROR] No valid values were provided in --cutoff-percentiles")
    return values


def pair_dir_parts(pair_dir: str):
    parts = str(pair_dir).split("/", 1)
    set_label = parts[0]
    pair_name = parts[1] if len(parts) > 1 else parts[0]
    match = PAIR_RE.match(pair_name)
    if not match:
        return None
    a = (
        match.group("chrom1"),
        int(match.group("s1")),
        int(match.group("e1")),
    )
    b = (
        match.group("chrom2"),
        int(match.group("s2")),
        int(match.group("e2")),
    )
    if a[0] == b[0] and a[1] > b[1]:
        a, b = b, a
    return set_label, a, b


def build_manifest_set_order(df: pd.DataFrame):
    set_order: Dict[str, Dict[str, List[Tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    seen = set()
    for pair_dir in df["pair_dir"].astype(str).drop_duplicates():
        parsed = pair_dir_parts(pair_dir)
        if not parsed:
            continue
        set_label, a, b = parsed
        for chrom, start, end in (a, b):
            key = (set_label, chrom, start, end)
            if key in seen:
                continue
            seen.add(key)
            set_order[set_label][chrom].append((start, end))
    for set_label in set_order:
        for chrom in set_order[set_label]:
            set_order[set_label][chrom].sort()
    return set_order


def read_blast_table_stats(out_path: Path, hit_evalue_max: float):
    empty = {
        "hits_count": 0,
        "best_evalue": math.nan,
        "best_bitscore": math.nan,
        "best_pident": math.nan,
        "best_length": 0,
        "best_qcovhsp": math.nan,
    }
    if (not out_path.exists()) or out_path.stat().st_size == 0:
        return empty

    df = pd.read_csv(out_path, sep="\t", names=BLAST_COLUMNS)
    if df.empty:
        return empty

    numeric_cols = [
        "pident", "length", "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore", "qcovhsp",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[df["evalue"] <= float(hit_evalue_max)].copy()
    if df.empty:
        return empty

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


def run_flank_blast_once(subject_fa: Path, query_fa: Path, out_path: Path, log_path: Path, args) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        args.blastn,
        "-subject", str(subject_fa),
        "-query", str(query_fa),
        "-evalue", str(args.flank_hit_evalue_max),
        "-outfmt", BLAST_OUTFMT,
    ] + shlex.split(args.flank_blast_args)
    with log_path.open("w") as logfh, out_path.open("w") as outfh:
        try:
            proc = sp.run(
                cmd,
                stdout=outfh,
                stderr=logfh,
                check=False,
                timeout=args.blast_timeout,
            )
            return proc.returncode
        except FileNotFoundError:
            logfh.write(f"[ERROR] blast executable not found: {args.blastn}\n")
            return 127


def flank_blast_worker(task, args):
    out_path = Path(task["blast_path"])
    log_path = Path(task["log_path"])
    if args.overwrite:
        out_path.unlink(missing_ok=True)
        log_path.unlink(missing_ok=True)

    if out_path.exists() and not args.overwrite:
        stats = read_blast_table_stats(out_path, args.flank_hit_evalue_max)
        status = "SKIP" if out_path.stat().st_size > 0 else "SKIP_EMPTY"
    else:
        rc = run_flank_blast_once(
            Path(task["subject_fa"]),
            Path(task["query_fa"]),
            out_path,
            log_path,
            args,
        )
        stats = read_blast_table_stats(out_path, args.flank_hit_evalue_max)
        status = "OK" if rc == 0 else f"WARN:rc={rc}"

    return {
        "pair_dir": task["pair_dir"],
        "comparison": task["comparison"],
        "status": status,
        "blast_path": str(out_path),
        "log_path": str(log_path),
        "is_close_gap": bool(task["is_close_gap"]),
        "is_adjacent_neighbors": bool(task["is_adjacent_neighbors"]),
        **stats,
    }


def safe_neg_log10(value: float, floor: float):
    if pd.isna(value):
        return math.nan
    return -math.log10(max(float(value), float(floor)))


def summarize_flank_blast_pairs(comp_df: pd.DataFrame, zero_evalue_floor: float):
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
            best_double_evalue = math.nan
            best_double_bitscore = math.nan
            best_double_comparison = ""
            selected_evalues: List[float] = []
            selected_bitscores: List[float] = []
            selected_neglog_evalues: List[float] = []
        else:
            positive_rows = positive_rows.sort_values(
                by=["best_evalue", "best_bitscore", "best_length"],
                ascending=[True, False, False],
                na_position="last",
                kind="mergesort",
            )
            best = positive_rows.iloc[0]
            best_double_evalue = float(best["best_evalue"])
            best_double_bitscore = float(best["best_bitscore"])
            best_double_comparison = str(best["comparison"])
            selected_evalues = [
                float(value)
                for value in positive_rows["best_evalue"].dropna().tolist()
            ]
            selected_bitscores = [
                float(value)
                for value in positive_rows["best_bitscore"].dropna().tolist()
            ]
            selected_neglog_evalues = [
                safe_neg_log10(value, zero_evalue_floor)
                for value in selected_evalues
            ]

        rows.append({
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
                best_double_evalue,
                zero_evalue_floor,
            ),
            "double_hit_comparisons": ",".join(selected),
            "double_hit_evalue_count": int(len(selected_evalues)),
            "double_hit_evalues": ",".join(f"{value:.16g}" for value in selected_evalues),
            "double_hit_bitscores": ",".join(f"{value:.16g}" for value in selected_bitscores),
            "double_hit_neg_log10_evalues": ",".join(
                f"{value:.16g}" for value in selected_neglog_evalues
            ),
        })
    return pd.DataFrame(rows)


def build_flank_cutoff_table(pair_df: pd.DataFrame, percentiles: Sequence[float]):
    positive = pair_df[pair_df["any_double_hit"] == 1].copy()
    if positive.empty:
        return pd.DataFrame(columns=[
            "metric",
            "percentile",
            "cutoff_value",
            "stronger_values",
            "n_positive_double_hits",
            "n_double_hit_distribution_values",
            "notes",
        ])

    rows = []
    n_positive = int(len(positive))
    neglog_values = []
    bitscore_values = []
    evalue_values = []
    for row in positive.itertuples(index=False):
        neglog_values.extend(
            float(raw)
            for raw in str(getattr(row, "double_hit_neg_log10_evalues", "")).split(",")
            if raw
        )
        bitscore_values.extend(
            float(raw)
            for raw in str(getattr(row, "double_hit_bitscores", "")).split(",")
            if raw
        )
        evalue_values.extend(
            float(raw)
            for raw in str(getattr(row, "double_hit_evalues", "")).split(",")
            if raw
        )
    neglog_values = np.asarray(neglog_values, dtype=float)
    bitscore_values = np.asarray(bitscore_values, dtype=float)
    evalue_values = np.asarray(evalue_values, dtype=float)
    n_distribution_values = int(len(evalue_values))

    for pct in percentiles:
        q = pct / 100.0
        if len(neglog_values):
            cutoff = float(np.quantile(neglog_values, q))
            rows.append({
                "metric": "neg_log10_evalue",
                "percentile": pct,
                "cutoff_value": cutoff,
                "stronger_values": "higher",
                "n_positive_double_hits": n_positive,
                "n_double_hit_distribution_values": n_distribution_values,
                "notes": f"Equivalent evalue ~= {10 ** (-cutoff):.6g}",
            })
        if len(bitscore_values):
            cutoff = float(np.quantile(bitscore_values, q))
            rows.append({
                "metric": "bitscore",
                "percentile": pct,
                "cutoff_value": cutoff,
                "stronger_values": "higher",
                "n_positive_double_hits": n_positive,
                "n_double_hit_distribution_values": n_distribution_values,
                "notes": "",
            })
        if len(evalue_values):
            cutoff = float(np.quantile(evalue_values, 1.0 - q))
            rows.append({
                "metric": "evalue",
                "percentile": pct,
                "cutoff_value": cutoff,
                "stronger_values": "lower",
                "n_positive_double_hits": n_positive,
                "n_double_hit_distribution_values": n_distribution_values,
                "notes": (
                    "Built from all flank evalues participating in the passing double-hit "
                    "comparisons for positive pairs."
                ),
            })
    return pd.DataFrame(rows)


def run_flank_null_blast_mode(args):
    flanks_root, manifest_path, outdir = resolve_flank_blast_inputs(args)
    if not flanks_root.exists():
        raise SystemExit(f"[ERROR] flanks root does not exist: {flanks_root}")
    if not manifest_path.exists():
        raise SystemExit(f"[ERROR] manifest does not exist: {manifest_path}")
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(manifest_path, sep="\t")
    required = {"gene_label", "pair_dir", "enhancer_id", "fasta_path"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"[ERROR] manifest missing columns: {sorted(missing)}")

    set_order = build_manifest_set_order(df)
    pair_groups = list(df.groupby("pair_dir", sort=True))
    if args.pair_limit and args.pair_limit > 0:
        pair_groups = pair_groups[: int(args.pair_limit)]

    tasks = []
    for pair_dir, group in pair_groups:
        parsed = pair_dir_parts(pair_dir)
        if not parsed:
            print(f"[WARN] Cannot parse pair_dir: {pair_dir}", flush=True)
            continue
        set_label, a, b = parsed

        dpaths = {}
        for _, row in group.iterrows():
            enh_id = str(row.get("enhancer_id") or row.get("flank_enhancer_id") or "")
            side = str(row.get("flank_side") or row.get("side") or "").upper()
            fasta_path = flanks_root / str(row["fasta_path"])
            if enh_id and side in {"LEFT", "RIGHT"} and fasta_path.exists():
                dpaths[(enh_id, side)] = fasta_path

        if not dpaths:
            print(f"[INFO] No flanks found for {pair_dir}", flush=True)
            continue

        planned_tags, is_adjacent, is_close = plan_comparisons(a, b, set_order, set_label, args)
        a_id = f"{a[0]}:{a[1]}-{a[2]}"
        b_id = f"{b[0]}:{b[1]}-{b[2]}"
        flank_lookup = {
            "LL": (dpaths.get((a_id, "LEFT")), dpaths.get((b_id, "LEFT"))),
            "RR": (dpaths.get((a_id, "RIGHT")), dpaths.get((b_id, "RIGHT"))),
            "LR": (dpaths.get((a_id, "LEFT")), dpaths.get((b_id, "RIGHT"))),
            "RL": (dpaths.get((a_id, "RIGHT")), dpaths.get((b_id, "LEFT"))),
        }

        for tag in planned_tags:
            subject_fa, query_fa = flank_lookup[tag]
            if not subject_fa or not query_fa:
                continue
            blast_dir = outdir / pair_dir / "blast" / tag
            tasks.append({
                "pair_dir": pair_dir,
                "comparison": tag,
                "subject_fa": str(subject_fa),
                "query_fa": str(query_fa),
                "blast_path": str(blast_dir / "blast_hits.tsv"),
                "log_path": str(blast_dir / "blast.log"),
                "is_close_gap": is_close,
                "is_adjacent_neighbors": is_adjacent,
            })

    if not tasks:
        print("[INFO] Planned comparisons: 0", flush=True)
        return

    print(
        f"[INFO] flank-null-blast: {len(pair_groups)} pair groups, {len(tasks)} comparisons, {args.threads} threads",
        flush=True,
    )
    results = []
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = [pool.submit(flank_blast_worker, task, args) for task in tasks]
        for idx, fut in enumerate(as_completed(futures), start=1):
            results.append(fut.result())
            if idx == len(futures) or idx % max(100, len(futures) // 20) == 0:
                print(f"[INFO] flank-null-blast: processed {idx}/{len(futures)} comparisons", flush=True)

    comp_df = pd.DataFrame(results).sort_values(by=["pair_dir", "comparison"], kind="mergesort")
    pair_df = summarize_flank_blast_pairs(comp_df, args.zero_evalue_floor)
    cutoff_df = build_flank_cutoff_table(pair_df, parse_percentiles(args.cutoff_percentiles))

    comp_path = outdir / "blast_comparisons.tsv"
    pair_path = outdir / "blast_pairs_summary.tsv"
    cutoff_path = outdir / "blast_null_cutoffs.tsv"
    comp_df.to_csv(comp_path, sep="\t", index=False)
    pair_df.to_csv(pair_path, sep="\t", index=False)
    cutoff_df.to_csv(cutoff_path, sep="\t", index=False)

    positive = int((pair_df["any_double_hit"] == 1).sum()) if not pair_df.empty else 0
    print(f"[DONE] Wrote per-comparison BLAST summary: {comp_path}", flush=True)
    print(f"[DONE] Wrote pair-level double-hit summary: {pair_path}", flush=True)
    print(f"[DONE] Wrote null cutoff table: {cutoff_path}", flush=True)
    print(f"[DONE] Positive double-hit pairs: {positive}", flush=True)


def run_validate_flank_manifest_mode(args):
    flanks_root, manifest_path, outdir = resolve_flank_blast_inputs(args)
    if not flanks_root.exists():
        raise SystemExit(f"[ERROR] flanks root does not exist: {flanks_root}")
    if not manifest_path.exists():
        raise SystemExit(f"[ERROR] manifest does not exist: {manifest_path}")
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(manifest_path, sep="\t")
    required = {"gene_label", "pair_dir", "enhancer_id", "fasta_path"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"[ERROR] manifest missing columns: {sorted(missing)}")

    side_col = "flank_side" if "flank_side" in df.columns else ("side" if "side" in df.columns else None)
    rows = []
    for idx, row in df.iterrows():
        rel_path = str(row["fasta_path"])
        abs_path = flanks_root / rel_path
        parent_dir = abs_path.parent
        side = str(row.get(side_col) or "").upper() if side_col else ""
        rows.append({
            "row_number": int(idx) + 2,  # header is line 1 in the TSV
            "gene_label": str(row.get("gene_label", "")),
            "pair_dir": str(row.get("pair_dir", "")),
            "enhancer_id": str(row.get("enhancer_id", "")),
            "flank_side": side,
            "fasta_path": rel_path,
            "resolved_path": str(abs_path),
            "parent_dir": str(parent_dir),
            "parent_dir_exists": int(parent_dir.is_dir()),
            "file_exists": int(abs_path.is_file()),
            "path_is_absolute_in_manifest": int(Path(rel_path).is_absolute()),
        })

    path_df = pd.DataFrame(rows)
    path_df.to_csv(outdir / "manifest_path_validation.tsv", sep="\t", index=False)

    summary_rows = [{
        "manifest_path": str(manifest_path),
        "flanks_root": str(flanks_root),
        "rows_total": int(len(path_df)),
        "rows_with_existing_parent_dir": int(path_df["parent_dir_exists"].sum()) if not path_df.empty else 0,
        "rows_with_existing_file": int(path_df["file_exists"].sum()) if not path_df.empty else 0,
        "rows_missing_parent_dir": int((path_df["parent_dir_exists"] == 0).sum()) if not path_df.empty else 0,
        "rows_missing_file": int((path_df["file_exists"] == 0).sum()) if not path_df.empty else 0,
        "rows_with_absolute_manifest_paths": int(path_df["path_is_absolute_in_manifest"].sum()) if not path_df.empty else 0,
        "unique_pair_dirs": int(df["pair_dir"].nunique()) if "pair_dir" in df.columns else 0,
    }]
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(outdir / "manifest_path_validation_summary.tsv", sep="\t", index=False)

    if side_col:
        side_presence_rows = []
        for pair_dir, group in df.groupby("pair_dir", sort=True):
            per_enhancer = defaultdict(set)
            for _, row in group.iterrows():
                enh_id = str(row.get("enhancer_id") or "")
                side = str(row.get(side_col) or "").upper()
                rel_path = str(row["fasta_path"])
                abs_path = flanks_root / rel_path
                if enh_id and side in {"LEFT", "RIGHT"} and abs_path.is_file():
                    per_enhancer[enh_id].add(side)
            for enh_id, present_sides in sorted(per_enhancer.items()):
                side_presence_rows.append({
                    "pair_dir": str(pair_dir),
                    "enhancer_id": enh_id,
                    "left_exists": int("LEFT" in present_sides),
                    "right_exists": int("RIGHT" in present_sides),
                    "both_sides_exist": int({"LEFT", "RIGHT"}.issubset(present_sides)),
                })
        side_df = pd.DataFrame(side_presence_rows)
        side_df.to_csv(outdir / "manifest_flank_side_validation.tsv", sep="\t", index=False)
        complete = int(side_df["both_sides_exist"].sum()) if not side_df.empty else 0
        total = int(len(side_df))
        print(f"[INFO] enhancers with both LEFT and RIGHT flank files: {complete}/{total}", flush=True)

    missing_file_count = int((path_df["file_exists"] == 0).sum()) if not path_df.empty else 0
    missing_dir_count = int((path_df["parent_dir_exists"] == 0).sum()) if not path_df.empty else 0
    print(f"[DONE] Wrote row-level path validation: {outdir / 'manifest_path_validation.tsv'}", flush=True)
    print(f"[DONE] Wrote summary: {outdir / 'manifest_path_validation_summary.tsv'}", flush=True)
    if side_col:
        print(f"[DONE] Wrote flank-side validation: {outdir / 'manifest_flank_side_validation.tsv'}", flush=True)
    print(f"[DONE] Missing parent directories: {missing_dir_count}", flush=True)
    print(f"[DONE] Missing files: {missing_file_count}", flush=True)


def any_flank_overlap(a: Tuple[str, str], b: Tuple[str, str], flank_coords) -> bool:
    ca, _, _ = parse_enhancer_id(a[1])
    cb, _, _ = parse_enhancer_id(b[1])
    if ca != cb:
        return False

    a_left, a_right = flank_coords[a]["LEFT"], flank_coords[a]["RIGHT"]
    b_left, b_right = flank_coords[b]["LEFT"], flank_coords[b]["RIGHT"]
    return (
        intervals_overlap(a_left[1], a_left[2], b_left[1], b_left[2]) or
        intervals_overlap(a_right[1], a_right[2], b_right[1], b_right[2]) or
        intervals_overlap(a_left[1], a_left[2], b_right[1], b_right[2]) or
        intervals_overlap(a_right[1], a_right[2], b_left[1], b_left[2])
    )


def blast_has_any_hit(seq_a: str, seq_b: str, args) -> bool:
    with tempfile.TemporaryDirectory(prefix="blast_pair_") as td:
        q = Path(td) / "q.fa"
        s = Path(td) / "s.fa"
        q.write_text(">q\n" + seq_a + "\n")
        s.write_text(">s\n" + seq_b + "\n")

        cmd = [
            args.blastn,
            "-query", str(q),
            "-subject", str(s),
            "-evalue", str(args.blast_evalue),
            "-word_size", str(args.blast_word_size),
            "-outfmt", "6",
        ]
        proc = sp.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            timeout=args.blast_timeout,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"BLAST failed for {q} vs {s}\n"
                f"cmd: {' '.join(cmd)}\n"
                f"stderr: {proc.stderr.strip()}"
            )
        return bool(proc.stdout.strip())


def identity_from_text(target_text: str, query_text: str) -> float:
    matches = 0
    aligned = 0
    for a, b in zip(target_text, query_text):
        if a == "-" or b == "-":
            continue
        aligned += 1
        if a.upper() == b.upper():
            matches += 1
    return (matches / aligned * 100.0) if aligned else 0.0


def maf_stats(maf_path: Path):
    if not maf_path.exists() or maf_path.stat().st_size == 0:
        return {
            "blocks": 0,
            "total_aligned_len": 0,
            "best_block_len": 0,
            "best_block_identity_pct": 0.0,
            "combined_identity_pct": 0.0,
        }

    blocks = 0
    total_len = 0
    best_len = 0
    best_ident = 0.0
    sum_matches = 0
    sum_aligned_nt = 0
    current = []

    with maf_path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                if len(current) >= 2:
                    target, query = current[:2]
                    blocks += 1
                    aln_len = min(target["size"], query["size"])
                    total_len += aln_len
                    ident = identity_from_text(target.get("text", ""), query.get("text", ""))
                    best_len = max(best_len, aln_len)
                    best_ident = max(best_ident, ident)

                    matches = 0
                    aligned = 0
                    for a, b in zip(target.get("text", ""), query.get("text", "")):
                        if a == "-" or b == "-":
                            continue
                        aligned += 1
                        if a.upper() == b.upper():
                            matches += 1
                    sum_matches += matches
                    sum_aligned_nt += aligned
                current = []
                continue

            if line.startswith("#") or line.startswith("a "):
                if line.startswith("a "):
                    current = []
                continue
            if line.startswith("s "):
                parts = line.split()
                if len(parts) >= 7:
                    _, src, start, size, strand, src_size, text = parts[:7]
                    current.append({
                        "src": src,
                        "start": int(start),
                        "size": int(size),
                        "strand": strand,
                        "srcSize": int(src_size),
                        "text": text,
                    })

        if len(current) >= 2:
            target, query = current[:2]
            blocks += 1
            aln_len = min(target["size"], query["size"])
            total_len += aln_len
            ident = identity_from_text(target.get("text", ""), query.get("text", ""))
            best_len = max(best_len, aln_len)
            best_ident = max(best_ident, ident)
            matches = 0
            aligned = 0
            for a, b in zip(target.get("text", ""), query.get("text", "")):
                if a == "-" or b == "-":
                    continue
                aligned += 1
                if a.upper() == b.upper():
                    matches += 1
            sum_matches += matches
            sum_aligned_nt += aligned

    combined_ident = (sum_matches / sum_aligned_nt * 100.0) if sum_aligned_nt else 0.0
    return {
        "blocks": blocks,
        "total_aligned_len": total_len,
        "best_block_len": best_len,
        "best_block_identity_pct": best_ident,
        "combined_identity_pct": combined_ident,
    }


def run_lastz_once(target: Path, query: Path, args) -> dict:
    with tempfile.TemporaryDirectory(prefix="lastz_pair_") as td:
        maf_path = Path(td) / "out.maf"
        cmd = [args.lastz, str(target), str(query)] + shlex.split(args.lastz_args)
        with maf_path.open("w") as mafh:
            proc = sp.run(
                cmd,
                stdout=mafh,
                stderr=sp.PIPE,
                text=True,
                check=False,
                timeout=args.lastz_timeout,
            )
        if proc.returncode != 0:
            raise RuntimeError(
                f"LASTZ failed for {target} vs {query}\n"
                f"cmd: {' '.join(cmd)}\n"
                f"stderr: {proc.stderr.strip()}"
            )
        return maf_stats(maf_path)


def gap_bp(a: Tuple[str, int, int], b: Tuple[str, int, int]) -> int:
    if a[0] != b[0]:
        return 10**12
    return max(0, b[1] - a[2])


def build_set_order(set_to_enhancers: Dict[str, Sequence[str]], enhancer_records: Dict[str, EnhancerRecord]):
    ordered: Dict[str, Dict[str, List[Tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    for set_label, enhancers in set_to_enhancers.items():
        for member_id in enhancers:
            chrom, start, end = parse_enhancer_id(enhancer_records[member_id].enhancer_id)
            ordered[set_label][chrom].append((start, end))
        for chrom in ordered[set_label]:
            ordered[set_label][chrom].sort()
    return ordered


def are_adjacent_neighbors(
    set_order: Dict[str, Dict[str, List[Tuple[int, int]]]],
    set_label: str,
    a: Tuple[str, int, int],
    b: Tuple[str, int, int],
) -> bool:
    if a[0] != b[0]:
        return False
    ordered = set_order.get(set_label, {}).get(a[0], [])
    if not ordered:
        return False
    try:
        ia = ordered.index((a[1], a[2]))
        ib = ordered.index((b[1], b[2]))
    except ValueError:
        return False
    return abs(ia - ib) == 1


def plan_comparisons(
    a: Tuple[str, int, int],
    b: Tuple[str, int, int],
    set_order,
    set_label: str,
    args,
):
    is_adjacent = are_adjacent_neighbors(set_order, set_label, a, b)
    is_close = gap_bp(a, b) < args.gap_threshold
    tags = ["LL", "RR", "LR"]
    if (not is_close) or (not is_adjacent):
        tags.append("RL")
    return tags, is_adjacent, is_close


class PairAnalyzer:
    def __init__(self, args, enhancer_records, enh_to_flanks, flank_coords, source_flanks_root):
        self.args = args
        self.enhancer_records = enhancer_records
        self.enh_to_flanks = enh_to_flanks
        self.flank_coords = flank_coords
        self.source_flanks_root = source_flanks_root
        self.seq_cache: Dict[str, str] = {}
        self.blast_cache: Dict[Tuple[str, str], bool] = {}
        self.lastz_cache: Dict[Tuple[Tuple[str, str], str], dict] = {}
        self.seq_lock = Lock()
        self.blast_lock = Lock()
        self.lastz_lock = Lock()

    def sequence_for(self, member_id: str) -> str:
        with self.seq_lock:
            if member_id in self.seq_cache:
                return self.seq_cache[member_id]
        _, seq = read_first_fasta_record(self.enhancer_records[member_id].seq_path)
        with self.seq_lock:
            self.seq_cache[member_id] = seq
        return seq

    def blast_pair(self, a: str, b: str) -> bool:
        pair_key = member_pair_ids(a, b)
        with self.blast_lock:
            if pair_key in self.blast_cache:
                return self.blast_cache[pair_key]
        hit = blast_has_any_hit(self.sequence_for(a), self.sequence_for(b), self.args)
        with self.blast_lock:
            self.blast_cache[pair_key] = hit
        return hit

    def lastz_stats(self, a: str, b: str, tag: str) -> dict:
        pair_key = member_pair_ids(a, b)
        cache_key = (pair_key, tag)
        with self.lastz_lock:
            if cache_key in self.lastz_cache:
                return self.lastz_cache[cache_key]

        left, right = pair_key
        flank_a = self.enh_to_flanks[member_manifest_key(left, self.enhancer_records)]
        flank_b = self.enh_to_flanks[member_manifest_key(right, self.enhancer_records)]

        if tag == "LL":
            target = self.source_flanks_root / flank_a["LEFT"]
            query = self.source_flanks_root / flank_b["LEFT"]
        elif tag == "RR":
            target = self.source_flanks_root / flank_a["RIGHT"]
            query = self.source_flanks_root / flank_b["RIGHT"]
        elif tag == "LR":
            target = self.source_flanks_root / flank_a["LEFT"]
            query = self.source_flanks_root / flank_b["RIGHT"]
        elif tag == "RL":
            target = self.source_flanks_root / flank_a["RIGHT"]
            query = self.source_flanks_root / flank_b["LEFT"]
        else:
            raise ValueError(tag)

        stats = run_lastz_once(target, query, self.args)
        with self.lastz_lock:
            self.lastz_cache[cache_key] = stats
        return stats


def choose_null_set(
    focal_enhancers: Sequence[str],
    valid_enhancers: Sequence[str],
    enh_to_source_set: Dict[str, str],
    rng: random.Random,
    used_globally: Set[str],
    args,
) -> Optional[List[str]]:
    n = len(focal_enhancers)
    forbid = set(focal_enhancers)
    forbid_biological_ids = {
        member_biological_id(member_id, args._enhancer_records)
        for member_id in focal_enhancers
    }
    pool = list(valid_enhancers)

    for _ in range(args.max_sample_attempts_per_set):
        rng.shuffle(pool)
        picked = []
        used_source_sets = set()
        used_biological_ids = set(forbid_biological_ids)

        for member_id in pool:
            if member_id in forbid:
                continue
            if args.without_replacement_across_sets and member_id in used_globally:
                continue
            biological_id = member_biological_id(member_id, args._enhancer_records)
            if biological_id in used_biological_ids:
                continue
            source_set = enh_to_source_set[member_id]
            if (not args.allow_same_source_set_within_null_set) and source_set in used_source_sets:
                continue
            picked.append(member_id)
            used_source_sets.add(source_set)
            used_biological_ids.add(biological_id)
            if len(picked) == n:
                return picked
    return None


def sample_null_sets(
    observed_sets: Dict[str, Sequence[str]],
    valid_enhancers: Sequence[str],
    enh_to_source_set: Dict[str, str],
    rng: random.Random,
    args,
):
    sampled = {}
    used_globally: Set[str] = set()
    ordered_sets = sorted(observed_sets.items(), key=lambda kv: (-len(kv[1]), kv[0]))

    for set_label, enhancers in ordered_sets:
        picked = choose_null_set(
            enhancers,
            valid_enhancers,
            enh_to_source_set,
            rng,
            used_globally,
            args,
        )
        if picked is None:
            raise RuntimeError(
                f"Could not sample a null set of size {len(enhancers)} for {set_label}. "
                f"Try relaxing sampling constraints."
            )
        sampled[set_label] = sorted(picked, key=lambda x: member_sort_key(x, args._enhancer_records))
        if args.without_replacement_across_sets:
            used_globally.update(picked)

    return sampled


def analyze_pair_task(
    task,
    pair_analyzer: PairAnalyzer,
    set_order,
    observed_pairs: Set[Tuple[str, str]],
    flank_coords,
    args,
):
    set_label, bucket, replicate, a, b = task
    pair_key = member_pair_ids(a, b)
    genomic_pair = genomic_pair_for_members(a, b, pair_analyzer.enhancer_records)
    pair_name = pair_name_for_members(a, b, pair_analyzer.enhancer_records)
    a_rec = pair_analyzer.enhancer_records[a]
    b_rec = pair_analyzer.enhancer_records[b]
    same_biological_enhancer = (a_rec.enhancer_file == b_rec.enhancer_file)

    if same_biological_enhancer:
        return {
            "set_label": set_label,
            "bucket": bucket,
            "replicate": replicate,
            "enhancer_a_member_id": pair_key[0],
            "enhancer_b_member_id": pair_key[1],
            "enhancer_a": a_rec.enhancer_id,
            "enhancer_b": b_rec.enhancer_id,
            "pair_name": pair_name,
            "analyzed": 0,
            "skip_reason": "same_enhancer_file",
            "blast_hit": 0,
            "qualifying_comparisons": "",
            "pair_is_duplication_hit": 0,
            "touched_member_ids": "",
        }

    if args.exclude_observed_pairs_globally and replicate is not None and genomic_pair in observed_pairs:
        return {
            "set_label": set_label,
            "bucket": bucket,
            "replicate": replicate,
            "enhancer_a_member_id": pair_key[0],
            "enhancer_b_member_id": pair_key[1],
            "enhancer_a": a_rec.enhancer_id,
            "enhancer_b": b_rec.enhancer_id,
            "pair_name": pair_name,
            "analyzed": 0,
            "skip_reason": "observed_pair",
            "blast_hit": 0,
            "qualifying_comparisons": "",
            "pair_is_duplication_hit": 0,
            "touched_member_ids": "",
        }

    if args.skip_overlapping_flanks and any_flank_overlap(
        member_manifest_key(pair_key[0], pair_analyzer.enhancer_records),
        member_manifest_key(pair_key[1], pair_analyzer.enhancer_records),
        flank_coords,
    ):
        return {
            "set_label": set_label,
            "bucket": bucket,
            "replicate": replicate,
            "enhancer_a_member_id": pair_key[0],
            "enhancer_b_member_id": pair_key[1],
            "enhancer_a": a_rec.enhancer_id,
            "enhancer_b": b_rec.enhancer_id,
            "pair_name": pair_name,
            "analyzed": 0,
            "skip_reason": "overlapping_flanks",
            "blast_hit": 0,
            "qualifying_comparisons": "",
            "pair_is_duplication_hit": 0,
            "touched_member_ids": "",
        }

    blast_hit = pair_analyzer.blast_pair(pair_key[0], pair_key[1])
    row = {
        "set_label": set_label,
        "bucket": bucket,
        "replicate": replicate,
        "enhancer_a_member_id": pair_key[0],
        "enhancer_b_member_id": pair_key[1],
        "enhancer_a": a_rec.enhancer_id,
        "enhancer_b": b_rec.enhancer_id,
        "enhancer_a_file": a_rec.enhancer_file,
        "enhancer_b_file": b_rec.enhancer_file,
        "pair_name": pair_name,
        "analyzed": 1,
        "skip_reason": "",
        "blast_hit": int(blast_hit),
        "qualifying_comparisons": "",
        "pair_is_duplication_hit": 0,
        "touched_member_ids": "",
    }

    a_tuple = parse_enhancer_id(a_rec.enhancer_id)
    b_tuple = parse_enhancer_id(b_rec.enhancer_id)
    planned_tags, is_adjacent, is_close = plan_comparisons(a_tuple, b_tuple, set_order, set_label, args)

    row["gap_bp"] = gap_bp(a_tuple, b_tuple)
    row["is_adjacent_neighbors"] = int(is_adjacent)
    row["is_close_pair"] = int(is_close)
    row["planned_comparisons"] = ",".join(planned_tags)

    qualifying_tags = []
    for tag in planned_tags:
        stats = pair_analyzer.lastz_stats(pair_key[0], pair_key[1], tag)
        row[f"{tag}_blocks"] = stats["blocks"]
        row[f"{tag}_best_block_len"] = stats["best_block_len"]
        row[f"{tag}_best_block_identity_pct"] = stats["best_block_identity_pct"]
        if (
            stats["blocks"] > 0 and
            stats["best_block_len"] > args.len_threshold and
            stats["best_block_identity_pct"] > args.identity_threshold
        ):
            qualifying_tags.append(tag)

    qualifies = (
        {"LL", "RR"}.issubset(set(qualifying_tags)) or
        {"LR", "RL"}.issubset(set(qualifying_tags))
    )
    row["qualifying_comparisons"] = ",".join(sorted(qualifying_tags))
    row["pair_is_duplication_hit"] = int(blast_hit and qualifies)
    if blast_hit and qualifies:
        row["touched_member_ids"] = ",".join(pair_key)
    return row


def analyze_sets(
    sets_by_label: Dict[str, Sequence[str]],
    pair_analyzer: PairAnalyzer,
    observed_pairs: Set[Tuple[str, str]],
    flank_coords,
    args,
    replicate: Optional[int],
    phase_label: str = "analysis",
):
    set_order = build_set_order(sets_by_label, pair_analyzer.enhancer_records)
    tasks = []
    sampled_rows = []

    for set_label, enhancers in sorted(sets_by_label.items()):
        enhancers_sorted = sorted(enhancers, key=lambda x: member_sort_key(x, pair_analyzer.enhancer_records))
        bucket = bucket_for_size(len(enhancers_sorted))
        for idx, member_id in enumerate(enhancers_sorted, start=1):
            rec = pair_analyzer.enhancer_records[member_id]
            sampled_rows.append({
                "replicate": replicate,
                "set_label": set_label,
                "bucket": bucket,
                "rank_within_set": idx,
                "member_id": member_id,
                "enhancer_id": rec.enhancer_id,
                "enhancer_file": rec.enhancer_file,
                "source_set_label": rec.set_label,
                "sequence_path": str(rec.seq_path),
            })
        for a, b in itertools.combinations(enhancers_sorted, 2):
            tasks.append((set_label, bucket, replicate, a, b))

    total_tasks = len(tasks)
    print(
        f"[INFO] {phase_label}: {len(sets_by_label)} sets, {total_tasks} pairs, {args.threads} threads",
        flush=True,
    )

    pair_rows = []
    if total_tasks:
        progress_every = max(100, total_tasks // 20)
        processed = 0
        last_report = time.time()
        with ThreadPoolExecutor(max_workers=args.threads) as pool:
            futures = [
                pool.submit(
                    analyze_pair_task,
                    task,
                    pair_analyzer,
                    set_order,
                    observed_pairs,
                    flank_coords,
                    args,
                )
                for task in tasks
            ]
            for fut in as_completed(futures):
                pair_rows.append(fut.result())
                processed += 1
                now = time.time()
                if processed == total_tasks or processed % progress_every == 0 or (now - last_report) >= 30:
                    print(
                        f"[INFO] {phase_label}: processed {processed}/{total_tasks} pairs",
                        flush=True,
                    )
                    last_report = now

    pair_df = pd.DataFrame(pair_rows)
    sampled_df = pd.DataFrame(sampled_rows)
    summary_rows = []
    for bucket in ["2", "3", ">=4"]:
        bucket_sets = {
            label: enhs for label, enhs in sets_by_label.items()
            if bucket_for_size(len(enhs)) == bucket
        }
        total_enhancers = sum(len(enhs) for enhs in bucket_sets.values())
        hit_enhancers: Set[str] = set()
        if not pair_df.empty:
            sub = pair_df[
                (pair_df["bucket"] == bucket) &
                (pair_df["pair_is_duplication_hit"] == 1)
            ]
            for touched in sub["touched_member_ids"].dropna():
                if not touched:
                    continue
                hit_enhancers.update(x for x in str(touched).split(",") if x)
            analyzed_pairs = int(pair_df[(pair_df["bucket"] == bucket) & (pair_df["analyzed"] == 1)].shape[0])
            blast_positive_pairs = int(pair_df[(pair_df["bucket"] == bucket) & (pair_df["blast_hit"] == 1)].shape[0])
            qualifying_pairs = int(sub.shape[0])
        else:
            analyzed_pairs = 0
            blast_positive_pairs = 0
            qualifying_pairs = 0

        summary_rows.append({
            "replicate": replicate,
            "bucket": bucket,
            "n_sets": len(bucket_sets),
            "total_enhancers": total_enhancers,
            "unique_hit_enhancers": len(hit_enhancers),
            "enhancer_hit_fraction": (
                len(hit_enhancers) / total_enhancers if total_enhancers else math.nan
            ),
            "analyzed_pairs": analyzed_pairs,
            "blast_positive_pairs": blast_positive_pairs,
            "qualifying_pairs": qualifying_pairs,
        })

    return sampled_df, pair_df, pd.DataFrame(summary_rows)


def summarize_null_vs_observed(null_df: pd.DataFrame, observed_df: Optional[pd.DataFrame]):
    observed_lookup = {}
    if observed_df is not None and not observed_df.empty:
        observed_lookup = (
            observed_df.set_index("bucket")["enhancer_hit_fraction"].to_dict()
        )

    rows = []
    for bucket in ["2", "3", ">=4"]:
        sub = null_df[null_df["bucket"] == bucket]["enhancer_hit_fraction"].dropna()
        obs = observed_lookup.get(bucket, math.nan)
        mu = float(sub.mean()) if len(sub) else math.nan
        sd = float(sub.std(ddof=1)) if len(sub) > 1 else math.nan
        z = ((obs - mu) / sd) if (pd.notna(obs) and pd.notna(sd) and sd > 0) else math.nan
        emp_p_ge = (
            float((sub >= obs).mean())
            if len(sub) and pd.notna(obs)
            else math.nan
        )
        rows.append({
            "bucket": bucket,
            "observed_enhancer_hit_fraction": obs,
            "null_mean_enhancer_hit_fraction": mu,
            "null_sd_enhancer_hit_fraction": sd,
            "z_score": z,
            "empirical_p_ge_observed": emp_p_ge,
            "n_reps": int(len(sub)),
        })
    return pd.DataFrame(rows)


def main():
    args = parse_args()

    if args.mode == "flank-null-blast":
        run_flank_null_blast_mode(args)
        return
    if args.mode == "validate-flank-manifest":
        run_validate_flank_manifest_mode(args)
        return

    enhancer_root, source_flanks_root, source_manifest, outdir = resolve_input_paths(args)
    if not enhancer_root.exists():
        raise SystemExit(f"[ERROR] enhancer root does not exist: {enhancer_root}")
    if not source_flanks_root.exists():
        raise SystemExit(f"[ERROR] source flanks root does not exist: {source_flanks_root}")
    if not source_manifest.exists():
        raise SystemExit(f"[ERROR] source manifest does not exist: {source_manifest}")

    if outdir.exists():
        if not args.overwrite:
            raise SystemExit(f"[ERROR] outdir exists. Use --overwrite to replace it: {outdir}")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(args.seed)

    enhancer_records, enhancer_sets, enhancer_catalog_df = read_enhancer_catalog(enhancer_root)
    args._enhancer_records = enhancer_records
    (
        enh_to_flanks,
        manifest_source_sets,
        _enh_to_coords,
        flank_coords,
        observed_pairs,
    ) = read_source_manifest(source_flanks_root, source_manifest)
    member_source_sets = {
        member_id: enhancer_records[member_id].set_label
        for member_id in enhancer_records
    }

    valid_enhancers = sorted(
        member_id for member_id in enhancer_records
        if member_manifest_key(member_id, enhancer_records) in enh_to_flanks
        and {"LEFT", "RIGHT"}.issubset(enh_to_flanks[member_manifest_key(member_id, enhancer_records)])
    )
    if len(valid_enhancers) < 2:
        raise SystemExit("[ERROR] Fewer than two enhancers have both sequence and flank files")

    observed_sets = {}
    set_rows = []
    for set_label, enhancers in sorted(enhancer_sets.items()):
        valid = sorted(
            [member_id for member_id in enhancers if member_manifest_key(member_id, enhancer_records) in enh_to_flanks],
            key=lambda x: member_sort_key(x, enhancer_records),
        )
        set_rows.append({
            "set_label": set_label,
            "total_enhancers_in_sequence_dir": len(enhancers),
            "valid_enhancers_with_flanks": len(valid),
            "bucket": bucket_for_size(len(valid)),
        })
        if len(valid) >= 2:
            observed_sets[set_label] = valid

    if not observed_sets:
        raise SystemExit("[ERROR] No observed sets with at least 2 enhancers and valid flanks were found")

    eligible_sets_df = pd.DataFrame(set_rows)
    eligible_sets_df.to_csv(outdir / "eligible_shadow_sets.tsv", sep="\t", index=False)
    enhancer_catalog_df.to_csv(outdir / "enhancer_catalog.tsv", sep="\t", index=False)

    pair_analyzer = PairAnalyzer(
        args=args,
        enhancer_records=enhancer_records,
        enh_to_flanks=enh_to_flanks,
        flank_coords=flank_coords,
        source_flanks_root=source_flanks_root,
    )

    observed_summary_df = None
    if args.compute_observed:
        print("[INFO] analyzing observed shadow sets")
        observed_sampled_df, observed_pair_df, observed_summary_df = analyze_sets(
            observed_sets,
            pair_analyzer,
            observed_pairs=observed_pairs,
            flank_coords=flank_coords,
            args=args,
            replicate=None,
            phase_label="observed",
        )
        observed_sampled_df.to_csv(outdir / "observed_sets_used.tsv", sep="\t", index=False)
        observed_pair_df.to_csv(outdir / "observed_pair_results.tsv", sep="\t", index=False)
        observed_summary_df.to_csv(outdir / "observed_bucket_summary.tsv", sep="\t", index=False)

    valid_source_pool = sorted(
        member_id for member_id in valid_enhancers
        if member_manifest_key(member_id, enhancer_records) in manifest_source_sets
    )
    if len(valid_source_pool) < 2:
        raise SystemExit("[ERROR] Null-sampling pool has fewer than two valid enhancers")

    null_sampled_frames = []
    null_summary_frames = []
    null_pair_frames = []

    for rep in range(1, args.reps + 1):
        sampled_null_sets = sample_null_sets(
            observed_sets,
            valid_source_pool,
            member_source_sets,
            rng,
            args,
        )
        sampled_df, pair_df, summary_df = analyze_sets(
            sampled_null_sets,
            pair_analyzer,
            observed_pairs=observed_pairs,
            flank_coords=flank_coords,
            args=args,
            replicate=rep,
            phase_label=f"null replicate {rep}",
        )

        null_sampled_frames.append(sampled_df)
        null_summary_frames.append(summary_df)
        if args.write_null_pair_results:
            null_pair_frames.append(pair_df)

        brief = summary_df[["bucket", "unique_hit_enhancers", "enhancer_hit_fraction"]]
        as_text = "; ".join(
            f"{row.bucket}: hits={row.unique_hit_enhancers}, frac={row.enhancer_hit_fraction:.4f}"
            for row in brief.itertuples(index=False)
        )
        print(f"[INFO] finished replicate {rep}/{args.reps}: {as_text}")

    null_sampled_df = pd.concat(null_sampled_frames, ignore_index=True)
    null_summary_df = pd.concat(null_summary_frames, ignore_index=True)
    null_sampled_df.to_csv(outdir / "null_sampled_sets.tsv", sep="\t", index=False)
    null_summary_df.to_csv(outdir / "null_replicate_bucket_summary.tsv", sep="\t", index=False)

    if args.write_null_pair_results:
        pd.concat(null_pair_frames, ignore_index=True).to_csv(
            outdir / "null_pair_results.tsv",
            sep="\t",
            index=False,
        )

    final_summary_df = summarize_null_vs_observed(null_summary_df, observed_summary_df)
    final_summary_df.to_csv(outdir / "null_distribution_summary.tsv", sep="\t", index=False)

    (outdir / "README.txt").write_text(
        "shadow_duplication_null_pipeline output\n"
        "Files:\n"
        "  eligible_shadow_sets.tsv          per observed set, how many enhancers had usable flank data\n"
        "  enhancer_catalog.tsv              enhancer sequence catalog built from FlilesOutput\n"
        "  observed_sets_used.tsv            observed sets analyzed when --compute-observed is used\n"
        "  observed_pair_results.tsv         pair-level observed BLAST/LASTZ results\n"
        "  observed_bucket_summary.tsv       observed per-bin enhancer hit fractions\n"
        "  null_sampled_sets.tsv             enhancers sampled in each null replicate\n"
        "  null_replicate_bucket_summary.tsv replicate-level null fractions by bin\n"
        "  null_pair_results.tsv             optional full null pair table\n"
        "  null_distribution_summary.tsv     mean/sd/z-score summary by bin\n"
        "\n"
        "A pair is counted as a duplication hit only if:\n"
        "  1. BLAST has at least one hit at the requested evalue threshold\n"
        "  2. LASTZ has flank comparisons passing best_block_len and best_block_identity_pct thresholds\n"
        "  3. The pair satisfies the tier rule: LL+RR or LR+RL\n"
    )

    print(f"[DONE] wrote pipeline output to {outdir}")


if __name__ == "__main__":
    main()
