#!/usr/bin/env python3
"""
make_shadow_flanks_rm.py

Builds flank FASTAs for pairwise "shadow enhancer" comparisons with:
  - RepeatMasker-based soft-masking (lowercase) preserved in output FASTAs
  - Flank restriction: flanks cannot include ANY other enhancers from the same set
  - Gap rule:
      * If same-chromosome pair's inner gap < --gap-threshold (default 5kb):
            BOTH inner flanks are the intervening region [up.end, down.start),
            with any other enhancers from the same set subtracted.
            The piece contiguous with each enhancer edge is kept (no overlapping flanks).
      * Else (>= threshold or different chromosomes):
            independent LEFT/RIGHT flanks up to --flank-bp (default 2500),
            with other enhancers subtracted and trimmed by chromosome bounds.

Soft-masking strategy
---------------------
Provide a **soft-masked genome FASTA** (lowercase = masked) via --masked-genome,
or let the script run RepeatMasker (-xsmall) on --genome-fasta (requires RepeatMasker).

LASTZ won’t seed in lowercase by default, so this preserves TE/low-complexity masking.

Outputs
-------
  out_root/
    <set_label>/<chr_start-end>__<chr_start-end>/
      chr_start-end_LEFT.fa
      chr_start-end_RIGHT.fa
      ...
  out_root/manifest.tsv   (one row per written flank FASTA)

Usage examples
--------------
# A) Use an already soft-masked dm6 (lowercase)
python make_shadow_flanks_rm.py \
  --shadow-bed 011925_all_shadowsets_DM6.bed \
  --masked-genome /path/to/dm6.softmasked.fa \
  --outdir shadow_pairs_flanks

# B) Let the script run RepeatMasker (requires RepeatMasker on PATH)
python make_shadow_flanks_rm.py \
  --shadow-bed 011925_all_shadowsets_DM6.bed \
  --genome-fasta /path/to/dm6.fasta \
  --outdir shadow_pairs_flanks \
  --flank-bp 2500 \
  --gap-threshold 5000 \
  --rm-run --rm-species drosophila --rm-threads 8
"""

import argparse
import csv
import sys
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict, Optional

# ---------- FASTA access ----------
try:
    from pyfaidx import Fasta
except ImportError:
    sys.stderr.write("[ERROR] pyfaidx is required. Install into your env with:\n"
                     "  conda install -c conda-forge pyfaidx\n"
                     "  # or\n"
                     "  python -m pip install -U pyfaidx\n")
    sys.exit(2)

@dataclass(frozen=True)
class Enh:
    chrom: str
    start: int
    end: int
    label: str

Interval = Tuple[int, int]  # [start, end)

def merge_intervals(ivls: List[Interval]) -> List[Interval]:
    if not ivls:
        return []
    ivls = sorted(ivls)
    out = [ivls[0]]
    for s,e in ivls[1:]:
        ls,le = out[-1]
        if s <= le:
            out[-1] = (ls, max(le, e))
        else:
            out.append((s,e))
    return out

def subtract_one(window: Interval, block: Interval) -> List[Interval]:
    a,b = window; x,y = block
    if y <= a or x >= b:
        return [window]
    pieces = []
    if x > a:
        pieces.append((a, x))
    if y < b:
        pieces.append((y, b))
    return pieces

def subtract_many(window: Interval, blocks: List[Interval]) -> List[Interval]:
    pieces = [window]
    for blk in blocks:
        next_pieces = []
        for p in pieces:
            next_pieces.extend(subtract_one(p, blk))
        pieces = next_pieces
        if not pieces:
            break
    return pieces

def pick_adjacent_piece(pieces: List[Interval], window: Interval, side: str) -> Optional[Interval]:
    """Choose the piece contiguous with the enhancer edge.
       For RIGHT: piece with start == window.start
       For LEFT:  piece with end   == window.end
    """
    if not pieces:
        return None
    a,b = window
    if side == "RIGHT":
        for s,e in sorted(pieces):
            if s == a:
                return (s,e)
        return None
    else:  # LEFT
        for s,e in sorted(pieces):
            if e == b:
                return (s,e)
        return None

def load_shadow_bed(bed_path: Path) -> List[Enh]:
    enhs: List[Enh] = []
    with bed_path.open() as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("#","track","browser")):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise ValueError("Shadow BED must have at least 4 columns: chrom start end label")
            chrom, s, e, label = parts[0], int(parts[1]), int(parts[2]), parts[3]
            enhs.append(Enh(chrom, s, e, label))
    return enhs

def group_by_label(enhs: List[Enh]) -> Dict[str, List[Enh]]:
    by = {}
    for e in enhs:
        by.setdefault(e.label, []).append(e)
    for k in by:
        by[k] = sorted(by[k], key=lambda x: (x.chrom, x.start, x.end))
    return by

def enh_id(e: Enh) -> str:
    return f"{e.chrom}:{e.start}-{e.end}"

def pair_dir_name(e1: Enh, e2: Enh) -> str:
    left = f"{e1.chrom}_{e1.start}-{e1.end}"
    right = f"{e2.chrom}_{e2.start}-{e2.end}"
    return f"{left}__{right}"

def build_union_of_other_enhancers(this_set: List[Enh], exclude: List[Enh]) -> Dict[str, List[Interval]]:
    ex_ids = { (x.chrom, x.start, x.end) for x in exclude }
    tmp: Dict[str, List[Interval]] = {}
    for e in this_set:
        if (e.chrom, e.start, e.end) in ex_ids:
            continue
        tmp.setdefault(e.chrom, []).append((e.start, e.end))
    return { c: merge_intervals(lst) for c,lst in tmp.items() }

def count_lowercase_bases(seq: str) -> int:
    return sum(1 for ch in seq if ch.islower())

def write_fa(path: Path, header: str, seq: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as out:
        out.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")

def make_outer_flank(e: Enh, side: str, flank_bp: int, chr_size: int,
                     others_union: Dict[str, List[Interval]]) -> Tuple[Optional[Interval], Dict]:
    meta = {"restricted_by_other_enhancers": False, "truncated_at_boundary": False}
    if side == "LEFT":
        raw = (max(0, e.start - flank_bp), e.start)
        if raw[0] == 0:
            meta["truncated_at_boundary"] = True
        blocks = others_union.get(e.chrom, [])
        allowed = subtract_many(raw, blocks)
        meta["restricted_by_other_enhancers"] = (allowed != [raw])
        chosen = pick_adjacent_piece(allowed, raw, "LEFT")
    else:  # RIGHT
        raw = (e.end, min(e.end + flank_bp, chr_size))
        blocks = others_union.get(e.chrom, [])
        allowed = subtract_many(raw, blocks)
        meta["restricted_by_other_enhancers"] = (allowed != [raw])
        chosen = pick_adjacent_piece(allowed, raw, "RIGHT")
    return chosen, meta

def make_inner_flanks_for_pair(up: Enh, down: Enh, flank_bp: int, gap_threshold: int,
                               up_chr_size: int, down_chr_size: int,
                               others_union_excl_pair: Dict[str, List[Interval]]
                               ) -> Tuple[Optional[Interval], Optional[Interval], Dict, Dict]:
    """Return up_RIGHT_interval, down_LEFT_interval (or None), with metadata."""
    metaR = {"is_intervening": False, "restricted_by_other_enhancers": False}
    metaL = {"is_intervening": False, "restricted_by_other_enhancers": False}

    gap = max(0, down.start - up.end) if up.chrom == down.chrom else None

    if (gap is not None) and (gap < gap_threshold):
        # Intervening region, subtract other enhancers, then keep pieces contiguous to edges
        raw = (up.end, down.start)
        blocks = others_union_excl_pair.get(up.chrom, [])
        allowed = subtract_many(raw, blocks)
        up_piece = pick_adjacent_piece(allowed, raw, "RIGHT")
        down_piece = pick_adjacent_piece(allowed, raw, "LEFT")
        metaR["is_intervening"] = True
        metaL["is_intervening"] = True
        metaR["restricted_by_other_enhancers"] = (allowed != [raw])
        metaL["restricted_by_other_enhancers"] = (allowed != [raw])
        return up_piece, down_piece, metaR, metaL
    else:
        # Independent inner flanks (pointing toward each other), restricted & clamped
        up_raw   = (up.end, min(up.end + flank_bp, up_chr_size))
        down_raw = (max(0, down.start - flank_bp), down.start)
        blocks_up   = others_union_excl_pair.get(up.chrom, [])
        blocks_down = others_union_excl_pair.get(down.chrom, [])
        up_allowed   = subtract_many(up_raw,   blocks_up)
        down_allowed = subtract_many(down_raw, blocks_down)
        up_piece   = pick_adjacent_piece(up_allowed,   up_raw,   "RIGHT")
        down_piece = pick_adjacent_piece(down_allowed, down_raw, "LEFT")
        metaR["restricted_by_other_enhancers"] = (up_allowed != [up_raw])
        metaL["restricted_by_other_enhancers"] = (down_allowed != [down_raw])
        return up_piece, down_piece, metaR, metaL

def run_cmd(cmd: List[str]) -> Tuple[int, str, str]:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return p.returncode, p.stdout, p.stderr

def maybe_run_repeatmasker(genome_fa: Path, rm_exec: str, rm_species: Optional[str],
                           rm_lib: Optional[Path], rm_threads: int, force: bool) -> Path:
    """
    Run RepeatMasker with -xsmall (lowercase) if needed.
    Returns path to the soft-masked fasta (typically <input>.masked).
    """
    masked = genome_fa.with_suffix(genome_fa.suffix + ".masked")
    if masked.exists() and not force:
        sys.stderr.write(f"[INFO] Using existing soft-masked genome: {masked}\n")
        return masked

    cmd = [rm_exec, "-xsmall", f"-pa={rm_threads}"]
    if rm_species:
        cmd += ["-species", rm_species]
    if rm_lib:
        cmd += ["-lib", str(rm_lib)]
    cmd += [str(genome_fa)]
    sys.stderr.write(f"[INFO] Running RepeatMasker: {' '.join(cmd)}\n")
    rc, out, err = run_cmd(cmd)
    if rc != 0:
        sys.stderr.write(f"[ERROR] RepeatMasker failed (rc={rc}). STDERR:\n{err}\n")
        sys.exit(2)
    # RepeatMasker writes <input>.masked (or sometimes <input>.fa.masked)
    if not masked.exists():
        alt = Path(str(genome_fa) + ".masked")
        if alt.exists():
            masked = alt
        else:
            sys.stderr.write("[ERROR] Could not find masked FASTA after RepeatMasker.\n")
            sys.exit(2)
    sys.stderr.write(f"[INFO] RepeatMasker output: {masked}\n")
    return masked

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shadow-bed", required=True, help="BED with shadow enhancers; col4 = set label")
    # One of these must be provided:
    ap.add_argument("--masked-genome", help="Soft-masked (lowercase) genome FASTA (preferred)")
    ap.add_argument("--genome-fasta", help="Unmasked genome FASTA; use with --rm-run to soft-mask automatically")
    # RepeatMasker options (used only if --rm-run and --genome-fasta)
    ap.add_argument("--rm-run", action="store_true", help="Run RepeatMasker (-xsmall) on --genome-fasta")
    ap.add_argument("--rm-exec", default="RepeatMasker", help="RepeatMasker executable (default: RepeatMasker)")
    ap.add_argument("--rm-species", default="drosophila", help="RepeatMasker species (default: drosophila)")
    ap.add_argument("--rm-lib", default=None, help="RepeatMasker custom library (optional)")
    ap.add_argument("--rm-threads", type=int, default=4, help="Threads for RepeatMasker (default: 4)")
    ap.add_argument("--rm-force", action="store_true", help="Re-run RepeatMasker even if masked file exists")
    # Output and flank rules
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--flank-bp", type=int, default=2500)
    ap.add_argument("--gap-threshold", type=int, default=5000)
    args = ap.parse_args()

    out_root = Path(args.outdir).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    # Prepare masked genome
    if args.masked_genome:
        masked_path = Path(args.masked_genome).resolve()
        if not masked_path.exists():
            sys.stderr.write(f"[ERROR] --masked-genome not found: {masked_path}\n")
            sys.exit(2)
    else:
        if not (args.genome_fasta and args.rm_run):
            sys.stderr.write("[ERROR] Provide either --masked-genome, or --genome-fasta with --rm-run to soft-mask.\n")
            sys.exit(2)
        genome_path = Path(args.genome_fasta).resolve()
        if not genome_path.exists():
            sys.stderr.write(f"[ERROR] --genome-fasta not found: {genome_path}\n")
            sys.exit(2)
        rm_lib = Path(args.rm_lib).resolve() if args.rm_lib else None
        masked_path = maybe_run_repeatmasker(genome_path, args.rm_exec, args.rm_species, rm_lib, args.rm_threads, args.rm_force)

    # Open masked genome with case preserved (lowercase retained)
    fa = Fasta(str(masked_path), as_raw=False, sequence_always_upper=False)

    # Load enhancers
    shadow_bed = Path(args.shadow_bed).resolve()
    enhancers = load_shadow_bed(shadow_bed)
    by_set = group_by_label(enhancers)

    manifest_rows = []

    # Iterate sets
    for label, enhs in by_set.items():
        # union of *all other* enhancers in this set (used for OUTER flank restriction)
        union_all = build_union_of_other_enhancers(enhs, exclude=[])

        n = len(enhs)
        for i in range(n):
            for j in range(i+1, n):
                eA, eB = enhs[i], enhs[j]
                if eA.chrom == eB.chrom:
                    up, down = (eA, eB) if eA.start <= eB.start else (eB, eA)
                else:
                    up, down = (eA, eB) if (eA.chrom, eA.start) <= (eB.chrom, eB.start) else (eB, eA)

                pair_dir = out_root / label / pair_dir_name(up, down)
                pair_dir.mkdir(parents=True, exist_ok=True)

                # For INNER flanks, exclude both members from the "other enhancers" union
                others_union_excl_pair = build_union_of_other_enhancers(enhs, exclude=[up, down])

                up_chr_len = len(fa[up.chrom])
                down_chr_len = len(fa[down.chrom])

                # OUTER flanks
                up_left_iv,  upL_meta  = make_outer_flank(up,  "LEFT",  args.flank_bp, up_chr_len,   union_all)
                down_right_iv,downR_meta=make_outer_flank(down,"RIGHT", args.flank_bp, down_chr_len, union_all)

                # INNER flanks (handles <5kb intervening rule)
                up_right_iv, down_left_iv, upR_meta, downL_meta = make_inner_flanks_for_pair(
                    up, down, args.flank_bp, args.gap_threshold, up_chr_len, down_chr_len, others_union_excl_pair
                )

                def write_if(iv: Optional[Interval], e: Enh, side: str, is_inner: bool, meta: Dict):
                    if iv is None:
                        return
                    s, eend = iv
                    if eend <= s:
                        return
                    # Extract preserving case (lowercase comes from RepeatMasker)
                    seq = fa[e.chrom][s:eend].seq
                    softmasked_bases = count_lowercase_bases(seq)
                    fname = f"{e.chrom}_{s}-{eend}_{side}.fa"
                    fpath = pair_dir / fname
                    header = f"{e.chrom}:{s}-{eend} side={side} set={label} pair={enh_id(up)}__{enh_id(down)}"
                    write_fa(fpath, header, seq)
                    manifest_rows.append({
                        "gene_label": label,
                        "pair_dir": str(pair_dir.relative_to(out_root)),
                        "enhancer_id": enh_id(e),
                        "flank_side": side,
                        "chrom": e.chrom,
                        "start": s,
                        "end": eend,
                        "length": eend - s,
                        "is_intervening_inner": "1" if meta.get("is_intervening", False) else "0",
                        "restricted_by_other_enhancers": "1" if meta.get("restricted_by_other_enhancers", False) else "0",
                        "truncated_at_boundary": "1" if meta.get("truncated_at_boundary", False) else "0",
                        "softmasked_bases": softmasked_bases,
                        "masked_source": "RepeatMasker",
                        "fasta_path": str((pair_dir.relative_to(out_root) / fname))  # <set>/<pair>/<file>
                    })

                # Up LEFT (outer)
                write_if(up_left_iv,  up,   "LEFT",  is_inner=False, meta=upL_meta)
                # Up RIGHT / Down LEFT (inner)
                write_if(up_right_iv, up,   "RIGHT", is_inner=True,  meta=upR_meta)
                write_if(down_left_iv, down,"LEFT",  is_inner=True,  meta=downL_meta)
                # Down RIGHT (outer)
                write_if(down_right_iv, down,"RIGHT", is_inner=False, meta=downR_meta)

    # Write manifest
    man_path = out_root / "manifest.tsv"
    with man_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "gene_label","pair_dir","enhancer_id","flank_side","chrom","start","end","length",
            "is_intervening_inner","restricted_by_other_enhancers","truncated_at_boundary",
            "softmasked_bases","masked_source","fasta_path"
        ])
        w.writeheader()
        for r in manifest_rows:
            w.writerow(r)

    print(f"[DONE] Wrote {len(manifest_rows)} flank FASTAs; manifest: {man_path}")

if __name__ == "__main__":
    main()
