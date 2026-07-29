#!/usr/bin/env python3
import argparse, os, re, csv, math
from collections import Counter, defaultdict
from typing import Dict, Tuple, Set, List, FrozenSet, Optional

# ---------------- basics ----------------
COORD_RE = re.compile(r'^(chr[^_:]+)[:_](\d+)-(\d+)$')

def to_eid(s: str) -> str:
    """Normalize:
       - chr2R:24854572-24855864
       - chr2R_24854572-24855864
       - filenames like chr2R_24854572-24855864.fa
       -> 'chr2R:24854572-24855864'
    """
    s = os.path.basename(s)
    s = re.sub(r'\.fa(sta)?$', '', s)
    m = COORD_RE.match(s)
    if m:
        chrom, a, b = m.groups()
        return f'{chrom}:{a}-{b}'
    parts = re.findall(r'(chr[^_:]+)[_:](\d+)-(\d+)', s)
    if parts:
        chrom, a, b = parts[-1]
        return f'{chrom}:{a}-{b}'
    raise ValueError(f"Cannot parse coordinates from: {s}")

def eid_to_basename(eid: str) -> str:
    chrom, span = eid.split(':', 1)
    return f'{chrom}_{span}'

def list_pair_dirs(enhancer_root: str) -> List[Dict]:
    """Find <shadowset>/<pair>/enhancers/ with >=2 FASTAs; derive the two EIDs."""
    out = []
    for shadowset in os.listdir(enhancer_root):
        ss_path = os.path.join(enhancer_root, shadowset)
        if not os.path.isdir(ss_path): continue
        for pairdir in os.listdir(ss_path):
            pd_path = os.path.join(ss_path, pairdir)
            if not os.path.isdir(pd_path): continue
            enh_dir = os.path.join(pd_path, 'enhancers')
            if not os.path.isdir(enh_dir): continue
            fastas = [f for f in os.listdir(enh_dir) if f.lower().endswith(('.fa','.fasta'))]
            if len(fastas) < 2: continue
            if '__' in pairdir:
                left, right = pairdir.split('__', 1)
                eida = to_eid(left); eidb = to_eid(right)
            else:
                eids = sorted({to_eid(x) for x in fastas})
                if len(eids) < 2: continue
                eida, eidb = eids[:2]
            out.append({'shadowset': shadowset, 'pair_dir': pd_path, 'eida': eida, 'eidb': eidb})
    if not out:
        raise SystemExit(f"No pair directories found under: {enhancer_root}")
    return out

def read_pairs_csv(csv_path: str) -> Set[FrozenSet[str]]:
    """Robust: header/no header, CSV/TSV, extra cols; collects first two coord-like tokens."""
    pairs: Set[FrozenSet[str]] = set()
    with open(csv_path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            parts = re.split(r'[,\t]+', line)
            coords = []
            for tok in parts:
                tok = tok.strip().strip('"').strip("'")
                m = COORD_RE.match(tok)
                if m:
                    chrom, a, b = m.groups()
                    coords.append(f'{chrom}:{a}-{b}')
                    if len(coords) == 2:
                        break
            if len(coords) == 2:
                pairs.add(frozenset((coords[0], coords[1])))
    return pairs

# --------- TF name canonicalization (Step 1) ----------
def canonical_tf(raw: str) -> str:
    """
    Turn MEME names like 'MA0247.1.tin', 'MA0535.1.Mad', 'dl', 'su(Hw)_v1'
    into a canonical TF symbol:
      - take substring after the LAST '.' if present (… .tin -> tin)
      - strip common suffix noise like _Dmel/_variantN
    """
    s = raw.strip()
    if '.' in s:
        s = s.split('.')[-1]
    s = re.sub(r'(_(Dmel|Dm|fly|variant\d+|var\d+))$', '', s, flags=re.IGNORECASE)
    return s

# --------- Parse MEME to get id -> raw name & canonical symbol (Step 2) ----------
def parse_meme_names(meme_path: str):
    """
    Returns (id2name_raw, id2symbol):
      - id2name_raw[id] = raw MEME name string (may include ID prefixes)
      - id2symbol[id]   = canonical TF symbol (e.g., 'tin', 'br', 'Mad')
    """
    id2name_raw = {}
    id2symbol   = {}
    with open(meme_path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            if not line.startswith('MOTIF '):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                mid = parts[1]
                raw_name = " ".join(parts[2:]) if len(parts) > 2 else parts[1]
                raw_name = raw_name.strip()
                if raw_name == mid or raw_name == "":
                    raw_name = mid
                id2name_raw[mid] = raw_name
                id2symbol[mid]   = canonical_tf(raw_name)
    return id2name_raw, id2symbol

def tf_root(name: str) -> str:
    """Coarse family/root from a TF symbol."""
    name = name.strip()
    name = re.sub(r'(_(Dmel|Dm|fly|variant\d+|var\d+))$', '', name, flags=re.IGNORECASE)
    m = re.match(r'[A-Za-z0-9()/+-]+', name)
    return m.group(0) if m else name

# ------------- FIMO readers -------------
def find_fimo_table(fimo_root: str, eid: str) -> Optional[str]:
    """Return path to fimo.tsv (preferred) or fimo.txt; else fimo.gff if present."""
    oc = os.path.join(fimo_root, eid_to_basename(eid))
    for cand in ("fimo.tsv", "fimo.txt"):
        p = os.path.join(oc, cand)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    gff = os.path.join(oc, "fimo.gff")
    return gff if os.path.exists(gff) and os.path.getsize(gff) > 0 else None

def _header_index_map(header_cells: List[str]) -> Dict[str,int]:
    """Map normalized header names to indices."""
    norm = {}
    for i, h in enumerate(header_cells):
        key = re.sub(r'[^a-z0-9]+', '_', h.strip().lower())
        norm[key] = i
    return norm

def load_fimo_hits_table(path: str, p_thresh: Optional[float], q_thresh: Optional[float]) -> List[Dict]:
    """Read fimo.tsv/txt or fimo.gff. Returns rows with motif_id, pvalue, qvalue, start, stop."""
    rows: List[Dict] = []
    if path.endswith(".gff"):
        with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'): continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9: continue
                attrs = cols[8]
                # motif id
                m = re.search(r'(?:Name|motif_id)=([^;]+)', attrs)
                mid = m.group(1) if m else None
                if not mid: continue
                # optional p & q in attributes like pvalue=1.2e-4; qvalue=0.05
                mp = re.search(r'pvalue=([0-9.eE+-]+)', attrs)
                mq = re.search(r'qvalue=([0-9.eE+-]+)', attrs)
                pval = float(mp.group(1)) if mp else None
                qval = float(mq.group(1)) if mq else None
                # optionally parse positions from columns 3/4 (1-based)
                try:
                    start = int(cols[3])
                    stop  = int(cols[4])
                except Exception:
                    start = stop = None
                keep = True
                if pval is not None and p_thresh is not None:
                    keep = (pval <= p_thresh)
                if keep and (q_thresh is not None) and (qval is not None):
                    keep = (qval <= q_thresh)
                if keep:
                    rows.append({'motif_id': mid, 'pvalue': pval, 'qvalue': qval, 'start': start, 'stop': stop})
        return rows

    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = None
        idx = {}
        for raw in reader:
            if not raw: continue
            if raw[0].startswith("#"): continue
            if header is None:
                header = raw
                hmap = _header_index_map(header)
                idx['motif'] = hmap.get('motif_id') or hmap.get('pattern_name') or 0
                idx['p'] = hmap.get('p_value') or hmap.get('pvalue')
                idx['q'] = hmap.get('q_value') or hmap.get('qvalue')
                idx['start'] = hmap.get('start')
                idx['stop']  = hmap.get('stop')
                continue
            try:
                motif_id = raw[idx['motif']].strip()
            except Exception:
                continue
            pval = float(raw[idx['p']]) if idx.get('p') is not None and raw[idx['p']] else None
            qval = float(raw[idx['q']]) if idx.get('q') is not None and raw[idx['q']] else None
            start = int(raw[idx['start']]) if idx.get('start') is not None and str(raw[idx['start']]).isdigit() else None
            stop  = int(raw[idx['stop']])  if idx.get('stop')  is not None and str(raw[idx['stop']]).isdigit()  else None
            keep = True
            if pval is not None and p_thresh is not None:
                keep = (pval <= p_thresh)
            if keep and (q_thresh is not None) and (qval is not None):
                keep = (qval <= q_thresh)
            if keep:
                rows.append({'motif_id': motif_id, 'pvalue': pval, 'qvalue': qval, 'start': start, 'stop': stop})
    return rows

# ------------- similarity metrics -------------
def jaccard(a: Set[str], b: Set[str]) -> float:
    u = a | b
    return (len(a & b) / float(len(u))) if u else 0.0

def cosine(counts_a: Counter, counts_b: Counter) -> float:
    if not counts_a and not counts_b: return 0.0
    keys = set(counts_a) | set(counts_b)
    dot = sum(counts_a.get(k, 0) * counts_b.get(k, 0) for k in keys)
    na = math.sqrt(sum(v*v for v in counts_a.values()))
    nb = math.sqrt(sum(v*v for v in counts_b.values()))
    return (dot / (na * nb)) if (na and nb) else 0.0

def overlap_coeff(a_set: Set[str], b_set: Set[str]) -> float:
    if not a_set or not b_set:
        return 0.0
    return len(a_set & b_set) / float(min(len(a_set), len(b_set)))

def bray_curtis_sim(counts_a: Counter, counts_b: Counter) -> float:
    keys = set(counts_a) | set(counts_b)
    num = sum(abs(counts_a.get(k,0) - counts_b.get(k,0)) for k in keys)
    den = sum(counts_a.get(k,0) + counts_b.get(k,0) for k in keys)
    return 1.0 - (num / den) if den else 0.0

def hellinger_sim(counts_a: Counter, counts_b: Counter) -> float:
    keys = set(counts_a) | set(counts_b)
    sa = sum(counts_a.get(k,0) for k in keys)
    sb = sum(counts_b.get(k,0) for k in keys)
    if sa == 0 or sb == 0:
        return 0.0
    pa = {k: counts_a.get(k,0)/sa for k in keys}
    pb = {k: counts_b.get(k,0)/sb for k in keys}
    s = sum((math.sqrt(pa[k]) - math.sqrt(pb[k]))**2 for k in keys)
    H = (1.0 / math.sqrt(2.0)) * math.sqrt(s)  # distance in [0,1]
    return max(0.0, 1.0 - H)

def weighted_jaccard(weights_a: Dict[str,float], weights_b: Dict[str,float]) -> float:
    keys = set(weights_a) | set(weights_b)
    num = sum(min(weights_a.get(k,0.0), weights_b.get(k,0.0)) for k in keys)
    den = sum(max(weights_a.get(k,0.0), weights_b.get(k,0.0)) for k in keys)
    return (num / den) if den else 0.0

def js_divergence(p: List[float], q: List[float]) -> float:
    """JSD with log base 2; p,q auto-normalized with epsilon."""
    eps = 1e-12
    sp, sq = sum(p), sum(q)
    p = [(x + eps) / (sp + eps*len(p)) for x in p]
    q = [(x + eps) / (sq + eps*len(q)) for x in q]
    m = [(pi + qi) / 2.0 for pi, qi in zip(p, q)]
    def kl(a, b):
        return sum(ai * math.log2(ai / bi) for ai, bi in zip(a, b))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)

def pos_jsd_sim(hist_a: List[int], hist_b: List[int]) -> float:
    """Similarity in [0,1], 1 - JSD(positional histograms)."""
    if not hist_a or not hist_b or len(hist_a) != len(hist_b):
        return 0.0
    jsd = js_divergence(hist_a, hist_b)
    return max(0.0, 1.0 - jsd)

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser(description="Compute TFBS similarity per shadow pair from EXISTING FIMO outputs; collapse motif IDs to canonical TF names.")
    ap.add_argument('enhancer_root', help="Root of enhancer_pairs tree (to find pair dirs)")
    ap.add_argument('--fimo_root', required=True, help="Root folder containing fimo_out/<enhancer_id>/fimo.tsv|fimo.txt|fimo.gff")
    ap.add_argument('--motifs', required=True, help="MEME PWM file used (for motif_id→TF name mapping)")
    ap.add_argument('--group1_csv', required=True, help="CSV of enhancer-hit pairs (col1, col2 in chr:start-end)")
    ap.add_argument('--group2_csv', required=True, help="CSV of flank-hit-only pairs (col1, col2)")
    ap.add_argument('--outdir', required=True, help="Where to write outputs")
    ap.add_argument('--collapse', choices=['id','name','root'], default='name',
                    help="Feature space: 'id' (raw motif IDs), 'name' (canonical TF symbols), 'root' (coarse family). Default: name")
    ap.add_argument('--p_thresh', type=float, default=1e-3, help="Keep hits with p-value ≤ this (default 1e-3). Use None to skip.")
    ap.add_argument('--q_thresh', type=float, default=None, help="(Optional) Also require q-value ≤ this.")
    ap.add_argument('--no-plot', action='store_true', help="Skip plotting (avoids matplotlib import).")
    args = ap.parse_args()

    enhancer_root = os.path.abspath(args.enhancer_root)
    fimo_root = os.path.abspath(args.fimo_root)
    outdir = os.path.abspath(args.outdir)
    motifs_path = os.path.abspath(args.motifs)
    os.makedirs(outdir, exist_ok=True)

    # list pairs & groups
    pairs = list_pair_dirs(enhancer_root)
    g1 = read_pairs_csv(args.group1_csv)
    g2 = read_pairs_csv(args.group2_csv)

    def grp_of(a: str, b: str) -> str:
        key = frozenset((a,b))
        if key in g1: return 'enhancer_hit'
        if key in g2: return 'flank_hit'
        return 'neither'

    # (Step 2) motif id -> raw name & canonical symbol
    id2name_raw, id2symbol = parse_meme_names(motifs_path)

    # cache per-enhancer features
    cache: Dict[str, Dict] = {}

    # (Step 3) extract features using canonical symbol for "name" collapse
    def features_for_eid(eid: str) -> Dict:
        if eid in cache: return cache[eid]
        table = find_fimo_table(fimo_root, eid)

        ids_set, ids_counts = set(), Counter()
        names_set, names_counts = set(), Counter()
        roots_set, roots_counts = set(), Counter()
        best_w_by_id, best_w_by_name, best_w_by_root = {}, {}, {}
        bins = 30
        pos_hist = [0]*bins

        if table is not None:
            rows = load_fimo_hits_table(table, p_thresh=args.p_thresh, q_thresh=args.q_thresh)
            # enhancer length from EID
            enh_len = None
            try:
                _, span = eid.split(':',1)
                a,b = map(int, span.split('-'))
                enh_len = max(1, b - a)
            except Exception:
                pass

            for r in rows:
                mid = r['motif_id']
                ids_set.add(mid); ids_counts[mid] += 1
                # weight: best -log10(p) per motif
                w = None
                if r['pvalue'] is not None and r['pvalue'] > 0:
                    w = -math.log10(r['pvalue'])
                    if w > best_w_by_id.get(mid, 0.0): best_w_by_id[mid] = w

                # canonical symbol for "name" collapse
                name = id2symbol.get(mid, canonical_tf(mid))
                names_set.add(name); names_counts[name] += 1
                if w is not None and w > best_w_by_name.get(name, 0.0): best_w_by_name[name] = w

                # root/family
                root = tf_root(name)
                roots_set.add(root); roots_counts[root] += 1
                if w is not None and w > best_w_by_root.get(root, 0.0): best_w_by_root[root] = w

                # positional histogram (midpoint)
                s, t = r.get('start'), r.get('stop')
                if s is not None and t is not None and enh_len:
                    midpos = ((s + t) / 2.0) / float(enh_len)
                    k = int(max(0, min(bins-1, math.floor(midpos * bins))))
                    pos_hist[k] += 1

        feat = {
            'id_set': ids_set, 'id_counts': ids_counts,
            'name_set': names_set, 'name_counts': names_counts,
            'root_set': roots_set, 'root_counts': roots_counts,
            'w_id': best_w_by_id, 'w_name': best_w_by_name, 'w_root': best_w_by_root,
            'pos_hist': pos_hist,
            'has_fimo': table is not None
        }
        cache[eid] = feat
        return feat

    # collect values for plots
    metrics_list = ['jaccard','cosine','overlap_coeff','weighted_jaccard','bray_curtis_sim','hellinger_sim','pos_jsd_sim']
    metrics_by_group = {m: {'enhancer_hit': [], 'flank_hit': [], 'neither': []} for m in metrics_list}
    labels = ['enhancer_hit', 'flank_hit', 'neither']

    # (Step 4) compute similarities
    out_tsv = os.path.join(outdir, f"tfbs_similarity_three_groups_{args.collapse}.tsv")
    with open(out_tsv, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['shadowset','pair_dir','enhancer_a','enhancer_b',
                    'feat_space','size_a','size_b','hits_a','hits_b',
                    'intersection','union',
                    'jaccard','cosine','overlap_coeff','weighted_jaccard',
                    'bray_curtis_sim','hellinger_sim','pos_jsd_sim',
                    'group','fimo_a','fimo_b'])
        for row in pairs:
            a, b = row['eida'], row['eidb']
            fa = features_for_eid(a); fb = features_for_eid(b)
            if args.collapse == 'id':
                Sa, Sb = fa['id_set'], fb['id_set']
                Ca, Cb = fa['id_counts'], fb['id_counts']
                Wa, Wb = fa['w_id'], fb['w_id']
            elif args.collapse == 'name':
                Sa, Sb = fa['name_set'], fb['name_set']
                Ca, Cb = fa['name_counts'], fb['name_counts']
                Wa, Wb = fa['w_name'], fb['w_name']
            else:
                Sa, Sb = fa['root_set'], fb['root_set']
                Ca, Cb = fa['root_counts'], fb['root_counts']
                Wa, Wb = fa['w_root'], fb['w_root']

            inter = len(Sa & Sb); uni = len(Sa | Sb)
            j = jaccard(Sa, Sb)
            c = cosine(Ca, Cb)
            oc = overlap_coeff(Sa, Sb)
            wj = weighted_jaccard(Wa, Wb)
            bc = bray_curtis_sim(Ca, Cb)
            hs = hellinger_sim(Ca, Cb)
            pjs = pos_jsd_sim(fa['pos_hist'], fb['pos_hist'])
            grp = grp_of(a, b)

            # collect for plots
            for m, val in [('jaccard', j),
                           ('cosine', c),
                           ('overlap_coeff', oc),
                           ('weighted_jaccard', wj),
                           ('bray_curtis_sim', bc),
                           ('hellinger_sim', hs),
                           ('pos_jsd_sim', pjs)]:
                metrics_by_group[m][grp].append(val)

            w.writerow([row['shadowset'], row['pair_dir'], a, b,
                        args.collapse, len(Sa), len(Sb), sum(Ca.values()), sum(Cb.values()),
                        inter, uni,
                        f"{j:.6f}", f"{c:.6f}", f"{oc:.6f}", f"{wj:.6f}",
                        f"{bc:.6f}", f"{hs:.6f}", f"{pjs:.6f}",
                        grp,
                        'yes' if fa['has_fimo'] else 'no',
                        'yes' if fb['has_fimo'] else 'no'])

    print(f"[done] wrote: {out_tsv}")

    # plotting
    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # ---- Jaccard histogram ----
        plt.figure(figsize=(7,5))
        bins = 25
        for lab in labels:
            vals = [x for x in metrics_by_group['jaccard'][lab] if x == x]
            if vals:
                plt.hist(vals, bins=bins, alpha=0.6, label=f'{lab} (n={len(vals)})', density=True)
        plt.xlabel(f"TFBS similarity (Jaccard, collapsed by {args.collapse}; allTFs)")
        plt.ylabel("Density")
        plt.legend()
        plt.tight_layout()
        out_png = os.path.join(outdir, f"tfbs_similarity_three_groups_{args.collapse}.png")
        plt.savefig(out_png, dpi=200)
        print(f"[done] wrote: {out_png}")

        # ---- Per-metric boxplots ----
        pretty = {
            'jaccard': 'Jaccard',
            'cosine': 'Cosine similarity',
            'overlap_coeff': 'Overlap coefficient',
            'weighted_jaccard': 'Weighted Jaccard (best −log10 p)',
            'bray_curtis_sim': 'Bray–Curtis similarity',
            'hellinger_sim': 'Hellinger similarity',
            'pos_jsd_sim': 'Positional JSD similarity',
        }
        for metric in metrics_list:
            box_labels, box_data = [], []
            for lab in labels:
                vals = [x for x in metrics_by_group[metric][lab] if x == x]
                if vals:
                    box_labels.append(f"{lab} (n={len(vals)})")
                    box_data.append(vals)
            if not box_data:
                print(f"[warn] No data for metric '{metric}', skipping boxplot.")
                continue
            plt.figure(figsize=(7,5))
            plt.boxplot(
                box_data,
                labels=box_labels,
                notch=True,
                showmeans=True,
                meanline=True
            )
            plt.ylabel(f"{pretty.get(metric, metric)}")
            plt.title(f"{pretty.get(metric, metric)} by group ({args.collapse}; allTFs)")
            plt.tight_layout()
            out_box = os.path.join(outdir, f"tfbs_similarity_{args.collapse}_{metric}_boxplot.png")
            plt.savefig(out_box, dpi=200)
            print(f"[done] wrote: {out_box}")

if __name__ == "__main__":
    main()
