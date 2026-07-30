# Shadow Birth Paper 2026

Comparative genomics project testing whether shadow enhancers (redundant enhancer pairs/sets driving the same expression pattern) arose through duplication of transposable elements (TEs), in *Drosophila* (fly) and mouse.

Pipeline: isolate shadow-enhancer sets and their flanking sequences → align/BLAST flanks and enhancer bodies against each other and against null/random controls to detect TE-mediated duplication signatures → quantify TE-class enrichment and "TE splitting" across enhancer boundaries genome-wide → compare TFBS similarity between shadow pairs → generate paper figures.

## Folder structure

### `TE_cooption/`
Fly and mouse subfolders (`fly/`, `mouse/`), each with `input/` (shadow vs. single enhancer BEDs, TE annotations), `scripts/` (TE-class enrichment / observed-vs-neighborhood stats), and output (`output/` for fly, `outputs/` for mouse — inconsistent naming, see notes). `fly/plot_figures/` holds figure-generation code for this section.

### `TE_splitting/`
Smaller companion analysis: whether TEs are split across enhancer boundaries. Same `fly/` and `mouse/` layout, each with `scripts/` and `output/` (hit tables + donut-plot scripts).

### `duplication_analysis/`
The core duplication-signature analysis

- **`FLY_flank_analysis/`** — `make_shadow_flanks_v2.py` builds flank FASTAs per shadow pair. `shadow_pairs_flanks2/<gene>/` holds one subfolder per gene/shadow-pair (e.g. `FBgn0023095_sm`, `Abd-B`) with BLAST results (`LL/RR/LR/RL` hit files). `null/` holds the null-model BLAST comparison and cutoff scripts/tables. `finaloutput_enhancerhitANDdoubleflank/` holds the final collated output.
- **`MOUSE_flank_analysis/`** — mouse equivalent: `shadow_duplication_null_pipeline.py`, `null/` and `final/` result tables.
- **`TFBS_similarity_analysis/`** — TFBS Jaccard similarity and EMBOSS-needle global alignment comparisons between shadow-pair sequences (`TFBS_similarity_v2_groupmotifs.py`, `run_Globalalignments_with_kmer.py`).
- **`blast_enhancerbody_shadow_sets_analysis/`** — earlier/simpler enhancer-body BLAST pipeline, largely superseded by the FLY/MOUSE flank analyses above (currently missing its expected `input/` folder — see notes).
- **`enhancer_hits_per_shadow_bin/`** — bins shadow sets by size and tests whether larger sets get more enhancer "hits." Has its own `Fly/` and `Mouse/` subpipelines and its own README (that README is stale — see notes).
- **`plots/`** — final figure-generation scripts (`Fly_plots_code/`, `mouse_plots_Figure5/`) reading from the analyses above, writing to `output_pngs/`.
