This folder collects the fly per-shadow-bin inputs and scripts that fed the hits-per-bin analysis.

Contents

- `original_scripts/final_NEWenhancer_fig2_tfbs_FIXED_violon_pct.py`
  - Original larger figure script containing `plot_enh_hit_by_shadow_category(...)` for the observed fly per-bin plot.
- `original_scripts/enhancer_newmethod_ineach_shadow_bin.py`
  - Original null/control script for the new flank-aware method using `null_distributions_enhancer_extra.csv` and `duplicated_enhancers_extra_copies_by_group.csv`.
- `original_scripts/redo_fig_size_bin_purple_corrected.py`
  - Later cleaned observed-only script using exact enhancer coordinates once per bin.

- `data/011925_all_shadowsets_DM6.bed`
  - Fly shadow-set BED used to define set size bins.
- `data/FINAL_breakdown_single_double_flanks_UPDATED_1e4_flankmethod.csv`
  - Final fly enhancer/flank breakdown used to define observed hit enhancers.
- `data/null_distributions_enhancer_extra.csv`
  - Null enhancer-per-bin distributions used by the new-method null/control plot.
- `data/duplicated_enhancers_extra_copies_by_group.csv`
  - Observed extra-copy summary used by the new-method null/control script.

- `scripts/run_fly_observed_hits_per_bin.py`
  - Clean local runner for the observed fly per-bin plot from the copied BED + final breakdown.
- `scripts/run_fly_null_control_extra_copies.py`
  - Clean local runner for the new-method null/control boxplot using only files in this folder.

Notes

- The original scripts still contain historical absolute paths and extra code.
- The two scripts in `scripts/` are the self-contained local versions.
- Output files are written to `outputs/` inside this folder.
