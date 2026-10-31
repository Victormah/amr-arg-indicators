# Verification: the cleaned code gives the same results as the original

Date: 2026-10-31. The original scripts (`00_run_analysis.R`, `01_data_models.R`,
`02_network_summary.R`, `03_network_figures.R`, `04_literature.R`, `05_make_figures.R`) and this
repository were run side by side on the same inputs and compared output by output
(`data.table::fread` + `all.equal`, tolerance 1e-12; text files compared after removing paths and
time stamps; figures compared pixel by pixel).

## What changed in the code

- The function files in `R/` are the original files with a header added; their code is
  unchanged except two default folder names (`"AMR_reanalysis_outputs_v5"` -> `"output"`).
- `run_all.R` is the original driver with: project root found by `here::here()` instead of
  `getwd()`; functions sourced from `R/`; output folder `output/`; `setwd()` replaced by
  `withr::with_dir()` for the optional literature step; a new final step that writes the figure
  data and draws the manuscript figures. Configuration values (seeds, folds, thresholds,
  feature caps, tuning grids) are unchanged.
- `scripts/build_card_mapping.R`: Windows line endings and trailing spaces removed (parsed code
  identical).

## Results

| Test | Outputs compared | Result |
|---|---|---|
| A. Network, tables and figures from the saved model results (`RUN_MODELS = FALSE`) | 16 CSV/MD files (network nodes, edges, metrics, sensitivity, prioritized genes, CARD audit, result tables, summaries) | 13 byte-identical; 3 equal: eigenvector centrality differs by < 1e-15 (numerical solver precision; descriptive only, not used for any result) |
| B. Full pipeline in quick mode (data preparation, regression, classification, project-blocked regression and classification; 3 folds x 1 repeat) | 48 CSV/MD files and 5 checkpoint objects (preparation, regression, classification, both project-blocked analyses) | all 48 files identical; all 5 objects equal |
| C. Figure data exported by `run_all.R` vs. archived `results/figure_data/` | sample flow, regression summary, sample-level predictions, network nodes and edges | equal |
| D. `scripts/make_manuscript_figures.R` from `results/figure_data/` | Figs 2, 3, 4, S1 (PNG) | pixel-identical to the submitted figures |

## Limits of this check

- The test environment (R 4.3.3) could not install Boruta or xgboost, so in test B the Boruta
  step fell back to the documented fallback and the XGBoost models were skipped, identically in
  both versions (the same two-line test-only patch was applied to both copies). The code of
  these steps is unchanged between the versions.
- The full-size run (5 x 5 folds, all packages) is done on the authors' computer with
  `source("run_all.R")`; its `output/sessionInfo.txt` and `renv.lock` record the exact package
  versions.
