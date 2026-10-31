# Changelog

## 1.0.0 (2026-10-31)

First public release, accompanying the submission to *Research in Veterinary Science*.

Changes from the authors' earlier internal (v5) analysis code:

1. **XGBoost** — the number of boosting rounds is taken from the inner cross-validation /
   validation evaluation log. In v5 early stopping never took effect with xgboost 3.x, so all
   v5 XGBoost models used 100 rounds.
2. **Elastic-net regression** — every mixing parameter (alpha) is compared on the same inner folds.
3. **Gene-sum reference model** — the unscaled sum of selected genes is linearly calibrated to
   the target on the training fold.
4. **ResFinder names** — reference names are parsed correctly for NC_/NZ_ accessions.
5. **CARD annotation** — matched by accession plus gene name, then unique exact name, then
   documented synonyms; never by accession alone (v5 merged aph(6)-Id with aph(3'')-Ib).
   `output/network/selected_ARG_CARD_mapping_audit.csv` lists the method for every gene.
6. **Prioritization** — betweenness is calculated on ecological (sample-level) edges; the
   all-layer ranking is reported as a sensitivity analysis.
7. **Project-blocked cross-validation** added as a robustness analysis.
8. **Reproducibility** — checkpointing, optional parallel regression (identical results),
   automatic result tables, figure-data export and manuscript figures; project root found
   with `here`, no `setwd()`.
