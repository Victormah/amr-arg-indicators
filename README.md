# Candidate antimicrobial resistance gene indicators for veterinary public health surveillance

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23832482.svg)](https://doi.org/10.5281/zenodo.23832482)

R code for the article:

> Mahawan T, Thomrongsuwannakij T, Turni C, Omaleki L. Candidate antimicrobial resistance gene
> indicators for veterinary public health surveillance: a leakage-controlled analysis of animal,
> food and environmental metagenomes. *Research in Veterinary Science* (submitted).

The analysis reuses a curated public resource of more than 214,000 metagenomes
(Martiny et al., 2022) and asks three questions for the five sources Avian, Compost, Feces,
Food and Wastewater:

1. Can a few antimicrobial resistance genes (ARGs) stand in for the total ARG burden of a sample?
2. Do ARG profiles identify the source of a sample, also for sequencing projects not used to build the model?
3. Which repeatedly selected genes occupy bridging positions in a network of genes that rise and fall together?

The result is a short list of seven candidate indicator genes (six loci) for targeted qPCR surveillance.

## Workflow

```
data/ (public inputs)
  │
  ├─ 1 Data preparation ........ reads per million bacterial reads, log(1 + x); CLR for networks
  ├─ 2 Burden reconstruction ... per source, 5-fold CV x 5 repeats; elastic net, random forest,
  │                               support vector regression, XGBoost; held-out R²
  ├─ 3 Source classification ... nested CV (outer 5-fold x 5, inner 3-fold); balanced accuracy, recall
  ├─ 4 Robustness .............. steps 2-3 again with project-blocked CV (whole projects held out)
  ├─ 5 Network ................. CARD annotation; ecological (within-source Spearman) edges;
  │                               betweenness-based prioritization
  └─ 6 Outputs ................. result tables, figure data, manuscript figures
```

Every data-driven step (gene filtering, screening, Boruta, tuning) uses the training
portion of each split only, so test samples never influence model building.

## Repository structure

```
run_all.R                       runs the whole pipeline (steps 1-6)
R/
  data_models.R                 data preparation, regression and classification functions
  network_summary.R             CARD annotation, network, prioritization, summary tables
  network_figures.R             draft network plots and Cytoscape export
  literature.R                  optional PubMed/Europe PMC search (not used for results)
scripts/
  make_manuscript_figures.R     manuscript Figs 2-4 and Fig. S1 (TIFF + PNG, 300 dpi)
  build_card_mapping.R          optional: rebuild the ResFinder-to-CARD mapping table
  install_packages.R            installs the required R packages
data/README.md                  where to download the input files
results/                        archived results reported in the manuscript
  figure_data/                  the numbers behind every figure
  network/                      network nodes, edges, metrics, prioritized genes, CARD audit
  tables/                       Supplementary Tables S1-S7 and model summaries
docs/                           Fig. 1 workflow diagram (drawn in Canva; PDF, SVG, PNG)
```

## Requirements

- R ≥ 4.3 (the manuscript analysis used R 4.5.1).
- Packages: data.table, glmnet, Boruta, ranger, e1071, xgboost, igraph, readxl, ggplot2,
  ggraph, ggrepel, cowplot, patchwork, ragg, scales, here, withr. Install them with
  `source("scripts/install_packages.R")`, or restore the exact versions with
  `renv::restore()` when `renv.lock` is present.
- Versions used for the manuscript: glmnet 5.0, Boruta 9.0.0, ranger 0.18.0, e1071 1.7-17,
  xgboost 3.2.1.1, igraph 2.3.1, ggplot2 4.0.3. The full list is written to
  `output/sessionInfo.txt` by every run.
- Run time: the full pipeline takes several hours on a desktop computer (the models are the
  slow part). Redrawing the figures from `results/` takes about one minute.

## How to run

1. Download or clone this repository and open `amr-arg-indicators.Rproj` in RStudio
   (or start R in the repository folder).
2. Install packages: `source("scripts/install_packages.R")`.
3. Put the input files in `data/` (see [data/README.md](data/README.md)).
4. Run the pipeline:
   ```r
   source("run_all.R")
   ```
   Finished steps are saved in `output/checkpoints/`. If R stops, run the same line again
   and it continues where it stopped. Optional settings (set before `source()`):

   | Variable | Default | Meaning |
   |---|---|---|
   | `RUN_MODELS` | `TRUE` | `FALSE` reuses checkpoints and only rebuilds network, tables and figures |
   | `RUN_PROJECT_CV` | `TRUE` | project-blocked cross-validation (robustness analysis) |
   | `CENTRALITY_LAYERS` | `"Ecological"` | edge layers used to rank genes; `c("Ecological", "Mechanism", "DrugClass")` gives the all-layer sensitivity ranking |
   | `N_CORES` | `1` | parallel workers for the regression step (results are identical) |
   | `RUN_LITERATURE` | `FALSE` | optional literature search (needs internet) |
   | `CONFIG_OVERRIDES` | none | e.g. `list(outer_repeats = 1, outer_folds = 3)` for a quick test (not manuscript results) |

5. To redraw only the manuscript figures from the archived results (no input data needed):
   ```r
   source("scripts/make_manuscript_figures.R")
   ```

## Scripts, figures and tables

| Manuscript item | Produced by | File(s) |
|---|---|---|
| Table 1 (samples, provenance) | `run_all.R` step 1 | `output/qc/sample_flow.csv`; `results/tables/Supplementary_Tables_S1-S7.xlsx` (S1) |
| Table 2, Fig. 2 (burden reconstruction) | `run_all.R` step 2; `scripts/make_manuscript_figures.R` | `output/models_regression/model_summary.csv`; `Fig2_burden_reconstruction.*` |
| Table 3, Fig. 3 (source classification) | `run_all.R` steps 3-4; `scripts/make_manuscript_figures.R` | `output/results_tables/`; `Fig3_source_classification.*` |
| Table 4, Fig. 4 (network, prioritized genes) | `run_all.R` step 5; `scripts/make_manuscript_figures.R` | `output/network/`; `Fig4_ecological_network.*` |
| Fig. 1 (workflow) | drawn in Canva | `docs/Fig1_workflow.*` |
| Fig. S1 (all-layer network) | `scripts/make_manuscript_figures.R` | `FigS1_all_layer_network.*` |
| Tables S1-S7 | `run_all.R` | `results/tables/Supplementary_Tables_S1-S7.xlsx` |

## Reproducibility note

The model performance values reported in the article (Tables 2-3, Figs 2-3) come from the
authors' analysis run made before the three model fixes listed in `CHANGELOG.md` (items 1-3:
XGBoost rounds, elastic-net tuning folds, calibrated gene-sum reference). Those results are
archived in `results/`, and `scripts/make_manuscript_figures.R` redraws the manuscript figures
from them exactly. Re-running the models with this version gives the same selected genes,
network and prioritized genes (feature selection is unaffected by the fixes) but slightly
different elastic-net and XGBoost performance values. The project-blocked classification
(elastic net and random forest) was run with this version of the code.

## Data

The input files are public but are not redistributed here; [data/README.md](data/README.md)
explains where to download them. `results/` contains only derived tables (sample accessions,
model performance, selected genes, network). CARD free-text descriptions are not included;
CARD accession numbers and names are given for identification.

## How to cite

Please cite the article above and this software:

> Mahawan T, Thomrongsuwannakij T, Turni C, Omaleki L. amr-arg-indicators: R code for
> leakage-controlled prioritization of antimicrobial resistance gene indicators (version 1.0.0)
> [software]. Zenodo; 2026. https://doi.org/10.5281/zenodo.23832482

Citation metadata are in `CITATION.cff`.

## Acknowledgement

Claude (Anthropic), an AI assistant, was used to help tidy and document the code and to check it
for errors. The authors reviewed all changes and take full responsibility for the code.

## License

Code: MIT License (see `LICENSE`). Derived results in `results/` and the figure in `docs/`:
CC BY 4.0. The input data keep the licenses of their original sources.

## Contact
Dr.Tanakamol Mahawan (analysis and code), m.tanakamol@gmail.com
