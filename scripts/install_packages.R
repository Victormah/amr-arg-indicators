# =============================================================================
# scripts/install_packages.R - install the R packages used by the pipeline
# Run once:  source("scripts/install_packages.R")
# Versions used for the manuscript are listed in README.md and renv.lock.
# =============================================================================
packages <- c(
  # analysis
  "data.table", "glmnet", "Boruta", "ranger", "e1071", "xgboost", "igraph", "readxl",
  # figures
  "ggplot2", "ggraph", "ggrepel", "cowplot", "patchwork", "ragg", "scales", "svglite",
  # project helpers
  "here", "withr", "renv",
  # optional: literature search (R/literature.R) and CARD mapping (scripts/build_card_mapping.R)
  "httr2", "jsonlite", "xml2", "digest", "dplyr", "tidyr", "readr", "writexl"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
message("All packages available.")
