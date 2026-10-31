# =============================================================================
# run_all.R - run the complete analysis pipeline
# Project : Candidate antimicrobial resistance gene indicators for veterinary
#           public health surveillance (Mahawan et al., Research in Veterinary
#           Science, submitted)
#
# How to run (from the project folder, e.g. in RStudio open the .Rproj file):
#   source("run_all.R")
#
# Inputs  (data/, see data/README.md): metadata.tsv, ARG.tsv (or the subset
#          ARG_eligible_5sources.tsv), aro.tsv, AMRgenelist_mapped_CARD_ALL.csv
# Outputs (output/): model results, network tables, result tables, figure data,
#          manuscript figures (output/figures_manuscript/) and sessionInfo.txt
#
# Steps : 1 data preparation -> 2 burden reconstruction (sample-level CV)
#         -> 3 nested source classification (sample-level CV)
#         -> 4 both again with project-blocked CV (robustness)
#         -> 5 CARD annotation and gene network -> 6 tables, figure data, figures
#
# Options (set before source()):
#   RUN_MODELS = TRUE   re-runs data preparation and all models (several hours).
#     Each finished step is saved in output/checkpoints/; if R stops, source
#     this script again and it resumes after the last finished step. Delete
#     output/checkpoints/ to start from scratch.
#   RUN_MODELS = FALSE  reuses finished checkpoints (or a saved results file,
#     see SAVED_MODELS below) and only rebuilds the network, tables and figures.
#   CONFIG_OVERRIDES = list(outer_repeats = 1, outer_folds = 3) gives a quick
#     test run (minutes); results are NOT those of the manuscript.
# =============================================================================

if (!exists("RUN_MODELS"))     RUN_MODELS     <- TRUE   # data preparation + all models
if (!exists("RUN_PROJECT_CV")) RUN_PROJECT_CV <- TRUE   # robustness: project-blocked CV
if (!exists("RUN_FIGURES"))    RUN_FIGURES    <- TRUE   # Figure 2 + Cytoscape exports
if (!exists("RUN_LITERATURE")) RUN_LITERATURE <- FALSE  # needs internet; see 04_literature.R
if (!exists("N_CORES"))        N_CORES        <- 1L     # >1 runs sources in parallel (regression)

# Edge layers used to rank nodes by betweenness. With complete CARD annotation
# the annotation layers dominate the all-layer network, so prioritization uses
# sample-level ecological associations; the all-layer ranking is reported as a
# sensitivity analysis (network_centrality_sensitivity.csv).
if (!exists("CENTRALITY_LAYERS")) CENTRALITY_LAYERS <- "Ecological"

# Project folder: located with here::here() (looks for the .here / .Rproj file),
# so the script works from any working directory inside the project.
project_dir <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
stopifnot(dir.exists(file.path(project_dir, "R")),
          dir.exists(file.path(project_dir, "data")))
input_dir   <- file.path(project_dir, "data")
output_dir  <- file.path(project_dir, if (identical(CENTRALITY_LAYERS, "Ecological"))
  "output" else "output_all_layers")
checkpoint_dir <- file.path(output_dir, "checkpoints")
# Optional saved model results used when RUN_MODELS = FALSE and no checkpoints exist.
if (!exists("SAVED_MODELS")) SAVED_MODELS <- file.path(project_dir, "data", "saved_model_results.rdata")
saved_models <- SAVED_MODELS

source(file.path(project_dir, "R", "data_models.R"))
if (RUN_MODELS && !isTRUE(get0("SKIP_PACKAGE_CHECK"))) check_required_packages()

# Input files. ARG.tsv / metadata.tsv are only read when RUN_MODELS = TRUE.
find_optional <- function(...) tryCatch(find_input_file(...), error = function(e) NA_character_)
metadata_path <- find_optional(input_dir, c("metadata.tsv", "metadata(1).tsv"), "^metadata.*\\.tsv$")
arg_path      <- find_optional(input_dir, c("ARG.tsv", "ARG_eligible_5sources.tsv", "ARG_sub.tsv.txt"),
                               "^(ARG|ARG_sub).*\\.(tsv|tsv\\.txt)$")
mapping_path  <- find_input_file(input_dir, c("AMRgenelist_mapped_CARD_ALL.csv", "AMRgenelist_mapped_CARD_ALL.xlsx"),
                                 "^AMRgenelist_mapped_CARD_ALL\\.(csv|xlsx)$")
aro_path      <- find_optional(input_dir, c("aro.tsv"), "^aro.*\\.tsv$")

# Analysis configuration (v4 values kept).
config <- list(
  metadata_path = metadata_path,
  arg_path = arg_path,
  mapping_path = mapping_path,
  aro_path = aro_path,
  output_dir = output_dir,

  source_map = data.frame(
    metadata_host = c("bird metagenome", "compost metagenome", "feces metagenome",
                      "food metagenome", "wastewater metagenome"),
    source_key = c("avian", "compost", "feces", "food", "wastewater"),
    display_name = c("Avian", "Compost", "Feces", "Food", "Wastewater"),
    stringsAsFactors = FALSE
  ),

  clr_pseudocount = 1e-6,

  # Leakage-controlled screening.
  min_prevalence_samples = 3,
  prevalence_fraction = 0.05,
  min_candidates = 10,
  max_candidates = 120,
  max_selected_features = 80,
  fallback_features = 20,
  screening_alpha = 0.50,
  use_boruta = TRUE,
  boruta_repeats = 1,
  boruta_max_runs = 50,
  stable_frequency = 0.50,

  # Outer validation.
  outer_folds = 5,
  outer_repeats = 5,
  seed = 20261013,
  num_threads = 1,

  alpha_grid = c(0, 0.25, 0.50, 0.75, 1.00),
  rf_trees = 500,
  rf_min_node_grid = c(3, 5, 10),
  svm_cost_grid = c(0.5, 1, 2),
  xgb_eta = 0.05,
  xgb_max_depth = 3L,
  xgb_min_child_weight = 3,
  xgb_subsample = 0.80,
  xgb_colsample_bytree = 0.80,
  xgb_nrounds = 300L,
  xgb_early_stopping_rounds = 20L,

  # Source classification feature selection.
  classification_prevalence_fraction = 0.03,
  classification_min_features = 30,
  classification_max_features = 150,

  # Network.
  correlation_threshold = 0.30,
  correlation_fdr = 0.05,
  min_sources_for_central = 2,
  centrality_quantile = 0.75,
  # CARD annotation: "tiered" (accession+gene, then unique exact name, then
  # curated synonyms) or "accession_and_gene" (conservative, tier 1 only).
  card_match = "tiered",
  # Layers used for betweenness-based prioritization. v4 used all three;
  # use "Ecological" to rank on sample-level associations only.
  centrality_layers = CENTRALITY_LAYERS
)
# Optional quick test: e.g. CONFIG_OVERRIDES <- list(outer_repeats = 1, outer_folds = 3)
if (exists("CONFIG_OVERRIDES")) config[names(CONFIG_OVERRIDES)] <- CONFIG_OVERRIDES
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_message("Project: ", project_dir)
log_message("CARD mapping: ", config$mapping_path)
log_message("Output directory: ", config$output_dir)

source(file.path(project_dir, "R", "network_summary.R"))

# ---------------------------------------------------------------------------
# Models (with checkpoints)
# ---------------------------------------------------------------------------
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint <- function(name, expr) {
  path <- file.path(checkpoint_dir, paste0(name, ".rds"))
  if (file.exists(path)) {
    log_message("Using checkpoint: ", name)
    return(readRDS(path))
  }
  log_message("Running: ", name)
  value <- force(expr)
  saveRDS(value, path)
  value
}

run_regression_parallel <- function(processed, config) {
  # Sources are independent and every fold has its own seed, so running them
  # in parallel gives the same results as running them one after another.
  if (N_CORES <= 1L || length(processed) < 2L) {
    return(run_repeated_nested_regression(processed, config))
  }
  cluster <- parallel::makeCluster(min(N_CORES, length(processed)))
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterExport(cluster, c("project_dir", "config"), envir = environment())
  parallel::clusterEvalQ(cluster, source(file.path(project_dir, "R", "data_models.R")))
  parts <- parallel::parLapply(cluster, names(processed), function(key) {
    single <- processed[key]
    cfg <- config
    cfg$source_order <- names(processed)
    cfg$output_dir <- file.path(config$output_dir, "parallel_parts", key)
    out <- run_repeated_nested_regression(single, cfg)
    out$key <- key
    out
  })
  # Recalculate seeds exactly as the sequential loop does (source position).
  combine_regression_parts(parts, processed, config)
}

models_available <- file.exists(file.path(checkpoint_dir, "classification.rds"))
if (RUN_MODELS || models_available) {
  if (RUN_MODELS && !isTRUE(get0("SKIP_PACKAGE_CHECK"))) check_required_packages()
  preparation    <- checkpoint("preparation", prepare_all_sources(config))
  regression     <- checkpoint("regression", run_regression_parallel(preparation$processed, config))
  classification <- checkpoint("classification",
                               run_source_classification_cv(preparation$processed, config))
  if (RUN_PROJECT_CV) {
    project_config <- config
    project_config$cv_grouping <- "project"
    project_config$output_dir <- file.path(output_dir, "project_blocked_cv")
    regression_project <- checkpoint("regression_project_blocked",
      run_regression_parallel(preparation$processed, project_config))
    classification_project <- checkpoint("classification_project_blocked",
      run_source_classification_cv(preparation$processed, project_config))
  }
} else {
  if (!file.exists(saved_models)) stop("Saved model results not found: ", saved_models)
  log_message("Loading saved model results: ", saved_models)
  saved <- new.env()
  suppressWarnings(load(saved_models, envir = saved))
  for (object in c("preparation", "regression", "classification")) {
    if (!exists(object, envir = saved, inherits = FALSE)) stop(object, " missing from ", saved_models)
  }
  preparation    <- saved$preparation
  regression     <- saved$regression
  classification <- saved$classification
  rm(saved)
  # Apply the corrected ResFinder name parsing to the saved feature tables.
  regression$stable_features <- repair_feature_labels(regression$stable_features)
  for (key in names(preparation$processed)) {
    preparation$processed[[key]]$feature_mapping <-
      repair_feature_labels(preparation$processed[[key]]$feature_mapping)
  }
}


# Network from the repeatedly selected (stable) regression indicators.
network <- build_corrected_network(preparation$processed, regression$stable_features, config)
summary_files <- write_analysis_summary(preparation, regression, classification, network, config)
write_csv_atomic(regression$stable_features, file.path(output_dir, "stable_indicator_features.csv"))
write_model_reports(regression, classification, output_dir,
  regression_project = if (exists("regression_project")) regression_project else NULL,
  classification_project = if (exists("classification_project")) classification_project else NULL)

log_message("Network: ", network$metrics$nodes, " nodes, ", network$metrics$edges, " edges; prioritized: ",
            paste(network$prioritized$node, collapse = ", "))

if (RUN_FIGURES) {
  source(file.path(project_dir, "R", "network_figures.R"))
  # Main figure: ecological (sample-level) edges used for prioritization.
  plot_arg_network(output_dir = output_dir, export = TRUE, show_layers = "Ecological",
                   file_stem = "ARG_network_ecological_two_panel")
  # Supplementary figure: all three edge layers.
  plot_arg_network(output_dir = output_dir, export = TRUE)
  export_network_for_cytoscape(output_dir = output_dir)
}

if (RUN_LITERATURE) {
  source(file.path(project_dir, "R", "literature.R"))
  LITERATURE_GENES <- network$prioritized$node
  # The literature functions read and write relative to the working directory.
  withr::with_dir(output_dir, run_literature_mapping())
}

# Figure data (the numbers behind every manuscript figure) and manuscript figures.
figure_data_dir <- file.path(output_dir, "figure_data")
dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)
write_csv_atomic(preparation$sample_flow, file.path(figure_data_dir, "sample_flow.csv"))
write_csv_atomic(regression$model_summary, file.path(figure_data_dir, "regression_model_summary.csv"))
write_csv_atomic(classification$predictions,
                 file.path(figure_data_dir, "classification_predictions_sample.csv"))
if (exists("classification_project")) {
  write_csv_atomic(classification_project$predictions,
                   file.path(figure_data_dir, "classification_predictions_project.csv"))
}
for (file in c("network_nodes.csv", "network_edges.csv")) {
  file.copy(file.path(output_dir, "network", file), file.path(figure_data_dir, file), overwrite = TRUE)
}
if (RUN_FIGURES && file.exists(file.path(figure_data_dir, "classification_predictions_project.csv"))) {
  FIG_DATA_DIR <- figure_data_dir
  OUT_DIR <- file.path(output_dir, "figures_manuscript")
  source(file.path(project_dir, "scripts", "make_manuscript_figures.R"))
} else if (RUN_FIGURES) {
  log_message("Manuscript figures skipped: project-blocked results not available ",
              "(run with RUN_PROJECT_CV = TRUE, or use results/figure_data).")
}

log_message("Done. Summary: ", summary_files$summary)
save(preparation, regression, classification, network, config,
     file = file.path(output_dir, "run_all_results.rdata"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
