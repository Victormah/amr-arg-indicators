# =============================================================================
# R/data_models.R - data preparation and leakage-controlled models
# Project : Candidate antimicrobial resistance gene indicators for veterinary
#           public health surveillance (Mahawan et al., Research in Veterinary
#           Science, submitted)
# Purpose : Functions to read and normalize the ARG counts (reads per million
#           bacterial reads), and to run (i) within-source reconstruction of total
#           ARG burden and (ii) nested source classification, with sample-level or
#           project-blocked cross-validation. All gene screening and tuning use
#           training samples only.
# Used by : run_all.R (sourced; defines functions only, runs nothing)
# Manuscript: Methods 2.2-2.5; Tables 2-3; Figs 2-3; Supplementary Tables S1-S4
# =============================================================================

# Data preparation and model functions. Original v4 calculation steps kept.

# This file intentionally does not install packages.

required_packages <- c(
  "data.table", "glmnet", "Boruta", "ranger", "e1071", "xgboost",
  "igraph", "readxl", "ggplot2"
)

check_required_packages <- function() {
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      ". Install them before running, for example: install.packages(c('",
      paste(missing_packages, collapse = "','"), "'))"
    )
  }
  invisible(TRUE)
}

log_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(...))
}

write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  data.table::fwrite(x, temporary)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) {
    stop("Could not move temporary output into place: ", path)
  }
  invisible(path)
}

rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

mae <- function(actual, predicted) {
  mean(abs(actual - predicted), na.rm = TRUE)
}

r_squared <- function(actual, predicted) {
  # Out-of-sample R-squared relative to the validation-set mean.
  valid <- is.finite(actual) & is.finite(predicted)
  actual <- actual[valid]
  predicted <- predicted[valid]
  denominator <- sum((actual - mean(actual))^2)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  1 - sum((actual - predicted)^2) / denominator
}

mean_ci_by_repetition <- function(values, conf = 0.95) {
  values <- values[is.finite(values)]
  n <- length(values)
  if (!n) return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_, n = 0))
  estimate <- mean(values)
  if (n == 1L) return(c(mean = estimate, lower = NA_real_, upper = NA_real_, n = 1))
  critical <- stats::qt((1 + conf) / 2, df = n - 1L)
  half_width <- critical * stats::sd(values) / sqrt(n)
  c(mean = estimate, lower = estimate - half_width, upper = estimate + half_width, n = n)
}

quantile_stratified_folds <- function(y, k, seed) {
  # Regression folds stratified by target quantiles.
  set.seed(seed)
  unique_breaks <- unique(stats::quantile(
    y,
    probs = seq(0, 1, length.out = min(10L, length(y)) + 1L),
    na.rm = TRUE,
    type = 8
  ))
  if (length(unique_breaks) < 3L) {
    strata <- rep(1L, length(y))
  } else {
    strata <- cut(y, breaks = unique_breaks, include.lowest = TRUE, labels = FALSE)
  }
  fold_id <- integer(length(y))
  for (stratum in sort(unique(strata))) {
    indices <- which(strata == stratum)
    indices <- sample(indices, length(indices))
    fold_id[indices] <- rep(seq_len(k), length.out = length(indices))
  }
  fold_id
}

# Project-blocked folds: every sample of a sequencing project goes to the same
# fold, so held-out samples come from projects unseen during training.
# Projects are assigned largest-first to the fold that currently holds the
# fewest samples of the project's class (or overall, for regression).
group_folds <- function(groups, k, seed, strata = NULL) {
  set.seed(seed)
  groups <- as.character(groups)
  if (is.null(strata)) strata <- rep("all", length(groups))
  strata <- as.character(strata)
  info <- data.frame(group = groups, stratum = strata, stringsAsFactors = FALSE)
  sizes <- stats::aggregate(rep(1L, nrow(info)), by = list(group = info$group), FUN = sum)
  names(sizes)[2] <- "n"
  sizes$stratum <- vapply(sizes$group, function(g) {
    tab <- table(info$stratum[info$group == g]); names(tab)[which.max(tab)]
  }, character(1))
  sizes <- sizes[sample(nrow(sizes)), , drop = FALSE]          # random tie-breaking
  sizes <- sizes[order(-sizes$n), , drop = FALSE]
  load <- matrix(0, nrow = length(unique(sizes$stratum)), ncol = k,
                 dimnames = list(unique(sizes$stratum), NULL))
  assignment <- integer(nrow(sizes))
  for (i in seq_len(nrow(sizes))) {
    row <- load[sizes$stratum[i], ]
    target <- which(row == min(row))
    target <- if (length(target) > 1L) sample(target, 1L) else target
    assignment[i] <- target
    load[sizes$stratum[i], target] <- load[sizes$stratum[i], target] + sizes$n[i]
  }
  assignment[match(groups, sizes$group)]
}

# Best boosting round from an xgboost evaluation log (works for xgboost 1.x-3.x).
xgb_best_round_from_log <- function(log, metric_pattern) {
  if (is.null(log)) return(NA_integer_)
  log <- as.data.frame(log)
  column <- grep(metric_pattern, names(log), value = TRUE)
  if (!length(column)) return(NA_integer_)
  values <- log[[column[1]]]
  if (!any(is.finite(values))) return(NA_integer_)
  as.integer(which.min(values))
}

xgb_model_log <- function(model) {
  log <- model$evaluation_log
  if (is.null(log)) log <- attributes(model)$evaluation_log
  log
}

class_stratified_folds <- function(y, k, seed) {
  # Classification folds stratified by source class.
  set.seed(seed)
  y <- as.factor(y)
  fold_id <- integer(length(y))
  for (class_value in levels(y)) {
    indices <- which(y == class_value)
    indices <- sample(indices, length(indices))
    fold_id[indices] <- rep(seq_len(k), length.out = length(indices))
  }
  fold_id
}

safe_gene_names <- function(original_names) {
  make.names(original_names, unique = TRUE)
}

# ResFinder reference names have the form <gene>_<variant>_<accession>.
# Accessions may themselves contain underscores (e.g. NC_010937, NZ_AUJS01000017),
# so the gene is everything before the FIRST "_<digits>_" and the accession is
# everything after it. (Fixed 2026-10-30: the earlier "last underscore" rule
# truncated NC_/NZ_ accessions and left the full string as the gene label.)
resfinder_ref_pattern <- "^(.+?)_([0-9]+)_(.+)$"

extract_accession <- function(ref_sequence) {
  ref_sequence <- as.character(ref_sequence)
  ok <- grepl(resfinder_ref_pattern, ref_sequence, perl = TRUE)
  out <- sub("^.*_", "", ref_sequence)
  out[ok] <- sub(resfinder_ref_pattern, "\\3", ref_sequence[ok], perl = TRUE)
  out
}

strip_allele_suffix <- function(ref_sequence) {
  # Remove the variant/accession suffix but keep the ARG gene label.
  ref_sequence <- as.character(ref_sequence)
  ok <- grepl(resfinder_ref_pattern, ref_sequence, perl = TRUE)
  out <- ref_sequence
  out[ok] <- sub(resfinder_ref_pattern, "\\1", ref_sequence[ok], perl = TRUE)
  out
}

# Re-derive accession and gene_label from ref_sequence in saved objects that
# were created before the parsing fix (used when loading saved model results).
repair_feature_labels <- function(x) {
  if (is.null(x) || !"ref_sequence" %in% names(x)) return(x)
  x$accession <- extract_accession(x$ref_sequence)
  x$gene_label <- strip_allele_suffix(x$ref_sequence)
  x
}

standardize_names <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "_", trimws(x)))
}

check_required_columns <- function(header, required, file_label) {
  # Stop early if an input file lacks required columns.
  missing <- setdiff(required, header)
  if (length(missing)) {
    stop(file_label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

find_input_file <- function(input_dir, preferred, pattern) {
  # Find an input file by exact preferred names first, then by regex pattern.
  candidates <- file.path(input_dir, preferred)
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) return(normalizePath(existing[1], winslash = "/", mustWork = TRUE))

  found <- list.files(
    input_dir,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(found)) return(normalizePath(found[1], winslash = "/", mustWork = TRUE))

  stop(
    "Could not find input file in ", input_dir,
    ". Preferred names: ", paste(preferred, collapse = ", "),
    ". Regex pattern: ", pattern
  )
}

format_interval <- function(estimate, lower, upper, digits = 3) {
  if (!is.finite(estimate)) return("not estimable")
  if (!is.finite(lower) || !is.finite(upper)) {
    return(format(round(estimate, digits), nsmall = digits))
  }
  paste0(
    format(round(estimate, digits), nsmall = digits),
    " (repeated-CV interval ",
    format(round(lower, digits), nsmall = digits),
    " to ",
    format(round(upper, digits), nsmall = digits),
    ")"
  )
}

scale_train_test <- function(x_train, x_test) {
  # Scale test data using training means/SDs only.
  center <- colMeans(x_train, na.rm = TRUE)
  scale <- apply(x_train, 2, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(
    train = sweep(sweep(x_train, 2, center, "-"), 2, scale, "/"),
    test = sweep(sweep(x_test, 2, center, "-"), 2, scale, "/")
  )
}

classification_metrics <- function(actual, predicted, class_levels) {
  actual <- factor(actual, levels = class_levels)
  predicted <- factor(predicted, levels = class_levels)
  confusion <- table(actual, predicted)
  support <- rowSums(confusion)
  recall <- diag(confusion) / support
  precision <- diag(confusion) / colSums(confusion)
  f1 <- 2 * precision * recall / (precision + recall)
  recall[!is.finite(recall)] <- NA_real_
  precision[!is.finite(precision)] <- NA_real_
  f1[!is.finite(f1)] <- NA_real_
  data.frame(
    n = sum(confusion),
    accuracy = sum(diag(confusion)) / sum(confusion),
    balanced_accuracy = mean(recall, na.rm = TRUE),
    macro_f1 = mean(f1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}


# Prepare source-specific ARG matrices from metadata.tsv and ARG.tsv.
# The full ARG file should be named ARG.tsv. ARG_sub.tsv.txt can be used only for testing.

prepare_all_sources <- function(config) {
  log_message("Reading metadata: ", config$metadata_path)

  # Step 1: Check metadata columns before reading the full metadata file.
  metadata_header <- names(data.table::fread(config$metadata_path, sep = "\t", nrows = 0))
  required_metadata <- c(
    "run_accession", "sample_accession", "project_accession", "country", "continent",
    "host", "instrument_platform", "trimmed_fragments"
  )
  check_required_columns(metadata_header, required_metadata, "metadata.tsv")

  metadata <- data.table::fread(
    config$metadata_path,
    sep = "\t",
    select = required_metadata,
    showProgress = TRUE
  )

  # Step 2: Keep only the five One Health source labels specified in config$source_map.
  metadata <- merge(
    metadata,
    data.table::as.data.table(config$source_map),
    by.x = "host",
    by.y = "metadata_host",
    all = FALSE
  )
  metadata <- unique(metadata)

  eligible_samples <- unique(metadata[, .(
    source_key, display_name, sample_accession
  )])
  eligible_counts <- eligible_samples[, .(
    metadata_eligible_samples = data.table::uniqueN(sample_accession)
  ), by = .(source_key, display_name)]

  if (!nrow(eligible_samples)) {
    stop("No metadata samples matched config$source_map. Check exact host labels.")
  }

  log_message("Reading ARG mapping table: ", config$arg_path)

  # Step 3: Check ARG.tsv columns. ARG_sub.tsv.txt has the same structure but fewer rows.
  arg_header <- names(data.table::fread(config$arg_path, sep = "\t", nrows = 0))
  required_arg <- c(
    "run_accession", "sample_accession", "project_accession", "trimmed_fragments",
    "refSequence", "refSequence_length", "fragmentCountAln", "bacterial_fragment"
  )
  check_required_columns(arg_header, required_arg, "ARG.tsv")

  # Step 4: Read only required ARG columns to reduce memory use with the full ARG.tsv.
  arg <- data.table::fread(
    config$arg_path,
    sep = "\t",
    select = required_arg,
    showProgress = TRUE
  )

  # Step 5: Restrict ARG rows to metadata-eligible samples.
  arg <- arg[sample_accession %chin% eligible_samples$sample_accession]
  if (!nrow(arg)) {
    stop("No ARG rows overlapped with the selected metadata samples.")
  }

  numeric_columns <- c(
    "refSequence_length", "fragmentCountAln", "bacterial_fragment", "trimmed_fragments"
  )
  arg[, (numeric_columns) := lapply(
    .SD,
    function(value) suppressWarnings(as.numeric(value))
  ), .SDcols = numeric_columns]

  arg <- merge(
    arg,
    eligible_samples,
    by = "sample_accession",
    all.x = TRUE,
    allow.cartesian = TRUE
  )

  # Step 6: Keep valid ARG hits only.
  arg <- arg[
    !is.na(source_key) &
      is.finite(fragmentCountAln) & fragmentCountAln >= 0 &
      is.finite(refSequence_length) & refSequence_length > 0 &
      is.finite(bacterial_fragment) & bacterial_fragment > 0
  ]

  overlap_counts <- arg[, .(
    arg_positive_samples = data.table::uniqueN(sample_accession)
  ), by = .(source_key, display_name)]

  # Step 7: Calculate sample-level bacterial denominators.
  run_denominators <- unique(arg[, .(
    source_key, sample_accession, run_accession, bacterial_fragment
  )])
  run_denominators <- run_denominators[
    is.finite(bacterial_fragment) & bacterial_fragment > 0
  ]
  run_denominators <- run_denominators[, .(
    bacterial_fragments = max(bacterial_fragment)
  ), by = .(source_key, sample_accession, run_accession)]
  sample_denominators <- run_denominators[, .(
    bacterial_fragments = sum(bacterial_fragments)
  ), by = .(source_key, sample_accession)]

  # Step 8: Calculate gene-level abundance using two scales.
  # gene_fpm is used for load reconstruction because it matches the target unit.
  gene_long <- arg[, .(
    arg_fragments = sum(fragmentCountAln),
    length_normalized_fragments = sum(fragmentCountAln / refSequence_length)
  ), by = .(source_key, sample_accession, refSequence)]
  gene_long <- merge(
    gene_long,
    sample_denominators,
    by = c("source_key", "sample_accession"),
    all = FALSE
  )
  gene_long[, gene_fpm := 1e6 * arg_fragments / bacterial_fragments]
  gene_long[, gene_fpkm_bacterial := 1e9 * length_normalized_fragments / bacterial_fragments]

  # Step 9: Define normalized total ARG load as the regression target.
  load_long <- gene_long[, .(
    arg_fragments = sum(arg_fragments),
    bacterial_fragments = unique(bacterial_fragments)[1]
  ), by = .(source_key, sample_accession)]
  load_long[, arg_fragments_per_million := 1e6 * arg_fragments / bacterial_fragments]
  load_long[, log_arg_load := log1p(arg_fragments_per_million)]

  processed_dir <- file.path(config$output_dir, "processed")
  qc_dir <- file.path(config$output_dir, "qc")
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

  processed <- list()
  matrix_qc <- list()

  # Step 10: Build one matrix per source.
  for (source_key in config$source_map$source_key) {
    # Use current_source to avoid data.table name masking.
    # In data.table expressions, a column named source_key can mask the loop variable source_key.
    current_source <- as.character(source_key)
    display_name <- config$source_map$display_name[
      match(current_source, config$source_map$source_key)
    ]
    log_message("Preparing ", display_name, " matrix")

    # Subset strictly to the current source.
    source_long <- gene_long[gene_long[["source_key"]] == current_source, ]
    source_load <- load_long[load_long[["source_key"]] == current_source, ]

    if (!nrow(source_long) || !nrow(source_load)) {
      log_message("Skipping ", display_name, ": no usable ARG-positive samples")
      next
    }

    wide_fpm <- data.table::dcast(
      source_long,
      sample_accession ~ refSequence,
      value.var = "gene_fpm",
      fill = 0,
      fun.aggregate = sum
    )
    wide_fpkm <- data.table::dcast(
      source_long,
      sample_accession ~ refSequence,
      value.var = "gene_fpkm_bacterial",
      fill = 0,
      fun.aggregate = sum
    )

    original_features <- setdiff(names(wide_fpm), "sample_accession")
    wide_fpkm <- wide_fpkm[, c("sample_accession", original_features), with = FALSE]

    wide_fpm <- merge(
      source_load[, .(sample_accession, log_arg_load, arg_fragments_per_million)],
      wide_fpm,
      by = "sample_accession",
      all = FALSE
    )
    wide_fpkm <- merge(
      source_load[, .(sample_accession)],
      wide_fpkm,
      by = "sample_accession",
      all = FALSE
    )

    abundance_fpm <- as.matrix(wide_fpm[, ..original_features])
    abundance_fpkm <- as.matrix(wide_fpkm[, ..original_features])
    storage.mode(abundance_fpm) <- "double"
    storage.mode(abundance_fpkm) <- "double"

    valid_rows <- is.finite(wide_fpm$log_arg_load) &
      is.finite(rowSums(abundance_fpm)) & rowSums(abundance_fpm) > 0 &
      is.finite(rowSums(abundance_fpkm)) & rowSums(abundance_fpkm) > 0

    wide_fpm <- wide_fpm[valid_rows]
    abundance_fpm <- abundance_fpm[valid_rows, , drop = FALSE]
    abundance_fpkm <- abundance_fpkm[valid_rows, , drop = FALSE]

    # Step 11: Model predictors use absolute log-normalized gene load.
    x_log_fpm <- log1p(abundance_fpm)

    # Step 12: Network predictors use compositional CLR-transformed abundance.
    row_total <- rowSums(abundance_fpkm)
    relative_abundance <- abundance_fpkm / row_total
    log_relative <- log(relative_abundance + config$clr_pseudocount)
    x_clr <- log_relative - rowMeans(log_relative)

    safe_names <- safe_gene_names(original_features)
    colnames(x_log_fpm) <- safe_names
    colnames(x_clr) <- safe_names
    colnames(abundance_fpm) <- safe_names
    colnames(abundance_fpkm) <- safe_names
    rownames(x_log_fpm) <- wide_fpm$sample_accession
    rownames(x_clr) <- wide_fpm$sample_accession
    rownames(abundance_fpm) <- wide_fpm$sample_accession
    rownames(abundance_fpkm) <- wide_fpm$sample_accession

    mapping <- data.frame(
      safe_name = safe_names,
      ref_sequence = original_features,
      accession = extract_accession(original_features),
      gene_label = strip_allele_suffix(original_features),
      stringsAsFactors = FALSE
    )

    sample_projects <- unique(metadata[source_key == current_source, .(sample_accession, project_accession)])
    sample_projects <- sample_projects[!duplicated(sample_accession)]

    object <- list(
      source_key = current_source,
      display_name = display_name,
      sample_accession = wide_fpm$sample_accession,
      project_accession = sample_projects$project_accession[
        match(wide_fpm$sample_accession, sample_projects$sample_accession)],
      x_log_fpm = x_log_fpm,
      x_clr = x_clr,
      abundance_fpm = abundance_fpm,
      abundance_fpkm_bacterial = abundance_fpkm,
      log_arg_load = wide_fpm$log_arg_load,
      arg_fragments_per_million = wide_fpm$arg_fragments_per_million,
      feature_mapping = mapping,
      normalization = list(
        model_features = "log1p(gene ARG fragments per million bacterial fragments)",
        network_features = "length-normalized abundance converted to relative abundance and CLR",
        target = "log1p(total ARG fragments per million bacterial fragments)",
        clr_pseudocount = config$clr_pseudocount
      )
    )

    saveRDS(object, file.path(processed_dir, paste0(current_source, "_matrix.rds")))
    processed[[current_source]] <- object

    matrix_qc[[current_source]] <- data.frame(
      source_key = current_source,
      display_name = display_name,
      modeled_samples = nrow(x_log_fpm),
      modeled_features = ncol(x_log_fpm),
      median_arg_fpm = stats::median(wide_fpm$arg_fragments_per_million),
      median_log_arg_load = stats::median(wide_fpm$log_arg_load),
      stringsAsFactors = FALSE
    )
  }

  if (!length(processed)) stop("No source matrices were created.")

  # Step 13: Write sample-flow/QC outputs.
  qc <- Reduce(
    function(x, y) merge(x, y, by = c("source_key", "display_name"), all = TRUE),
    list(
      as.data.frame(eligible_counts),
      as.data.frame(overlap_counts),
      do.call(rbind, matrix_qc)
    )
  )
  qc$arg_positive_samples[is.na(qc$arg_positive_samples)] <- 0
  qc$modeled_samples[is.na(qc$modeled_samples)] <- 0
  qc$modeled_features[is.na(qc$modeled_features)] <- 0
  qc$metadata_without_arg_hits <- with(
    qc,
    metadata_eligible_samples - arg_positive_samples
  )
  qc <- qc[match(config$source_map$source_key, qc$source_key), ]
  write_csv_atomic(qc, file.path(qc_dir, "sample_flow.csv"))

  metadata_qc <- metadata[, .(
    runs = data.table::uniqueN(run_accession),
    samples = data.table::uniqueN(sample_accession),
    projects = data.table::uniqueN(project_accession),
    countries_recorded = data.table::uniqueN(country[country != "NULL" & !is.na(country)]),
    missing_country_samples = data.table::uniqueN(
      sample_accession[country == "NULL" | is.na(country)]
    )
  ), by = .(source_key, display_name)]
  write_csv_atomic(metadata_qc, file.path(qc_dir, "metadata_completeness.csv"))

  list(processed = processed, sample_flow = qc, metadata_qc = metadata_qc)
}


# Leakage-controlled regression and source-classification models.
# not CLR predictors, because the target is total normalized ARG load.

screen_and_select_features <- function(
    x_train,
    abundance_train,
    y_train,
    config,
    seed) {
  # Step 1: Select features using training data only.
  set.seed(seed)

  minimum_present <- max(
    config$min_prevalence_samples,
    ceiling(config$prevalence_fraction * nrow(x_train))
  )
  prevalence <- colSums(abundance_train > 0)
  variances <- apply(x_train, 2, stats::var)
  eligible <- names(prevalence)[
    prevalence >= minimum_present &
      is.finite(variances) &
      variances > .Machine$double.eps
  ]
  if (!length(eligible)) {
    stop("No features passed the within-fold prevalence and variance filter.")
  }

  # Step 2: Rank features by univariate monotonic association with the target.
  correlations <- vapply(
    eligible,
    function(feature) {
      value <- suppressWarnings(stats::cor(
        x_train[, feature],
        y_train,
        method = "spearman",
        use = "complete.obs"
      ))
      ifelse(is.finite(value), abs(value), 0)
    },
    numeric(1)
  )
  correlation_rank <- names(sort(correlations, decreasing = TRUE))

  # Step 3: Use elastic-net screening to handle correlated ARG predictors.
  nfolds <- min(5L, max(3L, floor(nrow(x_train) / 8L)))
  elastic_features <- character()
  elastic_fit <- tryCatch(
    glmnet::cv.glmnet(
      x = x_train[, eligible, drop = FALSE],
      y = y_train,
      family = "gaussian",
      alpha = config$screening_alpha,
      nfolds = nfolds,
      standardize = TRUE,
      type.measure = "mse"
    ),
    error = function(error) NULL
  )
  if (!is.null(elastic_fit)) {
    coefficients <- stats::coef(elastic_fit, s = "lambda.1se")
    elastic_features <- rownames(coefficients)[as.vector(coefficients != 0)]
    elastic_features <- setdiff(elastic_features, "(Intercept)")
    if (!length(elastic_features)) {
      coefficients <- stats::coef(elastic_fit, s = "lambda.min")
      elastic_features <- rownames(coefficients)[as.vector(coefficients != 0)]
      elastic_features <- setdiff(elastic_features, "(Intercept)")
    }
  }

  # Step 4: Build a capped candidate set for Boruta.
  candidates <- unique(c(elastic_features, correlation_rank))
  candidates <- candidates[candidates %in% eligible]
  candidates <- head(candidates, config$max_candidates)
  if (length(candidates) < min(config$min_candidates, length(eligible))) {
    candidates <- head(correlation_rank, min(config$min_candidates, length(eligible)))
  }

  # Step 5: Run Boruta only within the outer training fold.
  # If runtime is high, set config$use_boruta <- FALSE; the workflow remains leakage-controlled.
  selected_runs <- list()
  if (isTRUE(config$use_boruta) && length(candidates) >= 2L) {
    boruta_data <- data.frame(
      x_train[, candidates, drop = FALSE],
      burden = y_train,
      check.names = FALSE
    )
    for (iteration in seq_len(config$boruta_repeats)) {
      set.seed(seed + iteration * 1009L)
      selected_runs[[iteration]] <- tryCatch(
        {
          fit <- Boruta::Boruta(
            burden ~ .,
            data = boruta_data,
            doTrace = 0,
            maxRuns = config$boruta_max_runs
          )
          selected <- Boruta::getSelectedAttributes(fit, withTentative = FALSE)
          if (!length(selected)) {
            fixed <- Boruta::TentativeRoughFix(fit)
            selected <- Boruta::getSelectedAttributes(fixed, withTentative = FALSE)
          }
          selected
        },
        error = function(error) {
          log_message("Boruta warning: ", conditionMessage(error))
          character()
        }
      )
    }
  }

  counts <- sort(table(unlist(selected_runs)), decreasing = TRUE)
  threshold <- max(1L, ceiling(config$boruta_repeats / 2))
  selected <- names(counts[counts >= threshold])

  # Step 6: Fallback to elastic-net/correlation candidates if Boruta selects none.
  if (!length(selected)) {
    selected <- head(candidates, min(config$fallback_features, length(candidates)))
  }

  selected <- selected[selected %in% colnames(x_train)]
  selected <- head(
    selected[order(
      -ifelse(selected %in% names(counts), as.numeric(counts[selected]), 0),
      -correlations[selected],
      na.last = TRUE
    )],
    config$max_selected_features
  )

  list(
    selected = selected,
    eligible = eligible,
    candidates = candidates,
    boruta_counts = counts,
    correlations = correlations
  )
}

fit_predict_elastic_net <- function(x_train, y_train, x_test, config, seed) {
  # Step 1: Tune alpha by inner CV and predict the outer test fold.
  set.seed(seed)
  nfolds <- min(5L, max(3L, floor(nrow(x_train) / 8L)))
  # The same inner folds are used for every alpha so alphas are compared fairly.
  foldid <- quantile_stratified_folds(y_train, nfolds, seed)
  best <- NULL
  for (alpha in config$alpha_grid) {
    fit <- tryCatch(
      glmnet::cv.glmnet(
        x = x_train,
        y = y_train,
        family = "gaussian",
        alpha = alpha,
        foldid = foldid,
        standardize = TRUE,
        type.measure = "mse"
      ),
      error = function(error) NULL
    )
    if (is.null(fit)) next
    score <- min(fit$cvm, na.rm = TRUE)
    if (is.null(best) || score < best$score) {
      best <- list(fit = fit, alpha = alpha, score = score)
    }
  }
  if (is.null(best)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(stats::predict(best$fit, newx = x_test, s = "lambda.min"))
}

fit_predict_random_forest <- function(x_train, y_train, x_test, config, seed) {
  # Step 1: Light OOB tuning of mtry and min.node.size within the training fold.
  set.seed(seed)
  p <- ncol(x_train)
  mtry_grid <- unique(pmax(1L, pmin(p, round(c(sqrt(p), p / 4, p / 3, p / 2)))))
  node_grid <- config$rf_min_node_grid
  best <- NULL
  for (mtry in mtry_grid) {
    for (node in node_grid) {
      fit <- tryCatch(
        ranger::ranger(
          x = data.frame(x_train, check.names = FALSE),
          y = y_train,
          num.trees = config$rf_trees,
          mtry = mtry,
          min.node.size = node,
          importance = "none",
          seed = seed,
          num.threads = config$num_threads
        ),
        error = function(error) NULL
      )
      if (is.null(fit)) next
      score <- fit$prediction.error
      if (is.null(best) || score < best$score) {
        best <- list(fit = fit, score = score)
      }
    }
  }
  if (is.null(best)) return(rep(NA_real_, nrow(x_test)))
  stats::predict(best$fit, data = data.frame(x_test, check.names = FALSE))$predictions
}

fit_predict_svm <- function(x_train, y_train, x_test, config, seed) {
  # Step 1: Radial SVM with training-fold scaling.
  set.seed(seed)
  scaled <- scale_train_test(x_train, x_test)
  p <- ncol(x_train)
  gamma_grid <- unique(c(1 / p, 2 / p))
  cost_grid <- config$svm_cost_grid
  best <- NULL

  # Step 2: Inner validation split from training data only.
  inner_folds <- quantile_stratified_folds(y_train, k = min(3L, length(y_train)), seed = seed)
  inner_test <- which(inner_folds == 1L)
  inner_train <- setdiff(seq_along(y_train), inner_test)

  for (cost in cost_grid) {
    for (gamma in gamma_grid) {
      fit <- tryCatch(
        e1071::svm(
          x = scaled$train[inner_train, , drop = FALSE],
          y = y_train[inner_train],
          type = "eps-regression",
          kernel = "radial",
          cost = cost,
          gamma = gamma,
          scale = FALSE
        ),
        error = function(error) NULL
      )
      if (is.null(fit)) next
      pred <- as.numeric(stats::predict(fit, scaled$train[inner_test, , drop = FALSE]))
      score <- rmse(y_train[inner_test], pred)
      if (is.null(best) || score < best$score) {
        best <- list(cost = cost, gamma = gamma, score = score)
      }
    }
  }
  if (is.null(best)) return(rep(NA_real_, nrow(x_test)))

  final_fit <- tryCatch(
    e1071::svm(
      x = scaled$train,
      y = y_train,
      type = "eps-regression",
      kernel = "radial",
      cost = best$cost,
      gamma = best$gamma,
      scale = FALSE
    ),
    error = function(error) NULL
  )
  if (is.null(final_fit)) return(rep(NA_real_, nrow(x_test)))
  as.numeric(stats::predict(final_fit, scaled$test))
}

fit_predict_xgboost <- function(x_train, y_train, x_test, config, seed) {
  # Step 1: Make sure XGBoost receives clean numeric matrices.
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  storage.mode(x_train) <- "double"
  storage.mode(x_test) <- "double"
  x_train[!is.finite(x_train)] <- 0
  x_test[!is.finite(x_test)] <- 0
  y_train <- as.numeric(y_train)

  # Step 2: If the response is almost constant in a fold, use a mean predictor.
  if (length(unique(round(y_train, 8))) < 2L || nrow(x_train) < 10L || ncol(x_train) < 1L) {
    return(rep(mean(y_train, na.rm = TRUE), nrow(x_test)))
  }

  set.seed(seed)
  train_matrix <- xgboost::xgb.DMatrix(data = x_train, label = y_train)

  params <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    eta = config$xgb_eta,
    max_depth = config$xgb_max_depth,
    min_child_weight = config$xgb_min_child_weight,
    subsample = config$xgb_subsample,
    colsample_bytree = config$xgb_colsample_bytree,
    nthread = config$num_threads
  )

  # Step 3: Tune nrounds only when the fold is large enough for stable inner CV.
  best_nrounds <- min(100L, config$xgb_nrounds)
  nfold <- min(5L, max(3L, floor(nrow(x_train) / 8L)))

  if (nrow(x_train) >= 30L && ncol(x_train) >= 2L) {
    cv <- tryCatch(
      xgboost::xgb.cv(
        params = params,
        data = train_matrix,
        nrounds = config$xgb_nrounds,
        nfold = nfold,
        early_stopping_rounds = config$xgb_early_stopping_rounds,
        verbose = 0
      ),
      error = function(error) NULL
    )

    # Step 4: Number of rounds = round with the lowest inner-CV RMSE
    # (read from the evaluation log, which works across xgboost versions).
    candidate <- if (is.null(cv)) NA_integer_ else
      xgb_best_round_from_log(cv$evaluation_log, "^test.*rmse.*mean$")
    if (length(candidate) == 1L && !is.na(candidate) && candidate >= 1L) {
      best_nrounds <- candidate
    }
  }

  if (length(best_nrounds) != 1L || !is.finite(best_nrounds) || is.na(best_nrounds) || best_nrounds < 1L) {
    best_nrounds <- min(100L, config$xgb_nrounds)
  }

  # Step 5: Fit final XGBoost model, falling back to the training mean if it fails.
  model <- tryCatch(
    xgboost::xgb.train(
      params = params,
      data = train_matrix,
      nrounds = best_nrounds,
      verbose = 0
    ),
    error = function(error) NULL
  )

  if (is.null(model)) {
    return(rep(mean(y_train, na.rm = TRUE), nrow(x_test)))
  }

  prediction <- as.numeric(stats::predict(model, xgboost::xgb.DMatrix(data = x_test)))
  if (any(!is.finite(prediction))) {
    prediction[!is.finite(prediction)] <- mean(y_train, na.rm = TRUE)
  }
  prediction
}

predict_regression_models <- function(
    x_train,
    y_train,
    x_test,
    abundance_test,
    selected,
    config,
    seed,
    abundance_train = NULL) {
  # Step 1: Restrict predictors to fold-selected indicators.
  x_train <- x_train[, selected, drop = FALSE]
  x_test <- x_test[, selected, drop = FALSE]
  abundance_test <- abundance_test[, selected, drop = FALSE]

  # Indicator Sum: log1p(summed normalized abundance of the selected genes),
  # linearly calibrated to the target on the TRAINING fold (a simple panel score).
  sum_test <- log1p(rowSums(abundance_test))
  indicator_sum <- if (!is.null(abundance_train)) {
    sum_train <- log1p(rowSums(abundance_train[, selected, drop = FALSE]))
    calibration <- stats::lm(y ~ s, data = data.frame(y = y_train, s = sum_train))
    as.numeric(stats::predict(calibration, newdata = data.frame(s = sum_test)))
  } else sum_test

  predictions <- list(
    Baseline = rep(mean(y_train), nrow(x_test)),
    `Indicator Sum` = indicator_sum
  )

  predictions[["Elastic Net"]] <- tryCatch(
    fit_predict_elastic_net(x_train, y_train, x_test, config, seed),
    error = function(error) {
      log_message("Elastic Net warning: ", conditionMessage(error))
      rep(NA_real_, nrow(x_test))
    }
  )

  predictions[["Random Forest"]] <- tryCatch(
    fit_predict_random_forest(x_train, y_train, x_test, config, seed),
    error = function(error) {
      log_message("Random Forest warning: ", conditionMessage(error))
      rep(NA_real_, nrow(x_test))
    }
  )

  predictions[["SVM"]] <- tryCatch(
    fit_predict_svm(x_train, y_train, x_test, config, seed),
    error = function(error) {
      log_message("SVM warning: ", conditionMessage(error))
      rep(NA_real_, nrow(x_test))
    }
  )

  predictions[["XGBoost"]] <- tryCatch(
    fit_predict_xgboost(x_train, y_train, x_test, config, seed),
    error = function(error) {
      log_message("XGBoost warning: ", conditionMessage(error))
      rep(NA_real_, nrow(x_test))
    }
  )

  predictions
}

run_repeated_nested_regression <- function(processed, config) {
  # Step 1: Run source-specific leakage-controlled regression.
  output_dir <- file.path(config$output_dir, "models_regression")
  figure_dir <- file.path(config$output_dir, "figures")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  all_predictions <- list()
  all_fold_metrics <- list()
  all_features <- list()
  cursor <- 0L

  for (source_key in names(processed)) {
    object <- processed[[source_key]]
    x <- object$x_log_fpm
    abundance <- object$abundance_fpm
    y <- object$log_arg_load
    n_folds <- min(config$outer_folds, max(2L, floor(length(y) / 8L)))
    project_blocked <- identical(config$cv_grouping, "project")
    if (project_blocked) {
      if (is.null(object$project_accession)) stop("project_accession missing from processed data")
      n_folds <- min(n_folds, length(unique(object$project_accession)))
    }

    log_message(
      "Regression modeling ", object$display_name,
      if (project_blocked) " [project-blocked CV]" else "",
      " (n=", length(y), ", p=", ncol(x), ", ",
      n_folds, " folds x ", config$outer_repeats, " repeats)"
    )

    for (repetition in seq_len(config$outer_repeats)) {
      # Seeds depend on the source's position in the full source list, so a
      # source gives identical results whether run alone or with the others.
      source_order <- if (is.null(config$source_order)) names(processed) else config$source_order
      repetition_seed <- config$seed + match(source_key, source_order) * 100000L + repetition * 1000L
      fold_ids <- if (project_blocked) {
        group_folds(object$project_accession, n_folds, repetition_seed)
      } else {
        quantile_stratified_folds(y, n_folds, repetition_seed)
      }

      for (fold in seq_len(n_folds)) {
        fold_seed <- repetition_seed + fold * 17L
        test_indices <- which(fold_ids == fold)
        train_indices <- setdiff(seq_along(y), test_indices)

        selection <- screen_and_select_features(
          x_train = x[train_indices, , drop = FALSE],
          abundance_train = abundance[train_indices, , drop = FALSE],
          y_train = y[train_indices],
          config = config,
          seed = fold_seed
        )

        model_predictions <- predict_regression_models(
          x_train = x[train_indices, , drop = FALSE],
          y_train = y[train_indices],
          x_test = x[test_indices, , drop = FALSE],
          abundance_test = abundance[test_indices, , drop = FALSE],
          selected = selection$selected,
          config = config,
          seed = fold_seed,
          abundance_train = abundance[train_indices, , drop = FALSE]
        )

        cursor <- cursor + 1L
        fold_prediction <- do.call(
          rbind,
          lapply(names(model_predictions), function(model_name) {
            data.frame(
              source_key = source_key,
              display_name = object$display_name,
              repetition = repetition,
              fold = fold,
              sample_accession = object$sample_accession[test_indices],
              model = model_name,
              actual = y[test_indices],
              predicted = model_predictions[[model_name]],
              n_selected_features = length(selection$selected),
              stringsAsFactors = FALSE
            )
          })
        )
        all_predictions[[cursor]] <- fold_prediction

        all_fold_metrics[[cursor]] <- do.call(
          rbind,
          lapply(split(fold_prediction, fold_prediction$model), function(rows) {
            valid <- is.finite(rows$actual) & is.finite(rows$predicted)
            data.frame(
              source_key = source_key,
              display_name = object$display_name,
              repetition = repetition,
              fold = fold,
              model = rows$model[1],
              n_validation = sum(valid),
              n_selected_features = length(selection$selected),
              RMSE = rmse(rows$actual[valid], rows$predicted[valid]),
              MAE = mae(rows$actual[valid], rows$predicted[valid]),
              R_squared = r_squared(rows$actual[valid], rows$predicted[valid]),
              stringsAsFactors = FALSE
            )
          })
        )

        mapping <- object$feature_mapping[
          match(selection$selected, object$feature_mapping$safe_name),
          ,
          drop = FALSE
        ]
        mapping$source_key <- source_key
        mapping$display_name <- object$display_name
        mapping$repetition <- repetition
        mapping$fold <- fold
        mapping$fold_seed <- fold_seed
        all_features[[cursor]] <- mapping

        log_message(
          object$display_name, " repetition ", repetition,
          " fold ", fold, ": selected ", length(selection$selected), " features"
        )
      }
    }
  }

  predictions <- data.table::rbindlist(all_predictions, fill = TRUE)
  fold_metrics <- data.table::rbindlist(all_fold_metrics, fill = TRUE)
  selected_features <- data.table::rbindlist(all_features, fill = TRUE)

  write_csv_atomic(predictions, file.path(output_dir, "heldout_predictions.csv"))
  write_csv_atomic(fold_metrics, file.path(output_dir, "fold_metrics.csv"))
  write_csv_atomic(selected_features, file.path(output_dir, "selected_features_by_fold.csv"))

  repetition_metrics <- predictions[, {
    valid <- is.finite(actual) & is.finite(predicted)
    .(
      n = sum(valid),
      RMSE = rmse(actual[valid], predicted[valid]),
      MAE = mae(actual[valid], predicted[valid]),
      R_squared = r_squared(actual[valid], predicted[valid])
    )
  }, by = .(source_key, display_name, repetition, model)]
  write_csv_atomic(repetition_metrics, file.path(output_dir, "repetition_metrics.csv"))

  summary_rows <- list()
  summary_cursor <- 0L
  for (metric in c("RMSE", "MAE", "R_squared")) {
    grouped <- split(
      repetition_metrics,
      interaction(repetition_metrics$source_key, repetition_metrics$model, drop = TRUE)
    )
    for (rows in grouped) {
      summary_cursor <- summary_cursor + 1L
      interval <- mean_ci_by_repetition(rows[[metric]])
      summary_rows[[summary_cursor]] <- data.frame(
        source_key = rows$source_key[1],
        display_name = rows$display_name[1],
        model = rows$model[1],
        metric = metric,
        estimate = interval[["mean"]],
        lower = interval[["lower"]],
        upper = interval[["upper"]],
        repetitions = interval[["n"]],
        stringsAsFactors = FALSE
      )
    }
  }
  model_summary <- data.table::rbindlist(summary_rows)
  write_csv_atomic(model_summary, file.path(output_dir, "model_summary.csv"))

  fold_denominators <- unique(selected_features[, .(
    source_key, repetition, fold
  )])[, .(outer_folds = .N), by = source_key]
  feature_frequency <- selected_features[, .(
    selection_count = .N
  ), by = .(
    source_key, display_name, safe_name, ref_sequence, accession, gene_label
  )]
  feature_frequency <- merge(feature_frequency, fold_denominators, by = "source_key")
  feature_frequency[, selection_frequency := selection_count / outer_folds]
  data.table::setorder(feature_frequency, source_key, -selection_frequency, gene_label)
  write_csv_atomic(feature_frequency, file.path(output_dir, "feature_selection_frequency.csv"))

  stable_features <- feature_frequency[selection_frequency >= config$stable_frequency]
  if (!nrow(stable_features)) {
    stable_features <- feature_frequency[, head(.SD, min(10L, .N)), by = source_key]
  }
  write_csv_atomic(stable_features, file.path(output_dir, "stable_indicator_features.csv"))

  # Step 2: Plot regression performance.
  plot_data <- model_summary[metric == "R_squared" & model != "Baseline"]
  performance_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = display_name, y = estimate, color = model, group = model)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.15,
      position = ggplot2::position_dodge(width = 0.6),
      na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 2.4, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::labs(
      x = NULL,
      y = expression(held-out ~ R^2),
      color = NULL,
      title = "Leakage-controlled reconstruction of normalized ARG load",
      subtitle = "Predictors: log gene FPM; feature selection repeated inside every training fold"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    )
  ggplot2::ggsave(
    file.path(figure_dir, "regression_model_performance.png"),
    performance_plot,
    width = 8.8,
    height = 5.2,
    dpi = 320,
    bg = "white"
  )

  list(
    predictions = predictions,
    fold_metrics = fold_metrics,
    repetition_metrics = repetition_metrics,
    model_summary = model_summary,
    feature_frequency = feature_frequency,
    stable_features = stable_features
  )
}

build_source_classification_matrix <- function(processed) {
  # Step 1: Combine source-specific matrices into a single union-feature matrix.
  all_features <- sort(unique(unlist(lapply(processed, function(object) colnames(object$x_log_fpm)))))
  total_samples <- sum(vapply(processed, function(object) nrow(object$x_log_fpm), integer(1)))
  x <- matrix(0, nrow = total_samples, ncol = length(all_features))
  colnames(x) <- all_features

  sample_accession <- character(total_samples)
  project_accession <- rep(NA_character_, total_samples)
  source_key <- character(total_samples)
  display_name <- character(total_samples)
  cursor <- 0L

  mapping_list <- list()
  for (key in names(processed)) {
    object <- processed[[key]]
    rows <- seq_len(nrow(object$x_log_fpm)) + cursor
    x[rows, colnames(object$x_log_fpm)] <- object$x_log_fpm
    sample_accession[rows] <- object$sample_accession
    if (!is.null(object$project_accession)) project_accession[rows] <- object$project_accession
    source_key[rows] <- object$source_key
    display_name[rows] <- object$display_name
    cursor <- cursor + nrow(object$x_log_fpm)
    mapping_list[[key]] <- object$feature_mapping
  }

  mapping <- unique(do.call(rbind, mapping_list))
  list(
    x = x,
    sample_accession = sample_accession,
    project_accession = project_accession,
    source_key = source_key,
    display_name = display_name,
    feature_mapping = mapping
  )
}

# ============================================================
# Nested, leakage-controlled five-class source classification.
#   Outer CV: unbiased, sample-level held-out performance estimates.
#   Inner CV: repeat feature screening INSIDE every inner training fold,
#             then tune elastic net, random forest and XGBoost.
# NOTE: This replaces ONLY the classification functions from v4.
#       The source-specific regression workflow above is unchanged.
# ============================================================

classification_config <- function(config, key, default) {
  value <- config[[key]]
  if (is.null(value) || length(value) == 0L) default else value
}

classification_balanced_accuracy <- function(actual, predicted, levels) {
  actual <- as.character(actual)
  predicted <- as.character(predicted)
  recalls <- vapply(levels, function(class_name) {
    present <- which(actual == class_name)
    if (!length(present) || anyNA(predicted[present])) return(NA_real_)
    mean(predicted[present] == class_name)
  }, numeric(1))
  if (anyNA(recalls)) return(NA_real_)
  mean(recalls)
}

classification_inner_folds <- function(y, config, seed) {
  y <- factor(y)
  max_folds <- min(as.integer(classification_config(config, "classification_inner_folds", 3L)),
                   min(table(y)))
  if (max_folds < 3L) {
    stop("Inner cross-validation requires at least three samples per source class.")
  }
  class_stratified_folds(y, max_folds, seed)
}

select_classification_features <- function(x_train, y_train, config, seed) {
  # Fold-specific, outcome-aware screening (never uses an inner/outer test fold).
  set.seed(seed)
  max_features <- max(1L, as.integer(classification_config(config, "classification_max_features", 150L)))
  fallback_features <- max(1L, as.integer(classification_config(config, "classification_fallback_features", 25L)))
  minimum_present <- max(
    as.integer(classification_config(config, "min_prevalence_samples", 2L)),
    ceiling(classification_config(config, "classification_prevalence_fraction", 0.05) * nrow(x_train))
  )
  prevalence <- colSums(x_train > 0)
  variances <- apply(x_train, 2, stats::var)
  eligible <- names(prevalence)[
    prevalence >= minimum_present & is.finite(variances) & variances > .Machine$double.eps
  ]
  if (!length(eligible)) stop("No classification ARG features passed prevalence/variance filtering.")

  # Multinomial elastic-net selection at lambda.1se: selected count is adaptive.
  # 'classification_max_features' is a CAP, not a target.
  y_train <- factor(y_train)
  nfolds <- min(5L, min(table(y_train)))
  if (nfolds < 3L) stop("Too few training samples in a source class for glmnet screening.")
  screening_alpha <- classification_config(config, "screening_alpha", 0.5)
  cv_fit <- tryCatch(
    glmnet::cv.glmnet(
      x = as.matrix(x_train[, eligible, drop = FALSE]),
      y = y_train,
      family = "multinomial",
      alpha = screening_alpha,
      nfolds = nfolds,
      standardize = TRUE,
      type.measure = "class"
    ),
    error = function(e) {
      log_message("Classification elastic-net screening warning: ", conditionMessage(e))
      NULL
    }
  )
  selected <- character()
  if (!is.null(cv_fit)) {
    coefs <- stats::coef(cv_fit, s = "lambda.1se")
    selected <- unique(unlist(lapply(coefs, function(m) {
      rownames(m)[as.vector(m != 0)]
    }), use.names = FALSE))
    selected <- intersect(setdiff(selected, "(Intercept)"), eligible)
  }

  # ANOVA ranks are ONLY used for fallback or to trim an oversize elastic-net
  # set. The original v4 appended the entire ANOVA ranking, filling to 150.
  n_elastic <- length(selected)
  did_fallback <- !n_elastic
  did_cap <- n_elastic > max_features
  if (!length(selected) || length(selected) > max_features) {
    ranking_pool <- if (length(selected)) selected else eligible
    f_scores <- vapply(ranking_pool, function(feature) {
      tab <- tryCatch(
        summary(stats::aov(x_train[, feature] ~ y_train))[[1L]],
        error = function(e) NULL
      )
      if (is.null(tab)) return(0)
      value <- tab$`F value`[1L]
      ifelse(is.finite(value), value, 0)
    }, numeric(1))
    ranking <- names(sort(f_scores, decreasing = TRUE))
    if (!length(selected)) {
      selected <- head(ranking, min(fallback_features, max_features, length(ranking)))
    } else {
      selected <- head(ranking, max_features)
    }
  }
  if (!length(selected)) stop("Classification screening returned no features.")
  # glmnet requires at least two predictors; add one training-ranked feature
  # if elastic net returns only one (never pad to the old 150-feature target).
  if (length(selected) == 1L && length(eligible) >= 2L && max_features >= 2L) {
    remaining <- setdiff(eligible, selected)
    extra_scores <- vapply(remaining, function(feature) {
      tab <- tryCatch(summary(stats::aov(x_train[, feature] ~ y_train))[[1L]],
                      error = function(e) NULL)
      if (is.null(tab)) return(0)
      score <- tab$`F value`[1L]
      ifelse(is.finite(score), score, 0)
    }, numeric(1))
    selected <- c(selected, names(extra_scores)[which.max(extra_scores)])
  }
  if (length(selected) < 2L) {
    stop("Classification needs at least two eligible ARG predictors for glmnet.")
  }
  attr(selected, "diagnostics") <- list(
    n_eligible = length(eligible), n_elastic = n_elastic,
    n_final = length(selected), feature_cap = max_features,
    fallback_used = did_fallback, capped = did_cap
  )
  selected
}

classification_rf_mtry <- function(mode, p) {
  p <- as.integer(p)
  if (mode == "sqrt") return(max(1L, min(p, round(sqrt(p)))))
  if (mode == "third") return(max(1L, min(p, round(p / 3))))
  if (mode == "half") return(max(1L, min(p, round(p / 2))))
  stop("Unknown random-forest mtry mode: ", mode)
}

classification_xgb_classes <- function(model, x_test, class_levels) {
  # xgboost R 3.x returns an [n_samples, n_classes] probability matrix.
  # Do NOT reshape an existing matrix with byrow=TRUE: it scrambles scores.
  prob <- stats::predict(model, xgboost::xgb.DMatrix(as.matrix(x_test)))
  if (is.null(dim(prob))) {
    prob <- matrix(prob, nrow = nrow(x_test), ncol = length(class_levels), byrow = TRUE)
  }
  if (!is.matrix(prob) || nrow(prob) != nrow(x_test) ||
      ncol(prob) != length(class_levels)) {
    stop("Unexpected XGBoost probability dimensions; expected samples x classes.")
  }
  class_levels[max.col(prob, ties.method = "first")]
}

classification_xgb_grid <- function(config) {
  custom <- config[["classification_xgb_grid"]]
  if (!is.null(custom)) {
    required <- c("eta", "max_depth", "min_child_weight", "subsample", "colsample_bytree")
    if (!is.data.frame(custom) || !all(required %in% names(custom)) || !nrow(custom)) {
      stop("classification_xgb_grid must be a nonempty data.frame with: ",
           paste(required, collapse = ", "))
    }
    return(custom)
  }
  # Compact, predeclared grid: avoids an impractically large Cartesian search.
  data.frame(
    eta = c(0.03, 0.10, 0.05, 0.15, 0.05, 0.10, 0.03, 0.15),
    max_depth = c(2L, 2L, 4L, 4L, 3L, 3L, 4L, 2L),
    min_child_weight = c(1, 5, 1, 5, 5, 1, 1, 5),
    subsample = c(0.8, 1.0, 0.8, 1.0, 1.0, 0.8, 1.0, 0.8),
    colsample_bytree = c(0.8, 1.0, 1.0, 0.8, 0.8, 1.0, 0.8, 1.0)
  )
}

classification_xgb_params <- function(params, class_levels, config) {
  list(
    objective = "multi:softprob", eval_metric = "mlogloss",
    num_class = length(class_levels),
    eta = as.numeric(params$eta),
    max_depth = as.integer(params$max_depth),
    min_child_weight = as.numeric(params$min_child_weight),
    subsample = as.numeric(params$subsample),
    colsample_bytree = as.numeric(params$colsample_bytree),
    nthread = as.integer(classification_config(config, "num_threads", 1L))
  )
}

classification_parameter_row <- function(model, inner_score, alpha = NA_real_,
                                         lambda_choice = NA_character_, mtry_mode = NA_character_,
                                         min_node = NA_integer_, eta = NA_real_, max_depth = NA_integer_,
                                         min_child_weight = NA_real_, subsample = NA_real_,
                                         colsample_bytree = NA_real_, nrounds = NA_integer_) {
  data.frame(
    model = model, inner_cv_balanced_accuracy = inner_score,
    alpha = alpha, lambda_choice = lambda_choice,
    mtry_mode = mtry_mode, min_node = min_node,
    eta = eta, max_depth = max_depth, min_child_weight = min_child_weight,
    subsample = subsample, colsample_bytree = colsample_bytree,
    nrounds = nrounds, stringsAsFactors = FALSE
  )
}

tune_classification_hyperparameters <- function(x_train, y_train, class_levels, config, seed) {
  # Nested design: feature selection is repeated in EACH inner training fold.
  y_train <- factor(y_train, levels = class_levels)
  inner_ids <- classification_inner_folds(y_train, config, seed)
  inner_data <- lapply(sort(unique(inner_ids)), function(fold) {
    validation_idx <- which(inner_ids == fold)
    training_idx <- which(inner_ids != fold)
    selected <- select_classification_features(
      x_train[training_idx, , drop = FALSE], y_train[training_idx],
      config, seed + fold * 101L
    )
    list(
      x_train = as.matrix(x_train[training_idx, selected, drop = FALSE]),
      y_train = factor(y_train[training_idx], levels = class_levels),
      x_valid = as.matrix(x_train[validation_idx, selected, drop = FALSE]),
      y_valid = factor(y_train[validation_idx], levels = class_levels)
    )
  })
  average_score <- function(values) {
    if (any(!is.finite(values))) return(-Inf)
    mean(values)
  }

  # Elastic net: inner CV chooses alpha AND lambda selection rule; glmnet's
  # own training-only CV chooses lambda along its regularization path.
  alphas <- unique(as.numeric(classification_config(config,
      "classification_alpha_grid", c(0.1, 0.5, 1.0))))
  lambda_rules <- c("lambda.1se", "lambda.min")
  elastic_candidates <- list()
  for (alpha in alphas) {
    score_by_rule <- lapply(lambda_rules, function(x) numeric(length(inner_data)))
    for (i in seq_along(inner_data)) {
      d <- inner_data[[i]]
      nfolds <- min(5L, min(table(d$y_train)))
      set.seed(seed + i * 127L + as.integer(1000 * alpha))
      fit <- tryCatch(glmnet::cv.glmnet(
        x = d$x_train, y = d$y_train, family = "multinomial",
        alpha = alpha, nfolds = nfolds, standardize = TRUE,
        type.measure = "class"
      ), error = function(e) NULL)
      for (j in seq_along(lambda_rules)) {
        score_by_rule[[j]][i] <- if (is.null(fit)) NA_real_ else {
          pred <- tryCatch(as.character(stats::predict(
            fit, newx = d$x_valid, s = lambda_rules[[j]], type = "class"
          )), error = function(e) rep(NA_character_, nrow(d$x_valid)))
          classification_balanced_accuracy(d$y_valid, pred, class_levels)
        }
      }
    }
    for (j in seq_along(lambda_rules)) {
      elastic_candidates[[length(elastic_candidates) + 1L]] <-
        classification_parameter_row("Elastic Net", average_score(score_by_rule[[j]]),
                                     alpha = alpha, lambda_choice = lambda_rules[[j]])
    }
  }
  elastic_table <- do.call(rbind, elastic_candidates)
  if (all(!is.finite(elastic_table$inner_cv_balanced_accuracy))) {
    stop("All elastic-net inner CV fits failed. Check class counts and R packages.")
  }
  best_elastic <- elastic_table[which.max(elastic_table$inner_cv_balanced_accuracy), , drop = FALSE]

  # Random forest: tune mtry rule and minimum terminal node size.
  rf_modes <- unique(as.character(classification_config(config,
      "classification_rf_mtry_modes", c("sqrt", "third", "half"))))
  rf_nodes <- unique(as.integer(classification_config(config,
      "classification_rf_min_node_grid", c(1L, 5L))))
  rf_grid <- expand.grid(mtry_mode = rf_modes, min_node = rf_nodes,
                         stringsAsFactors = FALSE)
  rf_candidates <- vector("list", nrow(rf_grid))
  for (g in seq_len(nrow(rf_grid))) {
    scores <- vapply(seq_along(inner_data), function(i) {
      d <- inner_data[[i]]
      fit <- tryCatch(ranger::ranger(
        x = data.frame(d$x_train, check.names = FALSE), y = d$y_train,
        num.trees = as.integer(classification_config(config, "rf_trees", 300L)),
        mtry = classification_rf_mtry(rf_grid$mtry_mode[g], ncol(d$x_train)),
        min.node.size = as.integer(rf_grid$min_node[g]),
        classification = TRUE,
        seed = seed + g * 101L + i,
        num.threads = as.integer(classification_config(config, "num_threads", 1L))
      ), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      pred <- tryCatch(as.character(stats::predict(fit, data =
        data.frame(d$x_valid, check.names = FALSE))$predictions),
        error = function(e) rep(NA_character_, nrow(d$x_valid)))
      classification_balanced_accuracy(d$y_valid, pred, class_levels)
    }, numeric(1))
    rf_candidates[[g]] <- classification_parameter_row(
      "Random Forest", average_score(scores),
      mtry_mode = rf_grid$mtry_mode[g], min_node = as.integer(rf_grid$min_node[g])
    )
  }
  rf_table <- do.call(rbind, rf_candidates)
  if (all(!is.finite(rf_table$inner_cv_balanced_accuracy))) {
    stop("All random-forest inner CV fits failed.")
  }
  best_rf <- rf_table[which.max(rf_table$inner_cv_balanced_accuracy), , drop = FALSE]

  # XGBoost: tune a joint parameter grid; early stopping on each INNER
  # validation fold determines boosting rounds, never using the outer test fold.
  xgb_grid <- classification_xgb_grid(config)
  max_rounds <- max(20L, as.integer(classification_config(config, "xgb_nrounds", 200L)))
  early_stop <- max(5L, as.integer(classification_config(config,
      "xgb_early_stopping_rounds", 20L)))
  xgb_candidates <- vector("list", nrow(xgb_grid))
  for (g in seq_len(nrow(xgb_grid))) {
    scores <- rep(NA_real_, length(inner_data))
    rounds <- rep(NA_integer_, length(inner_data))
    for (i in seq_along(inner_data)) {
      d <- inner_data[[i]]
      train_matrix <- xgboost::xgb.DMatrix(d$x_train,
        label = as.integer(d$y_train) - 1L)
      valid_matrix <- xgboost::xgb.DMatrix(d$x_valid,
        label = as.integer(d$y_valid) - 1L)
      set.seed(seed + g * 317L + i)
      fit <- tryCatch(xgboost::xgb.train(
        params = classification_xgb_params(xgb_grid[g, , drop = FALSE], class_levels, config),
        data = train_matrix, nrounds = max_rounds,
        evals = list(train = train_matrix, validation = valid_matrix),
        early_stopping_rounds = early_stop, verbose = 0
      ), error = function(e) NULL)
      if (is.null(fit)) next
      pred <- tryCatch(classification_xgb_classes(fit, d$x_valid, class_levels),
                       error = function(e) rep(NA_character_, nrow(d$x_valid)))
      scores[i] <- classification_balanced_accuracy(d$y_valid, pred, class_levels)
      # Best round = lowest validation mlogloss in the evaluation log
      # (works for xgboost 1.x-3.x; fit$best_iteration does not exist in 3.x).
      candidate_rounds <- xgb_best_round_from_log(xgb_model_log(fit), "^validation.*mlogloss")
      if (length(candidate_rounds) == 1L && !is.na(candidate_rounds) &&
          candidate_rounds >= 1L) rounds[i] <- candidate_rounds
    }
    chosen_rounds <- if (any(is.finite(rounds))) {
      as.integer(round(stats::median(rounds[is.finite(rounds)])))
    } else {
      min(100L, max_rounds)
    }
    xgb_candidates[[g]] <- classification_parameter_row(
      "XGBoost", average_score(scores),
      eta = as.numeric(xgb_grid$eta[g]),
      max_depth = as.integer(xgb_grid$max_depth[g]),
      min_child_weight = as.numeric(xgb_grid$min_child_weight[g]),
      subsample = as.numeric(xgb_grid$subsample[g]),
      colsample_bytree = as.numeric(xgb_grid$colsample_bytree[g]),
      nrounds = chosen_rounds
    )
  }
  xgb_table <- do.call(rbind, xgb_candidates)
  if (all(!is.finite(xgb_table$inner_cv_balanced_accuracy))) {
    stop("All XGBoost inner CV fits failed. Check installed xgboost version.")
  }
  best_xgb <- xgb_table[which.max(xgb_table$inner_cv_balanced_accuracy), , drop = FALSE]
  rbind(best_elastic, best_rf, best_xgb)
}

predict_classification_models <- function(x_train, y_train, x_test, selected,
                                          class_levels, config, seed, tuned) {
  # Refit tuned models using ALL outer-training samples and their selected genes.
  set.seed(seed)
  x_train <- as.matrix(x_train[, selected, drop = FALSE])
  x_test <- as.matrix(x_test[, selected, drop = FALSE])
  y_factor <- factor(y_train, levels = class_levels)
  get_tuned <- function(model_name) tuned[tuned$model == model_name, , drop = FALSE]

  predictions <- list(
    Baseline = rep(names(sort(table(y_factor), decreasing = TRUE))[1L], nrow(x_test))
  )

  en <- get_tuned("Elastic Net")
  predictions[["Elastic Net"]] <- tryCatch({
    nfolds <- min(5L, min(table(y_factor)))
    fit <- glmnet::cv.glmnet(
      x = x_train, y = y_factor, family = "multinomial",
      alpha = en$alpha[1L], nfolds = nfolds,
      standardize = TRUE, type.measure = "class"
    )
    as.character(stats::predict(fit, newx = x_test,
                                 s = en$lambda_choice[1L], type = "class"))
  }, error = function(e) {
    log_message("Classification Elastic Net warning: ", conditionMessage(e))
    rep(NA_character_, nrow(x_test))
  })

  rf <- get_tuned("Random Forest")
  predictions[["Random Forest"]] <- tryCatch({
    fit <- ranger::ranger(
      x = data.frame(x_train, check.names = FALSE), y = y_factor,
      num.trees = as.integer(classification_config(config, "rf_trees", 300L)),
      mtry = classification_rf_mtry(rf$mtry_mode[1L], ncol(x_train)),
      min.node.size = as.integer(rf$min_node[1L]), classification = TRUE,
      seed = seed, num.threads = as.integer(classification_config(config, "num_threads", 1L))
    )
    as.character(stats::predict(fit,
      data = data.frame(x_test, check.names = FALSE))$predictions)
  }, error = function(e) {
    log_message("Classification RF warning: ", conditionMessage(e))
    rep(NA_character_, nrow(x_test))
  })

  xgb <- get_tuned("XGBoost")
  predictions[["XGBoost"]] <- tryCatch({
    dtrain <- xgboost::xgb.DMatrix(x_train, label = as.integer(y_factor) - 1L)
    fit <- xgboost::xgb.train(
      params = classification_xgb_params(xgb, class_levels, config),
      data = dtrain, nrounds = as.integer(xgb$nrounds[1L]), verbose = 0
    )
    classification_xgb_classes(fit, x_test, class_levels)
  }, error = function(e) {
    log_message("Classification XGBoost warning: ", conditionMessage(e))
    rep(NA_character_, nrow(x_test))
  })
  predictions
}

run_source_classification_cv <- function(processed, config) {
  # Outer evaluation is unchanged: stratified, sample-level repeated CV.
  output_dir <- file.path(config$output_dir, "models_classification")
  figure_dir <- file.path(config$output_dir, "figures")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  data <- build_source_classification_matrix(processed)
  x <- data$x
  y <- data$display_name
  class_levels <- unique(config$source_map$display_name[
    config$source_map$display_name %in% unique(y)])
  y <- factor(y, levels = class_levels)

  n_folds <- min(config$outer_folds, max(2L, min(table(y))))
  project_blocked <- identical(config$cv_grouping, "project")
  if (project_blocked) {
    if (is.null(data$project_accession)) stop("project_accession missing from processed data")
    projects_per_class <- tapply(data$project_accession, y, function(v) length(unique(v)))
    n_folds <- min(n_folds, min(projects_per_class))
  }
  log_message("Nested source classification (n=", length(y), ", p=", ncol(x), ", ",
              n_folds, " outer folds x ", config$outer_repeats, " repeats",
              if (project_blocked) ", project-blocked" else "", ")")

  all_predictions <- list()
  all_metrics <- list()
  all_features <- list()
  all_tuning <- list()
  all_selection_audit <- list()
  cursor <- 0L

  for (repetition in seq_len(config$outer_repeats)) {
    repetition_seed <- config$seed + 999000L + repetition * 1000L
    fold_ids <- if (project_blocked) {
      group_folds(data$project_accession, n_folds, repetition_seed, strata = y)
    } else {
      class_stratified_folds(y, n_folds, repetition_seed)
    }
    for (fold in seq_len(n_folds)) {
      fold_seed <- repetition_seed + fold * 19L
      test_indices <- which(fold_ids == fold)
      train_indices <- setdiff(seq_along(y), test_indices)

      # IMPORTANT: inner tuning must happen BEFORE selecting features on the
      # full outer training set. Each inner validation fold is unseen by its
      # own feature selection step.
      tuned <- tune_classification_hyperparameters(
        x_train = x[train_indices, , drop = FALSE],
        y_train = y[train_indices], class_levels = class_levels,
        config = config, seed = fold_seed
      )
      selected <- select_classification_features(
        x_train = x[train_indices, , drop = FALSE],
        y_train = y[train_indices], config = config,
        seed = fold_seed + 200001L
      )
      selection_audit <- attr(selected, "diagnostics")
      all_selection_audit[[length(all_selection_audit) + 1L]] <- data.frame(
        repetition = repetition, fold = fold,
        n_eligible = selection_audit$n_eligible,
        n_elastic = selection_audit$n_elastic,
        n_final = selection_audit$n_final,
        feature_cap = selection_audit$feature_cap,
        fallback_used = selection_audit$fallback_used,
        capped = selection_audit$capped
      )
      model_predictions <- predict_classification_models(
        x_train = x[train_indices, , drop = FALSE],
        y_train = y[train_indices],
        x_test = x[test_indices, , drop = FALSE],
        selected = selected, class_levels = class_levels,
        config = config, seed = fold_seed, tuned = tuned
      )

      cursor <- cursor + 1L
      fold_prediction <- do.call(rbind,
        lapply(names(model_predictions), function(model_name) {
          data.frame(
            repetition = repetition, fold = fold,
            sample_accession = data$sample_accession[test_indices],
            actual = as.character(y[test_indices]),
            predicted = model_predictions[[model_name]],
            model = model_name, n_selected_features = length(selected),
            stringsAsFactors = FALSE
          )
        })
      )
      all_predictions[[cursor]] <- fold_prediction

      all_metrics[[cursor]] <- do.call(rbind,
        lapply(split(fold_prediction, fold_prediction$model), function(rows) {
          metrics <- classification_metrics(rows$actual, rows$predicted, class_levels)
          data.frame(repetition = repetition, fold = fold,
                     model = rows$model[1L], n_selected_features = length(selected),
                     metrics, stringsAsFactors = FALSE)
        })
      )
      tuned$repetition <- repetition
      tuned$fold <- fold
      tuned$n_selected_features <- length(selected)
      all_tuning[[cursor]] <- tuned

      mapping <- data$feature_mapping[
        match(selected, data$feature_mapping$safe_name), , drop = FALSE]
      mapping$repetition <- repetition
      mapping$fold <- fold
      all_features[[cursor]] <- mapping
      log_message("Nested classification repetition ", repetition, " fold ", fold,
                  ": selected ", length(selected), " ARGs")
    }
  }

  predictions <- data.table::rbindlist(all_predictions, fill = TRUE)
  fold_metrics <- data.table::rbindlist(all_metrics, fill = TRUE)
  selected_features <- data.table::rbindlist(all_features, fill = TRUE)
  tuned_hyperparameters <- data.table::rbindlist(all_tuning, fill = TRUE)
  selection_audit <- data.table::rbindlist(all_selection_audit, fill = TRUE)
  write_csv_atomic(selection_audit,
    file.path(output_dir, "source_classification_feature_selection_audit.csv"))
  write_csv_atomic(predictions,
    file.path(output_dir, "source_classification_predictions.csv"))
  write_csv_atomic(fold_metrics,
    file.path(output_dir, "source_classification_fold_metrics.csv"))
  write_csv_atomic(selected_features,
    file.path(output_dir, "source_classification_selected_features.csv"))
  write_csv_atomic(tuned_hyperparameters,
    file.path(output_dir, "source_classification_tuned_hyperparameters.csv"))

  repetition_metrics <- fold_metrics[, .(
    n = sum(n), accuracy = mean(accuracy, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    macro_f1 = mean(macro_f1, na.rm = TRUE)
  ), by = .(repetition, model)]
  write_csv_atomic(repetition_metrics,
    file.path(output_dir, "source_classification_repetition_metrics.csv"))

  summary_rows <- list()
  summary_cursor <- 0L
  for (metric in c("accuracy", "balanced_accuracy", "macro_f1")) {
    grouped <- split(repetition_metrics, repetition_metrics$model)
    for (rows in grouped) {
      summary_cursor <- summary_cursor + 1L
      interval <- mean_ci_by_repetition(rows[[metric]])
      summary_rows[[summary_cursor]] <- data.frame(
        model = rows$model[1], metric = metric, estimate = interval[["mean"]],
        lower = interval[["lower"]], upper = interval[["upper"]],
        repetitions = interval[["n"]], stringsAsFactors = FALSE
      )
    }
  }
  model_summary <- data.table::rbindlist(summary_rows)
  write_csv_atomic(model_summary,
    file.path(output_dir, "source_classification_summary.csv"))

  feature_frequency <- selected_features[, .(
    selection_count = .N
  ), by = .(safe_name, ref_sequence, accession, gene_label)]
  feature_frequency[, outer_folds := config$outer_repeats * n_folds]
  feature_frequency[, selection_frequency := selection_count / outer_folds]
  data.table::setorder(feature_frequency, -selection_frequency, gene_label)
  write_csv_atomic(feature_frequency,
    file.path(output_dir, "source_classification_feature_frequency.csv"))

  # A fixed 150-feature count is no longer expected. Report observed counts.
  feature_counts <- unique(fold_metrics[, .(repetition, fold, n_selected_features)])
  write_csv_atomic(feature_counts,
    file.path(output_dir, "source_classification_feature_counts.csv"))

  plot_data <- model_summary[metric == "balanced_accuracy" & model != "Baseline"]
  classification_plot <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = model, y = estimate)
  ) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                            width = 0.15, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(x = NULL, y = "Balanced accuracy",
      title = "Nested leakage-controlled source classification from ARG profiles",
      subtitle = "Feature selection repeated within every inner training fold") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
  ggplot2::ggsave(
    file.path(figure_dir, "source_classification_balanced_accuracy.png"),
    classification_plot, width = 7.2, height = 4.8, dpi = 320, bg = "white"
  )

  list(data = data, predictions = predictions, fold_metrics = fold_metrics,
       repetition_metrics = repetition_metrics, model_summary = model_summary,
       feature_frequency = feature_frequency,
       tuned_hyperparameters = tuned_hyperparameters,
       feature_counts = feature_counts, selection_audit = selection_audit)
}


# Combine per-source regression results run in parallel. All tables are
# per-source, so row-binding gives exactly the sequential result.
combine_regression_parts <- function(parts, processed, config) {
  order <- names(processed)
  parts <- parts[order(match(vapply(parts, `[[`, character(1), "key"), order))]
  bind <- function(name) data.table::rbindlist(lapply(parts, `[[`, name), fill = TRUE)
  out <- list(
    predictions = bind("predictions"), fold_metrics = bind("fold_metrics"),
    repetition_metrics = bind("repetition_metrics"), model_summary = bind("model_summary"),
    feature_frequency = bind("feature_frequency"), stable_features = bind("stable_features")
  )
  output_dir <- file.path(config$output_dir, "models_regression")
  for (name in names(out)) {
    file <- c(predictions = "heldout_predictions.csv", fold_metrics = "fold_metrics.csv",
              repetition_metrics = "repetition_metrics.csv", model_summary = "model_summary.csv",
              feature_frequency = "feature_selection_frequency.csv",
              stable_features = "stable_indicator_features.csv")[[name]]
    write_csv_atomic(out[[name]], file.path(output_dir, file))
  }
  out
}
