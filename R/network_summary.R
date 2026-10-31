# =============================================================================
# R/network_summary.R - CARD annotation, gene network and result tables
# Project : Candidate antimicrobial resistance gene indicators for veterinary
#           public health surveillance (Mahawan et al., Research in Veterinary
#           Science, submitted)
# Purpose : Annotate repeatedly selected ResFinder sequences with CARD terms
#           (accession+gene, unique name, documented synonyms; never accession
#           alone), build the multilayer network (ecological, mechanism and
#           drug-class edges), rank genes by betweenness, and write summary tables.
# Used by : run_all.R (sourced; defines functions only)
# Manuscript: Methods 2.6; Table 4; Supplementary Tables S5-S7
# =============================================================================

# Network and result summaries. v4 rules and values kept.

# Ecological edges are sample-level correlations; annotation edges are mechanism/drug-class links.

pairwise_edges <- function(values) {
  values <- sort(unique(values[!is.na(values) & nzchar(values)]))
  if (length(values) < 2L) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  pairs <- utils::combn(values, 2)
  data.frame(from = pairs[1, ], to = pairs[2, ], stringsAsFactors = FALSE)
}

# Conservative CARD crosswalk: match both reference accession and gene label.
# The full 6,484-entry ARO catalog is NOT a sample-feature mapping table.
card_normalize_accession <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x) | x %in% c("", "NA")] <- ""
  sub("\\.[0-9]+$", "", x)
}

card_normalize_gene <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[\u2018\u2019\u2032]", "'", x)
  x <- gsub("\u2033", "''", x, fixed = TRUE)
  x <- gsub("[[:space:]]+", "", x)
  x[is.na(x) | x == "na"] <- ""
  x
}

# Name key for matching ResFinder gene labels to CARD names. Case, spaces,
# brackets and hyphens are ignored, but PRIME marks are kept, so APH(3')-Ib and
# APH(3'')-Ib remain different genes.
card_name_key <- function(x) {
  x <- card_normalize_gene(x)
  gsub("[^a-z0-9']", "", x)
}

# Curated ResFinder -> CARD synonyms used only when no direct match exists.
# Each entry is a documented synonym; edit here if the CARD release changes.
default_card_aliases <- c(
  "ant(3'')-Ia" = "aadA",               # CARD names ANT(3'')-Ia as aadA
  "aadA1"       = "aadA",
  "aph(3')-III" = "APH(3')-IIIa",
  "aadD"        = "ANT(4')-Ia",
  "mdf(A)"      = "Escherichia coli mdfA"
)

# Reads the CARD table exported by build_full_CARD_mapping_online.R
# (data/AMRgenelist_mapped_CARD_ALL.csv). Columns used: Gene, ID (CARD DNA
# reference accession), ARO Accession, ARO Name, Drug Class, Resistance
# Mechanism, CARD Description, CARD Short Name.
read_annotation_tables <- function(config) {
  path <- config$mapping_path
  if (is.null(path) || !file.exists(path)) stop("CARD mapping file not found: ", path)
  card <- if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    sheets <- readxl::excel_sheets(path)
    as.data.frame(readxl::read_excel(path, sheet = if ("Mapping_Style" %in% sheets) "Mapping_Style" else 1))
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")
  }
  needed <- c("Gene", "ID", "ARO Accession", "ARO Name", "Drug Class", "Resistance Mechanism")
  if (!all(needed %in% names(card))) {
    stop("CARD mapping missing columns: ", paste(setdiff(needed, names(card)), collapse = ", "))
  }
  if (!"CARD Description" %in% names(card)) card[["CARD Description"]] <- NA_character_
  if (!"CARD Short Name" %in% names(card)) card[["CARD Short Name"]] <- NA_character_
  card <- card[!is.na(card[["ARO Accession"]]) & nzchar(card[["ARO Accession"]]), , drop = FALSE]
  data.frame(
    gene = as.character(card$Gene),
    accession_key = card_normalize_accession(card$ID),
    aro_accession = as.character(card[["ARO Accession"]]),
    aro_name = as.character(card[["ARO Name"]]),
    short_name = as.character(card[["CARD Short Name"]]),
    drug_class = as.character(card[["Drug Class"]]),
    resistance_mechanism = as.character(card[["Resistance Mechanism"]]),
    card_description = as.character(card[["CARD Description"]]),
    stringsAsFactors = FALSE
  )
}

# Annotate ResFinder features with CARD terms in three tiers:
#  1. accession + gene name agree with a CARD reference sequence (most specific);
#  2. the gene name (or name without the "bla" prefix) equals exactly ONE CARD
#     ARO Name/Short Name;
#  3. a curated synonym in config$card_aliases (default_card_aliases).
# Matching on accession alone is never used: ResFinder and CARD use different
# reference accessions, and some accessions carry more than one gene
# (e.g. M28829 carries aph(3'')-Ib and aph(6)-Id).
# config$card_match = "accession_and_gene" keeps only tier 1 (conservative).
annotate_features_with_card <- function(features, card, config) {
  mode <- if (is.null(config$card_match)) "tiered" else config$card_match
  aliases <- if (is.null(config$card_aliases)) default_card_aliases else config$card_aliases
  terms <- card[!duplicated(card$aro_accession), , drop = FALSE]
  cols <- c("aro_accession", "aro_name", "drug_class", "resistance_mechanism", "card_description")
  out <- features
  for (col in cols) out[[col]] <- NA_character_
  out$card_annotation_status <- "unmatched"

  # Tier 1
  card$tier1_key <- paste(card$accession_key, card_name_key(card$gene), sep = "\034")
  tier1 <- split(card$aro_accession, card$tier1_key)
  tier1 <- vapply(tier1, function(v) if (length(unique(v)) == 1L) v[1] else NA_character_, character(1))
  key1 <- paste(card_normalize_accession(out$accession), card_name_key(out$gene_label), sep = "\034")
  aro <- unname(tier1[key1])
  status <- ifelse(is.na(aro), "unmatched", "tier1_accession_and_gene")

  if (!identical(mode, "accession_and_gene")) {
    # Tier 2: unique exact name match against ARO Name or CARD Short Name
    name_tab <- rbind(
      data.frame(key = card_name_key(terms$aro_name), aro = terms$aro_accession, stringsAsFactors = FALSE),
      data.frame(key = card_name_key(terms$short_name), aro = terms$aro_accession, stringsAsFactors = FALSE)
    )
    name_tab <- unique(name_tab[nzchar(name_tab$key), , drop = FALSE])
    by_name <- split(name_tab$aro, name_tab$key)
    unique_name <- vapply(by_name, function(v) if (length(unique(v)) == 1L) v[1] else NA_character_, character(1))
    try_name <- function(labels) unname(unique_name[card_name_key(labels)])
    need <- is.na(aro)
    hit <- try_name(out$gene_label[need])
    aro[need] <- hit; status[need][!is.na(hit)] <- "tier2_unique_name"
    need <- is.na(aro)
    hit <- try_name(sub("^bla", "", out$gene_label[need]))
    aro[need] <- hit; status[need][!is.na(hit)] <- "tier2_unique_name_without_bla"
    # Tier 3: curated synonyms
    need <- is.na(aro) & out$gene_label %in% names(aliases)
    if (any(need)) {
      hit <- try_name(aliases[out$gene_label[need]])
      aro[need] <- hit; status[need][!is.na(hit)] <- "tier3_curated_synonym"
    }
  }
  idx <- match(aro, terms$aro_accession)
  for (col in cols) out[[col]] <- terms[[col]][idx]
  out$card_annotation_status <- status
  out
}

build_corrected_network <- function(processed, stable_features, config) {
  network_dir <- file.path(config$output_dir, "network")
  figure_dir <- file.path(config$output_dir, "figures")
  dir.create(network_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  # Step 1: Annotate every selected feature with CARD (see annotate_features_with_card).
  # Retain ALL selected features, including features with no CARD annotation.
  card <- read_annotation_tables(config)
  stable <- as.data.frame(stable_features)
  required <- c("source_key", "display_name", "safe_name", "ref_sequence", "selection_frequency")
  missing <- setdiff(required, names(stable))
  if (length(missing)) stop("stable_features missing: ", paste(missing, collapse = ", "))
  stable <- repair_feature_labels(stable)
  n_before <- nrow(stable)
  stable <- annotate_features_with_card(stable, card, config)
  if (nrow(stable) != n_before) stop("Annotation step unexpectedly changed feature count")
  write_csv_atomic(stable[, c("source_key", "ref_sequence", "gene_label", "accession",
    "selection_frequency", "aro_accession", "aro_name", "drug_class",
    "resistance_mechanism", "card_annotation_status")],
    file.path(network_dir, "selected_ARG_CARD_mapping_audit.csv"))
  message("CARD annotation: ", paste(names(table(stable$card_annotation_status)),
          table(stable$card_annotation_status), sep = " = ", collapse = "; "))

  # Step 2: Define nodes from the original gene labels, not unverified ARO names.
  stable$node <- trimws(as.character(stable$gene_label))
  if (any(is.na(stable$node) | !nzchar(stable$node))) {
    stop("Selected features have missing gene labels")
  }

  # Step 3: Summarize nodes.
  node_table <- data.table::as.data.table(stable)[, .(
    source_count = data.table::uniqueN(source_key),
    sources = paste(sort(unique(display_name)), collapse = "|"),
    accessions = paste(sort(unique(accession)), collapse = "|"),
    gene_labels = paste(sort(unique(gene_label)), collapse = "|"),
    aro_accessions = paste(sort(unique(stats::na.omit(aro_accession))), collapse = "|"),
    drug_classes = paste(sort(unique(stats::na.omit(drug_class))), collapse = "|"),
    mechanisms = paste(sort(unique(stats::na.omit(resistance_mechanism))), collapse = "|"),
    card_descriptions = paste(sort(unique(stats::na.omit(card_description))), collapse = "|"),
    selection_frequency_max = max(selection_frequency, na.rm = TRUE)
  ), by = node]

  # Step 4: Build ecological edges from sample-level positive Spearman associations.
  ecological_edges <- list()
  ecological_cursor <- 0L
  for (current_source in intersect(names(processed), unique(stable$source_key))) {
    object <- processed[[current_source]]
    selected <- stable[
      stable$source_key == current_source & stable$safe_name %in% colnames(object$x_clr),
      c("safe_name", "node"),
      drop = FALSE
    ]
    selected <- unique(selected)
    if (nrow(selected) < 2L) next

    # Collapse reference features sharing the same original gene label (for plotting network).
    node_matrix <- sapply(
      split(selected$safe_name, selected$node),
      function(features) {
        if (length(features) == 1L) {
          object$x_clr[, features]
        } else {
          rowMeans(object$x_clr[, features, drop = FALSE])
        }
      }
    )
    node_matrix <- as.matrix(node_matrix)
    if (ncol(node_matrix) < 2L) next

    pairs <- utils::combn(colnames(node_matrix), 2)
    source_edges <- lapply(seq_len(ncol(pairs)), function(index) {
      from <- pairs[1, index]
      to <- pairs[2, index]
      test <- suppressWarnings(stats::cor.test(
        node_matrix[, from],
        node_matrix[, to],
        method = "spearman",
        exact = FALSE
      ))
      data.frame(
        from = min(from, to),
        to = max(from, to),
        layer = "Ecological",
        source_key = current_source,
        source = object$display_name,
        weight = unname(test$estimate),
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    })
    source_edges <- data.table::rbindlist(source_edges)
    source_edges[, q_value := stats::p.adjust(p_value, method = "BH")]
    source_edges <- source_edges[
      is.finite(weight) &
        weight >= config$correlation_threshold &
        q_value <= config$correlation_fdr
    ]
    if (nrow(source_edges)) {
      ecological_cursor <- ecological_cursor + 1L
      ecological_edges[[ecological_cursor]] <- source_edges
    }
  }
  if (length(ecological_edges)) {
    ecological_edges <- data.table::rbindlist(ecological_edges, fill = TRUE)
  } else {
    ecological_edges <- data.table::data.table(
      from = character(), to = character(), layer = character(), source_key = character(),
      source = character(), weight = numeric(), p_value = numeric(), q_value = numeric()
    )
  }

  # Step 5: Annotation edges. Split multi-valued CARD categories into
  # individual terms: "class A; class B" is NOT a new combined class.
  card_category_edges <- function(rows, column, layer) {
    cells <- as.character(rows[[column]])
    valid <- which(!is.na(cells) & nzchar(trimws(cells)))
    if (!length(valid)) return(data.table::data.table())
    expanded <- lapply(valid, function(i) {
      categories <- trimws(strsplit(cells[i], ";", fixed = TRUE)[[1]])
      categories <- unique(categories[nzchar(categories)])
      if (!length(categories)) return(NULL)
      data.frame(node = rows$node[i], category = categories,
                 stringsAsFactors = FALSE)
    })
    expanded <- Filter(Negate(is.null), expanded)
    if (!length(expanded)) return(data.table::data.table())
    expanded <- unique(do.call(rbind, expanded))
    groups <- split(expanded$node, expanded$category)
    pairs <- lapply(names(groups), function(category) {
      vals <- pairwise_edges(groups[[category]])
      if (!nrow(vals)) return(NULL)
      vals$source <- category
      vals
    })
    pairs <- Filter(Negate(is.null), pairs)
    if (!length(pairs)) return(data.table::data.table())
    edges <- data.table::rbindlist(pairs)
    edges[, `:=`(layer = layer, source_key = NA_character_,
                 weight = 1, p_value = NA_real_, q_value = NA_real_)]
    unique(edges)
  }
  mechanism_edges <- card_category_edges(stable, "resistance_mechanism", "Mechanism")
  drug_edges <- card_category_edges(stable, "drug_class", "DrugClass")

  # Step 6: Combine edges and calculate network centrality.
  edges <- data.table::rbindlist(list(ecological_edges, mechanism_edges, drug_edges), fill = TRUE)
  if (nrow(edges)) {
    edges[, `:=`(from = pmin(from, to), to = pmax(from, to))]
    edges <- unique(edges)
  }
  write_csv_atomic(edges, file.path(network_dir, "network_edges.csv"))

  graph <- igraph::graph_from_data_frame(
    edges,
    directed = FALSE,
    vertices = as.data.frame(node_table)
  )

  # Centrality is calculated on the layers in config$centrality_layers
  # (default: all three, as in v4). The full graph is still saved and plotted.
  centrality_layers <- if (is.null(config$centrality_layers)) {
    c("Ecological", "Mechanism", "DrugClass")
  } else config$centrality_layers
  centrality_graph <- igraph::subgraph.edges(
    graph, which(igraph::E(graph)$layer %in% centrality_layers), delete.vertices = FALSE)
  centrality_scores <- function(g) {
    if (igraph::ecount(g)) {
      edge_strength <- pmax(abs(igraph::E(g)$weight), 1e-6)
      distances <- 1 / edge_strength
      list(
        betweenness = igraph::betweenness(g, directed = FALSE, weights = distances, normalized = TRUE),
        closeness = igraph::closeness(g, weights = distances, normalized = TRUE),
        eigenvector = igraph::eigen_centrality(g, directed = FALSE, weights = edge_strength)$vector
      )
    } else {
      zero <- stats::setNames(rep(0, igraph::vcount(g)), igraph::V(g)$name)
      list(betweenness = zero, closeness = zero, eigenvector = zero)
    }
  }
  scores <- centrality_scores(centrality_graph)
  betweenness <- scores$betweenness
  closeness <- scores$closeness
  eigenvector <- scores$eigenvector

  centrality <- data.frame(
    node = igraph::V(graph)$name,
    degree = igraph::degree(graph),
    degree_centrality_layers = igraph::degree(centrality_graph),
    strength = igraph::strength(centrality_graph, weights = abs(igraph::E(centrality_graph)$weight)),
    betweenness = betweenness,
    closeness = closeness,
    eigenvector = eigenvector,
    stringsAsFactors = FALSE
  )
  centrality <- merge(centrality, as.data.frame(node_table), by = "node", all.x = TRUE)

  positive_betweenness <- centrality$betweenness[
    is.finite(centrality$betweenness) & centrality$betweenness > 0
  ]
  threshold <- if (length(positive_betweenness)) {
    stats::quantile(positive_betweenness, config$centrality_quantile, names = FALSE)
  } else {
    0
  }
  centrality$prioritized_central <-
    centrality$source_count >= config$min_sources_for_central & centrality$betweenness >= threshold
  centrality <- centrality[
    order(-centrality$prioritized_central, -centrality$betweenness, centrality$node),
  ]

  write_csv_atomic(centrality, file.path(network_dir, "network_nodes.csv"))
  write_csv_atomic(
    centrality[centrality$prioritized_central, ],
    file.path(network_dir, "prioritized_central_args.csv")
  )

  network_metrics <- data.frame(
    nodes = igraph::vcount(graph),
    edges = igraph::ecount(graph),
    ecological_edges = sum(edges$layer == "Ecological"),
    mechanism_edges = sum(edges$layer == "Mechanism"),
    drug_class_edges = sum(edges$layer == "DrugClass"),
    components = igraph::components(graph)$no,
    density = igraph::edge_density(graph, loops = FALSE),
    prioritized_nodes = sum(centrality$prioritized_central),
    stringsAsFactors = FALSE
  )
  network_metrics$centrality_layers <- paste(centrality_layers, collapse = "+")
  network_metrics$betweenness_threshold <- threshold
  network_metrics$annotated_features <- sum(!is.na(stable$aro_accession))
  network_metrics$selected_features <- nrow(stable)
  write_csv_atomic(network_metrics, file.path(network_dir, "network_metrics.csv"))

  # Sensitivity: prioritization under alternative layer sets (same rule).
  layer_sets <- list(
    "Ecological only" = "Ecological",
    "Ecological + Mechanism + DrugClass" = c("Ecological", "Mechanism", "DrugClass")
  )
  sensitivity <- data.table::rbindlist(lapply(names(layer_sets), function(label) {
    g <- igraph::subgraph.edges(graph, which(igraph::E(graph)$layer %in% layer_sets[[label]]),
                                delete.vertices = FALSE)
    b <- centrality_scores(g)$betweenness
    pb <- b[is.finite(b) & b > 0]
    thr <- if (length(pb)) stats::quantile(pb, config$centrality_quantile, names = FALSE) else 0
    sc <- node_table$source_count[match(names(b), node_table$node)]
    data.table::data.table(layers = label, node = names(b), degree = igraph::degree(g),
      betweenness = b, threshold = thr, source_count = sc,
      prioritized = sc >= config$min_sources_for_central & b >= thr)
  }))
  data.table::setorder(sensitivity, layers, -prioritized, -betweenness)
  write_csv_atomic(sensitivity, file.path(network_dir, "network_centrality_sensitivity.csv"))

  # Step 7: Save a simple network overview figure.
  grDevices::png(
    file.path(figure_dir, "network_overview.png"),
    width = 2400,
    height = 1800,
    res = 240,
    bg = "white"
  )
  graphics::par(mar = c(1, 1, 4, 1))
  set.seed(config$seed)
  layout <- igraph::layout_with_fr(graph)
  node_values <- centrality$source_count[match(igraph::V(graph)$name, centrality$node)]
  node_values[!is.finite(node_values)] <- 1
  node_colors <- grDevices::colorRampPalette(c("#D9F0D3", "#31A354"))(
    max(2L, max(node_values, na.rm = TRUE))
  )[pmax(1L, node_values)]
  label_nodes <- centrality$node[seq_len(min(12L, nrow(centrality)))]
  node_labels <- ifelse(igraph::V(graph)$name %in% label_nodes, igraph::V(graph)$name, NA_character_)
  edge_colors <- c(Ecological = "#238B45AA", Mechanism = "#2B8CBE80", DrugClass = "#F1691380")

  igraph::plot.igraph(
    graph,
    layout = layout,
    vertex.size = 7 + 2 * sqrt(igraph::degree(graph)),
    vertex.label = node_labels,
    vertex.label.cex = 0.65,
    vertex.label.color = "#222222",
    vertex.color = node_colors,
    vertex.frame.color = ifelse(
      igraph::V(graph)$name %in% centrality$node[centrality$prioritized_central],
      "#B2182B",
      "#FFFFFF"
    ),
    edge.color = edge_colors[igraph::E(graph)$layer],
    edge.width = 0.6 + 1.4 * pmin(abs(igraph::E(graph)$weight), 1)
  )
  graphics::title(
    main = "Stable ARG indicator network with sample-level ecological edges",
    line = 1,
    cex.main = 1.2
  )
  graphics::legend(
    "topleft",
    inset = 0.01,
    legend = names(edge_colors),
    col = edge_colors,
    lty = 1,
    lwd = 3,
    bty = "n",
    cex = 0.85
  )
  grDevices::dev.off()

  list(
    graph = graph,
    nodes = centrality,
    edges = edges,
    prioritized = centrality[centrality$prioritized_central, ],
    metrics = network_metrics,
    sensitivity = sensitivity,
    annotations = stable
  )
}


write_analysis_summary <- function(preparation, regression, classification, network, config) {
  output_path <- file.path(config$output_dir, "analysis_summary.md")
  snippet_path <- file.path(config$output_dir, "manuscript_result_snippets.md")

  lines <- c(
    "# Revised AMR analysis summary",
    "",
    paste0("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "."),
    paste0("Input metadata: ", config$metadata_path),
    paste0("Input ARG table: ", config$arg_path),
    paste0("Input annotation workbook: ", config$mapping_path),
    paste0("Input ARO table: ", config$aro_path),
    "",
    "## Key interpretation boundary",
    "",
    "Regression models reconstruct normalized ARG load from selected ARG indicators. The target is still derived from ARG mapping data, so these models should be interpreted as internal dimension-reduction/reconstruction, not prediction of clinical resistance, transmission, host linkage, or plasmid mobility.",
    "",
    "The source-classification task is included as a stronger complementary analysis because the outcome is an independent metadata label rather than a mathematical sum of ARG features.",
    "",
    "## Sample flow",
    "",
    "| Source | Metadata-eligible | ARG-positive | Modeled | Features |",
    "|---|---:|---:|---:|---:|"
  )

  for (index in seq_len(nrow(preparation$sample_flow))) {
    row <- preparation$sample_flow[index, ]
    lines <- c(
      lines,
      paste0(
        "| ", row$display_name,
        " | ", row$metadata_eligible_samples,
        " | ", row$arg_positive_samples,
        " | ", row$modeled_samples,
        " | ", row$modeled_features,
        " |"
      )
    )
  }

  lines <- c(
    lines,
    "",
    paste0(
      "Across the five requested metadata labels, ",
      format(sum(preparation$sample_flow$metadata_eligible_samples, na.rm = TRUE), big.mark = ","),
      " samples were metadata-eligible and ",
      format(sum(preparation$sample_flow$modeled_samples, na.rm = TRUE), big.mark = ","),
      " had ARG mappings plus valid bacterial-fragment denominators for modeling."
    )
  )

  if (!is.null(regression)) {
    best <- regression$model_summary[
      regression$model_summary$metric == "R_squared" & regression$model_summary$model != "Baseline",
    ]
    best <- best[order(best$source_key, -best$estimate), ]
    best <- best[!duplicated(best$source_key), ]
    lines <- c(lines, "", "## Leakage-controlled ARG-load reconstruction", "")
    for (index in seq_len(nrow(best))) {
      row <- best[index, ]
      lines <- c(
        lines,
        paste0(
          "- ", row$display_name, ": ", row$model,
          " held-out R-squared ", format_interval(row$estimate, row$lower, row$upper), "."
        )
      )
    }
  }

  if (!is.null(classification)) {
    best_class <- classification$model_summary[
      classification$model_summary$metric == "balanced_accuracy" & classification$model_summary$model != "Baseline",
    ]
    best_class <- best_class[order(-best_class$estimate), ]
    lines <- c(lines, "", "## Source-classification performance", "")
    for (index in seq_len(nrow(best_class))) {
      row <- best_class[index, ]
      lines <- c(
        lines,
        paste0(
          "- ", row$model,
          ": balanced accuracy ", format_interval(row$estimate, row$lower, row$upper), "."
        )
      )
    }
  }

  if (!is.null(network)) {
    prioritized <- network$prioritized$node
    lines <- c(
      lines,
      "",
      "## Corrected network",
      "",
      paste0(
        "The corrected network contained ", network$metrics$nodes,
        " nodes and ", network$metrics$edges,
        " edges, including ", network$metrics$ecological_edges,
        " sample-level ecological correlation edges. ",
        length(prioritized),
        " nodes met the prespecified centrality and multi-source criteria."
      )
    )
    if (length(prioritized)) {
      lines <- c(lines, paste0("Prioritized nodes: ", paste(prioritized, collapse = ", "), "."))
    }
  }

  writeLines(lines, output_path, useBytes = TRUE)

  snippet <- c(
    "# Manuscript-ready result wording",
    "",
    "## Methods wording",
    "",
    "To avoid feature-selection leakage, prevalence screening, elastic-net screening, Boruta selection, model fitting, and hyperparameter tuning were performed within each outer training fold. Validation samples were held out until final prediction. Regression used log1p gene-level ARG fragments per million bacterial fragments as predictors and log1p total ARG fragments per million bacterial fragments as the reconstruction target. CLR-transformed length-normalized abundances were retained only for sample-level ecological network analysis.",
    "",
    "## Results wording",
    "",
    "After leakage control and bacterial-depth normalization, source-specific ARG-load reconstruction showed modest-to-moderate held-out performance. These values should be interpreted as reconstruction of normalized ARG load from selected indicator genes, not as independent prediction of clinical resistance or transmission.",
    "",
    "A complementary source-classification analysis was included because the outcome was an independent metadata label rather than a mathematical sum of ARG abundances. This analysis evaluates whether ARG profiles contain source-discriminatory surveillance signal across avian, compost, feces, food, and wastewater metagenomes.",
    "",
    "## Limitation wording",
    "",
    "Although log gene-level normalized abundance improved alignment between predictors and the ARG-load target, the regression task remains an internal reconstruction problem because both target and predictors are derived from the same ARG mapping table. External validation against phenotypic antimicrobial resistance, curated surveillance outcomes, assembly-based host linkage, or plasmid-resolved evidence is still required."
  )
  writeLines(snippet, snippet_path, useBytes = TRUE)

  invisible(list(summary = output_path, snippets = snippet_path))
}

# ---------------------------------------------------------------------------
# Model result tables for the manuscript (Table 2/3 and supplementary tables).
# ---------------------------------------------------------------------------
summarise_regression <- function(regression, validation) {
  if (is.null(regression)) return(NULL)
  ms <- data.table::as.data.table(regression$model_summary)[metric == "R_squared"]
  ms[, validation := validation]
  ms[, .(validation, source = display_name, model, R2 = estimate, lower, upper)]
}

summarise_classification <- function(classification, validation) {
  if (is.null(classification)) return(NULL)
  p <- data.table::as.data.table(classification$predictions)
  reps <- length(unique(p$repetition))
  recall <- p[, .(recall = mean(actual == predicted, na.rm = TRUE)), by = .(model, class = actual)]
  recall[, validation := validation]
  confusion <- p[model != "Baseline", .(mean_per_repetition = .N / reps),
                 by = .(model, actual, predicted)]
  confusion[, validation := validation]
  overall <- data.table::as.data.table(classification$model_summary)
  overall[, validation := validation]
  list(overall = overall, recall = recall, confusion = confusion)
}

write_model_reports <- function(regression, classification, output_dir,
                                regression_project = NULL, classification_project = NULL) {
  dir <- file.path(output_dir, "results_tables")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  reg <- data.table::rbindlist(list(
    summarise_regression(regression, "sample-level CV"),
    summarise_regression(regression_project, "project-blocked CV")), fill = TRUE)
  write_csv_atomic(reg, file.path(dir, "regression_R2_all_models.csv"))
  cls <- Filter(Negate(is.null), list(
    summarise_classification(classification, "sample-level CV"),
    summarise_classification(classification_project, "project-blocked CV")))
  overall <- data.table::rbindlist(lapply(cls, `[[`, "overall"), fill = TRUE)
  recall <- data.table::rbindlist(lapply(cls, `[[`, "recall"), fill = TRUE)
  confusion <- data.table::rbindlist(lapply(cls, `[[`, "confusion"), fill = TRUE)
  write_csv_atomic(overall, file.path(dir, "classification_summary.csv"))
  write_csv_atomic(recall, file.path(dir, "classification_class_recall.csv"))
  write_csv_atomic(confusion, file.path(dir, "classification_confusion_matrix.csv"))

  # Short human-readable digest.
  f3 <- function(x) formatC(x, digits = 3, format = "f")
  lines <- c("# Key model results", "")
  best <- reg[model != "Baseline", .SD[which.max(R2)], by = .(validation, source)]
  lines <- c(lines, "## Best held-out R2 by source", "",
             "| Validation | Source | Model | R2 (interval) |", "|---|---|---|---|",
             best[, sprintf("| %s | %s | %s | %s (%s to %s) |", validation, source, model,
                            f3(R2), f3(lower), f3(upper))])
  ind <- reg[model == "Indicator Sum"]
  if (nrow(ind)) lines <- c(lines, "", "Calibrated indicator-sum baseline R2: ",
    ind[, paste0("- ", validation, ", ", source, ": ", f3(R2))])
  ba <- overall[metric == "balanced_accuracy"]
  lines <- c(lines, "", "## Source classification (balanced accuracy)", "",
             "| Validation | Model | Balanced accuracy (interval) |", "|---|---|---|",
             ba[, sprintf("| %s | %s | %s (%s to %s) |", validation, model, f3(estimate),
                          f3(lower), f3(upper))])
  rc <- data.table::dcast(recall[model != "Baseline"], validation + model ~ class, value.var = "recall")
  lines <- c(lines, "", "## Class-specific recall", "",
             paste0("| ", paste(names(rc), collapse = " | "), " |"),
             paste0("|", paste(rep("---", ncol(rc)), collapse = "|"), "|"),
             apply(rc, 1, function(r) paste0("| ", paste(ifelse(is.na(suppressWarnings(as.numeric(r))), r,
               f3(as.numeric(r))), collapse = " | "), " |")))
  writeLines(lines, file.path(dir, "key_model_results.md"), useBytes = TRUE)
  invisible(dir)
}
