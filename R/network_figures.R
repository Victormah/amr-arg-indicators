# =============================================================================
# R/network_figures.R - network plots and Cytoscape export
# Project : Candidate antimicrobial resistance gene indicators for veterinary
#           public health surveillance (Mahawan et al., Research in Veterinary
#           Science, submitted)
# Purpose : Two-panel network plots (ecological and all-layer) and Cytoscape
#           node/edge/SIF/GraphML exports from output/network/.
# Used by : run_all.R (sourced; defines functions only)
# Manuscript: draft network figures; final Fig. 4 and Fig. S1 are drawn by
#           scripts/make_manuscript_figures.R
# =============================================================================

# Figure 2 from the saved network tables. No centrality is recalculated.
# Prioritized genes are read from <output_dir>/network/prioritized_central_args.csv
# unless supplied explicitly, so the figure always matches the network output.
read_prioritized_genes <- function(output_dir) {
  path <- file.path(output_dir, "network", "prioritized_central_args.csv")
  if (!file.exists(path)) stop("Missing ", path, ". Run the network step first.")
  genes <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)$node
  if (!length(genes)) stop("No prioritized genes in ", path)
  trimws(genes)
}

check_against_metrics <- function(output_dir, g) {
  path <- file.path(output_dir, "network", "network_metrics.csv")
  if (!file.exists(path)) return(invisible(NULL))
  m <- utils::read.csv(path)
  if (igraph::vcount(g) != m$nodes[1] || igraph::ecount(g) != m$edges[1]) {
    warning("Figure graph (", igraph::vcount(g), " nodes, ", igraph::ecount(g),
            " edges) differs from network_metrics.csv (", m$nodes[1], ", ", m$edges[1], ").")
  }
}

plot_arg_network <- function(output_dir = "output",
                             export = FALSE, neighbours_per_gene = 5L,
                             key_genes = NULL, show_layers = NULL,
                             file_stem = "ARG_network_two_panel_FINAL") {
  # Panel A: complete multilayer network (all original nodes / edges)
  # Exactly ONE shared legend (cowplot); no figure title; no label overlap by design.
  # Required CSVs (from the ORIGINAL v4 network pipeline):
  #   AMR_reanalysis_outputs_v4/network/network_nodes.csv
  #   AMR_reanalysis_outputs_v4/network/network_edges.csv
  # Preview with plot_arg_network(export = FALSE).
  # Save with plot_arg_network(export = TRUE).
  # This script changes only the DISPLAY, never original network data/centralities.
  # Existing sequence-to-ARO annotation discrepancies are not corrected here.

  OUTPUT_DIR <- output_dir
  EXPORT <- export
  SEED <- 1234
  NEIGHBOURS_PER_GENE <- neighbours_per_gene        # Usually 4-6 gives readable Panel B
  SHOW_NEIGHBOUR_TO_NEIGHBOUR_EDGES <- FALSE

  KEY_GENES <- if (is.null(key_genes)) read_prioritized_genes(OUTPUT_DIR) else key_genes

  EDGE_PALETTE <- c(
    Ecological = "#13A995",     # teal
    Mechanism  = "#4E8CCC",     # blue
    DrugClass  = "#E8A050"      # orange
  )
  if (!is.null(show_layers)) EDGE_PALETTE <- EDGE_PALETTE[names(EDGE_PALETTE) %in% show_layers]
  NODE_PALETTE <- c(
    "Other ARGs"       = "#B6E3EC",
    "Prioritized ARGs" = "#E9507B"
  )
  BORDER_PALETTE <- c(
    "Other ARGs"       = "#4A7180",
    "Prioritized ARGs" = "#812746"
  )

  needed <- c("igraph", "ggraph", "ggplot2", "ggrepel", "cowplot")
  not_installed <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(not_installed)) {
    stop("Please install these packages first: install.packages(c(",
         paste(sprintf('"%s"', not_installed), collapse = ", "), "))")
  }

  # 3. READ ORIGINAL NETWORK TABLES --------------------------------------------
  network_dir <- file.path(OUTPUT_DIR, "network")
  edge_file <- file.path(network_dir, "network_edges.csv")
  node_file <- file.path(network_dir, "network_nodes.csv")

  if (!file.exists(edge_file) || !file.exists(node_file)) {
    stop("Cannot find the existing network CSV files. Set OUTPUT_DIR correctly.\n",
         "Expected: \n  ", edge_file, "\n  ", node_file)
  }

  edges <- utils::read.csv(edge_file, check.names = FALSE, stringsAsFactors = FALSE)
  nodes <- utils::read.csv(node_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("from", "to", "layer") %in% names(edges))) {
    stop("network_edges.csv must contain columns: from, to, layer")
  }
  if (!"node" %in% names(nodes)) {
    stop("network_nodes.csv must contain a 'node' column")
  }

  if (!is.null(show_layers)) edges <- edges[edges$layer %in% show_layers, , drop = FALSE]
  edges$from <- trimws(as.character(edges$from))
  edges$to <- trimws(as.character(edges$to))
  edges$layer <- trimws(as.character(edges$layer))
  nodes$node <- trimws(as.character(nodes$node))

  if (anyNA(nodes$node) || any(!nzchar(nodes$node)) || anyDuplicated(nodes$node)) {
    stop("Node identifiers must be nonmissing, nonempty and unique")
  }
  if (anyNA(edges$from) || anyNA(edges$to) || anyNA(edges$layer) ||
      any(!nzchar(edges$from)) || any(!nzchar(edges$to))) {
    stop("Edge endpoints and layer values cannot be empty or NA")
  }
  unknown_endpoints <- setdiff(unique(c(edges$from, edges$to)), nodes$node)
  if (length(unknown_endpoints)) {
    stop("These endpoints are missing from network_nodes.csv: ",
         paste(unknown_endpoints, collapse = ", "))
  }
  unknown_layers <- setdiff(unique(edges$layer), names(EDGE_PALETTE))
  if (length(unknown_layers)) {
    stop("Unrecognized edge types: ", paste(unknown_layers, collapse = ", "))
  }
  if (!"weight" %in% names(edges)) edges$weight <- 1
  edges$weight <- suppressWarnings(as.numeric(edges$weight))

  vertices <- nodes[, c("node", setdiff(names(nodes), "node")), drop = FALSE]
  names(vertices)[1] <- "name"
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
  # Important: do not simplify 'g'; separate evidence layers may form parallel edges.

  missing_focal <- setdiff(KEY_GENES, igraph::V(g)$name)
  if (length(missing_focal)) {
    stop("Priority gene label(s) not found: ", paste(missing_focal, collapse = ", "),
         ". Check spelling in network_nodes.csv.")
  }

  message("Panel A — full network: ", igraph::vcount(g), " nodes, ",
          igraph::ecount(g), " edges.")
  print(table(igraph::E(g)$layer))
  if (is.null(show_layers)) check_against_metrics(OUTPUT_DIR, g)

  # 4. BUILD A TRANSPARENT, READABLE FOCUS VIEW (DISPLAY ONLY) ------------------
  # then by stronger observed ecological rho. This is only a plotting decision;
  edge_tbl <- igraph::as_data_frame(g, what = "edges")
  visible <- KEY_GENES

  for (key in KEY_GENES) {
    hit <- edge_tbl[edge_tbl$from == key | edge_tbl$to == key, , drop = FALSE]
    if (!nrow(hit)) next
    hit$neighbour <- ifelse(hit$from == key, hit$to, hit$from)
    hit <- hit[!(hit$neighbour %in% KEY_GENES), , drop = FALSE]
    if (!nrow(hit)) next

    rho <- suppressWarnings(as.numeric(hit$weight))
    rho[!is.finite(rho)] <- 0
    hit$display_score <- ifelse(hit$layer == "Ecological",
                                2 + pmax(0, pmin(1, rho)), 1)
    neighbour_rank <- stats::aggregate(
      display_score ~ neighbour, data = hit, FUN = max
    )
    neighbour_rank <- neighbour_rank[
      order(-neighbour_rank$display_score, neighbour_rank$neighbour),
      , drop = FALSE
    ]
    visible <- c(visible, utils::head(neighbour_rank$neighbour, NEIGHBOURS_PER_GENE))
  }
  visible <- unique(visible)

  zoom_edges_keep <- edge_tbl$from %in% visible & edge_tbl$to %in% visible
  if (!SHOW_NEIGHBOUR_TO_NEIGHBOUR_EDGES) {
    zoom_edges_keep <- zoom_edges_keep &
      (edge_tbl$from %in% KEY_GENES | edge_tbl$to %in% KEY_GENES)
  }
  zoom_edges <- edge_tbl[zoom_edges_keep, , drop = FALSE]
  zoom_vertices <- igraph::as_data_frame(g, what = "vertices")
  zoom_vertices <- zoom_vertices[zoom_vertices$name %in% visible, , drop = FALSE]
  g_zoom <- igraph::graph_from_data_frame(
    zoom_edges, directed = FALSE, vertices = zoom_vertices
  )
  message("Panel B — selected focus view: ", igraph::vcount(g_zoom),
          " nodes, ", igraph::ecount(g_zoom), " shown edges.")

  # Use a simple undirected backbone only to position nodes. The plotted graph
  # still contains the original layer-specific parallel edges.
  manual_layout <- function(graph, seed = 1234, iterations = 2700L) {
    layout_graph <- igraph::simplify(graph, remove.multiple = TRUE,
                                    remove.loops = TRUE)
    set.seed(seed)
    xy <- igraph::layout_with_fr(
      layout_graph,
      niter = iterations,
      weights = rep(1, igraph::ecount(layout_graph))
    )
    rownames(xy) <- igraph::V(layout_graph)$name
    # Isolated nodes would otherwise be placed far from the network and squeeze
    # the plot; put them in a row just below the connected nodes.
    isolated <- igraph::degree(layout_graph) == 0
    if (any(isolated) && any(!isolated)) {
      rng <- apply(xy[!isolated, , drop = FALSE], 2, range)
      xy[isolated, 1] <- seq(rng[1, 1], rng[2, 1], length.out = sum(isolated) + 2)[-c(1, sum(isolated) + 2)]
      xy[isolated, 2] <- rng[1, 2] - 0.08 * diff(rng[, 2])
    }
    xy <- xy[igraph::V(graph)$name, , drop = FALSE]
    ggraph::create_layout(graph, layout = "manual", x = xy[, 1], y = xy[, 2])
  }

  make_panel <- function(graph, seed, full = FALSE, show_legends = FALSE) {
    focal <- igraph::V(graph)$name %in% KEY_GENES
    igraph::V(graph)$node_type <- factor(
      ifelse(focal, "Prioritized ARGs", "Other ARGs"),
      levels = c("Other ARGs", "Prioritized ARGs")
    )
    # These are visual sizes, not measurements of biological importance.
    degree_display <- igraph::degree(graph)
    igraph::V(graph)$dot_size <- ifelse(
      focal, 6.5, pmin(4.65, 3.1 + 0.28 * sqrt(degree_display))
    )

    xy <- manual_layout(graph, seed = seed,
                        iterations = if (full) 3300L else 2700L)
    only_labels <- xy[xy$name %in% KEY_GENES, , drop = FALSE]

    p <- ggraph::ggraph(xy) +
      # Edge widths are fixed: do not treat correlation rho and binary annotations
      ggraph::geom_edge_fan(
        ggplot2::aes(edge_colour = layer),
        edge_width = if (full) 0.38 else 0.52,
        alpha = if (full) 0.43 else 0.62,
        show.legend = show_legends
      ) +
      ggraph::scale_edge_colour_manual(
        values = EDGE_PALETTE,
        breaks = names(EDGE_PALETTE),
        drop = FALSE,
        name = "Edge type"
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(size = dot_size, fill = node_type, colour = node_type),
        shape = 21, stroke = if (full) 0.68 else 0.80,
        alpha = 1,
        show.legend = show_legends
      ) +
      ggplot2::scale_size_identity(guide = "none") +
      ggplot2::scale_fill_manual(
        values = NODE_PALETTE,
        breaks = names(NODE_PALETTE),
        drop = FALSE,
        name = "Node type"
      ) +
      ggplot2::scale_colour_manual(values = BORDER_PALETTE, guide = "none") +
      ggraph::geom_node_label(
        data = only_labels,
        ggplot2::aes(label = name),
        repel = TRUE,
        seed = seed,
        size = if (full) 3.25 else 3.65,
        fontface = "bold",
        colour = "#21374A",
        fill = "white",
        label.size = 0.20,
        label.padding = grid::unit(0.16, "lines"),
        box.padding = grid::unit(0.50, "lines"),
        point.padding = grid::unit(0.25, "lines"),
        max.overlaps = Inf,
        max.time = 5,
        force = 4,
        segment.color = "#748493",
        segment.size = 0.35,
        show.legend = FALSE
      ) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.13)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.13)) +
      ggplot2::coord_equal(clip = "off") +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(16, 18, 16, 18),
        legend.title = ggplot2::element_text(face = "bold", size = 10.5),
        legend.text = ggplot2::element_text(size = 9.5),
        legend.key.width = grid::unit(1.5, "cm"),
        legend.key.height = grid::unit(0.55, "cm"),
        legend.spacing.y = grid::unit(0.04, "cm"),
        legend.box = "vertical",
        legend.box.just = "center",
        legend.direction = "horizontal"
      )

    if (show_legends) {
      p <- p + ggplot2::guides(
        edge_colour = ggplot2::guide_legend(
          order = 1, nrow = 1,
          override.aes = list(edge_width = 2.5, edge_alpha = 1, alpha = 1)
        ),
        fill = ggplot2::guide_legend(
          order = 2, nrow = 1,
          override.aes = list(
            shape = 21, size = 3.4,
            colour = unname(BORDER_PALETTE), stroke = 0.75
          )
        )
      ) + ggplot2::theme(legend.position = "bottom")
    } else {
      p <- p + ggplot2::theme(legend.position = "none")
    }
    p
  }

  plot_A <- make_panel(g, seed = SEED, full = TRUE, show_legends = FALSE)
  plot_B <- make_panel(g_zoom, seed = SEED + 1L,
                       full = FALSE, show_legends = FALSE)
  legend_source <- make_panel(g, seed = SEED, full = TRUE, show_legends = TRUE)
  shared_legend <- cowplot::get_legend(legend_source)
  if (is.null(shared_legend)) stop("Could not build shared legend")

  combined_panels <- cowplot::plot_grid(
    plot_A, plot_B,
    nrow = 1,
    rel_widths = c(1.00, 1.10),
    labels = c("A", "B"),
    label_size = 18,
    label_fontface = "bold",
    label_colour = "#1E354C",
    hjust = -0.25,
    vjust = 1.45,
    align = "hv"
  )

  figure <- cowplot::plot_grid(
    combined_panels,
    shared_legend,
    ncol = 1,
    rel_heights = c(1, 0.16)
  )

  print(figure)

  # 9. EXPORT ONLY AFTER SETTING EXPORT <- TRUE --------------------------------
  if (EXPORT) {
    figure_dir <- file.path(OUTPUT_DIR, "figures")
    dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

    ggplot2::ggsave(
      filename = file.path(figure_dir, paste0(file_stem, ".pdf")),
      plot = figure, width = 15.2, height = 8.4, units = "in",
      device = "pdf", useDingbats = FALSE, bg = "white"
    )
    ggplot2::ggsave(
      filename = file.path(figure_dir, paste0(file_stem, "_450dpi.png")),
      plot = figure, width = 15.2, height = 8.4, units = "in",
      dpi = 450, bg = "white"
    )
    ggplot2::ggsave(
      filename = file.path(figure_dir, paste0(file_stem, "_600dpi.tiff")),
      plot = figure, width = 15.2, height = 8.4, units = "in",
      dpi = 600, compression = "lzw", bg = "white"
    )
    message("Saved final figure to: ", normalizePath(figure_dir))
  }

  # SUGGESTED FIGURE 2 CAPTION is supplied in accompanying text file.

  invisible(figure)
}

# Optional single-panel figure and Cytoscape files.
# Reads the saved edges and nodes; no network fitting here.
export_network_for_cytoscape <- function(output_dir = "output",
                                         key_genes = NULL) {
  # 05_export_ARG_network_publication_revised.R
  # Publication-ready multilayer ARG network + Cytoscape exports
  # Input: <output_dir>/network/network_edges.csv and network_nodes.csv
  # Exports: PNG (300 dpi), LZW TIFF (300 dpi), GraphML, SIF, node/edge CSV
  # IMPORTANT: Preserves all parallel ecological/mechanism/drug-class edges.

  needed <- c("data.table", "igraph", "ggplot2", "ggraph", "ggrepel")
  for (pkg in needed) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Please install package first: ", pkg)
    }
  }

  # Edit ONLY this path to match your v4 analysis output directory.

  network_dir <- file.path(output_dir, "network")
  figure_dir  <- file.path(output_dir, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  edge_file <- file.path(network_dir, "network_edges.csv")
  node_file <- file.path(network_dir, "network_nodes.csv")
  if (!file.exists(edge_file) || !file.exists(node_file)) {
    stop(
      "Missing required v4 network files.\n",
      "Expected:\n  ", edge_file, "\n  ", node_file, "\n",
      "Run the v4 network analysis first, or change output_dir."
    )
  }

  # 2. Import and validate the ORIGINAL multilayer network --------------------
  edges <- data.table::fread(edge_file, na.strings = c("", "NA"))
  nodes <- data.table::fread(node_file, na.strings = c("", "NA"))

  stopifnot(all(c("from", "to", "layer") %in% names(edges)))
  stopifnot("node" %in% names(nodes))

  edges[, `:=`(
    from = trimws(as.character(from)),
    to = trimws(as.character(to)),
    layer = trimws(as.character(layer))
  )]
  nodes[, node := trimws(as.character(node))]

  if (anyNA(nodes$node) || anyDuplicated(nodes$node)) {
    stop("Node identifiers must be non-missing and unique (nodes$node).")
  }
  if (anyNA(edges$from) || anyNA(edges$to) || anyNA(edges$layer)) {
    stop("The edge table contains missing from/to/layer values.")
  }
  unknown_nodes <- setdiff(unique(c(edges$from, edges$to)), nodes$node)
  if (length(unknown_nodes)) {
    stop("Edge endpoints absent from node table: ", paste(unknown_nodes, collapse = ", "))
  }

  edge_layers <- c("Ecological", "Mechanism", "DrugClass")
  unknown_layers <- setdiff(unique(edges$layer), edge_layers)
  if (length(unknown_layers)) {
    stop("Unknown edge layer(s): ", paste(unknown_layers, collapse = ", "))
  }

  if (!"weight" %in% names(edges)) edges[, weight := 1]
  edges[, weight := suppressWarnings(as.numeric(weight))]
  # Original weight and significance fields are preserved for Cytoscape.
  # abs_weight is only a convenience attribute; figure edge widths are uniform
  # so annotation weights (=1) are not mistaken for comparable correlations.
  edges[, abs_weight := abs(weight)]

  node_cols <- c("node", setdiff(names(nodes), "node"))
  vertices <- as.data.frame(nodes[, ..node_cols])

  g <- igraph::graph_from_data_frame(
    d = as.data.frame(edges),
    directed = FALSE,
    vertices = vertices
  )
  # DO NOT simplify(g): it would merge distinct layer-specific parallel edges.

  cat("Network:", igraph::vcount(g), "nodes,", igraph::ecount(g), "edges\n")
  print(table(igraph::E(g)$layer))
  check_against_metrics(output_dir, g)

  priority_genes <- if (is.null(key_genes)) read_prioritized_genes(output_dir) else key_genes
  priority_label <- paste0("Prioritized ARGs (n = ", length(priority_genes), ")")
  missing_priority <- setdiff(priority_genes, igraph::V(g)$name)
  if (length(missing_priority)) {
    stop(
      "These prioritized gene names are missing from the node table: ",
      paste(missing_priority, collapse = ", "),
      "\nCheck the ARO node names in network_nodes.csv."
    )
  }

  igraph::V(g)$is_priority <- igraph::V(g)$name %in% priority_genes
  igraph::V(g)$node_group <- ifelse(
    igraph::V(g)$is_priority, priority_label, "Other ARGs"
  )
  igraph::V(g)$degree_calc <- igraph::degree(g)
  igraph::V(g)$plot_size <- ifelse(
    igraph::V(g)$is_priority,
    7.0 + 0.23 * sqrt(igraph::V(g)$degree_calc),
    2.4 + 0.46 * sqrt(igraph::V(g)$degree_calc)
  )
  igraph::V(g)$label_plot <- ifelse(
    igraph::V(g)$is_priority, igraph::V(g)$name, ""
  )

  # One reproducible layout is used for all figure exports. Isolated nodes are
  set.seed(1234)
  layout_tbl <- ggraph::create_layout(g, layout = "fr", niter = 2500)

  edge_palette <- c(
    Ecological = "#149B70",  # emerald/teal
    Mechanism  = "#2775D9",  # bright blue
    DrugClass  = "#E98635"   # warm orange
  )
  node_palette <- stats::setNames(c("#A7D4DD", "#ED4661"), c("Other ARGs", priority_label))

  p <- ggraph::ggraph(layout_tbl) +
    # All original edges retained, including parallel multilayer connections.
    ggraph::geom_edge_link(
      ggplot2::aes(edge_colour = layer),
      edge_width = 0.53,
      alpha = 0.70,
      show.legend = TRUE
    ) +
    ggraph::scale_edge_colour_manual(
      values = edge_palette,
      breaks = c("Ecological", "Mechanism", "DrugClass"),
      name = "Edge layer"
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(size = plot_size, fill = node_group),
      shape = 21,
      colour = "#24475A",
      stroke = 0.45,
      alpha = 0.98,
      show.legend = TRUE
    ) +
    ggplot2::scale_size_identity(guide = "none") +
    ggplot2::scale_fill_manual(
      values = node_palette,
      breaks = c(priority_label, "Other ARGs"),
      name = "Nodes"
    ) +
    # Repelled white labels, ONLY for the prioritized ARGs.
    ggraph::geom_node_label(
      ggplot2::aes(label = label_plot),
      repel = TRUE,
      size = 4.2,
      fontface = "bold",
      colour = "#92233A",
      fill = "#FFFFFF",
      label.size = 0.25,
      label.padding = grid::unit(0.22, "lines"),
      box.padding = grid::unit(0.75, "lines"),
      point.padding = grid::unit(0.40, "lines"),
      force = 3,
      max.overlaps = Inf,
      max.time = 4,
      seed = 1234,
      segment.color = "#687C8A",
      segment.size = 0.45,
      show.legend = FALSE
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1, override.aes = list(size = 5)),
      edge_colour = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title = ggplot2::element_text(
        family = "sans", face = "bold", size = 16, colour = "#1C3445",
        hjust = 0.5, margin = ggplot2::margin(b = 5)
      ),
      plot.subtitle = ggplot2::element_text(
        family = "sans", size = 11, colour = "#506777", hjust = 0.5
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(face = "bold", colour = "#254358"),
      legend.text = ggplot2::element_text(size = 10, colour = "#2F4150"),
      legend.key.width = grid::unit(1.2, "cm"),
      plot.margin = ggplot2::margin(t = 20, r = 65, b = 16, l = 65)
    )

  if (interactive()) print(p)

  # 7. Export 300-dpi figures -------------------------------------------------
  png_file  <- file.path(figure_dir, "ARG_network_publication_300dpi.png")
  tiff_file <- file.path(figure_dir, "ARG_network_publication_300dpi.tiff")

  ggplot2::ggsave(
    filename = png_file,
    plot = p,
    width = 13, height = 9, units = "in", dpi = 300,
    bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = tiff_file,
    plot = p,
    width = 13, height = 9, units = "in", dpi = 300,
    compression = "lzw", bg = "white", limitsize = FALSE
  )

  # 8. Export Cytoscape-ready multilayer data ---------------------------------
  graphml_file <- file.path(network_dir, "ARG_network_for_Cytoscape.graphml")
  igraph::write_graph(g, graphml_file, format = "graphml")

  # Plain tables help modify styling and audit annotations in Cytoscape.
  node_export <- data.table::as.data.table(
    igraph::as_data_frame(g, what = "vertices")
  )
  edge_export <- data.table::as.data.table(
    igraph::as_data_frame(g, what = "edges")
  )
  data.table::fwrite(
    node_export,
    file.path(network_dir, "ARG_network_Cytoscape_nodes.csv")
  )
  data.table::fwrite(
    edge_export,
    file.path(network_dir, "ARG_network_Cytoscape_edges.csv")
  )

  # SIF is offered for basic compatibility, but it cannot retain edge weights,
  # source-specific q values, or the full attribute structure like GraphML.
  sif_export <- edge_export[, .(from, layer, to)]
  data.table::fwrite(
    sif_export,
    file.path(network_dir, "ARG_network_Cytoscape.sif"),
    sep = "\t", col.names = FALSE
  )

  cat("\nExport complete:\n")
  cat("PNG:     ", normalizePath(png_file), "\n")
  cat("TIFF:    ", normalizePath(tiff_file), "\n")
  cat("GraphML: ", normalizePath(graphml_file), "\n")
  cat("Cytoscape node, edge and SIF files:", network_dir, "\n")
  cat("Prioritized nodes:", paste(priority_genes, collapse = "; "), "\n")

  invisible(graphml_file)
}
