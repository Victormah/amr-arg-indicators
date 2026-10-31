# =============================================================================
# scripts/make_manuscript_figures.R - manuscript Figs 2-4 and Supplementary Fig. S1
# Project : Candidate antimicrobial resistance gene indicators for veterinary
#           public health surveillance (Mahawan et al., Research in Veterinary
#           Science, submitted)
#
# Purpose : Draw the data figures from the figure-data tables and export each as
#           TIFF (LZW) and PNG at 300 dpi, at Elsevier column widths.
#           (Fig. 1, the workflow diagram, was drawn in Canva: docs/Fig1_workflow.*)
# Run     : called at the end of run_all.R, or on its own from the project folder:
#             source("scripts/make_manuscript_figures.R")
#           On its own it uses the archived tables in results/figure_data/, so
#           the manuscript figures can be redrawn without re-running the models.
# Inputs  : FIG_DATA_DIR (default results/figure_data/)
#   regression_model_summary.csv           -> Fig. 2
#   classification_predictions_sample.csv  -> Fig. 3 (sample-level CV)
#   classification_predictions_project.csv -> Fig. 3 (project-blocked CV)
#   sample_flow.csv                        -> Fig. 2 (sample numbers)
#   network_nodes.csv, network_edges.csv   -> Fig. 4, Fig. S1
# Outputs : OUT_DIR (default output/figures_manuscript/)
#   Fig2_burden_reconstruction, Fig3_source_classification,
#   Fig4_ecological_network, FigS1_all_layer_network (.tiff and .png)
# Packages: data.table, ggplot2, patchwork, igraph, ggraph, ggrepel, ragg, here
# =============================================================================

## ---- settings you may change ----------------------------------------------
FIG_DATA_DIR <- if (exists("FIG_DATA_DIR")) FIG_DATA_DIR else here::here("results", "figure_data")
NETWORK_DIR  <- if (exists("NETWORK_DIR"))  NETWORK_DIR  else FIG_DATA_DIR
OUT_DIR      <- if (exists("OUT_DIR"))      OUT_DIR      else here::here("output", "figures_manuscript")
DPI          <- 300
FONT         <- "sans"   # Arial on Windows, Helvetica on macOS
SEED         <- 1234     # fixes the network layout

## Elsevier figure widths (mm): single column 90, 1.5 column 140, full width 190

## ---- packages ---------------------------------------------------------------
pkgs <- c("data.table", "ggplot2", "patchwork", "igraph", "ggraph", "ggrepel", "ragg", "here")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                          ". Install them first (see README).")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(igraph); library(ggraph)
})
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- shared style -----------------------------------------------------------
SOURCES  <- c("Avian", "Compost", "Feces", "Food", "Wastewater")
INK      <- "#1f2a33"; MUTED <- "#5b6670"; GRID <- "#e3e6e8"
BLUE     <- "#0072B2"; LIGHT <- "#e8f1f8"
# Okabe-Ito palette (colour-blind safe)
SRC_COL  <- c(Avian = "#E69F00", Compost = "#56B4E9", Feces = "#009E73",
              Food = "#CC79A7", Wastewater = "#0072B2")
MODEL_COL <- c("Elastic net" = "#0072B2", "Random forest" = "#D55E00",
               "Support vector regression" = "#009E73", "XGBoost" = "#CC79A7")

theme_ms <- function(base = 8) {
  theme_minimal(base_size = base, base_family = FONT) +
    theme(text = element_text(colour = INK),
          axis.text = element_text(colour = INK),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
          axis.line.x = element_line(colour = MUTED, linewidth = 0.3),
          plot.title = element_text(face = "bold", size = base + 1, hjust = 0),
          plot.title.position = "plot",
          legend.position = "bottom",
          plot.background = element_rect(fill = "white", colour = NA))
}

# Export one figure as TIFF (LZW) and PNG at the same size and resolution
save_fig <- function(plot, name, width_mm, height_mm, dpi = DPI) {
  tif <- file.path(OUT_DIR, paste0(name, ".tiff"))
  png <- file.path(OUT_DIR, paste0(name, ".png"))
  ragg::agg_tiff(tif, width = width_mm, height = height_mm, units = "mm",
                 res = dpi, compression = "lzw", background = "white")
  print(plot); invisible(dev.off())
  ragg::agg_png(png, width = width_mm, height = height_mm, units = "mm",
                res = dpi, background = "white")
  print(plot); invisible(dev.off())
  message(sprintf("saved %s (.tiff, .png): %d x %d mm, %d dpi",
                  name, width_mm, height_mm, dpi))
}

# ResFinder names use ASCII primes; show proper prime symbols in figures
pretty_gene <- function(x) gsub("'", "\u2032", gsub("''", "\u2033", x, fixed = TRUE), fixed = TRUE)

rd <- function(f, dir = FIG_DATA_DIR) {
  p <- file.path(dir, f)
  if (!file.exists(p)) stop("File not found: ", p, call. = FALSE)
  fread(p)
}

flow <- rd("sample_flow.csv")
n_mod <- setNames(flow$modeled_samples, flow$display_name)

## =============================================================================
## Figure 1. Study workflow
## =============================================================================
make_fig1 <- function() {
  box <- function(x, y, w, h, title, body, fill = "white", txt = INK, line = BLUE) {
    list(
      annotate("rect", xmin = x, xmax = x + w, ymin = y, ymax = y + h,
               fill = fill, colour = line, linewidth = 0.4),
      annotate("text", x = x + w / 2, y = y + h - 2, label = title, family = FONT,
               fontface = "bold", size = 2.9, colour = txt, vjust = 1),
      annotate("text", x = x + w / 2, y = y + h - 6.2, label = body, family = FONT,
               size = 2.4, colour = txt, vjust = 1, lineheight = 1.05))
  }
  arr <- function(x, y, xend, yend)
    annotate("segment", x = x, y = y, xend = xend, yend = yend, colour = MUTED,
             linewidth = 0.4, arrow = arrow(length = unit(1.6, "mm"), type = "closed"))
  seg <- function(x, y, xend, yend)
    annotate("segment", x = x, y = y, xend = xend, yend = yend, colour = MUTED, linewidth = 0.4)

  n_elig <- format(sum(flow$metadata_eligible_samples), big.mark = ",")
  n_used <- format(sum(flow$modeled_samples), big.mark = ",")

  ggplot() + xlim(0, 100) + ylim(0, 76) + theme_void(base_family = FONT) +
    theme(plot.background = element_rect(fill = "white", colour = NA)) +
    box(1, 58, 30, 16, "1  Public metagenomes",
        "Curated 214K-metagenome\nresource; five sources:\nAvian, Compost, Feces,\nFood, Wastewater") +
    box(35, 58, 30, 16, "2  Sample selection",
        paste0(n_elig, " eligible samples\n", n_used, " with resistance-gene\nreads and a bacterial\nread count")) +
    box(69, 58, 30, 16, "3  Normalization",
        "Gene reads per million\nbacterial reads; log scale\nfor models, CLR for\nnetworks") +
    arr(31, 66, 35, 66) + arr(65, 66, 69, 66) +
    seg(84, 58, 84, 53) + seg(25, 53, 84, 53) + arr(25, 53, 25, 50) + arr(75, 53, 75, 50) +
    box(4, 31, 42, 19, "4a  Burden reconstruction",
        "Within each source: can a few genes\nstand in for the total resistance signal?\n5-fold cross-validation x 5 repeats\nOutput: held-out R\u00b2; genes selected\nin \u2265 50% of training splits", fill = LIGHT) +
    box(54, 31, 42, 19, "4b  Source classification",
        "Across sources: do gene profiles show\nwhere a sample came from?\nNested cross-validation x 5 repeats\nOutput: balanced accuracy, recall\nCheck: whole projects held out", fill = LIGHT) +
    arr(25, 31, 25, 27) +
    box(4, 6, 42, 21, "5  Gene network",
        "Nodes: 74 repeatedly selected genes\nLinks: genes that rise and fall together\nwithin a source (\u03c1 \u2265 0.30, FDR \u2264 0.05)\nCARD mechanism and drug-class\nlinks shown for description only") +
    box(54, 6, 42, 21, "6  Candidate indicator genes",
        "Selected in \u2265 2 sources and in the top\nquartile of network betweenness\n\u2192 7 genes (6 loci)\n\u2192 candidate qPCR panel for\nveterinary surveillance",
        fill = BLUE, txt = "white") +
    arr(46, 16.5, 54, 16.5)
}
# Fig. 1 (final version) was drawn in Canva (docs/). An R draft is kept here:
if (isTRUE(get0("DRAW_FIG1_DRAFT"))) save_fig(make_fig1(), "Fig1_workflow_draft", 190, 140)

## =============================================================================
## Figure 2. Within-source reconstruction of total ARG burden (held-out R²)
## =============================================================================
make_fig2 <- function() {
  reg <- rd("regression_model_summary.csv")
  reg <- reg[metric == "R_squared" & model %in% c("Elastic Net", "Random Forest", "SVM", "XGBoost")]
  reg[, model := factor(c("Elastic Net" = "Elastic net", "Random Forest" = "Random forest",
                          "SVM" = "Support vector regression", "XGBoost" = "XGBoost")[model],
                        levels = names(MODEL_COL))]
  ord <- c("Compost", "Avian", "Food", "Feces", "Wastewater")
  reg[, src := factor(sprintf("%s (n = %s)", display_name, format(n_mod[display_name], big.mark = ",", trim = TRUE)),
                      levels = sprintf("%s (n = %s)", ord, format(n_mod[ord], big.mark = ",", trim = TRUE)))]
  pd <- position_dodge(width = 0.65)
  ggplot(reg, aes(x = estimate, y = src, colour = model, shape = model)) +
    geom_vline(xintercept = 0, colour = "#9aa3ab", linewidth = 0.3) +
    geom_linerange(aes(xmin = lower, xmax = upper), position = pd, linewidth = 0.5) +
    geom_point(position = pd, size = 2.1, stroke = 0.3, fill = "white") +
    scale_colour_manual(values = MODEL_COL, name = NULL) +
    scale_shape_manual(values = c(16, 15, 18, 17), name = NULL) +
    scale_x_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = "Held-out R\u00b2 (mean of 5 repeated cross-validation runs; lines, repeated-CV interval)",
         y = NULL) +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
    theme_ms()
}
save_fig(make_fig2(), "Fig2_burden_reconstruction", 140, 85)

## =============================================================================
## Figure 3. Source classification: confusion matrix and recall by source
## =============================================================================
recall_tab <- function(file, label) {
  p <- rd(file)[model %in% c("Elastic Net", "Random Forest") & !is.na(predicted)]
  p[, .(recall = mean(actual == predicted)), by = .(model, actual)][, validation := label][]
}
make_fig3 <- function() {
  ps <- rd("classification_predictions_sample.csv")[model == "Elastic Net" & !is.na(predicted)]
  cm <- ps[, .N, by = .(actual, predicted)]
  cm <- merge(CJ(actual = SOURCES, predicted = SOURCES), cm, all.x = TRUE)
  cm[is.na(N), N := 0][, prop := N / sum(N), by = actual]
  cm[, actual := factor(actual, levels = rev(SOURCES))][, predicted := factor(predicted, levels = SOURCES)]
  pA <- ggplot(cm, aes(predicted, actual, fill = prop)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(prop >= 0.005, sprintf("%.2f", prop), ""),
                  colour = prop > 0.55), size = 2.4, family = FONT) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
    scale_fill_gradient(low = "#f4f8fb", high = "#08306b", limits = c(0, 1),
                        name = "Proportion of\nrecorded source") +
    labs(title = "A  Elastic net, sample-level CV", x = "Predicted source", y = "Recorded source") +
    coord_fixed() + theme_ms() +
    theme(panel.grid = element_blank(), axis.line.x = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "right", legend.key.height = unit(5, "mm"),
          legend.title = element_text(size = 6.5))

  rc <- rbind(recall_tab("classification_predictions_sample.csv", "sample-level CV"),
              recall_tab("classification_predictions_project.csv", "project-blocked CV"))
  rc[, model := c("Elastic Net" = "Elastic net", "Random Forest" = "Random forest")[model]]
  rc[, grp := factor(paste(model, validation, sep = ", "),
                     levels = c("Elastic net, sample-level CV", "Elastic net, project-blocked CV",
                                "Random forest, sample-level CV", "Random forest, project-blocked CV"))]
  rc[, y := as.numeric(factor(actual, levels = rev(SOURCES))) + ifelse(model == "Elastic net", 0.15, -0.15)]
  pB <- ggplot(rc, aes(recall, y)) +
    geom_vline(xintercept = 0.2, linetype = 2, colour = "#9aa3ab", linewidth = 0.3) +
    annotate("text", x = 0.215, y = 0.55, label = "chance", hjust = 0, size = 2.2,
             colour = MUTED, family = FONT) +
    geom_line(aes(group = interaction(model, actual), colour = model), alpha = 0.6, linewidth = 0.35) +
    geom_point(aes(colour = model, shape = grp), size = 1.9, stroke = 0.6, fill = "white") +
    scale_colour_manual(values = c("Elastic net" = BLUE, "Random forest" = "#D55E00"), guide = "none") +
    scale_shape_manual(values = c(16, 21, 15, 22), name = NULL) +
    guides(shape = guide_legend(ncol = 1, override.aes = list(
      colour = c(BLUE, BLUE, "#D55E00", "#D55E00")))) +
    scale_y_continuous(breaks = 1:5, labels = rev(SOURCES), limits = c(0.5, 5.5)) +
    scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.2)) +
    labs(title = "B  Recall by source", x = "Recall (proportion classified correctly)", y = NULL) +
    theme_ms() + theme(legend.text = element_text(size = 6.5), legend.key.height = unit(3.2, "mm"))
  pA + pB + plot_layout(widths = c(1, 1.05))
}
save_fig(make_fig3(), "Fig3_source_classification", 190, 105)

## =============================================================================
## Figure 4. Ecological gene network (A) and prioritized genes with neighbours (B)
## Supplementary Fig. S1. Network with all three edge layers
## =============================================================================
nodes <- rd("network_nodes.csv", NETWORK_DIR)
edges <- rd("network_edges.csv", NETWORK_DIR)
nodes[, label := pretty_gene(node)]
nodes[, eco_degree := degree_centrality_layers]
nodes[, type := factor(ifelse(prioritized_central, "Prioritized gene", "Other selected gene"),
                       levels = c("Other selected gene", "Prioritized gene"))]

layout_xy <- function(g, seed) {
  set.seed(seed)
  # stronger correlations pull genes closer; isolated genes placed in a row below
  xy <- igraph::layout_with_fr(g, weights = igraph::E(g)$weight, niter = 3000)
  iso <- igraph::degree(g) == 0
  if (any(iso) && any(!iso)) {
    rng <- range(xy[!iso, 1]); y0 <- min(xy[!iso, 2]) - 0.12 * diff(range(xy[!iso, 2]))
    xy[iso, 1] <- seq(rng[1], rng[2], length.out = sum(iso) + 2)[-c(1, sum(iso) + 2)]
    xy[iso, 2] <- y0
  }
  xy
}

net_plot <- function(g, seed, title, edge_colour_by = "source", edge_alpha = 0.55,
                     label_all = FALSE) {
  xy <- layout_xy(g, seed)
  lay <- ggraph::create_layout(g, layout = "manual", x = xy[, 1], y = xy[, 2])
  p <- ggraph(lay) +
    geom_edge_link(aes(colour = .data[[edge_colour_by]]), alpha = edge_alpha, width = 0.35) +
    geom_node_point(aes(fill = type, size = type), shape = 21, colour = "#2b3a45", stroke = 0.3) +
    scale_fill_manual(values = c("Other selected gene" = "#bfe3ec", "Prioritized gene" = "#e8496f"),
                      name = NULL, drop = FALSE) +
    scale_size_manual(values = c("Other selected gene" = 1.8, "Prioritized gene" = 3.6),
                      name = NULL, drop = FALSE) +
    ggrepel::geom_label_repel(
      data = function(d) if (label_all) d else d[d$type == "Prioritized gene", ],
      aes(x = x, y = y, label = label), family = FONT, fontface = "italic",
      size = 2.3, label.size = 0.15, label.padding = unit(0.12, "lines"),
      box.padding = 0.35, min.segment.length = 0.2, seed = seed, max.overlaps = Inf) +
    labs(title = title) + theme_void(base_family = FONT) +
    theme(plot.title = element_text(face = "bold", size = 9, colour = INK),
          legend.position = "bottom", legend.text = element_text(size = 6.8),
          legend.title = element_text(size = 7, face = "bold"),
          plot.background = element_rect(fill = "white", colour = NA))
  p
}

make_fig4 <- function() {
  eco <- edges[layer == "Ecological"]
  g <- igraph::graph_from_data_frame(eco[, .(from, to, source, weight)], directed = FALSE,
                                     vertices = nodes[, .(node, label, type, eco_degree)])
  # Panel B: prioritized genes and their direct ecological neighbours
  pri <- nodes[prioritized_central == TRUE, node]
  keep <- unique(c(pri, unlist(igraph::adjacent_vertices(g, pri) |> lapply(names))))
  gB <- igraph::induced_subgraph(g, keep)
  edge_scale <- ggraph::scale_edge_colour_manual(values = SRC_COL, name = "Correlation observed in",
                                                 drop = FALSE, limits = SOURCES[SOURCES != "Compost"])
  pA <- net_plot(g, SEED, "A  Ecological network (74 genes)") + edge_scale
  pB <- net_plot(gB, SEED + 1, "B  Prioritized genes and direct neighbours") + edge_scale
  (pA | pB) + plot_layout(guides = "collect", widths = c(1.15, 1)) &
    theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))
}
save_fig(make_fig4(), "Fig4_ecological_network", 190, 130)

make_figS1 <- function() {
  e <- copy(edges)
  e[, layer_lab := c(Ecological = "Ecological (sample-level correlation)",
                     Mechanism = "Shared CARD resistance mechanism",
                     DrugClass = "Shared CARD drug class")[layer]]
  g <- igraph::graph_from_data_frame(e[, .(from, to, layer_lab, weight)], directed = FALSE,
                                     vertices = nodes[, .(node, label, type, eco_degree)])
  igraph::E(g)$weight <- ifelse(igraph::E(g)$layer_lab == "Ecological (sample-level correlation)",
                                igraph::E(g)$weight, 0.15)
  net_plot(g, SEED, "All-layer network (sensitivity analysis)", edge_colour_by = "layer_lab",
           edge_alpha = 0.25) +
    ggraph::scale_edge_colour_manual(values = c("Ecological (sample-level correlation)" = "#009E73",
                                                "Shared CARD resistance mechanism" = "#9aa3ab",
                                                "Shared CARD drug class" = "#E69F00"),
                                     name = NULL) +
    theme(legend.box = "vertical")
}
save_fig(make_figS1(), "FigS1_all_layer_network", 190, 160)

message("All figures written to: ", normalizePath(OUT_DIR))
