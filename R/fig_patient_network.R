#!/usr/bin/env Rscript
# ============================================================================
# fig_patient_network.R
# ----------------------------------------------------------------------------
# Patient similarity network coloured by metabolic ecotype. Each node is a
# patient; edges link patients whose 12-dimensional metabolic-cluster
# composition vectors correlate at Pearson r >= 0.40. Ecotypes are the
# Louvain communities of that network.
#
# Design choices worth noting:
#   * The graph is taken directly from the input object rather than rebuilt,
#     so the figure reflects exactly the analysed network.
#   * Layout uses graphlayouts::layout_as_backbone, designed for dense
#     networks where stress majorisation fails to separate communities.
#     Stress treats edge weights as LENGTHS (higher weight -> farther apart),
#     which is backwards for a similarity network; backbone avoids this by
#     laying out the network's structural skeleton.
#   * The outlier is excluded from the layout computation and placed manually
#     just outside the main cluster, so a disconnected singleton doesn't get
#     flung to a far corner and open up dead whitespace.
#   * Communities are mapped to ecotype labels programmatically, by MC3/MC8
#     dominance, so the figure does not depend on Louvain's arbitrary numbering.
#
# Reads:   data/example_input.rds
# Writes:  output/patient_network_example.pdf / .png
# Run from the repository root.
# ============================================================================

required_packages <- c("dplyr", "tidyr", "tibble", "ggplot2",
                       "igraph", "ggraph", "graphlayouts",
                       "ggforce", "concaveman")

suppressPackageStartupMessages({
  missing_pkgs <- required_packages[
    !vapply(required_packages,
            function(p) requireNamespace(p, quietly = TRUE),
            logical(1))
  ]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ",
         paste(missing_pkgs, collapse = ", "),
         "\n  Install with: install.packages(c(",
         paste(sprintf("'%s'", missing_pkgs), collapse = ", "), "))",
         call. = FALSE)
  }
  for (p in required_packages) library(p, character.only = TRUE)
})

# ----------------------------------------------------------------------------
# 1. Load
# ----------------------------------------------------------------------------
cat("Loading network object...\n")
typeA <- readRDS("data/example_input.rds")

required <- c("graph", "communities", "patient_communities",
              "mc_matrix", "modularity")
missing  <- required[!required %in% names(typeA)]
if (length(missing) > 0) {
  stop("Input is missing components: ", paste(missing, collapse = ", "))
}

g            <- typeA$graph
comm         <- typeA$communities
patient_meta <- typeA$patient_communities
mod_val      <- typeA$modularity

cat("Network diagnostics:\n")
cat(sprintf("  Nodes:      %d\n", vcount(g)))
cat(sprintf("  Edges:      %d\n", ecount(g)))
cat(sprintf("  Modularity: %.3f\n", mod_val))
if (!is.null(E(g)$weight)) {
  w_range <- range(E(g)$weight)
  cat(sprintf("  Edge weight range: [%.3f, %.3f]\n", w_range[1], w_range[2]))
}

mem        <- comm$membership
comm_sizes <- sort(table(mem), decreasing = TRUE)
cat("\nCommunity sizes:\n")
print(comm_sizes)

# ----------------------------------------------------------------------------
# 2. Map communities -> ecotype labels (by MC3 / MC8 dominance)
# ----------------------------------------------------------------------------
cat("\nMapping communities to ecotype labels...\n")

if (!all(V(g)$name == patient_meta$Study.Patient.ID)) {
  patient_meta <- patient_meta[match(V(g)$name,
                                     patient_meta$Study.Patient.ID), ]
}

community_ids <- as.integer(names(comm_sizes))
top_two       <- community_ids[1:2]

mc3_means <- tapply(patient_meta$MC3, mem, mean)
mc8_means <- tapply(patient_meta$MC8, mem, mean)

if (mc3_means[as.character(top_two[1])] >
    mc3_means[as.character(top_two[2])]) {
  mc3_dom_id <- top_two[1]; mc8_dom_id <- top_two[2]
} else {
  mc3_dom_id <- top_two[2]; mc8_dom_id <- top_two[1]
}

ecotype_labels <- rep("Outlier", vcount(g))
ecotype_labels[mem == mc3_dom_id] <- "Ecotype 1"
ecotype_labels[mem == mc8_dom_id] <- "Ecotype 2"

V(g)$ecotype <- factor(ecotype_labels,
                       levels = c("Ecotype 1", "Ecotype 2", "Outlier"))

n_e1  <- sum(V(g)$ecotype == "Ecotype 1")
n_e2  <- sum(V(g)$ecotype == "Ecotype 2")
n_out <- sum(V(g)$ecotype == "Outlier")
cat(sprintf("  Ecotype 1: n = %d  |  Ecotype 2: n = %d  |  Outlier: n = %d\n",
            n_e1, n_e2, n_out))

# ----------------------------------------------------------------------------
# 3. Layout — backbone on the main subgraph, outlier placed manually
# ----------------------------------------------------------------------------
cat("\nComputing layout...\n")

outlier_idx <- which(V(g)$ecotype == "Outlier")

if (length(outlier_idx) > 0) {
  g_main     <- delete_vertices(g, outlier_idx)
  main_names <- V(g_main)$name
} else {
  g_main     <- g
  main_names <- V(g)$name
}

set.seed(42)
layout_result <- tryCatch(
  graphlayouts::layout_as_backbone(g_main, keep = 0.4),
  error = function(e) {
    cat("  Backbone failed; falling back to stress (inverted weights)\n")
    inv_w <- 1 - E(g_main)$weight
    list(xy = graphlayouts::layout_with_stress(g_main, weights = inv_w))
  }
)

main_layout <- layout_result$xy
rownames(main_layout) <- main_names

full_layout <- matrix(NA_real_, nrow = vcount(g), ncol = 2)
rownames(full_layout) <- V(g)$name
full_layout[main_names, ] <- main_layout

if (length(outlier_idx) > 0) {
  x_range <- range(main_layout[, 1]); y_range <- range(main_layout[, 2])
  out_x <- x_range[1] - 0.10 * diff(x_range)
  out_y <- y_range[2] - 0.05 * diff(y_range)
  for (i in outlier_idx) full_layout[V(g)$name[i], ] <- c(out_x, out_y)
}

V(g)$x <- full_layout[, 1]
V(g)$y <- full_layout[, 2]

# ----------------------------------------------------------------------------
# 4. Styling
# ----------------------------------------------------------------------------
palette_ecotype <- c("Ecotype 1" = "#2C7BB6",
                     "Ecotype 2" = "#D7301F",
                     "Outlier"   = "#999999")

plot_subtitle <- bquote(
  italic(n) ~ "=" ~ .(vcount(g)) ~ "patients;" ~
    .(format(ecount(g), big.mark = ",")) ~ "edges (Pearson" ~
    italic(r) ~ "\u2265" ~ "0.40);" ~
    "Louvain modularity" ~ italic(Q) ~ "=" ~ .(sprintf("%.2f", mod_val))
)

ecotype_legend_labels <- c(
  "Ecotype 1" = sprintf("Ecotype 1 (MC3-dominant, n = %d)", n_e1),
  "Ecotype 2" = sprintf("Ecotype 2 (MC8-dominant, n = %d)", n_e2),
  "Outlier"   = sprintf("Outlier (n = %d)", n_out)
)

# ----------------------------------------------------------------------------
# 5. Plot
# ----------------------------------------------------------------------------
cat("Building plot...\n")

outlier_coords <- if (length(outlier_idx) > 0) {
  data.frame(x = V(g)$x[outlier_idx], y = V(g)$y[outlier_idx])
} else NULL

p <- ggraph(g, layout = "manual", x = V(g)$x, y = V(g)$y) +
  geom_edge_link0(aes(alpha = weight),
                  edge_width = 0.10, edge_colour = "grey45") +
  scale_edge_alpha_continuous(range = c(0.03, 0.15), guide = "none") +
  geom_mark_hull(
    aes(x = x, y = y, group = ecotype, fill = ecotype,
        filter = ecotype %in% c("Ecotype 1", "Ecotype 2")),
    concavity = 4, expand = unit(1.5, "mm"),
    radius = unit(2.5, "mm"), alpha = 0.10, color = NA
  ) +
  geom_node_point(aes(fill = ecotype), shape = 21,
                  size = 2.9, stroke = 0.4, color = "black") +
  scale_fill_manual(
    values = palette_ecotype, labels = ecotype_legend_labels, name = NULL,
    guide = guide_legend(override.aes = list(size = 3.2, alpha = 1), nrow = 1)
  ) +
  labs(title = NULL, subtitle = plot_subtitle) +
  theme_void(base_size = 7, base_family = "Helvetica") +
  theme(
    text            = element_text(family = "Helvetica", colour = "black"),
    plot.subtitle   = element_text(colour = "grey25", size = 7, hjust = 0.5,
                                   margin = margin(t = 2, b = 8)),
    legend.text     = element_text(colour = "black", size = 6.5),
    legend.position = "bottom",
    legend.key.size = unit(7, "pt"),
    plot.margin     = margin(8, 10, 6, 10, "pt")
  ) +
  coord_fixed(clip = "off")

if (!is.null(outlier_coords)) {
  p <- p + annotate("text",
                    x = outlier_coords$x, y = outlier_coords$y,
                    label = "outlier", hjust = -0.30, vjust = 0.5,
                    size = 1.9, colour = "grey35", fontface = "italic")
}

# ----------------------------------------------------------------------------
# 6. Save
# ----------------------------------------------------------------------------
dir.create("output", recursive = TRUE, showWarnings = FALSE)

ggsave("output/patient_network_example.pdf", p,
       width = 130, height = 125, units = "mm", dpi = 600, device = cairo_pdf)
ggsave("output/patient_network_example.png", p,
       width = 130, height = 125, units = "mm", dpi = 600)

cat("\nSaved:\n")
cat("  output/patient_network_example.pdf\n")
cat("  output/patient_network_example.png\n")
