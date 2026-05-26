#!/usr/bin/env Rscript
# ============================================================================
# fig_mc3_mc8_scatter.R
# ----------------------------------------------------------------------------
# Patient-level scatter of MC3 vs MC8 abundance, coloured by metabolic
# ecotype, with median quadrant lines and association statistics (Spearman
# correlation and a Fisher test on the median-split quadrants).
#
# Ecotypes are assigned programmatically from the community labels, by MC3 /
# MC8 dominance, so the script does not depend on Louvain's arbitrary numbering.
#
# Reads:   data/example_input.rds
# Writes:  output/mc3_mc8_scatter_example.pdf / .png
# Run from the repository root.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ----------------------------------------------------------------------------
# 1. Load
# ----------------------------------------------------------------------------
cat("Loading patient composition data...\n")
typeA <- readRDS("data/example_input.rds")

if (is.null(typeA$patient_communities) ||
    !all(c("MC3", "MC8", "Community") %in%
         colnames(typeA$patient_communities))) {
  stop("Expected patient_communities with MC3, MC8 and Community columns.")
}
patient_df <- as.data.frame(typeA$patient_communities)
cat(sprintf("  n = %d patients\n", nrow(patient_df)))

# ----------------------------------------------------------------------------
# 2. Assign ecotypes from community labels (by MC3 / MC8 dominance)
# ----------------------------------------------------------------------------
comm_sizes  <- sort(table(patient_df$Community), decreasing = TRUE)
main_two    <- as.integer(names(comm_sizes))[1:2]
mc3_by_comm <- tapply(patient_df$MC3, patient_df$Community, mean)

if (mc3_by_comm[as.character(main_two[1])] >=
    mc3_by_comm[as.character(main_two[2])]) {
  e1_comm <- main_two[1]; e2_comm <- main_two[2]
} else {
  e1_comm <- main_two[2]; e2_comm <- main_two[1]
}

patient_df <- patient_df %>%
  mutate(
    Ecotype = case_when(
      Community == e1_comm ~ "Ecotype 1",   # MC3-dominant
      Community == e2_comm ~ "Ecotype 2",   # MC8-dominant
      TRUE                 ~ "Outlier"
    ),
    Ecotype = factor(Ecotype, levels = c("Ecotype 1", "Ecotype 2", "Outlier"))
  )

cat("  Ecotype distribution:\n"); print(table(patient_df$Ecotype))

# ----------------------------------------------------------------------------
# 3. Statistics
# ----------------------------------------------------------------------------
rho_test <- suppressWarnings(
  cor.test(patient_df$MC3, patient_df$MC8, method = "spearman")
)
rho_val <- as.numeric(rho_test$estimate)

mc3_med <- median(patient_df$MC3)
mc8_med <- median(patient_df$MC8)

patient_df <- patient_df %>%
  mutate(MC3_high = MC3 >= mc3_med,
         MC8_high = MC8 >= mc8_med)

q_n <- c(
  "Both high"          = sum(patient_df$MC3_high  & patient_df$MC8_high),
  "Both low"           = sum(!patient_df$MC3_high & !patient_df$MC8_high),
  "MC3 high / MC8 low" = sum(patient_df$MC3_high  & !patient_df$MC8_high),
  "MC3 low / MC8 high" = sum(!patient_df$MC3_high & patient_df$MC8_high)
)

fisher_mat <- matrix(c(
  q_n["Both low"], q_n["MC3 low / MC8 high"],
  q_n["MC3 high / MC8 low"], q_n["Both high"]
), nrow = 2, byrow = TRUE)
fisher_res <- fisher.test(fisher_mat)
or_val <- as.numeric(fisher_res$estimate)

cat(sprintf("\n  Spearman rho = %.3f\n", rho_val))
cat(sprintf("  Fisher OR    = %.3f\n", or_val))

# ----------------------------------------------------------------------------
# 4. Styling
# ----------------------------------------------------------------------------
palette_ecotype <- c("Ecotype 1" = "#2C7BB6",
                     "Ecotype 2" = "#D7301F",
                     "Outlier"   = "#999999")

theme_nature <- function() {
  theme_classic(base_size = 7, base_family = "Helvetica") +
    theme(
      text          = element_text(family = "Helvetica", colour = "black"),
      axis.text     = element_text(colour = "black", size = 7),
      axis.title    = element_text(colour = "black", size = 8),
      plot.subtitle = element_text(colour = "grey30", size = 7, hjust = 0,
                                   margin = margin(b = 6)),
      axis.line     = element_line(colour = "black", linewidth = 0.4),
      axis.ticks    = element_line(colour = "black", linewidth = 0.4),
      panel.grid    = element_blank(),
      legend.key    = element_blank(),
      legend.text   = element_text(colour = "black", size = 6.5),
      plot.margin   = margin(8, 12, 8, 8, "pt")
    )
}

# ----------------------------------------------------------------------------
# 5. Plot
# ----------------------------------------------------------------------------
x_max <- max(patient_df$MC3) * 1.04
y_max <- max(patient_df$MC8) * 1.04
x_pad <- 0.8
y_pad <- 1.2

plot_subtitle <- bquote(
  "Spearman" ~ rho ~ "=" ~ .(sprintf("%.2f", rho_val)) *
    ", " ~ italic(P) ~ "<" ~ "0.001" ~ "    " ~
    "Fisher OR" ~ "=" ~ .(sprintf("%.2f", or_val)) *
    ", " ~ italic(P) ~ "<" ~ "0.001"
)

p_main <- ggplot(patient_df, aes(x = MC3, y = MC8)) +
  geom_hline(yintercept = mc8_med, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = mc3_med, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  annotate("text", x = x_max, y = y_max,
           label = sprintf("Both high\nn = %d", q_n["Both high"]),
           hjust = 1, vjust = 1, size = 2.0, lineheight = 0.95,
           colour = "grey25", fontface = "italic") +
  annotate("text", x = mc3_med - x_pad, y = y_max,
           label = sprintf("MC3 low / MC8 high\nn = %d",
                           q_n["MC3 low / MC8 high"]),
           hjust = 1, vjust = 1, size = 2.0, lineheight = 0.95,
           colour = "grey25", fontface = "italic") +
  annotate("text", x = x_max, y = mc8_med - y_pad,
           label = sprintf("MC3 high / MC8 low\nn = %d",
                           q_n["MC3 high / MC8 low"]),
           hjust = 1, vjust = 1, size = 2.0, lineheight = 0.95,
           colour = "grey25", fontface = "italic") +
  annotate("text", x = 0, y = mc8_med - y_pad,
           label = sprintf("Both low\nn = %d", q_n["Both low"]),
           hjust = 0, vjust = 1, size = 2.0, lineheight = 0.95,
           colour = "grey25", fontface = "italic") +
  geom_point(aes(fill = Ecotype), shape = 21, size = 2.4,
             stroke = 0.3, colour = "black", alpha = 0.88) +
  scale_fill_manual(values = palette_ecotype, name = NULL,
                    guide = guide_legend(override.aes = list(size = 3))) +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.01, 0.04))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.01, 0.04))) +
  labs(title = NULL, subtitle = plot_subtitle,
       x = "MC3 abundance (% of patient cells)",
       y = "MC8 abundance (% of patient cells)") +
  theme_nature() +
  theme(
    legend.position   = c(0.86, 0.80),
    legend.background = element_rect(fill = "white", colour = "grey85",
                                     linewidth = 0.3),
    legend.margin     = margin(2, 4, 2, 4),
    legend.key.size   = unit(7, "pt")
  )

# ----------------------------------------------------------------------------
# 6. Save
# ----------------------------------------------------------------------------
dir.create("output", recursive = TRUE, showWarnings = FALSE)

ggsave("output/mc3_mc8_scatter_example.pdf", p_main,
       width = 183, height = 95, units = "mm", dpi = 600, device = cairo_pdf)
ggsave("output/mc3_mc8_scatter_example.png", p_main,
       width = 183, height = 95, units = "mm", dpi = 600)

cat("\nSaved:\n")
cat("  output/mc3_mc8_scatter_example.pdf\n")
cat("  output/mc3_mc8_scatter_example.png\n")
