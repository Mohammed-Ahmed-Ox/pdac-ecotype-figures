#!/usr/bin/env Rscript
# ============================================================================
# generate_example_data.R
# ----------------------------------------------------------------------------
# Builds a small SYNTHETIC patient-similarity dataset that mirrors the
# structure of the real network object used by the figure scripts, so the
# repository can be cloned and run end-to-end without any patient data.
#
# The simulation produces 95 patients scored across 12 metabolic clusters
# (MC1-MC12) and organised into two ecotypes plus a single outlier:
#
#   Ecotype 1 (n = 48)  MC3-dominant
#   Ecotype 2 (n = 46)  MC8-dominant
#   Outlier   (n =  1)  idiosyncratic profile, weakly correlated to the rest
#
# Variation is layered on top of a PRONOUNCED shared baseline. Several
# metabolic clusters (notably MC12, MC1, MC5, MC9) are abundant in every
# patient; this strong shared structure is what makes patients resemble one
# another enough to stay connected ACROSS ecotype lines at r >= 0.40, giving
# a dense, organic network rather than two isolated cliques. The ecotype axis
# (MC3 vs MC8) then carves the two soft communities on top of it, and a small
# shared background axis adds within-ecotype heterogeneity.
#
# Patients are linked into a similarity network (Pearson r >= 0.40 between
# their composition vectors) and Louvain community detection recovers the
# ecotypes. Detected communities are renamed to a fixed scheme (2 = MC3-
# dominant, 3 = MC8-dominant, 1 = anything else) so the figure scripts can
# rely on consistent labels.
#
# Writes:  data/example_input.rds
# Run from the repository root:  Rscript data/generate_example_data.R
# ============================================================================

suppressPackageStartupMessages(library(igraph))

set.seed(42)

# ----------------------------------------------------------------------------
# Parameters worth tuning
# ----------------------------------------------------------------------------
# ecotype_strength  how far MC3 / MC8 pull the two ecotypes apart; higher =
#                   stronger MC3/MC8 anti-correlation and sharper communities
#                   (higher modularity, lower density)
# bg_strength       within-ecotype heterogeneity along the shared background
#                   axis; modest values add realistic spread, but pushing it
#                   too high begins to fragment the communities
n_e1             <- 48
n_e2             <- 46
n_out            <- 1
r_threshold      <- 0.40
ecotype_strength <- 12
bg_strength      <- 1.0
noise_sd         <- 3

n_patients <- n_e1 + n_e2 + n_out
mc_names   <- sprintf("MC%d", 1:12)

# ----------------------------------------------------------------------------
# 1. Shared baseline composition
# ----------------------------------------------------------------------------
# A pronounced shared "shape" across the 12 MCs. The high, shared clusters
# (MC12, MC1, MC5, MC9) dominate every patient's profile, so any two patients
# correlate strongly on that common structure before the ecotype signal pulls
# MC3 / MC8 apart. This is what keeps the network dense and connected.
baseline <- c(MC1 = 18, MC2 = 6, MC3 = 8,  MC4 = 8,  MC5 = 14, MC6 = 3,
              MC7 = 3,  MC8 = 8, MC9 = 12, MC10 = 6, MC11 = 6, MC12 = 28)

# Loadings for the shared background axis (deliberately NOT on MC3 / MC8, so
# this axis is orthogonal to the ecotype signal).
bg_load <- c(MC1 = 4, MC5 = 5, MC9 = 4, MC10 = 4, MC11 = 4, MC12 = 5) *
           bg_strength

# ----------------------------------------------------------------------------
# 2. Simulate per-patient compositions
# ----------------------------------------------------------------------------
ecotype <- c(rep("E1", n_e1), rep("E2", n_e2), rep("Outlier", n_out))
latent  <- rnorm(n_patients, 0, 1)   # each patient's background-axis score

mc_matrix <- matrix(
  0, nrow = n_patients, ncol = 12,
  dimnames = list(sprintf("P%03d", seq_len(n_patients)), mc_names)
)

for (i in seq_len(n_patients)) {
  comp <- baseline + rnorm(12, 0, noise_sd)

  if (ecotype[i] == "E1") {
    comp["MC3"] <- comp["MC3"] + ecotype_strength       + rnorm(1, 0, 5)
    comp["MC8"] <- comp["MC8"] - ecotype_strength * 0.6 + rnorm(1, 0, 3)
    comp[names(bg_load)] <- comp[names(bg_load)] + latent[i] * bg_load
  } else if (ecotype[i] == "E2") {
    comp["MC8"] <- comp["MC8"] + ecotype_strength * 1.4 + rnorm(1, 0, 6)
    comp["MC3"] <- comp["MC3"] - ecotype_strength * 0.5 + rnorm(1, 0, 3)
    comp[names(bg_load)] <- comp[names(bg_load)] + latent[i] * bg_load
  } else {
    # Outlier: dominated by MCs that are neither high in the baseline nor on
    # the background axis, so it stays weakly correlated with everyone.
    comp[]                               <- 2
    comp[c("MC2", "MC4", "MC6", "MC7")]  <- c(24, 22, 30, 28)
  }

  comp[comp < 0] <- 0
  mc_matrix[i, ] <- comp / sum(comp) * 100   # % of each patient's cells
}

# ----------------------------------------------------------------------------
# 3. Patient similarity network
# ----------------------------------------------------------------------------
correlation_matrix <- cor(t(mc_matrix))      # 95 x 95, on MC profiles

adj <- correlation_matrix
adj[adj < r_threshold] <- 0
diag(adj) <- 0

g <- graph_from_adjacency_matrix(adj, mode = "undirected",
                                 weighted = TRUE, diag = FALSE)

# ----------------------------------------------------------------------------
# 4. Community detection
# ----------------------------------------------------------------------------
comm    <- cluster_louvain(g, weights = E(g)$weight)
mod_val <- modularity(comm)
raw_mem <- membership(comm)

# ----------------------------------------------------------------------------
# 5. Rename communities to a fixed scheme
#    The two largest communities are the ecotypes; the one with higher mean
#    MC3 becomes 2, the other becomes 3. Everything else (the singleton) is 1.
# ----------------------------------------------------------------------------
sizes       <- sort(table(raw_mem), decreasing = TRUE)
two_biggest <- as.integer(names(sizes))[1:2]
mc3_means   <- tapply(mc_matrix[names(raw_mem), "MC3"], raw_mem, mean)

if (mc3_means[as.character(two_biggest[1])] >=
    mc3_means[as.character(two_biggest[2])]) {
  mc3_dom <- two_biggest[1]; mc8_dom <- two_biggest[2]
} else {
  mc3_dom <- two_biggest[2]; mc8_dom <- two_biggest[1]
}

label_map <- setNames(rep(1L, length(unique(raw_mem))),
                      as.character(sort(unique(raw_mem))))
label_map[as.character(mc3_dom)] <- 2L
label_map[as.character(mc8_dom)] <- 3L

comm$membership <- as.integer(label_map[as.character(raw_mem)])

# ----------------------------------------------------------------------------
# 6. Patient metadata table
# ----------------------------------------------------------------------------
studies <- sprintf("Study_%s", LETTERS[1:9])

patient_communities <- data.frame(
  Study.Patient.ID = V(g)$name,
  Community        = comm$membership,
  Study            = sample(studies, n_patients, replace = TRUE),
  Treatment        = sample(c("Naive", "Neo-CTx"), n_patients,
                            replace = TRUE, prob = c(0.84, 0.16)),
  n_cells          = sample(500:8000, n_patients, replace = TRUE),
  stringsAsFactors = FALSE
)
patient_communities <- cbind(patient_communities,
                             as.data.frame(mc_matrix[V(g)$name, ]))

# ----------------------------------------------------------------------------
# 7. Save
# ----------------------------------------------------------------------------
typeA <- list(
  graph               = g,
  communities         = comm,
  patient_communities = patient_communities,
  mc_matrix           = mc_matrix,
  correlation_matrix  = correlation_matrix,
  threshold           = r_threshold,
  modularity          = mod_val,
  timestamp           = Sys.time()
)

dir.create("data", showWarnings = FALSE)
saveRDS(typeA, "data/example_input.rds")

# ----------------------------------------------------------------------------
# 8. Report
# ----------------------------------------------------------------------------
rho <- cor(mc_matrix[, "MC3"], mc_matrix[, "MC8"], method = "spearman")

cat("Wrote data/example_input.rds\n")
cat(sprintf("  Patients:    %d\n", vcount(g)))
cat(sprintf("  Edges:       %d  (density %.2f)\n", ecount(g), edge_density(g)))
cat(sprintf("  Communities: %d  (sizes: %s)\n",
            length(unique(comm$membership)),
            paste(sort(table(comm$membership), decreasing = TRUE),
                  collapse = " / ")))
cat(sprintf("  Modularity:  %.3f\n", mod_val))
cat(sprintf("  MC3-MC8 Spearman rho: %.2f\n", rho))
