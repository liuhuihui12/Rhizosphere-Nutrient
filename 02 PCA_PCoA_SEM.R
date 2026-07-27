## PCA, phylogenetic PCoA and SEM
## Run this script from R/RStudio. Key tables are printed to the console.
## Please download packages: install.packages(c("FactoMineR", "factoextra", "ape", "piecewiseSEM"))

rm(list = ls())
source("myplot.r")
path <- "data/"

suppressPackageStartupMessages({
  library(FactoMineR)
  library(factoextra)
  library(ape)
  library(piecewiseSEM)
})

normalize_species <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\s+", "_", x)
  gsub("_+", "_", x)
}

dat <- read.csv(paste0(path, "67_species_data.csv"), check.names = FALSE)
dat$Species <- normalize_species(dat$genus_species)

if (nrow(dat) != 67) stop("Expected 67 species rows.")
if (any(duplicated(dat$Species))) stop("Duplicated species names.")

pca_sets <- list(
  Leaf = list(
    title = "Leaf traits PCA loading plot",
    cols = c("leaf_organic_carbon", "leaf_total_nitrogen", "leaf_total_phosphor",
             "leaf_area", "leaf_dry_matter_content", "specific_leaf_area"),
    labels = c("LOC", "LTN", "LTP", "LA", "LDMC", "SLA")
  ),
  Root = list(
    title = "Root traits PCA loading plot",
    cols = c("root_organic_carbon", "root_total_nitrogen", "root_total_phosphor"),
    labels = c("ROC", "RTN", "RTP")
  ),
  Soil = list(
    title = "Nrh soil nutrient PCA loading plot",
    cols = c("non_rhizosphere_soil_organic_carbon",
             "non_rhizosphere_soil_total_nitrogen",
             "log_non_rhizosphere_soil_total_phosphor",
             "non_rhizosphere_soil_available_phosphor",
             "non_rhizosphere_soil_nitrate_nitrogen",
             "non_rhizosphere_soil_ammonia_nitrogen"),
    labels = c("SOC", "STN", "STP", "SAP", "SNN", "SAN")
  ),
  Strategy = list(
    title = "Strategy PCA loading plot",
    cols = c("strategy_log_rh_soc_nrh_soc", "strategy_log_rh_stn_nrh_stn",
             "strategy_log_rh_stp_nrh_stp", "strategy_log_rh_snn_nrh_snn",
             "strategy_log_rh_san_nrh_san", "strategy_log_rh_sap_nrh_sap"),
    labels = c("strategy_soc", "strategy_stn", "strategy_stp",
               "strategy_snn", "strategy_san", "strategy_sap")
  )
)

run_pca <- function(set) {
  x <- dat[, set$cols, drop = FALSE]
  x[] <- lapply(x, as.numeric)
  x <- x[complete.cases(x), , drop = FALSE]
  names(x) <- set$labels
  PCA(x, scale.unit = TRUE, ncp = min(6, ncol(x)), graph = FALSE)
}

pca <- lapply(pca_sets, run_pca)

pca_eigen <- do.call(rbind, lapply(names(pca), function(nm) {
  eig <- get_eigenvalue(pca[[nm]])
  data.frame(
    Data_set = paste(nm, "PCA"),
    Axis = rownames(eig),
    Eigenvalue = eig[, "eigenvalue"],
    Proportion_of_variance_percent = eig[, "variance.percent"],
    Cumulative_variance_percent = eig[, "cumulative.variance.percent"],
    row.names = NULL
  )
}))

cat("\nPCA eigenvalues and variance explained:\n")
print(pca_eigen, row.names = FALSE)

pca_plots <- lapply(names(pca), function(nm) {
  fviz_pca_var(
    pca[[nm]],
    col.var = "contrib",
    gradient.cols = c("#FFA9D2", "#A0C9FF", "#76C235"),
    repel = TRUE,
    title = pca_sets[[nm]]$title
  ) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8))
})
pca_loading_figure <- wrap_plots(pca_plots, ncol = 2) +
  plot_annotation(tag_levels = "a")
print(pca_loading_figure)

scores <- dat[, c("X", "genus_species", "Species", "fixation_states", "n_individuals")]
for (nm in names(pca)) {
  coord <- as.data.frame(pca[[nm]]$ind$coord)
  names(coord) <- paste0(nm, "_PC", seq_len(ncol(coord)))
  scores <- cbind(scores, coord)
}

tree <- read.tree(paste0(path, "phylogeny.treefile"))
phylo_dist <- cophenetic.phylo(tree)
phylo_pcoa <- cmdscale(phylo_dist, k = 5, eig = TRUE)

phylo_eig <- phylo_pcoa$eig
pcoa_summary <- data.frame(
  PC = paste0("PC", seq_along(phylo_eig)),
  Eigenvalue = phylo_eig,
  Proportion_of_variance_percent = phylo_eig / sum(phylo_eig) * 100,
  Cumulative_variance_percent = cumsum(phylo_eig / sum(phylo_eig) * 100)
)

cat("\nPhylogenetic PCoA eigenvalues and variance explained:\n")
print(head(pcoa_summary, 5), row.names = FALSE)

pcoa_scores <- as.data.frame(phylo_pcoa$points)
names(pcoa_scores) <- paste0("PCoA", seq_len(ncol(pcoa_scores)))
pcoa_scores$Species <- normalize_species(rownames(pcoa_scores))

if (!setequal(scores$Species, pcoa_scores$Species)) {
  stop("Species names do not match between PCA and PCoA outputs.")
}

sem_data <- merge(
  scores[, c("X", "genus_species", "Species", "fixation_states", "n_individuals",
             "Leaf_PC1", "Leaf_PC2", "Root_PC1", "Soil_PC1", "Soil_PC2", "Strategy_PC1")],
  pcoa_scores[, c("Species", "PCoA1")],
  by = "Species",
  sort = FALSE
)
sem_data <- sem_data[match(scores$Species, sem_data$Species), ]
names(sem_data)[names(sem_data) == "PCoA1"] <- "Phylo_PC1"
sem_data$Fixation_binary <- ifelse(sem_data$fixation_states == "Fixer", 1L, 0L)
sem_data$Strategy_binary <- ifelse(sem_data$Strategy_PC1 > 0, 1L, 0L)

model_equations <- list(
  M1 = c(
    "Fixation_binary ~ Phylo_PC1",
    "Root_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1",
    "Leaf_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1 + Root_PC1",
    "Strategy_binary ~ Soil_PC1 + Fixation_binary + Leaf_PC1 + Root_PC1 + Phylo_PC1"
  ),
  M1a = c(
    "Fixation_binary ~ Phylo_PC1",
    "Leaf_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1",
    "Root_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1",
    "Leaf_PC1 %~~% Root_PC1",
    "Strategy_binary ~ Soil_PC1 + Fixation_binary + Leaf_PC1 + Root_PC1 + Phylo_PC1"
  ),
  M2_root_to_leaf = c(
    "Fixation_binary ~ Phylo_PC1",
    "Root_PC1 ~ Soil_PC1 + Soil_PC2 + Fixation_binary + Phylo_PC1",
    "Leaf_PC1 ~ Soil_PC1 + Soil_PC2 + Fixation_binary + Phylo_PC1 + Root_PC1",
    "Leaf_PC2 ~ Soil_PC1 + Soil_PC2 + Fixation_binary + Phylo_PC1 + Root_PC1",
    "Strategy_binary ~ Soil_PC1 + Soil_PC2 + Leaf_PC1 + Leaf_PC2 + Root_PC1 + Fixation_binary + Phylo_PC1"
  ),
  M3_root_to_leaf = c(
    "Fixation_binary ~ Phylo_PC1",
    "Root_PC1 ~ Soil_PC1 + Soil_PC2 + Fixation_binary + Phylo_PC1",
    "Leaf_PC1 ~ Soil_PC1 + Soil_PC2 + Fixation_binary + Phylo_PC1 + Root_PC1",
    "Strategy_binary ~ Soil_PC1 + Soil_PC2 + Leaf_PC1 + Root_PC1 + Fixation_binary + Phylo_PC1"
  ),
  M4_root_to_leaf = c(
    "Fixation_binary ~ Phylo_PC1",
    "Root_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1",
    "Leaf_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1 + Root_PC1",
    "Leaf_PC2 ~ Soil_PC1 + Fixation_binary + Phylo_PC1 + Root_PC1",
    "Strategy_binary ~ Soil_PC1 + Leaf_PC1 + Leaf_PC2 + Root_PC1 + Fixation_binary + Phylo_PC1"
  )
)

sem_models <- list(
  M1 = psem(
    glm(Fixation_binary ~ Phylo_PC1, data = sem_data, family = binomial),
    lm(Root_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1, data = sem_data),
    lm(Leaf_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1 + Root_PC1, data = sem_data),
    glm(Strategy_binary ~ Soil_PC1 + Fixation_binary + Leaf_PC1 + Root_PC1 + Phylo_PC1,
        data = sem_data, family = binomial),
    data = sem_data
  ),
  M1a = psem(
    glm(Fixation_binary ~ Phylo_PC1, data = sem_data, family = binomial),
    lm(Leaf_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1, data = sem_data),
    lm(Root_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1, data = sem_data),
    glm(Strategy_binary ~ Soil_PC1 + Fixation_binary + Leaf_PC1 + Root_PC1 + Phylo_PC1,
        data = sem_data, family = binomial),
    Leaf_PC1 %~~% Root_PC1,
    data = sem_data
  ),
  M2_root_to_leaf = psem(
    glm(Fixation_binary ~ Phylo_PC1, data = sem_data, family = binomial),
    lm(Root_PC1 ~ Soil_PC1 + Soil_PC2 + Fixation_binary + Phylo_PC1, data = sem_data),
    lm(Leaf_PC1 ~ Soil_PC1 + Soil_PC2 + Fixation_binary + Phylo_PC1 + Root_PC1, data = sem_data),
    lm(Leaf_PC2 ~ Soil_PC1 + Soil_PC2 + Fixation_binary + Phylo_PC1 + Root_PC1, data = sem_data),
    glm(Strategy_binary ~ Soil_PC1 + Soil_PC2 + Leaf_PC1 + Leaf_PC2 + Root_PC1 +
          Fixation_binary + Phylo_PC1, data = sem_data, family = binomial),
    data = sem_data
  ),
  M3_root_to_leaf = psem(
    glm(Fixation_binary ~ Phylo_PC1, data = sem_data, family = binomial),
    lm(Root_PC1 ~ Soil_PC1 + Soil_PC2 + Fixation_binary + Phylo_PC1, data = sem_data),
    lm(Leaf_PC1 ~ Soil_PC1 + Soil_PC2 + Fixation_binary + Phylo_PC1 + Root_PC1, data = sem_data),
    glm(Strategy_binary ~ Soil_PC1 + Soil_PC2 + Leaf_PC1 + Root_PC1 +
          Fixation_binary + Phylo_PC1, data = sem_data, family = binomial),
    data = sem_data
  ),
  M4_root_to_leaf = psem(
    glm(Fixation_binary ~ Phylo_PC1, data = sem_data, family = binomial),
    lm(Root_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1, data = sem_data),
    lm(Leaf_PC1 ~ Soil_PC1 + Fixation_binary + Phylo_PC1 + Root_PC1, data = sem_data),
    lm(Leaf_PC2 ~ Soil_PC1 + Fixation_binary + Phylo_PC1 + Root_PC1, data = sem_data),
    glm(Strategy_binary ~ Soil_PC1 + Leaf_PC1 + Leaf_PC2 + Root_PC1 +
          Fixation_binary + Phylo_PC1, data = sem_data, family = binomial),
    data = sem_data
  )
)

sem_summaries <- lapply(sem_models, summary)
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
fmt_test <- function(value, df, p) sprintf("%.3f (%s, p = %s)", value, df, fmt_p(p))

sem_fit <- do.call(rbind, lapply(names(sem_summaries), function(nm) {
  s <- sem_summaries[[nm]]
  data.frame(
    Model = nm,
    Structural_equations = paste(model_equations[[nm]], collapse = "; "),
    Chi_square_df_p = fmt_test(s$ChiSq$Chisq, s$ChiSq$df, s$ChiSq$P.Value),
    Fisher_C_df_p = fmt_test(s$Cstat$Fisher.C, s$Cstat$df, s$Cstat$P.Value),
    AIC = s$AIC$AIC,
    row.names = NULL
  )
}))

cat("\nSEM model fit statistics:\n")
print(sem_fit[order(sem_fit$AIC), ], row.names = FALSE)

cat("\nFinal M1 model summary:\n")
print(sem_summaries$M1)
