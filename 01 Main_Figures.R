## Main Figures 2, 3 and 5
## Run this script from R/RStudio. The figures are printed to the Plots pane.

rm(list = ls())

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_file <- gsub("~\\+~", " ", script_file)
script_dir <- getwd()
if (length(script_file) > 0) {
  script_dir <- dirname(normalizePath(script_file[1]))
} else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  active_file <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(active_file)) script_dir <- dirname(normalizePath(active_file))
}

rcode_dir <- file.path(dirname(script_dir), "20260624", "20260705 NAT manuscript LHH", "03 Rcode")
if (!file.exists(file.path(rcode_dir, "myplot.r"))) rcode_dir <- getwd()

source(file.path(rcode_dir, "myplot.r"))
path <- file.path(rcode_dir, "data")

## Figure 2: species-level LRR distribution and NES proportion ----
dat_0603 <- read.csv(file.path(path, "data_all.csv"), check.names = FALSE)
dat_0603 <- dat_0603[, nzchar(names(dat_0603)), drop = FALSE]

species_lrr <- dat_0603 %>%
  group_by(Species_IDs) %>%
  summarise(across(all_of(unname(ratio_cols)), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(all_of(unname(ratio_cols)), names_to = "variable", values_to = "LRR") %>%
  mutate(
    nutrient = factor(names(ratio_cols)[match(variable, ratio_cols)], levels = names(ratio_cols)),
    strategy = factor(if_else(LRR > 0, "NES", "NDS"), levels = c("NES", "NDS"))
  ) %>%
  filter(!is.na(LRR), LRR != 0)

fig2_balance <- species_lrr %>%
  group_by(nutrient) %>%
  summarise(
    NES = sum(strategy == "NES"),
    NDS = sum(strategy == "NDS"),
    n = NES + NDS,
    NES_proportion = NES / n,
    p = binom.test(NES, n, p = 0.5)$p.value,
    CI_low = binom.test(NES, n, p = 0.5)$conf.int[1],
    CI_high = binom.test(NES, n, p = 0.5)$conf.int[2],
    .groups = "drop"
  ) %>%
  mutate(p_label = sprintf("italic(p) == \"%.3f\"", p))

fig2a <- ggplot(species_lrr, aes(nutrient, LRR)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.28, color = "grey35") +
  geom_violin(fill = "grey94", color = "grey55", linewidth = 0.25, width = 0.82,
              trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white",
               color = "black", linewidth = 0.25) +
  geom_point(aes(color = strategy),
             position = position_jitter(width = 0.16, height = 0, seed = 123),
             shape = 1, size = 0.9, stroke = 0.35) +
  scale_color_manual(values = strategy_cols) +
  labs(tag = "a", x = "Log10(rh/nrh) nutrient", y = "Rhizosphere nutrient effect size") +
  theme_ms +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.03, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "horizontal"
  )

fig2b <- ggplot(fig2_balance, aes(y = nutrient, x = NES_proportion)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.28, color = "grey35") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0,
                 linewidth = 0.42, color = "black") +
  geom_point(size = 1.85, shape = 21, stroke = 0.45, fill = "black", color = "black") +
  geom_text(aes(x = 0.705, label = p_label), hjust = 0, size = 1.8,
            family = "Helvetica", parse = TRUE) +
  scale_x_continuous(limits = c(0.30, 0.86), breaks = seq(0.3, 0.7, 0.1),
                     labels = sprintf("%.1f", seq(0.3, 0.7, 0.1)), expand = c(0, 0)) +
  scale_y_discrete(limits = rev(names(ratio_cols))) +
  coord_cartesian(clip = "off") +
  labs(tag = "b", x = "Proportion as NES", y = "Log10(rh/nrh) nutrient") +
  theme_ms +
  theme(legend.position = "none", plot.margin = margin(3, 20, 3, 3, "mm"))

fig2 <- fig2a + fig2b + plot_layout(widths = c(1.38, 1))
print(fig2)

## Figure 3: fixation status and rhizosphere nutrient strategy ----
strategy_vars <- paste0("strategy_", ratio_cols)
names(strategy_vars) <- names(ratio_cols)

fig3_dat <- dat_0603 %>%
  group_by(Species_IDs, fixation_states) %>%
  summarise(across(all_of(unname(ratio_cols)), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(all_of(unname(ratio_cols)), names_to = "variable", values_to = "LRR") %>%
  mutate(
    nutrient = factor(names(ratio_cols)[match(variable, ratio_cols)], levels = names(ratio_cols)),
    fixation_states = factor(fixation_states, levels = c("Fixer", "Non-fixer")),
    NES = as.integer(LRR > 0),
    strategy = factor(if_else(NES == 1, "NES", "NDS"), levels = c("NES", "NDS"))
  ) %>%
  filter(!is.na(LRR), LRR != 0)

fig3_models <- bind_rows(lapply(names(ratio_cols), function(nutrient_i) {
  d <- filter(fig3_dat, nutrient == nutrient_i)
  m <- glm(NES ~ fixation_states, data = d, family = binomial)
  pred <- predict(m, newdata = data.frame(
    fixation_states = factor(c("Fixer", "Non-fixer"), levels = c("Fixer", "Non-fixer"))
  ), type = "response")
  n_fix <- sum(d$fixation_states == "Fixer")
  n_non <- sum(d$fixation_states == "Non-fixer")
  diff <- pred[2] - pred[1]
  se <- sqrt(pred[1] * (1 - pred[1]) / n_fix + pred[2] * (1 - pred[2]) / n_non)
  data.frame(
    nutrient = nutrient_i,
    diff = diff,
    low = diff - 1.96 * se,
    high = diff + 1.96 * se,
    p = coef(summary(m))["fixation_statesNon-fixer", "Pr(>|z|)"]
  )
})) %>%
  mutate(
    nutrient = factor(nutrient, levels = names(ratio_cols)),
    p_label = ifelse(p < 0.001, "italic(p) < 0.001", sprintf("italic(p) == \"%.3f\"", p))
  )

fig3a <- ggplot(fig3_dat, aes(fixation_states, LRR)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, color = "grey40") +
    geom_violin(aes(group = fixation_states), width = 0.86, fill = "grey95",
                color = "grey55", linewidth = 0.22, trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.22, fill = NA, color = "black",
                 linewidth = 0.25, outlier.shape = NA) +
    geom_point(aes(color = strategy), position = position_jitter(width = 0.14, seed = 20260603),
               shape = 1, size = 0.95, stroke = 0.22) +
    facet_wrap(~ nutrient, nrow = 2) +
    scale_color_manual(values = strategy_cols) +
    labs(tag = "a", x = "Nitrogen fixation status",
         y = "Rhizosphere nutrient effect size, log10(rh/nrh)") +
    theme_ms +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "inside",
      legend.position.inside = c(0.50, 0.98),
      legend.justification = c(0.5, 1),
      legend.direction = "horizontal",
      legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
    )
  
fig3b <- ggplot(fig3_models, aes(y = nutrient, x = diff)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.28, color = "grey40") +
    geom_errorbarh(aes(xmin = low, xmax = high), height = 0, linewidth = 0.45, color = "black") +
    geom_point(shape = 21, size = 2.0, stroke = 0.4, fill = "black", color = "black") +
    geom_text(aes(x = 0.50, label = p_label), hjust = 0, family = "Helvetica",
              size = 2.0, parse = TRUE) +
    scale_y_discrete(limits = rev(names(ratio_cols))) +
    scale_x_continuous(limits = c(-0.34, 0.86), breaks = c(-0.3, 0, 0.3),
                       labels = sprintf("%.1f", c(-0.3, 0, 0.3))) +
    labs(tag = "b", x = "Proportion as NES\nNon-fixer vs. Fixer",
         y = "Log10(rh/nrh) nutrient") +
    theme_ms +
    theme(plot.margin = margin(3, 5, 3, 3, "mm"))
  
fig3 <- fig3a + fig3b + plot_layout(widths = c(1.7, 1)) &
    theme(aspect.ratio = 2.5)
print(fig3)

## Figure 5: species-level associations between plant traits/bulk soil and strategy ----
font_family <- "Helvetica"

fig5_dat <- read.csv(file.path(path, "fig5_plot_data_species_level.csv"), check.names = FALSE) %>%
  mutate(
    nutrient = factor(nutrient, levels = names(ratio_cols)),
    row_label = factor(row_label, levels = c("NRH", "RTP", "RTN", "ROC",
                                             "SLA", "LDMC", "LA", "LTP", "LTN", "LOC")),
    row_group = factor(row_group, levels = c("Leaf traits", "Root traits", "Soil")),
    nominal_significant = as.logical(nominal_significant),
    p_label = ifelse(P_value < 0.05,
                     ifelse(P_value < 0.001, "***", ifelse(P_value < 0.01, "**", "*")),
                     "")
  )

beta_limit <- max(0.5, ceiling(max(abs(fig5_dat$beta_standardized), na.rm = TRUE) * 10) / 10)
beta_breaks <- seq(-beta_limit, beta_limit, length.out = 5)
r2_max <- max(fig5_dat$pseudo_R2_McFadden, na.rm = TRUE)

fig5_theme <- theme_classic(base_family = font_family, base_size = 6) +
  theme(
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 6, color = "black"),
    strip.background = element_blank(),
    legend.title = element_markdown(size = 6, color = "black", lineheight = 0.9),
    legend.text = element_text(size = 6, color = "black"),
    panel.spacing.y = unit(2.5, "mm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(3, 3, 3, 3, "mm")
  )

fig5_main <- ggplot(fig5_dat, aes(nutrient, row_label)) +
  geom_point(aes(fill = beta_standardized, size = pseudo_R2_McFadden),
             shape = 21, color = "grey55", stroke = 0.22, alpha = 0.95) +
  geom_point(data = subset(fig5_dat, nominal_significant),
             aes(fill = beta_standardized, size = pseudo_R2_McFadden),
             shape = 21, color = "black", stroke = 0.55) +
  geom_text(data = subset(fig5_dat, P_value < 0.05), aes(label = p_label),
            family = font_family, size = 2.2, color = "black", vjust = 0.45) +
  facet_grid(row_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#E69F00",
                       midpoint = 0, limits = c(-beta_limit, beta_limit), guide = "none") +
  scale_size_continuous(range = c(1.4, 5.2), limits = c(0, max(0.05, r2_max)), guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = "Nutrient strategy", y = NULL) +
  fig5_theme +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 6),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 7, face = "bold"),
    legend.position = "none",
    plot.margin = margin(3, 1, 3, 3, "mm")
  )

r2_scale_max <- max(0.05, r2_max)
r2_breaks <- pretty(c(0, r2_scale_max), n = 4)
r2_breaks <- r2_breaks[r2_breaks <= r2_scale_max]

bar_x <- c(0.08, 0.30)
bar_y <- c(0.59, 0.8)
legend_text_x <- bar_x[1]
legend_y <- function(value) {
  bar_y[1] + (value + beta_limit) / (2 * beta_limit) * diff(bar_y)
}

color_legend_data <- expand.grid(
  x = seq(bar_x[1], bar_x[2], length.out = 20),
  value = seq(-beta_limit, beta_limit, length.out = 300)
) %>%
  mutate(y = legend_y(value))

beta_label_data <- data.frame(
  value = beta_breaks,
  y = legend_y(beta_breaks),
  label = sprintf("%.1f", beta_breaks)
)

r2_legend_data <- data.frame(
  r2 = r2_breaks,
  y = seq(0.41, 0.23, length.out = length(r2_breaks))
)

fig5_legend <- ggplot() +
  geom_raster(data = color_legend_data, aes(x = x, y = y, fill = value), interpolate = TRUE) +
  geom_segment(
    data = subset(beta_label_data, value != min(value) & value != max(value)),
    aes(x = bar_x[1], xend = bar_x[2], y = y, yend = y),
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(data = beta_label_data, aes(x = 0.40, y = y, label = label),
            hjust = 0, family = font_family, size = 2.1) +
  annotate("text", x = legend_text_x, y = bar_y[2] + 0.023, label = "NES",
           hjust = 0, family = font_family, size = 2.1) +
  annotate("text", x = legend_text_x, y = bar_y[1] - 0.023, label = "NDS",
           hjust = 0, family = font_family, size = 2.1) +
  geom_richtext(aes(x = legend_text_x, y = 0.47),
                label = "McFadden<br>pseudo <i>R</i><sup>2</sup>",
                hjust = 0, vjust = 1, family = font_family, size = 2.1,
                lineheight = 0.9, fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  geom_point(data = r2_legend_data, aes(size = r2), x = 0.22, y = r2_legend_data$y,
             shape = 21, fill = NA, color = "black", stroke = 0.45) +
  geom_text(data = r2_legend_data,
            aes(x = 0.48, y = y, label = sprintf("%.2f", r2)),
            hjust = 0, family = font_family, size = 2.1) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#E69F00",
                       midpoint = 0, limits = c(-beta_limit, beta_limit), guide = "none") +
  scale_size_continuous(range = c(1.4, 5.2), limits = c(0, r2_scale_max), guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void(base_family = font_family, base_size = 6) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 1, 1, 0, "mm")
  )

fig5 <- fig5_main + fig5_legend +
  plot_layout(widths = c(1, 0.15)) &
  theme(plot.background = element_rect(fill = "white", color = NA))
print(fig5)
