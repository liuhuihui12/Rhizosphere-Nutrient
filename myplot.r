## Shared settings for reviewer reproducibility scripts.
## Please download packages：install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "ggtext"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})

theme_ms <- theme_classic(base_family = "Helvetica", base_size = 6) +
  theme(
    axis.title = element_text(size = 7, color = "black"),
    axis.text = element_text(size = 6, color = "black"),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    plot.tag = element_text(size = 7, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 7, face = "bold")
  )

strategy_cols <- c(NES = "#E69F00", NDS = "#0072B2")

ratio_cols <- c(
  SOC = "log_rh_soc_nrh_soc",
  STN = "log_rh_stn_nrh_stn",
  STP = "log_rh_stp_nrh_stp",
  SNN = "log_rh_snn_nrh_snn",
  SAN = "log_rh_san_nrh_san",
  SAP = "log_rh_sap_nrh_sap"
)

