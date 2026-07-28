# Rhizosphere Nutrient Enhancement and Depletion Strategies in Tropical Legumes

This repository contains the R code and processed data required to reproduce the statistical analyses and figures on rhizosphere nutrient enrichment and depletion strategies in tropical legumes. 

The scripts were tested using R 4.4.2. No non-standard hardware is required. Install the required R packages before running the scripts. Package installation typically requires less than 10 minutes on a standard desktop computer. Typical run time is approximately 5–10 minutes.

## System requirements
- macOS / Windows / Linux
- R (version 4.4.2)
- R packages: 
  Required R packages are described within the scripts and include, but are not limited to, ggplot2, dplyr, tidyr, FactoMineR, factoextra, ape, phytools and piecewiseSEM.

## Run
Run "01_Main_Figures.R" to reproduce the figures.

Run "02_PCA_PCoA_SEM.R" to reproduce the PCA, PCoA and SEM analyses.

**Note**

Species identifiers were standardized as anonymous labels to ensure consistent matching among individual-level data, species-level data, and the phylogenetic tree. 

The individual-level dataset uses sp01, sp02, ... as species identifiers, whereas the 67-species dataset and phylogenetic tree use Sp_001, Sp_002, .... 

These identifiers are used only for data matching and do not affect any numerical trait, soil, PCA, PCoA or SEM results.
