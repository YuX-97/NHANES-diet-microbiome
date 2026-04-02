# NHANES-diet-microbiome
This repository contains code for the analysis of the manuscript based on NHANES data

## Project structure
NHANES-diet-microbiome/

├─ README.md

├─ code/ 

│  ├─ 00_dietary index calculation.R 

│  ├─ 01_dbRDA analysis.R 

│  ├─ 02_PERMANOVA.R 

│  ├─ 03_MaAsLin2 analysis.R 

│  ├─ 04_Mediation analysis.R 

│  ├─ 05_Unsupervised clustering.R 

│  ├─ 06_Differential abundance analysis.R 

## Data source
NHANES data are publicly available from the CDC NHANES website (https://wwwn.cdc.gov/nchs/nhanes/).

FEPD 37 food types information aligned with NHANES participants are publicly available at https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/.

## Notes
The code is provided to facilitate reproducibility of the study results. Analyses was done on R ver 4.4.1

## Required packageds may include:
readr, readxl, dplyr, openxlsx, parallel, doParallel, vegan, lavann, stringr, dietaryindex, foreach, tibble, ggplot2, grid, janitor, Maaslin2, FNN, cluster, mclust, fpc, ANCOMBC, phyloseq
