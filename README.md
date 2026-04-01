# NHANES-diet-microbiome
This repository contains code for the analysis of the manuscript based on NHANES data.

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

## Notes
The code is provided to facilitate reproducibility of the study results. Analyses was done on R ver 4.4.1

## Required packageds may include:
readr 

dplyr 

openxlsx 

parallel 

vegan 

lavann 

stringr 

dietaryindex 

foreach 

purrr 

tibble 

cluster 

mclust 

ANCOMBC 
