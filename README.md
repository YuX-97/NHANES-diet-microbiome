# NHANES Diet–Oral Microbiome Analysis

This repository provides the key R scripts for the downstream statistical analyses conducted for the manuscript **“Oral health status-contextualized associations between dietary intake and oral-rinse microbiome in US adults.”** The study used data from adults participating in the 2011–2014 National Health and Nutrition Examination Survey (NHANES).

The shared scripts focus on the major statistical procedures used to investigate associations among dietary intake, oral health conditions, and the oral-rinse microbiome. They have been curated for code availability: local file paths, dataset-specific column positions, and visualization routines have been removed, while the principal analytical steps and tabular outputs have been retained.

## Repository structure

```text
NHANES-diet-microbiome/
├── README.md
└── code/
    ├── ANCOM-BC2_adjustment.R
    ├── CLR transformation.R
    ├── Dietary intake Liner regression_adjustment.R
    ├── Dietary intake cluster_variable standardized differences.R
    ├── Dietary intake_PCA_kmeans clustering.R
    ├── MaAsLin2 analysis_adjustment.R
    ├── dbRDA analysis.R
    └── multivariable PERMANOVA analysis.R
```

## Description of analysis scripts

| Script | Main analysis |
| --- | --- |
| `CLR transformation.R` | Replaces zero values in the filtered genus-level read-count matrix with a fixed pseudocount (0.005), applies total-sum scaling, and performs centered log-ratio (CLR) transformation. Exports a tab-separated file. |
| `multivariable PERMANOVA analysis.R` | Evaluates the marginal contributions of dietary intake and other covariates to oral-rinse microbiome variation using Aitchison distance and multivariable PERMANOVA (`vegan::adonis2`, `by = "margin"`; 1,001 permutations). |
| `dbRDA analysis.R` | Performs distance-based redundancy analysis (dbRDA) of dietary predictors and the CLR-transformed microbiome in the overall population and separately within four oral-condition groups. Includes sequential and marginal permutation tests, variable-specific explained variance, and ordination scores. |
| `ANCOM-BC2_adjustment.R` | Examines genus-level differential abundance across oral conditions using ANCOM-BC2, with both unadjusted and covariate-adjusted models. |
| `Dietary intake Liner regression_adjustment.R` | Assesses differences in individual dietary variables across oral-condition groups using unadjusted and covariate-adjusted linear regression models. |
| `MaAsLin2 analysis_adjustment.R` | Tests associations between individual dietary variables and CLR-transformed microbial genera using MaAsLin2. Model 1 includes demographic, socioeconomic, smoking, and tooth-number covariates; Model 2 additionally adjusts for oral condition. |
| `Dietary intake_PCA_kmeans clustering.R` | Preprocesses dietary variables, performs principal component analysis (PCA), assesses cluster tendency and clustering solutions for `k = 2–5`, and conducts the final k-means clustering with `k = 3` using PC1–PC16. |
| `Dietary intake cluster_variable standardized differences.R` | Calculates cluster-specific mean dietary intakes and their standardized deviations from overall means. Exports the numerical data used for the dietary-pattern heatmap. |

## Data sources

- **NHANES:** Publicly available survey, examination, dietary, and laboratory datasets for the 2011–2012 and 2013–2014 cycles can be accessed through the [CDC NHANES website](https://wwwn.cdc.gov/nchs/nhanes/).
- **Food Patterns Equivalents Database (FPED):** Food-pattern component intake data are available through the [USDA Agricultural Research Service FPED databases](https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/).

The combined analysis-ready datasets used by the scripts are not included in this repository. Data assembly, variable harmonization, and upstream NHANES processing should be distinguished from the downstream analyses documented here.

## Software requirements

The analyses were conducted using **R version 4.4.1**. Depending on the script, required R packages may include:

`vegan`, `ANCOMBC`, `phyloseq`, `Maaslin2`, `readr`, `readxl`, `dplyr`, `janitor`, `openxlsx`, `FNN`, `cluster`, `mclust`, and `fpc`.

Base R functions and the `parallel` package are also used for data processing and, where applicable, parallel computation. Specific package requirements can be found at the beginning of each script.

## Usage and reproducibility notes

1. Download the relevant NHANES and FPED data and prepare the analysis-ready datasets appropriate to each statistical procedure.
2. Set the input/output filenames and identify the corresponding variables by **column name** as indicated in each script. Some dataset-specific variable lists are intentionally left as placeholders.
3. Run the relevant preprocessing step(s), followed by the statistical analyses of interest. The scripts are not intended to be executed as an entirely automated end-to-end pipeline.
4. Consult the manuscript Methods and supplementary information for the analytical population, exclusions, and interpretation of the reported findings.

**Note:** These files are curated representations of the key downstream analytical procedures, intended to support methodological transparency. Reproducing numerical results requires the corresponding analysis-ready data, variable definitions, sample exclusions, and software environment.
