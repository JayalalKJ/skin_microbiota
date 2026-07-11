# Skin Microbiota Analysis Pipeline

This repository contains R scripts used for bioinformatics and microbial ecology analysis of Atlantic salmon skin microbiota.

## Pipeline scripts

| Script                                            | Purpose                                                                  |
| ------------------------------------------------- | ------------------------------------------------------------------------ |
| `00_project_setup.R`                              | Loads packages, defines folders, and prepares the analysis environment.  |
| `01_import_mothur_to_microeco.R`                  | Imports mothur output and creates `microeco`/`phyloseq` objects.         |
| `02A_Blanks_visual_proofpack.R`                   | Visualises extraction blanks and negative controls.                      |
| `02B_Negative_control_samples_Composition.R`      | Examines the taxonomic composition of negative controls.                 |
| `02C_contaminant_removal_from_blanks.R`           | Identifies and removes potential contaminants.                           |
| `03_rarefaction_Skin_sample_for_alpha.R`          | Generates rarefaction curves and applies the sequencing-depth threshold. |
| `04_alpha_diversity_Skin_samples.R`               | Calculates and compares alpha-diversity metrics.                         |
| `05_Beta_diversity_Skin_samples.R`                | Performs beta-diversity analysis of skin microbiota.                     |
| `06_Beta_diversity_analysis_Skin_Water_samples.R` | Compares skin and water microbial communities.                           |
| `07_composition_analysis_Skin_water_samples.R`    | Analyses and visualises taxonomic composition.                           |
| `Run_all_pipeline.R`                              | Runs all scripts in the correct order.                                   |

## Run the pipeline

Open R or RStudio in the repository folder and run:

```r
source("Run_all_pipeline.R")
```

The scripts should be executed in numerical order because later steps depend on files and objects created in earlier steps.

## Main R packages

```r
library(phyloseq)
library(microeco)
library(file2meco)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ape)
```
