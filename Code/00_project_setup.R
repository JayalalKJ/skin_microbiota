#!/usr/bin/env Rscript
# =============================================================================
# 00_project_setup.R
# Jayalal K Jayanthan
# PhD candidate (2021– 2026)
# Research group: Seafood Science
# Institute: The Norwegian College of Fishery Science
# Faculty: Faculty of Biosciences, Fisheries, and Economics
# Campus: UiT Campus Tromsø
# UiT The Arctic University of Norway
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 1000)

dir.create("scripts", recursive = TRUE, showWarnings = FALSE)
dir.create("data",    recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

out_dir <- file.path("results", "Step_00_setup")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(requireNamespace("utils", quietly = TRUE))

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
suppressPackageStartupMessages(library(renv))

if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE)
} else {
  renv::activate()
  renv::restore(prompt = FALSE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
suppressPackageStartupMessages(library(BiocManager))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
suppressPackageStartupMessages(library(remotes))

cran_pkgs <- c(
  "microeco", "file2meco",
  "tidyverse", "readr", "dplyr", "tibble", "tidyr", "stringr", "purrr", "glue",
  "ggplot2", "ggpubr", "cowplot",
  "vegan", "ape", "sessioninfo",
  "svglite"
)

bioc_pkgs <- c(
  "phyloseq", "decontam", "DESeq2", "edgeR"
)

github_repos <- c(
  "ChiLiubio/microeco",
  "gmteunisse/ggnested",
  "taowenmicro/ggClusterNet"
)

install_cran <- function(pkgs) {
  ip <- rownames(installed.packages())
  miss <- setdiff(pkgs, ip)
  if (length(miss)) install.packages(miss, dependencies = TRUE)
  invisible(miss)
}

install_bioc <- function(pkgs) {
  ip <- rownames(installed.packages())
  miss <- setdiff(pkgs, ip)
  if (length(miss)) BiocManager::install(miss, ask = FALSE, update = FALSE)
  invisible(miss)
}

install_github <- function(repos) {
  for (repo in repos) {
    pkg <- sub(".*/", "", repo)
    if (!requireNamespace(pkg, quietly = TRUE)) {
      remotes::install_github(repo, dependencies = TRUE, upgrade = "never")
    }
  }
  invisible(TRUE)
}

install_cran(cran_pkgs)
install_bioc(bioc_pkgs)
install_github(github_repos)

renv::snapshot(prompt = FALSE)

if (!requireNamespace("sessioninfo", quietly = TRUE)) install.packages("sessioninfo")
si <- capture.output(sessioninfo::session_info())
writeLines(si, con = file.path(out_dir, "session_info.txt"))

pk <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
pk <- pk[order(pk$Package), , drop = FALSE]
readr::write_tsv(pk, file.path(out_dir, "package_versions.tsv"))

message("✅ Step 00 complete")
message("• renv.lock: ", normalizePath("renv.lock"))
message("• session:   ", normalizePath(file.path(out_dir, "session_info.txt")))
message("• versions:  ", normalizePath(file.path(out_dir, "package_versions.tsv")))
