#!/usr/bin/env Rscript
# =============================================================================
# 03_contaminant_removal_from_blanks.R
# Jayalal K Jayanthan
# PhD candidate (2021– 2026)
# Research group: Seafood Science
# Institute: The Norwegian College of Fishery Science
# Faculty: Faculty of Biosciences, Fisheries, and Economics
# Campus: UiT Campus Tromsø
# UiT The Arctic University of Norway
# =============================================================================

suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

# ---- I/O -------------------------------------------------------------------
input_rds  <- "results/Step_02_remove_spikein_taxa/microeco_spikein_removed.rds"
output_dir <- "results/Step_03C_contaminant_removal_from_blanks"
fig_dir    <- file.path(output_dir, "figures")
rep_dir    <- file.path(output_dir, "reports")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(rep_dir,    recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds))

mt <- readRDS(input_rds)
if (!inherits(mt, "microtable")) stop("❌ Input is not a microeco::microtable")

# ---- Config ----------------------------------------------------------------
# 1) Negative controls (sample IDs or rownames in sample_table)
negative_ids <- c("2N1","2N3","3N1","3N2","3N3","EB4","N2","N3","NGHF")

# 2) Optional manual sample removals (known problems)
manual_drop_ids <- character(0)  # e.g., c("HM2D")

# 3) Contaminant labels to silence (must match the genus label column exactly)
#    Recommended: paste your final list here (from your blanks analysis).
contaminant_taxa <- c(
  "g__Truepera_Otu00002",
  "g__Bacillus_Otu00017",
  "g__Salmonella_Otu00018",
  "g__Enterococcus_Otu00071",
  "g__Phyllobacterium_Otu00011",
  "g__Sphingomonas_Otu00003",
  "g__Sphingomonas_Otu00008",
  "g__Sphingomonas_Otu00197",
  "g__Novosphingobium_Otu00006",
  "g__Novosphingobium_Otu00045",
  "g__Methylobacterium_Otu00013",
  "g__Methylobacterium_Otu00024",
  "g__Methylobacterium_Otu00052",
  "g__Methylobacterium_Otu00068",
  "g__Methylobacterium_Otu00085",
  "g__Caulobacter_Otu00036",
  "g__Bradyrhizobium_Otu00077",
  "g__Sphingobium_Otu00026",
  "g__Brevundimonas_Otu00094",
  "g__Brevundimonas_Otu00086",
  "g__Escherichia-Shigella_Otu00019",
  "g__Methyloversatilis_Otu00004",
  "g__Acinetobacter_Otu00009",
  "g__Acinetobacter_Otu00010",
  "g__Acinetobacter_Otu00084",
  "g__Pseudomonas_Otu00039",
  "g__Pseudomonas_Otu00035",
  "g__Pseudomonas_Otu00145",
  "g__Pseudomonas_Otu00103",
  "g__Pseudomonas_Otu00107",
  "g__Acidovorax_Otu00029",
  "g__Acidovorax_Otu00048",
  "g__Acidovorax_Otu00289",
  "g__Burkholderia-Caballeronia-Paraburkholderia_Otu00092",
  "g__Ralstonia_Otu00173",
  "g__Hydrogenophaga_Otu00138",
  "g__Curvibacter_Otu00517",
  "g__Microbacterium_Otu00015",
  "g__Cutibacterium_Otu00025",
  "g__Micrococcus_Otu00043",
  "g__Meiothermus_Otu00033",
  "g__Listeria_Otu00076",
  "g__Staphylococcus_Otu00047",
  "g__Cloacibacterium_Otu00250"
)

# Plot quick before/after top taxa bars (safe, simple)
make_qc_plots <- TRUE
top_n <- 30

# ---- Helpers ----------------------------------------------------------------
ensure_dir <- function(p) { dir.create(p, recursive = TRUE, showWarnings = FALSE); invisible(p) }

safe_write_tsv <- function(df, path) {
  ensure_dir(dirname(path))
  readr::write_tsv(as.data.frame(df), path)
}

pick_tax_col <- function(tt_cols) {
  if ("Genera_level_OTUs" %in% tt_cols) return("Genera_level_OTUs")
  if ("Genus_level_OTUs" %in% tt_cols) return("Genus_level_OTUs")
  if ("Genus" %in% tt_cols) return("Genus")
  if ("Genera" %in% tt_cols) return("Genera")
  NA_character_
}

drop_samples <- function(m, ids, label, out_log) {
  st <- as.data.frame(m$sample_table)
  id_col <- if ("Sample_ID" %in% names(st)) "Sample_ID" else NA_character_

  if (!is.na(id_col)) {
    rm_mask <- st[[id_col]] %in% ids | rownames(st) %in% ids
  } else {
    rm_mask <- rownames(st) %in% ids
  }

  removed <- rownames(st)[rm_mask]
  if (!length(removed)) {
    return(list(mt = m, removed = character(0)))
  }

  safe_write_tsv(tibble(Removed_Rowname = removed,
                        Sample_ID = if (!is.na(id_col)) st[[id_col]][rm_mask] else removed,
                        Reason = label),
                 out_log)

  m$sample_table <- m$sample_table[!rm_mask, , drop = FALSE]
  m$otu_table    <- m$otu_table[, rownames(m$sample_table), drop = FALSE]
  stopifnot(identical(rownames(m$sample_table), colnames(m$otu_table)))

  list(mt = m, removed = removed)
}

prune_zero_otus <- function(m) {
  keep <- rownames(m$otu_table)[rowSums(m$otu_table) > 0]
  m$otu_table <- m$otu_table[keep, , drop = FALSE]
  m$tax_table <- m$tax_table[keep, , drop = FALSE]
  m
}

plot_top_taxa <- function(m, tax_col, out_stub) {
  tt <- as.data.frame(m$tax_table)
  ok <- !is.na(tt[[tax_col]]) & tt[[tax_col]] != ""
  tt <- tt[ok, , drop = FALSE]

  counts <- rowSums(m$otu_table[rownames(tt), , drop = FALSE])
  df <- tibble(Taxon = as.character(tt[[tax_col]]), Reads = as.numeric(counts)) %>%
    group_by(Taxon) %>%
    summarise(Reads = sum(Reads), .groups = "drop") %>%
    arrange(desc(Reads)) %>%
    slice_head(n = top_n)

  p <- ggplot(df, aes(x = reorder(Taxon, Reads), y = Reads)) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 13) +
    labs(x = NULL, y = "Total reads", title = NULL)

  ggsave(paste0(out_stub, ".png"), p, width = 10, height = 7, dpi = 450, bg = "white")
  ggsave(paste0(out_stub, ".pdf"), p, width = 10, height = 7, device = "pdf")
  invisible(p)
}

# ---- Resolve taxonomy column ------------------------------------------------
tax_col <- pick_tax_col(colnames(mt$tax_table))
if (is.na(tax_col)) stop("❌ No usable genus label column in mt$tax_table")

# ---- Snapshot QC (before) --------------------------------------------------
qc_before <- tibble(
  n_samples = nrow(mt$sample_table),
  n_otus    = nrow(mt$otu_table),
  min_depth = min(colSums(mt$otu_table)),
  median_depth = as.numeric(stats::median(colSums(mt$otu_table))),
  max_depth = max(colSums(mt$otu_table))
)
safe_write_tsv(qc_before, file.path(rep_dir, "qc_before.tsv"))

if (isTRUE(make_qc_plots)) {
  plot_top_taxa(mt, tax_col, file.path(fig_dir, "top_taxa_before"))
}

# ---- 1) Manual drop(s) ------------------------------------------------------
if (length(manual_drop_ids)) {
  res <- drop_samples(
    mt,
    ids = manual_drop_ids,
    label = "manual_drop",
    out_log = file.path(rep_dir, "samples_dropped_manual.tsv")
  )
  mt <- res$mt
}

# ---- 2) Silence contaminants by tax label ----------------------------------
tax_df <- as.data.frame(mt$tax_table)
hit_idx <- which(tax_df[[tax_col]] %in% contaminant_taxa)
hit_otus <- rownames(tax_df)[hit_idx]
hit_otus <- intersect(hit_otus, rownames(mt$otu_table))

silence_log <- tibble(
  tax_col = tax_col,
  n_taxa_requested = length(contaminant_taxa),
  n_tax_rows_matched = length(hit_idx),
  n_otus_silenced = length(hit_otus)
)
safe_write_tsv(silence_log, file.path(rep_dir, "contaminant_silencing_summary.tsv"))

if (length(hit_idx)) {
  safe_write_tsv(
    tibble(OTU = rownames(tax_df)[hit_idx], TaxLabel = tax_df[[tax_col]][hit_idx]),
    file.path(rep_dir, "otus_silenced_by_label.tsv")
  )
}

if (length(hit_otus)) {
  mt$otu_table[hit_otus, ] <- 0
}

# ---- 3) Prune zeros + refresh ------------------------------------------------
mt <- prune_zero_otus(mt)

# ---- 4) Remove negative controls -------------------------------------------
if (length(negative_ids)) {
  res2 <- drop_samples(
    mt,
    ids = negative_ids,
    label = "negative_control",
    out_log = file.path(rep_dir, "samples_dropped_negative_controls.tsv")
  )
  mt <- res2$mt
  mt <- prune_zero_otus(mt)
}

mt$tidy_dataset()
mt$cal_abund()

# ---- Snapshot QC (after) ----------------------------------------------------
qc_after <- tibble(
  n_samples = nrow(mt$sample_table),
  n_otus    = nrow(mt$otu_table),
  min_depth = min(colSums(mt$otu_table)),
  median_depth = as.numeric(stats::median(colSums(mt$otu_table))),
  max_depth = max(colSums(mt$otu_table))
)
safe_write_tsv(qc_after, file.path(rep_dir, "qc_after.tsv"))

if (isTRUE(make_qc_plots)) {
  plot_top_taxa(mt, tax_col, file.path(fig_dir, "top_taxa_after"))
}

# ---- Save ------------------------------------------------------------------
out_rds <- file.path(output_dir, "microeco_cleaned_after_blank_contaminants.rds")
saveRDS(mt, out_rds)

message("✅ Step 03 complete")
message("• Input : ", normalizePath(input_rds))
message("• Output: ", normalizePath(out_rds))
message("• Reports: ", normalizePath(rep_dir))
message("• Figures: ", normalizePath(fig_dir))
