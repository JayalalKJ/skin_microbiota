#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# Script : Step_11_taxonomic_composition_bar_hindgut.R
# Purpose: Taxonomic composition barplots (multiple ranks)
#          - Facet by Tanks × Hindgut_compartments
#          - Drop unwanted tank labels (controls/extra)
#          - OTU prevalence filter (optionally + abundance)
#          - Save TSS-normalized microtable for composition plots
# Inputs : microeco object (.RData/.rds)
# Outputs:
#   - QC/otu_prevalence_report.tsv (+ optional abundance columns)
#   - Composition_TSS.RData
#   - PNG per rank (top-N taxa + Others) using FINAL TSS object
# Author : Jayalal K Jayanthan
# Updated: 2026-02-05
# ─────────────────────────────────────────────────────────────

options(stringsAsFactors = FALSE)

## 0 │ Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(ggplot2)
  library(ggh4x)
  library(readr)
  library(tibble)
  library(stringr)
  library(grid)
  library(tools)
})

## 0.1 │ Logging --------------------------------------------------------------
stamp <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
}

## 1 │ I/O -------------------------------------------------------------------
dataset_path <- "results/Step_03C_contaminant_removal_from_blanks/microeco_cleaned_after_blank_contaminants.rds"
# Note -> ConQuR alters marginal distributions, don't need to use for compositinal analysis
#dataset_path  <- "results/Step_03D_Batchcorrection_ConQuR_2/Set_6groups/microeco_conqur_penalized_Genus_level_OTUs.rds"
output_root  <- "results/Step_08_Taxonomic_composition_HD_HM_Env_feed_water"
output_dir   <- output_root  # (kept for your requested TSS save snippet)

## 2 │ Config ----------------------------------------------------------------
# Tanks to drop (controls/extras)
drop_tanks <- c("EBs","PS_1","PS_2","PS_4","PS_6","PS_7")

# Ranks to plot (tax_table columns must exist)
rank_cfg <- list(
  Phyla             = list(taxrank = "Phyla",             ntaxa = 8),
  Classes           = list(taxrank = "Classes",           ntaxa = 10),
  Orders            = list(taxrank = "Orders",            ntaxa = 12),
  Families          = list(taxrank = "Families",          ntaxa = 15),
  Genera            = list(taxrank = "Genera",            ntaxa = 25),
  Genera_level_OTUs = list(taxrank = "Genera_level_OTUs", ntaxa = 35)
)

# Facets and required metadata
facet_vars <- c("Tanks", "Hindgut_compartments")

# Plot controls
BAR_WIDTH     <- 14.9
BAR_HEIGHT    <- 8
DPI           <- 300
X_TEXT_ANGLE  <- 65
LEGEND_NCOL   <- 1
OTHERS_COLOR  <- "grey70"
BARWIDTH      <- 0.92
BASE_FAMILY   <- "Times New Roman"

# ── OTU prevalence filtering (requested) ───────────────────────────────────
ENABLE_PREVALENCE_FILTER <- TRUE
PREVALENCE_FRAC          <- 0.10   # keep OTUs present in ≥10% of samples in THIS analysis set

# Optional extra guard (set FALSE if you only want prevalence)
ENABLE_ABUNDANCE_FILTER  <- TRUE
ABUND_MIN_READS          <- 10L    # keep OTUs with ≥10 total reads in THIS analysis set

QC_DIR <- file.path(output_root, "QC")

## 3 │ Helpers ---------------------------------------------------------------

# Robust loader (.RData/.rds) → microeco::microtable
load_microeco_object <- function(path){
  if (!file.exists(path)) stop("File not found: ", path)
  ext <- tolower(file_ext(path))
  
  if (ext == "rds") {
    obj <- readRDS(path)
    if (!inherits(obj, "microtable"))
      stop("RDS does not contain a microtable. Got: ", paste(class(obj), collapse = ", "))
    return(obj)
  }
  
  if (ext %in% c("rdata","rda")) {
    env <- new.env(parent = emptyenv())
    objs <- load(path, envir = env)
    for (nm in objs) {
      if (inherits(env[[nm]], "microtable")) return(env[[nm]])
    }
    stop("No microtable found in ", path, ". Objects: ", paste(objs, collapse = ", "))
  }
  
  stop("Unsupported extension: ", ext)
}

# PNG-only save (NO PDF) → avoids Windows path/cairo errors
save_plot_png <- function(out_dir, name, plot, w, h, dpi = 300, bg = "white") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(out_dir, paste0(name, ".png"))
  ggsave(png_path, plot = plot, width = w, height = h, dpi = dpi,
         limitsize = FALSE, bg = bg)
  stamp("Saved: ", png_path)
}

# Hard checks for required columns
assert_cols <- function(df, cols, where = "sample_table"){
  miss <- setdiff(cols, colnames(df))
  if (length(miss)) {
    stop("Missing columns in ", where, ": ", paste(miss, collapse = ", "),
         "\nAvailable: ", paste(colnames(df), collapse = ", "))
  }
  invisible(TRUE)
}

## 4 │ Load data -------------------------------------------------------------
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

stamp("Loading: ", dataset_path)
mt <- load_microeco_object(dataset_path)
if (!inherits(mt, "microtable")) stop("Input is not a microeco::microtable")

stamp("MANUSCRIPT COUNTS (input object):")
stamp("  • Total OTUs   : ", nrow(mt$otu_table))
stamp("  • Total samples: ", ncol(mt$otu_table))

st <- as.data.frame(mt$sample_table)
stamp("SAMPLE TABLE: ", nrow(st), " rows × ", ncol(st), " cols")
assert_cols(st, facet_vars, where = "sample_table")

## 5 │ Filter tanks + align OTU table ----------------------------------------
mt$sample_table <- mt$sample_table %>%
  mutate(across(where(is.character), ~trimws(.x))) %>%
  filter(!.data[["Tanks"]] %in% drop_tanks)

mt$sample_table$Tanks <- droplevels(as.factor(mt$sample_table$Tanks))
mt$sample_table$Hindgut_compartments <- droplevels(as.factor(mt$sample_table$Hindgut_compartments))

keep_samp <- rownames(mt$sample_table)
mt$otu_table <- mt$otu_table[, keep_samp, drop = FALSE]

if (ncol(mt$otu_table) == 0) stop("No samples left after filtering Tanks. Check drop_tanks names.")
stopifnot(identical(rownames(mt$sample_table), colnames(mt$otu_table)))

mt$tidy_dataset()

stamp("AFTER dropping unwanted tank labels:")
stamp("  • OTUs   : ", nrow(mt$otu_table))
stamp("  • Samples: ", ncol(mt$otu_table))

write_tsv(tibble(Sample_ID = rownames(mt$sample_table)),
          file.path(QC_DIR, "samples_used_after_tank_filter.tsv"))

## 5b │ OTU prevalence (and optional abundance) filtering --------------------
# NOTE: mt$otu_table is OTUs (rows) × Samples (cols)
if (isTRUE(ENABLE_PREVALENCE_FILTER)) {
  
  otu_tab <- mt$otu_table
  n_samp  <- ncol(otu_tab)
  if (n_samp < 2) stop("Not enough samples for prevalence filtering.")
  
  prev_threshold <- max(1L, ceiling(PREVALENCE_FRAC * n_samp))
  prev_counts    <- rowSums(otu_tab > 0)
  tot_reads      <- rowSums(otu_tab)
  
  keep_prev <- prev_counts >= prev_threshold
  
  if (isTRUE(ENABLE_ABUNDANCE_FILTER)) {
    keep_abun <- tot_reads >= ABUND_MIN_READS
  } else {
    keep_abun <- rep(TRUE, length(tot_reads))
  }
  
  keep_otus <- keep_prev & keep_abun
  
  report <- tibble(
    OTU                = rownames(otu_tab),
    prevalence_samples = prev_counts,
    prevalence_pass    = keep_prev,
    total_reads        = tot_reads,
    abundance_pass     = keep_abun,
    kept               = keep_otus
  )
  write_tsv(report, file.path(QC_DIR, "otu_prevalence_report.tsv"))
  
  stamp("🔍 OTU filtering summary:")
  stamp("  • Samples used: ", n_samp)
  stamp("  • Prevalence threshold: ceiling(", PREVALENCE_FRAC, " × ", n_samp, ") = ", prev_threshold)
  stamp("  • Original OTUs: ", nrow(otu_tab))
  stamp("  • Retained OTUs: ", sum(keep_otus))
  
  mt$otu_table <- otu_tab[keep_otus, , drop = FALSE]
  
  if (!is.null(mt$tax_table) && nrow(mt$tax_table) > 0) {
    common <- intersect(rownames(mt$tax_table), rownames(mt$otu_table))
    mt$tax_table <- mt$tax_table[common, , drop = FALSE]
    mt$otu_table <- mt$otu_table[common, , drop = FALSE]
  }
  
  zero_samp <- colSums(mt$otu_table) == 0
  if (any(zero_samp)) {
    write_tsv(tibble(Sample_ID = colnames(mt$otu_table)[zero_samp]),
              file.path(QC_DIR, "samples_dropped_zero_after_otu_filter.tsv"))
    keep_ids <- colnames(mt$otu_table)[!zero_samp]
    mt$otu_table    <- mt$otu_table[, keep_ids, drop = FALSE]
    mt$sample_table <- mt$sample_table[keep_ids, , drop = FALSE]
    stamp("⚠️ Samples removed (zero after OTU filter): ", sum(zero_samp))
  }
  
  mt$tidy_dataset()
  
  stamp("AFTER OTU filtering:")
  stamp("  • OTUs   : ", nrow(mt$otu_table))
  stamp("  • Samples: ", ncol(mt$otu_table))
} else {
  stamp("OTU prevalence filtering is DISABLED.")
}

## 5c │ TSS normalization + save (requested) ---------------------------------
tn  <- trans_norm$new(dataset = mt)
TSS <- tn$norm(method = "TSS")

save(TSS,
     file = file.path(output_dir, "Composition_TSS.RData"),
     compress = TRUE)

message("✓ Saved Composition_TSS.RData (sample-only) for compositional barplots.")

mt_plot <- TSS
mt_plot$tidy_dataset()
mt_plot$cal_abund()

## 6 │ Plot all ranks (PNG only) ---------------------------------------------
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

for (nm in names(rank_cfg)) {
  
  taxrank <- rank_cfg[[nm]]$taxrank
  ntaxa   <- rank_cfg[[nm]]$ntaxa
  
  if (!taxrank %in% colnames(mt_plot$tax_table)) {
    stop("taxrank not found in tax_table: ", taxrank,
         "\nAvailable: ", paste(colnames(mt_plot$tax_table), collapse = ", "))
  }
  
  stamp("Plotting: ", nm, " | taxrank=", taxrank, " | top=", ntaxa)
  
  t1 <- trans_abund$new(dataset = mt_plot, taxrank = taxrank, ntaxa = ntaxa)
  
  p <- t1$plot_bar(
    others_color = OTHERS_COLOR,
    facet = facet_vars,
    xtext_keep = FALSE,
    legend_text_italic = FALSE,
    barwidth = BARWIDTH
  ) +
    theme_classic(base_family = BASE_FAMILY) +
    theme(
      strip.text       = element_text(size = 12, face = "bold"),
      axis.title.x     = element_text(size = 10, face = "bold"),
      axis.title.y     = element_text(size = 10, face = "bold"),
      axis.text.y      = element_text(size = 10),
      axis.text.x      = element_text(angle = X_TEXT_ANGLE, hjust = 1, vjust = 1, size = 9),
      legend.position  = "right",
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 9),
      legend.key.size  = unit(0.40, "cm"),
      legend.spacing   = unit(0.05, "cm"),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin      = margin(3, 3, 3, 3)
    ) +
    labs(
      title = paste0(nm, " (Top ", ntaxa, ")"),
      x     = "Samples",
      y     = "Relative abundance (%)",
      fill  = nm
    ) +
    guides(fill = guide_legend(ncol = LEGEND_NCOL, byrow = TRUE), colour = "none")
  
  out_dir_rank <- file.path(output_root, nm)
  
  # shorter filename helps Windows too
  out_name <- paste0("BAR_", nm, "_top", ntaxa, "_TSS")
  
  save_plot_png(
    out_dir = out_dir_rank,
    name    = out_name,
    plot    = p,
    w       = BAR_WIDTH,
    h       = BAR_HEIGHT,
    dpi     = DPI
  )
}

stamp("✅ Done. Output folder: ", normalizePath(output_root))
