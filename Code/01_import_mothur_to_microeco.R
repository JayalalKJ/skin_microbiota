#!/usr/bin/env Rscript
# =============================================================================
# Jayalal K Jayanthan
# PhD candidate (2021– 2026)
# Research group: Seafood Science
# Institute: The Norwegian College of Fishery Science
# Faculty: Faculty of Biosciences, Fisheries, and Economics
# Campus: UiT Campus Tromsø
# UiT The Arctic University of Norway
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(microeco)
  library(file2meco)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ape)
})

# ---- I/O -------------------------------------------------------------------
INPUT_DIR  <- "data/97"
OUTPUT_DIR <- "results/Step_01_import_mothur_to_microeco"
REPORT_DIR <- file.path(OUTPUT_DIR, "reports")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

FILES <- list(
  metadata  = "Gut_updated_metadata_modified_V1.csv",
  shared    = "AlgOpti1.trim.contigs.renamed.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.shared",
  taxonomy  = "AlgOpti1.trim.contigs.renamed.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.cons.taxonomy",
  tree      = "AlgOpti1.trim.contigs.renamed.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.rep.phylip.tre"
)

meta_path   <- file.path(INPUT_DIR, FILES$metadata)
shared_path <- file.path(INPUT_DIR, FILES$shared)
tax_path    <- file.path(INPUT_DIR, FILES$taxonomy)
tree_path   <- file.path(INPUT_DIR, FILES$tree)

stopifnot(file.exists(meta_path), file.exists(shared_path), file.exists(tax_path), file.exists(tree_path))

# ---- Helpers ---------------------------------------------------------------
trim_chr <- function(x) trimws(as.character(x))

safe_write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(as.data.frame(df), path)
}

# Robust genus cleaning: removes rank prefixes and flags common "unknown" tokens
clean_rank_prefix <- function(x) {
  x <- trimws(as.character(x))
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_replace_all(x, "^g__", "")
  x <- str_replace_all(x, "^g_", "")
  trimws(x)
}

is_unclassified_token <- function(x_clean) {
  x0 <- tolower(trimws(as.character(x_clean)))
  x0 %in% c("", "na", "n/a", "unknown", "unclassified", "uncultured", "unidentified", "metagenome") |
    grepl("uncultured|unidentified|unclassified|unknown|metagenome", x0)
}

# ---- 1) Metadata -----------------------------------------------------------
metadata <- read_csv(meta_path, show_col_types = FALSE) %>%
  mutate(Sample_ID = trim_chr(Sample_ID)) %>%
  distinct(Sample_ID, .keep_all = TRUE)

stopifnot("Sample_ID" %in% names(metadata))
message("✅ metadata: ", nrow(metadata), " rows × ", ncol(metadata), " cols")

# ---- 2) Shared -> OTU table ------------------------------------------------
otu_df <- read_tsv(shared_path, show_col_types = FALSE)

req  <- c("label", "Group", "numOtus")
miss <- setdiff(req, names(otu_df))
if (length(miss)) stop("❌ mothur shared missing columns: ", paste(miss, collapse = ", "))

otu_counts <- otu_df %>%
  dplyr::select(-label, -numOtus) %>%   # qualify to avoid MASS conflicts
  dplyr::rename(Sample_ID = Group) %>%
  dplyr::mutate(Sample_ID = trim_chr(Sample_ID)) %>%
  tibble::column_to_rownames("Sample_ID")

otu_mat <- as.matrix(otu_counts)
storage.mode(otu_mat) <- "integer"
message("✅ shared OTU table: ", nrow(otu_mat), " samples × ", ncol(otu_mat), " OTUs")

# ---- 3) Taxonomy (parse + enforce explicit unclassified genus labels) ------
tax_raw <- read_tsv(
  tax_path,
  col_types = cols(
    OTU      = col_character(),
    Taxonomy = col_character(),
    Size     = col_double()
  )
) %>%
  rename(otu = OTU, taxonomy = Taxonomy, size = Size) %>%
  mutate(
    taxonomy = taxonomy %>%
      str_remove_all("\\(\\d+\\)") %>%
      str_remove_all(";$")
  ) %>%
  separate(
    taxonomy,
    into  = c("Kingdoms","Phyla","Classes","Orders","Families","Genera", "Genera_level_OTUs"),
    sep   = ";",
    fill  = "right",
    extra = "drop"
  ) %>%
  mutate(across(c(Kingdoms, Phyla, Classes, Orders, Families, Genera), ~ na_if(str_squish(.), ""))) %>%
  distinct(otu, .keep_all = TRUE) %>%
  column_to_rownames("otu")

message("✅ taxonomy rows: ", nrow(tax_raw), " OTUs")

# 3b) Diagnostics: raw genus frequency (helps confirm presence of "g__" etc.)
genus_freq <- tibble(Genera_raw = as.character(tax_raw$Genera)) %>%
  mutate(Genera_raw = ifelse(is.na(Genera_raw), "NA", trimws(Genera_raw))) %>%
  count(Genera_raw, sort = TRUE)

safe_write_tsv(genus_freq, file.path(REPORT_DIR, "genus_raw_frequency.tsv"))

# 3c) Robust unclassified-genus labeling
genus_raw   <- as.character(tax_raw$Genera)
genus_clean <- clean_rank_prefix(genus_raw)

unclassified_mask <- is.na(genus_raw) | is_unclassified_token(genus_clean)

tax_raw <- tax_raw %>%
  mutate(
    Genus_clean = ifelse(unclassified_mask, "Unclassified_genus", genus_clean),
    Genera_level_OTUs = paste0(Genus_clean, "_", rownames(.))
  )

unclassified_log <- tibble(
  OTU         = rownames(tax_raw),
  Genera_raw  = genus_raw,
  Genus_clean = tax_raw$Genus_clean,
  Label_used  = tax_raw$Genera_level_OTUs
) %>% filter(Genus_clean == "Unclassified_genus")

safe_write_tsv(unclassified_log, file.path(REPORT_DIR, "unclassified_genus_OTUs.tsv"))
message("✅ unclassified genus OTUs: ", nrow(unclassified_log), " (logged)")

# ---- 4) Harmonize OTUs -----------------------------------------------------
otus_common <- intersect(colnames(otu_mat), rownames(tax_raw))
if (!length(otus_common)) stop("❌ No overlapping OTU IDs between shared and taxonomy")

dropped_from_otu <- setdiff(colnames(otu_mat), otus_common)
dropped_from_tax <- setdiff(rownames(tax_raw), otus_common)

if (length(dropped_from_otu)) safe_write_tsv(tibble(OTU = dropped_from_otu), file.path(REPORT_DIR, "otus_dropped_missing_in_taxonomy.tsv"))
if (length(dropped_from_tax)) safe_write_tsv(tibble(OTU = dropped_from_tax), file.path(REPORT_DIR, "otus_dropped_missing_in_shared.tsv"))

otu_mat <- otu_mat[, otus_common, drop = FALSE]
tax_df  <- tax_raw[otus_common, , drop = FALSE]

# ---- 5) Harmonize Samples --------------------------------------------------
samples_before <- rownames(otu_mat)
meta_ids       <- metadata$Sample_ID

common_samples <- samples_before[samples_before %in% meta_ids]
lost_no_meta   <- setdiff(samples_before, meta_ids)
meta_only      <- setdiff(meta_ids, samples_before)

if (length(lost_no_meta)) safe_write_tsv(tibble(Sample_ID = lost_no_meta), file.path(REPORT_DIR, "samples_dropped_no_metadata.tsv"))
if (length(meta_only))    safe_write_tsv(tibble(Sample_ID = meta_only),    file.path(REPORT_DIR, "metadata_rows_without_counts.tsv"))

otu_mat <- otu_mat[common_samples, , drop = FALSE]

metadata_ps <- metadata %>%
  filter(Sample_ID %in% common_samples) %>%
  arrange(match(Sample_ID, rownames(otu_mat))) %>%
  column_to_rownames("Sample_ID")

stopifnot(identical(rownames(metadata_ps), rownames(otu_mat)))
message("✅ working dataset: ", nrow(otu_mat), " samples × ", ncol(otu_mat), " OTUs")

# ---- 6) Tree ---------------------------------------------------------------
tr <- read.tree(tree_path)
if (is.null(tr$tip.label)) stop("❌ tree has no tip labels")

keep_tips <- intersect(tr$tip.label, colnames(otu_mat))
drop_tips <- setdiff(tr$tip.label, keep_tips)
if (length(drop_tips)) tr <- drop.tip(tr, drop_tips)

stopifnot(all(colnames(otu_mat) %in% tr$tip.label))

# ---- 7) Build phyloseq -----------------------------------------------------
ps <- phyloseq(
  otu_table(t(otu_mat), taxa_are_rows = TRUE),
  tax_table(as.matrix(tax_df)),
  sample_data(metadata_ps),
  phy_tree(tr)
)

saveRDS(ps, file.path(OUTPUT_DIR, "phyloseq_raw.rds"))
message("💾 saved: ", file.path(OUTPUT_DIR, "phyloseq_raw.rds"))

# ---- 8) Convert to microeco ------------------------------------------------
mt <- phyloseq2meco(ps)
mt$tidy_dataset()
mt$cal_abund()

saveRDS(mt, file.path(OUTPUT_DIR, "microeco_raw.rds"))
message("💾 saved: ", file.path(OUTPUT_DIR, "microeco_raw.rds"))

# ---- 9) Minimal QC exports -------------------------------------------------
qc <- tibble(
  n_samples    = nsamples(ps),
  n_otus       = ntaxa(ps),
  min_depth    = min(sample_sums(ps)),
  median_depth = as.numeric(stats::median(sample_sums(ps))),
  max_depth    = max(sample_sums(ps))
)
safe_write_tsv(qc, file.path(REPORT_DIR, "import_qc_summary.tsv"))

depths <- tibble(Sample_ID = names(sample_sums(ps)), LibSize = as.integer(sample_sums(ps)))
safe_write_tsv(depths, file.path(REPORT_DIR, "library_sizes.tsv"))

message("🎉 Step 01 complete")
message("• Output: ", normalizePath(OUTPUT_DIR))
