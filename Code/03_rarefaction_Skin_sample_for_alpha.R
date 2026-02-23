#!/usr/bin/env Rscript
# =============================================================================
# 04_rarefaction_HMHD_for_alpha.R
# Jayalal K Jayanthan
# =============================================================================

set.seed(123)

suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(vegan)
})

# ---- I/O -------------------------------------------------------------------
#input_rds  <- "results/Step_03C_contaminant_removal_from_blanks/microeco_cleaned_after_blank_contaminants.rds"
input_rds  <- "results/Step_03D_Batchcorrection_ConQuR_2/Set_6groups/microeco_conqur_penalized_Genus_level_OTUs.rds"

output_dir <- "results/Step_04_rarefaction_HMHD_for_alpha"
qc_dir     <- file.path(output_dir, "QC")
fig_dir    <- file.path(qc_dir, "figures")
tab_dir    <- file.path(qc_dir, "tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir,    recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds))

mt <- readRDS(input_rds)
if (!inherits(mt, "microtable")) stop("❌ Input is not a microeco::microtable")

# ---- Config ----------------------------------------------------------------
cap_depth_plot      <- 60000L
cap_depth_analysis  <- 60000L
retain_prop_target  <- 0.80     # keep >= 80% HM/HD samples
min_per_group       <- 8L
min_floor           <- 1000L
n_iter_sim          <- 30L

exclude_ids <- c(
  "Swabcontrol",
  "T1","T2","T3",
  "T1FD1","T1FD2","T1FD3",
  "T2FD1","T2FD2","T2FD3",
  "T3FD1","T3FD2","T3FD3"
)

compartment_candidates <- c("Hindgut_compartments","Hindgut_compartments_1","Compartment")
tank_candidates        <- c("Tank","Tanks","Tanks_1")
fish_candidates        <- c("Fish_ID","Fish_number","FishID")

# ---- Helpers ----------------------------------------------------------------
trim_chr <- function(x) trimws(as.character(x))
ensure_dir <- function(p) { dir.create(p, recursive = TRUE, showWarnings = FALSE); invisible(p) }

safe_write_tsv <- function(df, file) {
  ensure_dir(dirname(file))
  readr::write_tsv(as.data.frame(df), file)
}

detect_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)][1]
  if (length(hit) == 0) NA_character_ else hit
}

subset_microtable_by_samples <- function(m, sample_ids) {
  mm <- m$clone(deep = TRUE)
  sample_ids <- intersect(sample_ids, colnames(mm$otu_table))
  if (!length(sample_ids)) stop("❌ subset: 0 samples after intersect")
  mm$otu_table    <- mm$otu_table[, sample_ids, drop = FALSE]
  mm$sample_table <- mm$sample_table[sample_ids, , drop = FALSE]
  mm$tidy_dataset()
  mm
}

safe_table_to_tsv <- function(tab, file) {
  df <- as.data.frame(tab, stringsAsFactors = FALSE)
  if ("Freq" %in% names(df)) names(df)[names(df) == "Freq"] <- "Count"
  safe_write_tsv(df, file)
}

get_depths <- function(m) {
  d <- as.integer(colSums(m$otu_table))
  names(d) <- colnames(m$otu_table)
  d
}

# ---- 1) Harmonize sample IDs -----------------------------------------------
st <- as.data.frame(mt$sample_table)
otu_cols <- colnames(mt$otu_table)

overlap <- if (!is.null(rownames(st))) mean(rownames(st) %in% otu_cols) else 0

if (!is.finite(overlap) || overlap < 0.6) {
  id_candidates <- intersect(c("Sample_ID","Sample_id","SampleID","Sample","Sample2"), names(st))
  if (!length(id_candidates)) stop("❌ Cannot align sample_table to otu_table (no ID column found)")
  score <- vapply(id_candidates, function(cc) mean(trim_chr(st[[cc]]) %in% otu_cols), numeric(1))
  id_col <- names(which.max(score))
  if (!is.finite(score[id_col]) || score[id_col] < 0.2) stop("❌ Cannot align sample IDs (low overlap)")
  st[[id_col]] <- trim_chr(st[[id_col]])
  rownames(st) <- st[[id_col]]
}

common <- intersect(colnames(mt$otu_table), rownames(st))
if (!length(common)) stop("❌ No overlap between otu_table and sample_table IDs")
common <- sort(common)

mt$otu_table    <- mt$otu_table[, common, drop = FALSE]
mt$sample_table <- st[common, , drop = FALSE]
stopifnot(identical(rownames(mt$sample_table), colnames(mt$otu_table)))

# ---- 2) Drop explicit non-hindgut IDs --------------------------------------
exclude_present <- intersect(colnames(mt$otu_table), exclude_ids)
if (length(exclude_present)) {
  safe_write_tsv(tibble(Sample_ID = exclude_present), file.path(tab_dir, "dropped_excluded_ids.tsv"))
  mt <- subset_microtable_by_samples(mt, setdiff(colnames(mt$otu_table), exclude_present))
}

# ---- 3) Subset to HM/HD -----------------------------------------------------
st2 <- as.data.frame(mt$sample_table)

col_comp <- detect_col(st2, compartment_candidates)
col_tank <- detect_col(st2, tank_candidates)
col_fish <- detect_col(st2, fish_candidates)

if (is.na(col_comp)) stop("❌ Missing hindgut compartment column (expected Hindgut_compartments or similar)")
if (is.na(col_tank)) stop("❌ Missing Tank column")

st2[[col_comp]] <- trim_chr(st2[[col_comp]])
st2[[col_tank]] <- trim_chr(st2[[col_tank]])
if (!is.na(col_fish)) st2[[col_fish]] <- trim_chr(st2[[col_fish]])

hindgut_ids <- rownames(st2)[st2[[col_comp]] %in% c("HM","HD")]
if (!length(hindgut_ids)) stop("❌ No HM/HD samples found. Check metadata values (HM/HD).")

hg <- subset_microtable_by_samples(mt, hindgut_ids)

# Minimal derived fields used downstream (+ persist FishUID for pairing)
st_hg <- as.data.frame(hg$sample_table) %>%
  mutate(
    Compartment = trim_chr(.data[[col_comp]]),
    Tank        = trim_chr(.data[[col_tank]]),
    Diet = case_when(
      grepl("Tank_1", Tank) ~ "0% (Control)",
      grepl("Tank_2", Tank) ~ "1% Poro",
      grepl("Tank_3", Tank) ~ "2.2% Poro",
      TRUE ~ Tank
    ),
    Diet = factor(Diet, levels = c("0% (Control)","1% Poro","2.2% Poro"))
  )

if (!is.na(col_fish) && col_fish %in% names(st_hg)) {
  st_hg <- st_hg %>%
    mutate(
      Fish_ID = trim_chr(.data[[col_fish]]),
      FishUID = paste(Tank, Fish_ID, sep = "__")
    )
}

hg$sample_table <- st_hg
hg$tidy_dataset()

# ---- 4) Deduplicate resequenced libraries ----------------------------------
depths0 <- get_depths(hg)
hg$sample_table$LibSize <- as.integer(depths0[rownames(hg$sample_table)])

if ("FishUID" %in% names(hg$sample_table)) {
  st3 <- as.data.frame(hg$sample_table) %>%
    rownames_to_column("Sample_ID")
  
  keep_best <- st3 %>%
    group_by(FishUID, Compartment) %>%
    slice_max(order_by = LibSize, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    pull(Sample_ID)
  
  dropped <- setdiff(rownames(hg$sample_table), keep_best)
  if (length(dropped)) {
    safe_write_tsv(tibble(Sample_ID = dropped), file.path(tab_dir, "dropped_resequenced_duplicates.tsv"))
  }
  
  hg <- subset_microtable_by_samples(hg, keep_best)
}

# ---- 5) QC tables -----------------------------------------------------------
st_qc <- as.data.frame(hg$sample_table)

safe_table_to_tsv(table(st_qc$Tank, useNA = "ifany"), file.path(tab_dir, "QC_Tank_counts.tsv"))
safe_table_to_tsv(table(st_qc$Compartment, useNA = "ifany"), file.path(tab_dir, "QC_Compartment_counts.tsv"))
safe_table_to_tsv(table(st_qc$Diet, st_qc$Compartment, useNA = "ifany"), file.path(tab_dir, "QC_Diet_by_Compartment.tsv"))

depths <- get_depths(hg)
depth_df <- tibble(Sample_ID = names(depths), LibSize = as.integer(depths)) %>%
  left_join(
    tibble(Sample_ID = rownames(st_qc),
           Diet = st_qc$Diet,
           Tank = st_qc$Tank,
           Compartment = st_qc$Compartment),
    by = "Sample_ID"
  )
safe_write_tsv(depth_df, file.path(tab_dir, "QC_library_sizes_HMHD.tsv"))

# ---- 6) Rarefaction curves (QC) --------------------------------------------
max_depth_plot <- min(cap_depth_plot, max(depths, na.rm = TRUE))

depth_grid <- sort(unique(c(
  10L, 50L, 200L, 500L, 1000L, 2000L, 4000L, 6000L,
  10000L, 15000L, 20000L, 30000L, 40000L, 50000L, 60000L
)))
depth_grid <- depth_grid[depth_grid <= max_depth_plot]

mat <- t(as.matrix(hg$otu_table))  # samples × features

sim_one_depth <- function(d) {
  eligible <- rowSums(mat) >= d
  if (!any(eligible)) return(NULL)
  mat_e <- mat[eligible, , drop = FALSE]
  
  obs_acc <- matrix(NA_real_, nrow = nrow(mat_e), ncol = n_iter_sim, dimnames = list(rownames(mat_e), NULL))
  sha_acc <- matrix(NA_real_, nrow = nrow(mat_e), ncol = n_iter_sim, dimnames = list(rownames(mat_e), NULL))
  
  for (i in seq_len(n_iter_sim)) {
    rare_e <- vegan::rrarefy(mat_e, sample = d)
    obs_acc[, i] <- rowSums(rare_e > 0)
    sha_acc[, i] <- vegan::diversity(rare_e, index = "shannon")
  }
  
  tibble(Sample_ID = rownames(mat_e), Depth = d,
         Observed = rowMeans(obs_acc, na.rm = TRUE),
         Shannon  = rowMeans(sha_acc, na.rm = TRUE))
}

df_sim <- bind_rows(lapply(depth_grid, sim_one_depth)) %>%
  left_join(depth_df %>% select(Sample_ID, Diet, Compartment), by = "Sample_ID") %>%
  tidyr::pivot_longer(cols = c("Observed","Shannon"), names_to = "Metric", values_to = "Value")

safe_write_tsv(df_sim, file.path(tab_dir, "rarefaction_curve_simulated_values.tsv"))

df_sum <- df_sim %>%
  group_by(Compartment, Diet, Metric, Depth) %>%
  summarize(n = sum(!is.na(Value)),
            mean = mean(Value, na.rm = TRUE),
            se = sd(Value, na.rm = TRUE) / sqrt(n),
            .groups = "drop")

plot_metric <- function(metric_name, out_stub) {
  p <- ggplot(df_sum %>% filter(Metric == metric_name),
              aes(x = Depth, y = mean, color = Diet, fill = Diet)) +
    geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), alpha = 0.15, linewidth = 0) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~Compartment, scales = "free_y") +
    theme_classic(base_size = 14) +
    labs(x = "Sequencing depth (reads)", y = metric_name) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(fig_dir, paste0(out_stub, ".png")), p, width = 10, height = 6, dpi = 350, bg = "white")
  ggsave(file.path(fig_dir, paste0(out_stub, ".pdf")), p, width = 10, height = 6, device = "pdf")
}

plot_metric("Observed", "RarefactionCurve_Observed_HMHD")
plot_metric("Shannon",  "RarefactionCurve_Shannon_HMHD")

# ---- 7) Choose rarefaction depth -------------------------------------------
max_depth_analysis <- min(cap_depth_analysis, max(depths, na.rm = TRUE))

q05 <- as.integer(stats::quantile(depths, 0.05, na.rm = TRUE))
q10 <- as.integer(stats::quantile(depths, 0.10, na.rm = TRUE))
q20 <- as.integer(stats::quantile(depths, 0.20, na.rm = TRUE))

depth_grid_choice <- sort(unique(c(
  min_floor,
  1500L, 2000L, 2400L, 3000L, 4000L, 5000L, 8000L, 10000L, 15000L, 20000L,
  30000L, 40000L, 50000L, 60000L,
  q05, q10, q20,
  seq(2500L, max_depth_analysis, by = 2500L)
)))
depth_grid_choice <- depth_grid_choice[depth_grid_choice <= max_depth_analysis]

ret_tbl <- purrr::map_dfr(depth_grid_choice, function(d) {
  kept <- depth_df %>% filter(LibSize >= d)
  
  n_total <- nrow(depth_df)
  n_kept  <- nrow(kept)
  prop_kept <- if (n_total > 0) n_kept / n_total else NA_real_
  
  group_counts <- kept %>% count(Diet, Compartment, name = "n") %>% pull(n)
  min_group_n <- if (length(group_counts)) min(group_counts) else 0L
  
  tibble(
    rare_depth_candidate = as.integer(d),
    n_total = n_total,
    n_retained = n_kept,
    prop_retained = prop_kept,
    min_group_n = as.integer(min_group_n)
  )
})

safe_write_tsv(ret_tbl, file.path(tab_dir, "rarefaction_depth_candidates_HMHD.tsv"))

eligible <- ret_tbl %>% filter(prop_retained >= retain_prop_target, min_group_n >= min_per_group)

rare_depth <- if (nrow(eligible) > 0) {
  max(eligible$rare_depth_candidate)
} else {
  ret_tbl %>% arrange(desc(prop_retained), desc(rare_depth_candidate)) %>% slice(1) %>% pull(rare_depth_candidate)
}

rare_depth <- max(min_floor, min(cap_depth_analysis, as.integer(rare_depth)))

depth_summary <- tibble(
  n_hindgut_samples     = length(depths),
  min_depth             = as.integer(min(depths, na.rm = TRUE)),
  q05_depth             = q05,
  q10_depth             = q10,
  q20_depth             = q20,
  max_depth             = as.integer(max(depths, na.rm = TRUE)),
  retain_prop_target    = retain_prop_target,
  min_per_group         = min_per_group,
  min_floor             = min_floor,
  cap_depth_analysis    = cap_depth_analysis,
  rare_depth_used       = rare_depth,
  n_retained            = sum(depths >= rare_depth),
  prop_retained         = mean(depths >= rare_depth)
)
safe_write_tsv(depth_summary, file.path(tab_dir, "rarefaction_depth_choice_HMHD.tsv"))

below <- names(depths)[depths < rare_depth]
if (length(below)) {
  safe_write_tsv(tibble(Sample_ID = below, LibSize = as.integer(depths[below])),
                 file.path(tab_dir, "dropped_below_rarefaction_depth_HMHD.tsv"))
}

# ---- 8) Rarefy + alpha diversity -------------------------------------------
keep_analysis <- names(depths)[depths >= rare_depth]
hg2 <- subset_microtable_by_samples(hg, keep_analysis)

tn <- trans_norm$new(dataset = hg2)
hg_rarefy <- tn$norm(method = "rarefy", sample.size = rare_depth)
hg_rarefy$cal_alphadiv()

# ---- 8b) Post-rarefaction QC (NEW; recommended) ----------------------------
# 1) Verify all retained samples are exactly rare_depth reads
lib_post <- as.integer(colSums(hg_rarefy$otu_table))
names(lib_post) <- colnames(hg_rarefy$otu_table)

if (!all(lib_post == as.integer(rare_depth))) {
  bad <- names(lib_post)[lib_post != as.integer(rare_depth)]
  # Save mismatches for debugging
  safe_write_tsv(
    tibble(Sample_ID = bad, LibSize_post = lib_post[bad]),
    file.path(tab_dir, "QC_rarefaction_mismatched_library_sizes.tsv")
  )
  stop("❌ Rarefaction check failed: not all samples have ", rare_depth,
       " reads. See: ", file.path(tab_dir, "QC_rarefaction_mismatched_library_sizes.tsv"))
}

# 2) Save post-rarefaction Diet × Compartment counts (so you can report final n)
st_post <- as.data.frame(hg_rarefy$sample_table)
post_counts <- st_post %>%
  count(Diet, Compartment, name = "Count") %>%
  arrange(Diet, Compartment)

safe_write_tsv(post_counts, file.path(tab_dir, "QC_Diet_by_Compartment_after_rarefaction.tsv"))

# 3) Enforce min_per_group AFTER rarefaction (optional but consistent with your criterion)
min_post <- suppressWarnings(min(post_counts$Count, na.rm = TRUE))
if (is.finite(min_post) && min_post < min_per_group) {
  stop("❌ Post-rarefaction group size < min_per_group (", min_per_group, "). ",
       "Smallest group = ", min_post, ". ",
       "Lower rarefaction depth or adjust min_per_group. ",
       "See: ", file.path(tab_dir, "QC_Diet_by_Compartment_after_rarefaction.tsv"))
}

# Save post-rarefaction library sizes (should all be rare_depth)
safe_write_tsv(
  tibble(Sample_ID = names(lib_post), LibSize_post = lib_post),
  file.path(tab_dir, "QC_library_sizes_HMHD_post_rarefaction.tsv")
)

# ---- 8c) Export alpha table -------------------------------------------------
alpha_tbl <- hg_rarefy$alpha_diversity %>%
  rownames_to_column("Sample_ID") %>%
  left_join(
    as.data.frame(hg_rarefy$sample_table) %>%
      rownames_to_column("Sample_ID") %>%
      select(Sample_ID, Diet, Tank, Compartment, dplyr::any_of(c("FishUID","Fish_ID"))),
    by = "Sample_ID"
  ) %>%
  mutate(LibSize = as.integer(lib_post[Sample_ID]))

safe_write_tsv(alpha_tbl, file.path(output_dir, "AlphaDiversity_metrics_HMHD_rarefied.tsv"))
saveRDS(hg_rarefy, file.path(output_dir, "microeco_HMHD_rarefied_for_alpha.rds"))

message("✅ Step 04 complete")
message("• Rarefaction depth used: ", rare_depth)
message("• Output (microtable): ", normalizePath(file.path(output_dir, "microeco_HMHD_rarefied_for_alpha.rds")))
message("• Output (alpha table): ", normalizePath(file.path(output_dir, "AlphaDiversity_metrics_HMHD_rarefied.tsv")))
message("• Post-QC (counts): ", normalizePath(file.path(tab_dir, "QC_Diet_by_Compartment_after_rarefaction.tsv")))
message("• Post-QC (libsizes): ", normalizePath(file.path(tab_dir, "QC_library_sizes_HMHD_post_rarefaction.tsv")))
