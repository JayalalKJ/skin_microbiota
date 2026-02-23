#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Script: Step_06_beta_diversity_Aitchison_ConQuR_Set6groups_NOLOCK_PUBLICATION.R
# Purpose:
#   β-diversity using Aitchison distance (Euclidean on CLR/rCLR) directly on the
#   ConQuR batch-corrected Set_6groups microeco object.
#   NO rarefaction lock / NO QC files.
#
# Output highlights:
#   - figures/PCoA_rCLR_PUBLICATION.(png/pdf)   main paper-ready PCoA (microeco)
#   - figures/PERMDISP_rCLR_distances.(png/pdf) optional dispersion plot
#   - tables/PERMANOVA_*.tsv, PERMDISP_*.tsv
#   - data/*.rds (transformed microtables + distance matrices)
# ─────────────────────────────────────────────────────────────────────────────

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(microeco)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(glue)
  library(readr)
  library(stringr)
  library(tidyr)
  library(vegan)
  library(tools)
  library(rlang)
  library(tibble)
})

# ───────────────────────────── Reproducibility ─────────────────────────────
SEED_NMDS <- 1234
set.seed(20260215)

# ───────────────────────────── I/O ──────────────────────────────────────────
input_rds   <- "results/Step_03D_Batchcorrection_ConQuR_2/Set_6groups/microeco_conqur_penalized_Genus_level_OTUs.rds"
output_root <- "results/Step_06_Beta_diversity"

# ───────────────────────────── Config ───────────────────────────────────────
# Use your real grouping column name in sample_table:
group_var <- "Sample_types"   # <- common in your pipeline; change to "SampleTypes" if needed

manual_colors <- c(
  "Tank_1_HD_Contrl_0%" = "#999999",
  "Tank_1_HM_Contrl_0%" = "#333333",
  "Tank_2_HD_Poro_1%"   = "#A6CEE3",
  "Tank_2_HM_Poro_1%"   = "#1F78B4",
  "Tank_3_HD_Poro_2_2%" = "#B2DF8A",
  "Tank_3_HM_Poro_2_2%" = "#33A02C"
)

keep_types_HM  <- c("Tank_1_HM_Contrl_0%", "Tank_2_HM_Poro_1%", "Tank_3_HM_Poro_2_2%")
keep_types_HD  <- c("Tank_1_HD_Contrl_0%", "Tank_2_HD_Poro_1%", "Tank_3_HD_Poro_2_2%")
keep_types_all <- c(keep_types_HM, keep_types_HD)

MALACO_PATTERNS <- c("(?i)\\b(?:candidatus\\s+)?malacoplasma\\b", "(?i)\\bmycoplasma\\b")
SILENCE_FAMILY  <- TRUE
FAMILY_PATTERNS <- c("(?i)\\bmycoplasmataceae\\b")

MAKE_FULL_CLR_PLOTS <- TRUE
N_PERM <- 999
PSEUDOCOUNT <- 0.5

# ───────────────────────────── Helpers ──────────────────────────────────────
ensure_dir <- function(p) { dir.create(p, recursive = TRUE, showWarnings = FALSE); invisible(p) }
trim_chr   <- function(x) trimws(as.character(x))

save_plot <- function(p, stub, w = 8, h = 6, dpi = 600) {
  ggsave(paste0(stub, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(paste0(stub, ".pdf"), p, width = w, height = h, device = "pdf")
}

save_tsv <- function(df, file) {
  if (is.null(df)) return(invisible(NULL))
  ensure_dir(dirname(file))
  readr::write_tsv(as.data.frame(df), file)
}

load_microeco_object <- function(path){
  if (!file.exists(path)) stop("❌ File not found: ", path)
  obj <- readRDS(path)
  if (!inherits(obj, "microtable"))
    stop("❌ .rds does not contain a microeco::microtable. Got: ", paste(class(obj), collapse = ", "))
  obj
}

preserve_col <- function(df, col, new_base){
  if (!col %in% colnames(df)) return(df)
  new_name <- new_base
  if (new_name %in% colnames(df)) {
    i <- 1L
    while (paste0(new_base, "_", i) %in% colnames(df)) i <- i + 1L
    new_name <- paste0(new_base, "_", i)
  }
  df[[new_name]] <- df[[col]]
  df
}

derive_from_sampletypes <- function(st, group_var = "Sample_types") {
  if (!group_var %in% colnames(st)) return(st)
  x <- as.character(st[[group_var]])
  
  if (!"Tank" %in% colnames(st)) st$Tank <- str_extract(x, "Tank_[0-9]+")
  if (!"Compartment" %in% colnames(st)) {
    st$Compartment <- ifelse(str_detect(x, "_HM_"), "HM",
                             ifelse(str_detect(x, "_HD_"), "HD", NA_character_))
  }
  if (!"Diet" %in% colnames(st)) {
    st$Diet <- dplyr::case_when(
      str_detect(x, "Contrl_0%") ~ "0% (Control)",
      str_detect(x, "Poro_1%")   ~ "1% Poro",
      str_detect(x, "Poro_2_2%") ~ "2.2% Poro",
      TRUE ~ NA_character_
    )
    st$Diet <- factor(st$Diet, levels = c("0% (Control)", "1% Poro", "2.2% Poro"))
  }
  st
}

sanitize_microtable <- function(mt, group_var = "Sample_types") {
  m <- mt$clone(deep = TRUE)
  
  st <- as.data.frame(m$sample_table)
  colnames(st) <- make.unique(colnames(st))
  otu_ids <- colnames(m$otu_table)
  
  # Align rownames(sample_table) to colnames(otu_table)
  if (!identical(rownames(st), otu_ids)) {
    common <- intersect(rownames(st), otu_ids)
    if (!length(common)) {
      id_candidates <- intersect(c("Sample_ID","Sample_id","SampleID","Sample","Sample2"), colnames(st))
      if (!length(id_candidates)) stop("❌ Cannot align sample_table to otu_table (no overlap + no ID columns).")
      score <- vapply(id_candidates, function(cc) mean(trim_chr(st[[cc]]) %in% otu_ids), numeric(1))
      id_col <- names(which.max(score))
      st[[id_col]] <- trim_chr(st[[id_col]])
      rownames(st) <- st[[id_col]]
      common <- intersect(rownames(st), otu_ids)
    }
    if (!length(common)) stop("❌ Cannot align sample_table and otu_table after ID repair.")
    common <- sort(common)
    m$otu_table    <- m$otu_table[, common, drop = FALSE]
    m$sample_table <- st[common, , drop = FALSE]
  }
  
  st <- as.data.frame(m$sample_table)
  colnames(st) <- make.unique(colnames(st))
  
  st <- preserve_col(st, "Sample_ID", "Sample_ID_in_table")
  st$Sample_ID <- rownames(st)
  
  # Ensure group_var exists (try common alternatives)
  if (!group_var %in% colnames(st)) {
    alt <- intersect(c("SampleTypes","Sample_types","Sample_Type","Groups","Group"), colnames(st))
    if (length(alt)) st[[group_var]] <- st[[alt[1]]]
  }
  
  st <- derive_from_sampletypes(st, group_var = group_var)
  st[[group_var]] <- trim_chr(st[[group_var]])
  
  m$sample_table <- st
  m$tidy_dataset()
  m
}

drop_allzero_otus <- function(mt_in) {
  m <- mt_in$clone(deep = TRUE)
  keep <- rowSums(m$otu_table) > 0
  m$otu_table <- m$otu_table[keep, , drop = FALSE]
  if (!is.null(m$tax_table) && nrow(m$tax_table)) {
    common <- intersect(rownames(m$tax_table), rownames(m$otu_table))
    m$tax_table <- m$tax_table[common, , drop = FALSE]
    m$otu_table <- m$otu_table[common, , drop = FALSE]
  }
  m$tidy_dataset()
  m
}

rank_groups <- list(Family = c("Family","Families"), Genus = c("Genus","Genera"))
pick_rank_label <- function(group_key, available_cols){
  opts <- rank_groups[[group_key]]
  hit  <- opts[opts %in% available_cols]
  if (length(hit)) hit[1] else NA_character_
}
normalize_tax <- function(x){
  x %>%
    as.character() %>%
    str_replace_all("^.?__", "") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("[\\[\\]\\(\\)\"'`]", "") %>%
    str_squish()
}

silence_taxa <- function(mt_in, genus_patterns, use_family = FALSE, family_patterns = character()){
  mt_new <- mt_in$clone(deep = TRUE)
  tt <- mt_new$tax_table
  if (is.null(tt) || !nrow(tt)) return(mt_new)
  
  genus_col  <- pick_rank_label("Genus",  colnames(tt))
  family_col <- pick_rank_label("Family", colnames(tt))
  if (is.na(genus_col)) return(mt_new)
  
  g_txt <- normalize_tax(tt[[genus_col]])
  g_hit <- rep(FALSE, length(g_txt))
  for (pat in genus_patterns) g_hit <- g_hit | str_detect(g_txt, regex(pat))
  g_hit[is.na(g_hit)] <- FALSE
  
  f_hit <- rep(FALSE, nrow(tt))
  if (use_family && !is.na(family_col) && length(family_patterns)) {
    f_txt <- normalize_tax(tt[[family_col]])
    for (pat in family_patterns) f_hit <- f_hit | str_detect(f_txt, regex(pat))
    f_hit[is.na(f_hit)] <- FALSE
  }
  
  hit <- g_hit | f_hit
  otu_ids <- intersect(rownames(tt)[hit], rownames(mt_new$otu_table))
  if (!length(otu_ids)) return(mt_new)
  
  mt_new$otu_table[otu_ids, ] <- 0L
  drop_allzero_otus(mt_new)
}

subset_by_levels <- function(mt_in, group_var, keep_levels) {
  m <- mt_in$clone(deep = TRUE)
  st <- as.data.frame(m$sample_table)
  colnames(st) <- make.unique(colnames(st))
  
  if (!group_var %in% colnames(st)) stop("❌ Missing group_var: ", group_var)
  
  st[[group_var]] <- trim_chr(st[[group_var]])
  st <- st[st[[group_var]] %in% keep_levels, , drop = FALSE]
  if (!nrow(st)) stop("❌ No samples after filtering keep_levels: ", paste(keep_levels, collapse = ", "))
  
  st[[group_var]] <- factor(st[[group_var]], levels = keep_levels)
  
  keep_ids <- intersect(rownames(st), colnames(m$otu_table))
  st <- st[keep_ids, , drop = FALSE]
  
  m$otu_table    <- m$otu_table[, keep_ids, drop = FALSE]
  m$sample_table <- st
  m$tidy_dataset()
  m
}

make_clr_objects <- function(mt_in) {
  tn <- trans_norm$new(dataset = mt_in)
  list(
    CLR  = tn$norm(method = "CLR"),
    rCLR = tn$norm(method = "rclr")
  )
}

aitchison_distance <- function(mt_transformed) {
  X <- t(as.matrix(mt_transformed$otu_table))
  vegan::vegdist(X, method = "euclidean")
}

run_permanova <- function(D, meta_df, group_var, permutations = 999) {
  labs <- attr(D, "Labels")
  meta_df <- meta_df[labs, , drop = FALSE]
  dat <- data.frame(group = droplevels(meta_df[[group_var]]))
  as.data.frame(vegan::adonis2(D ~ group, data = dat, permutations = permutations))
}

run_permdisp <- function(D, meta_df, group_var, permutations = 999) {
  labs <- attr(D, "Labels")
  meta_df <- meta_df[labs, , drop = FALSE]
  
  grp <- droplevels(meta_df[[group_var]])
  bd  <- vegan::betadisper(D, group = grp)
  pt  <- vegan::permutest(bd, permutations = permutations)
  
  tuk <- stats::TukeyHSD(bd)$group
  tuk_df <- as.data.frame(tuk) %>% tibble::rownames_to_column("Comparison")
  
  list(
    bd = bd,
    tab = as.data.frame(pt$tab),
    tukey = tuk_df
  )
}

annot_line <- function(perma_df) {
  if (is.null(perma_df) || !nrow(perma_df)) return("PERMANOVA: NA")
  idx <- which(rownames(perma_df) == "group")
  if (!length(idx)) idx <- 1
  glue("PERMANOVA: R²={sprintf('%.3f', perma_df$R2[idx])}  F={sprintf('%.2f', perma_df$F[idx])}  p={format.pval(perma_df$`Pr(>F)`[idx], digits = 2)}")
}

permdisp_line <- function(disp_tab) {
  if (is.null(disp_tab) || !nrow(disp_tab)) return("PERMDISP: NA")
  idx <- which(rownames(disp_tab) == "Groups")
  if (!length(idx)) idx <- 1
  p <- disp_tab[idx, "Pr(>F)"]
  glue("PERMDISP: p={format.pval(p, digits = 2)}")
}

microeco_pcoa_plot <- function(m_sub, D_dist, group_var, colour_vec,
                               title = "PCoA (Aitchison from rCLR)",
                               annot_text = NULL) {
  
  Dm <- as.matrix(D_dist)
  stopifnot(identical(rownames(Dm), colnames(Dm)))
  
  labs <- rownames(Dm)
  if (!all(labs %in% rownames(m_sub$sample_table))) {
    stop("❌ Distance labels do not match sample_table rownames.")
  }
  
  tb <- trans_beta$new(
    dataset = m_sub,
    group   = group_var,
    measure = Dm
  )
  tb$cal_ordination(method = "PCoA", ncomp = 2)
  
  p <- tb$plot_ordination(
    plot_color = group_var,
    plot_type  = c("point", "ellipse")
  )
  
  fit <- cmdscale(as.dist(Dm), k = 2, eig = TRUE)
  eig <- fit$eig
  eig_pos <- eig[eig > 0]
  var1 <- if (length(eig_pos)) 100 * eig_pos[1] / sum(eig_pos) else NA_real_
  var2 <- if (length(eig_pos) > 1) 100 * eig_pos[2] / sum(eig_pos) else NA_real_
  
  p <- p +
    scale_color_manual(values = colour_vec, breaks = names(colour_vec), drop = FALSE) +
    scale_fill_manual(values  = colour_vec, breaks = names(colour_vec), drop = FALSE) +
    labs(
      title = title,
      x = if (is.finite(var1)) sprintf("PCoA1 (%.1f%%)", var1) else "PCoA1",
      y = if (is.finite(var2)) sprintf("PCoA2 (%.1f%%)", var2) else "PCoA2"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  if (!is.null(annot_text) && nzchar(annot_text)) {
    p <- p + annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.15,
                      label = annot_text, size = 4, fontface = "bold") +
      coord_cartesian(clip = "off")
  }
  
  list(p = p, tb = tb)
}

plot_dispersion_distances <- function(disp_obj, meta_df, group_var, colour_vec, title) {
  dd <- data.frame(
    Sample_ID = names(disp_obj$bd$distances),
    Distance  = as.numeric(disp_obj$bd$distances),
    stringsAsFactors = FALSE
  ) %>% left_join(meta_df, by = "Sample_ID")
  
  dd[[group_var]] <- factor(trim_chr(dd[[group_var]]), levels = names(colour_vec))
  
  ggplot(dd, aes(x = .data[[group_var]], y = Distance, color = .data[[group_var]])) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 1.6) +
    scale_color_manual(values = colour_vec, breaks = names(colour_vec), drop = FALSE) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold")
    ) +
    labs(title = title, x = NULL, y = "Distance to group centroid")
}

# ───────────────────────────── STEP 1: Load + sanitize ──────────────────────
ensure_dir(output_root)
mt <- load_microeco_object(input_rds)
mt <- sanitize_microtable(mt, group_var = group_var)

# Defensive check: are all 6 groups present?
st0 <- as.data.frame(mt$sample_table)
if (!group_var %in% colnames(st0)) stop("❌ group_var not found in sample_table: ", group_var)
present_levels <- sort(unique(trim_chr(st0[[group_var]])))
message("✅ group_var levels present: ", paste(present_levels, collapse = " | "))

# ───────────────────────────── Core runner ──────────────────────────────────
run_beta_aitchison <- function(microtable_obj, out_dir, keep_levels) {
  message(glue("▶ β-diversity Aitchison: {out_dir}"))
  figs <- ensure_dir(file.path(out_dir, "figures"))
  tabs <- ensure_dir(file.path(out_dir, "tables"))
  dats <- ensure_dir(file.path(out_dir, "data"))
  
  # subset
  m_sub <- subset_by_levels(microtable_obj, group_var, keep_levels)
  st <- as.data.frame(m_sub$sample_table)
  colnames(st) <- make.unique(colnames(st))
  
  meta <- st
  meta <- preserve_col(meta, "Sample_ID", "Sample_ID_in_table")
  meta$Sample_ID <- rownames(meta)
  meta[[group_var]] <- factor(trim_chr(meta[[group_var]]), levels = keep_levels)
  colnames(meta) <- make.unique(colnames(meta))
  
  # colors
  colour_vec <- manual_colors[keep_levels]
  if (any(is.na(colour_vec))) colour_vec[is.na(colour_vec)] <- "#999999"
  names(colour_vec) <- keep_levels
  
  # CLR/rCLR
  tr <- make_clr_objects(m_sub)
  CLR  <- sanitize_microtable(tr$CLR,  group_var = group_var)
  rCLR <- sanitize_microtable(tr$rCLR, group_var = group_var)
  
  saveRDS(CLR,  file.path(dats, "microeco_CLR.rds"))
  saveRDS(rCLR, file.path(dats, "microeco_rCLR.rds"))
  
  # distances
  D_clr  <- aitchison_distance(CLR)
  D_rclr <- aitchison_distance(rCLR)
  saveRDS(D_clr,  file.path(dats, "Distance_Aitchison_from_CLR.rds"))
  saveRDS(D_rclr, file.path(dats, "Distance_Aitchison_from_rCLR.rds"))
  
  # PERMANOVA
  perm_clr  <- run_permanova(D_clr,  meta, group_var, permutations = N_PERM)
  perm_rclr <- run_permanova(D_rclr, meta, group_var, permutations = N_PERM)
  save_tsv(perm_clr,  file.path(tabs, "PERMANOVA_Aitchison_CLR.tsv"))
  save_tsv(perm_rclr, file.path(tabs, "PERMANOVA_Aitchison_rCLR.tsv"))
  
  # PERMDISP
  disp_clr  <- run_permdisp(D_clr,  meta, group_var, permutations = N_PERM)
  disp_rclr <- run_permdisp(D_rclr, meta, group_var, permutations = N_PERM)
  save_tsv(disp_clr$tab,    file.path(tabs, "PERMDISP_Aitchison_CLR.tsv"))
  save_tsv(disp_rclr$tab,   file.path(tabs, "PERMDISP_Aitchison_rCLR.tsv"))
  save_tsv(disp_clr$tukey,  file.path(tabs, "PERMDISP_TukeyHSD_CLR.tsv"))
  save_tsv(disp_rclr$tukey, file.path(tabs, "PERMDISP_TukeyHSD_rCLR.tsv"))
  
  # Publication plot (rCLR)
  ann_pub <- paste(annot_line(perm_rclr), permdisp_line(disp_rclr$tab), sep = "\n")
  pub <- microeco_pcoa_plot(
    m_sub      = rCLR,
    D_dist     = D_rclr,
    group_var  = group_var,
    colour_vec = colour_vec,
    title      = "PCoA (Aitchison from rCLR)",
    annot_text = ann_pub
  )
  save_plot(pub$p, file.path(figs, "PCoA_rCLR_PUBLICATION"), w = 8, h = 6, dpi = 600)
  
  p_disp <- plot_dispersion_distances(
    disp_rclr, meta, group_var, colour_vec,
    title = "PERMDISP (rCLR Aitchison): distance to centroid"
  )
  save_plot(p_disp, file.path(figs, "PERMDISP_rCLR_distances"), w = 9, h = 5.5, dpi = 600)
  
  # Optional extra plots
  if (isTRUE(MAKE_FULL_CLR_PLOTS)) {
    p_clr <- microeco_pcoa_plot(
      m_sub      = CLR,
      D_dist     = D_clr,
      group_var  = group_var,
      colour_vec = colour_vec,
      title      = "PCoA (Aitchison from CLR)",
      annot_text = annot_line(perm_clr)
    )
    save_plot(p_clr$p, file.path(figs, "PCoA_CLR"), w = 8, h = 6, dpi = 600)
  }
  
  # Samples used
  keep_cols <- intersect(c("Sample_ID", group_var, "Diet", "Tank", "Compartment"), colnames(meta))
  save_tsv(meta[, keep_cols, drop = FALSE], file.path(tabs, "Samples_used.tsv"))
  
  invisible(TRUE)
}

# ───────────────────────────── Branches ─────────────────────────────────────
branches <- list(
  With_Malacoplasma = mt,
  Without_Malacoplasma = silence_taxa(mt, MALACO_PATTERNS,
                                      use_family = SILENCE_FAMILY,
                                      family_patterns = FAMILY_PATTERNS)
)

# ───────────────────────────── Run All / HM / HD ────────────────────────────
for (bname in names(branches)) {
  bmt   <- sanitize_microtable(branches[[bname]], group_var = group_var)
  broot <- ensure_dir(file.path(output_root, bname))
  
  run_beta_aitchison(bmt, file.path(broot, "All"), keep_types_all)
  run_beta_aitchison(bmt, file.path(broot, "HM"),  keep_types_HM)
  run_beta_aitchison(bmt, file.path(broot, "HD"),  keep_types_HD)
}

message(glue("✅ DONE: {normalizePath(output_root)}"))

# Main-paper recommended figures:
#   results/Step_06_Beta_diversity_ConQuR_Set6groups_NOLOCK/With_Malacoplasma/HM/figures/PCoA_rCLR_PUBLICATION.png
#   results/Step_06_Beta_diversity_ConQuR_Set6groups_NOLOCK/With_Malacoplasma/HD/figures/PCoA_rCLR_PUBLICATION.png
