#!/usr/bin/env Rscript
# =============================================================================
# Script : Step_05_alpha_diversity_UPGRADED_redMedian.R
# Purpose: Publication-safe α-diversity for HM/HD rarefied microeco object
#          - Descriptive α-plots (All / HD / HM) with fixed palette order
#          - Export α-values tables
#          - Publication-safe inference ONLY: paired HM vs HD (overall + by tank)
#          - Effect sizes (median/mean HM-HD) by tank
#          - NEW: draw an explicit RED median line on each boxplot (matches legend)
#
# Notes (reviewer-safe):
#   - 1 tank per diet => diet is confounded with tank -> NO diet-level inference here
#   - Optional KW/Dunn remains EXPLORATORY ONLY (disabled by default)
#
# Outputs:
#   results/Step_05_alpha_diversity/
#     All_HMHD/figures/alpha_*.png|pdf + alpha_composite_ABCD.png|pdf
#     All_HMHD/tables/alpha_values.tsv
#     HD_only/... (same)
#     HM_only/... (same)
#     inference_HM_vs_HD_byTank/tables/
#       HM_vs_HD_paired_Wilcoxon_byTank.tsv
#       HM_vs_HD_paired_Wilcoxon_overall.tsv
#       HM_minus_HD_effectsizes_byTank.tsv
#     sessionInfo.txt
#
# Author : Jayalal K Jayanthan
# =============================================================================

set.seed(123)

suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(stringr)
  library(tidyr)
})

# ---- I/O -------------------------------------------------------------------
input_rds   <- "results/Step_04_rarefaction_HMHD_for_alpha/microeco_HMHD_rarefied_for_alpha.rds"
output_dir  <- "results/Step_05_alpha_diversity"
fig_dir     <- file.path(output_dir, "figures")
tab_dir     <- file.path(output_dir, "tables")
infer_dir   <- file.path(output_dir, "inference_HM_vs_HD_byTank")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(infer_dir,   recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds))

mt <- readRDS(input_rds)
if (!inherits(mt, "microtable")) stop("❌ Input is not a microeco::microtable")

# ---- Config ----------------------------------------------------------------
alpha_measures <- c("Observed", "Coverage", "Shannon", "InvSimpson")
group_var <- "SampleTypes"

# IMPORTANT: with 1 tank per diet, diet-level inference is not replicated.
# Keep FALSE for publication-safe pipeline; TRUE only for exploratory plots/tables.
do_exploratory_kw_dunn <- FALSE

manual_colors <- c(
  "Tank_1_HD_Contrl_0%" = "#4D4D4D",
  "Tank_1_HM_Contrl_0%" = "#000000",
  "Tank_2_HD_Poro_1%"   = "#FFB74D",
  "Tank_2_HM_Poro_1%"   = "#E69F00",
  "Tank_3_HD_Poro_2_2%" = "#66C2A5",
  "Tank_3_HM_Poro_2_2%" = "#009E73"
)

keep_levels_all <- names(manual_colors)
keep_levels_hm  <- grep("_HM_", keep_levels_all, value = TRUE)
keep_levels_hd  <- grep("_HD_", keep_levels_all, value = TRUE)

# NEW: Median line style
median_line_colour <- "red"
median_line_shape  <- 95   # horizontal line
median_line_size   <- 10   # thickness/length (tune 6–12)

# ---- Helpers ----------------------------------------------------------------
ensure_dir <- function(p) { dir.create(p, recursive = TRUE, showWarnings = FALSE); invisible(p) }

safe_write_tsv <- function(df, file) {
  ensure_dir(dirname(file))
  readr::write_tsv(as.data.frame(df), file)
}

save_plot <- function(p, base, w = 8, h = 6, dpi = 600) {
  ensure_dir(dirname(base))
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, device = "pdf")
}

theme_pub <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position   = "bottom",
      legend.title      = element_blank(),
      axis.title.x      = element_blank(),
      plot.title        = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold")
    )
}

# Normalize Compartment labels to exactly HM/HD
normalize_compartment <- function(x) {
  xx <- tolower(trimws(as.character(x)))
  ifelse(grepl("muc|hm", xx), "HM",
         ifelse(grepl("dig|hd", xx), "HD", NA_character_))
}

make_sampletypes <- function(st) {
  req <- c("Tank", "Compartment", "Diet")
  miss <- setdiff(req, names(st))
  if (length(miss)) stop("❌ sample_table missing: ", paste(miss, collapse = ", "))
  
  diet_tag <- as.character(st$Diet)
  diet_tag <- gsub("^0% \\(Control\\)$", "Contrl_0%", diet_tag)
  diet_tag <- gsub("^1% Poro$",         "Poro_1%",   diet_tag)
  diet_tag <- gsub("^2\\.2% Poro$",     "Poro_2_2%", diet_tag)
  diet_tag <- gsub("[^A-Za-z0-9_%]+", "_", diet_tag)
  
  tank_tag <- gsub("[^A-Za-z0-9_]+", "_", as.character(st$Tank))
  comp_tag <- gsub("[^A-Za-z0-9_]+", "_", normalize_compartment(st$Compartment))
  
  paste0(tank_tag, "_", comp_tag, "_", diet_tag)
}

# Rowname-safe subset (with alignment check)
subset_mt <- function(m, keep_levels) {
  mm <- m$clone(deep = TRUE)
  
  st <- as.data.frame(mm$sample_table) %>%
    rownames_to_column("Sample_ID") %>%
    mutate(.grp = as.character(.data[[group_var]])) %>%
    filter(.grp %in% keep_levels) %>%
    mutate(!!group_var := factor(.grp, levels = keep_levels)) %>%
    select(-.grp) %>%
    column_to_rownames("Sample_ID")
  
  if (!nrow(st)) stop("❌ Subset has 0 samples")
  if (!all(rownames(st) %in% colnames(mm$otu_table))) {
    stop("❌ sample_table rownames do not match otu_table column names after subsetting.")
  }
  
  mm$sample_table <- st
  mm$otu_table    <- mm$otu_table[, rownames(st), drop = FALSE]
  mm$tidy_dataset()
  mm
}

# ---- Ensure metadata is safe ------------------------------------------------
mt$sample_table <- as.data.frame(mt$sample_table)
colnames(mt$sample_table) <- make.unique(colnames(mt$sample_table))

# Ensure group var exists
if (!group_var %in% names(mt$sample_table)) {
  mt$sample_table[[group_var]] <- make_sampletypes(mt$sample_table)
}
mt$sample_table[[group_var]] <- trimws(as.character(mt$sample_table[[group_var]]))
mt$sample_table[[group_var]] <- factor(mt$sample_table[[group_var]], levels = keep_levels_all)

# Normalize Compartment if present
if ("Compartment" %in% names(mt$sample_table)) {
  mt$sample_table$Compartment <- normalize_compartment(mt$sample_table$Compartment)
}

# ---- Plotting (descriptive by default) -------------------------------------
plot_alpha_set <- function(microtable_obj, label, keep_levels) {
  out  <- ensure_dir(file.path(output_dir, label))
  figs <- ensure_dir(file.path(out, "figures"))
  tabs <- ensure_dir(file.path(out, "tables"))
  
  m <- subset_mt(microtable_obj, keep_levels)
  
  # Compute ONLY requested measures (consistent & fast)
  m$cal_alphadiv(measures = alpha_measures, PD = FALSE)
  
  alpha_df <- as.data.frame(m$alpha_diversity) %>%
    rownames_to_column("Sample_ID") %>%
    left_join(
      as.data.frame(m$sample_table) %>%
        rownames_to_column("Sample_ID") %>%
        select(Sample_ID, all_of(group_var), any_of(c("Diet", "Tank", "Compartment", "FishUID", "Fish_ID"))),
      by = "Sample_ID"
    )
  
  safe_write_tsv(alpha_df, file.path(tabs, "alpha_values.tsv"))
  
  # Optional exploratory KW/Dunn (NOT recommended for 1 tank per diet)
  if (isTRUE(do_exploratory_kw_dunn)) {
    aobj <- trans_alpha$new(dataset = m, group = group_var)
    aobj$cal_diff(measure = alpha_measures, method = "KW_dunn",
                  p_adjust_method = "holm", KW_dunn_letter = TRUE)
    if (!is.null(aobj$res_diff)) {
      safe_write_tsv(aobj$res_diff, file.path(tabs, "alpha_stats_KW_Dunn_EXPLORATORY.tsv"))
    }
  }
  
  cols <- manual_colors[keep_levels]
  names(cols) <- keep_levels
  
  plots <- list()
  for (mname in alpha_measures) {
    dfp <- alpha_df %>%
      mutate(SampleTypes = factor(.data[[group_var]], levels = keep_levels)) %>%
      filter(!is.na(SampleTypes))
    
    g <- ggplot(dfp, aes(x = SampleTypes, y = .data[[mname]], fill = SampleTypes, color = SampleTypes)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.7) +
      # NEW: explicit red median line (matches legend)
      stat_summary(
        fun = median,
        geom = "point",
        shape = median_line_shape,
        size  = median_line_size,
        colour = median_line_colour,
        show.legend = FALSE
      ) +
      geom_point(position = position_jitter(width = 0.12), size = 1.6, alpha = 0.75) +
      scale_fill_manual(values = cols, drop = FALSE) +
      scale_color_manual(values = cols, drop = FALSE) +
      theme_pub(14) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      labs(y = mname)
    
    save_plot(g, file.path(figs, paste0("alpha_", mname)), w = 8, h = 6, dpi = 600)
    plots[[mname]] <- g
  }
  
  comp <- ggpubr::ggarrange(
    plotlist      = plots[alpha_measures],
    ncol          = 2, nrow = 2,
    labels        = LETTERS[seq_along(alpha_measures)],
    label.x       = 0.01, label.y = 0.99,
    hjust         = 0, vjust = 1,
    font.label    = list(size = 16, face = "bold"),
    align         = "hv",
    common.legend = TRUE,
    legend        = "bottom"
  )
  save_plot(comp, file.path(figs, "alpha_composite_ABCD"), w = 12, h = 10, dpi = 600)
  
  invisible(alpha_df)
}

# ---- Publication-safe inference: paired HM vs HD (overall + by tank) --------
run_hm_hd_paired <- function(microtable_obj) {
  idir <- ensure_dir(infer_dir)
  tabs <- ensure_dir(file.path(idir, "tables"))
  
  m <- subset_mt(microtable_obj, keep_levels_all)
  m$cal_alphadiv(measures = alpha_measures, PD = FALSE)
  
  st <- as.data.frame(m$sample_table) %>%
    rownames_to_column("Sample_ID")
  
  # Ensure Compartment exists
  if (!"Compartment" %in% names(st)) {
    st$Compartment <- ifelse(grepl("_HM_", st[[group_var]]), "HM", "HD")
  }
  st$Compartment <- normalize_compartment(st$Compartment)
  
  # Ensure Tank exists
  if (!"Tank" %in% names(st)) {
    st$Tank <- sub("^(Tank_[0-9]+)_.*$", "\\1", st[[group_var]])
  } else {
    st$Tank <- as.character(st$Tank)
  }
  
  # Ensure FishUID exists (pairing key)
  if (!"FishUID" %in% names(st) || all(is.na(st$FishUID))) {
    if ("Fish_ID" %in% names(st)) {
      st$FishUID <- paste(st$Tank, st$Fish_ID, sep = "__")
    } else {
      warning("FishUID/Fish_ID not found. Paired HM–HD tests will be skipped.")
      return(invisible(NULL))
    }
  }
  
  alpha <- as.data.frame(m$alpha_diversity) %>%
    rownames_to_column("Sample_ID") %>%
    left_join(st %>% select(Sample_ID, Tank, Compartment, FishUID), by = "Sample_ID") %>%
    mutate(Compartment = normalize_compartment(Compartment)) %>%
    filter(!is.na(Compartment))
  
  # (A) Paired HM vs HD WITHIN each tank
  by_tank <- lapply(alpha_measures, function(mname) {
    alpha %>%
      select(Tank, FishUID, Compartment, value = all_of(mname)) %>%
      pivot_wider(names_from = Compartment, values_from = value) %>%
      filter(!is.na(HM) & !is.na(HD)) %>%
      group_by(Tank) %>%
      summarize(
        Measure = mname,
        n_pairs = n(),
        median_diff_HM_minus_HD = median(HM - HD),
        mean_diff_HM_minus_HD   = mean(HM - HD),
        p_value = suppressWarnings(wilcox.test(HM, HD, paired = TRUE)$p.value),
        .groups = "drop"
      )
  }) %>% bind_rows() %>%
    group_by(Measure) %>%
    mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
    ungroup()
  
  safe_write_tsv(by_tank, file.path(tabs, "HM_vs_HD_paired_Wilcoxon_byTank.tsv"))
  
  # (B) Paired HM vs HD OVERALL (all tanks pooled; paired within fish)
  overall <- lapply(alpha_measures, function(mname) {
    df <- alpha %>%
      select(FishUID, Compartment, value = all_of(mname)) %>%
      pivot_wider(names_from = Compartment, values_from = value) %>%
      filter(!is.na(HM) & !is.na(HD))
    
    tibble(
      Measure = mname,
      n_pairs = nrow(df),
      median_diff_HM_minus_HD = median(df$HM - df$HD),
      mean_diff_HM_minus_HD   = mean(df$HM - df$HD),
      p_value = suppressWarnings(wilcox.test(df$HM, df$HD, paired = TRUE)$p.value)
    )
  }) %>% bind_rows() %>%
    mutate(FDR = p.adjust(p_value, method = "fdr"))
  
  safe_write_tsv(overall, file.path(tabs, "HM_vs_HD_paired_Wilcoxon_overall.tsv"))
  
  # (C) Effect size summary only (nice for supplementary)
  eff <- by_tank %>%
    select(Tank, Measure, n_pairs, median_diff_HM_minus_HD, mean_diff_HM_minus_HD)
  safe_write_tsv(eff, file.path(tabs, "HM_minus_HD_effectsizes_byTank.tsv"))
  
  invisible(list(by_tank = by_tank, overall = overall, effects = eff))
}

# ---- Run -------------------------------------------------------------------
plot_alpha_set(mt, "All_HMHD", keep_levels_all)
plot_alpha_set(mt, "HD_only",  keep_levels_hd)
plot_alpha_set(mt, "HM_only",  keep_levels_hm)

run_hm_hd_paired(mt)

# Reproducibility log
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message("✅ Step 05 complete")
message("• Outputs: ", normalizePath(output_dir))
