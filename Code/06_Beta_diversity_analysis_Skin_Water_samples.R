#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# Script: Step_07_Beta_diversity_HMHD_ENV_feed_PUBLICATION_Sample_types.R
# Purpose: β-diversity (With/Without Malacoplasma) for:
#          - All (HM + HD + Water + Pellets)
#          - HM only
#          - HD only
#          - ENV only (Water + Pellets)
#          using unified colors, clean legends, and A/B/C/D panel labels.
#
# Notes:
#  - Uses sample_table column: Sample_types  ✅ (as you requested)
#  - Silences genus Malacoplasma/Mycoplasma (and optionally Mycoplasmataceae)
#  - Outputs .png + .pdf for each ordination + composites
# ─────────────────────────────────────────────────────────────

options(stringsAsFactors = FALSE)

## 0 │ Libraries
suppressPackageStartupMessages({
  library(microeco)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(glue)
  library(readr)
  library(tidyr)
  library(stringr)
  library(grid)
  library(tools)
})

## 0.1 │ Reproducibility
SEED_NMDS <- 1234
set.seed(20260215)

## 1 │ I/O
input_rdata  <- "results/Step_03D_Batchcorrection_ConQuR_2/Set_6groups/microeco_conqur_penalized_Genus_level_OTUs.rds"
output_root  <- "results/Step_07_Beta_diversity_HMHD_ENv_feed"

## 2 │ Config
group_var       <- "Sample_types"           # ✅ USE THIS COLUMN
beta_measures   <- c("bray", "jaccard")
use_shapes      <- FALSE

# Manual colors for each Sample_types level
manual_colors <- c(
  # Control (neutral greys)
  "Tank_1_HD_Contrl_0%" = "#999999",
  "Tank_1_HM_Contrl_0%" = "#333333",
  
  # 1% inclusion – Blues
  "Tank_2_HD_Poro_1%"   = "#A6CEE3",
  "Tank_2_HM_Poro_1%"   = "#1F78B4",
  
  # 2.2% inclusion – Greens
  "Tank_3_HD_Poro_2_2%" = "#B2DF8A",
  "Tank_3_HM_Poro_2_2%" = "#33A02C",
  
  # ENV: Water
  "Tank_1_Water"        = "#BDBDBD",
  "Tank_2_Water"        = "#6BAED6",
  "Tank_3_Water"        = "#74C476",
  
  # ENV: Pellets
  "Tank_1_Pellets"      = "#636363",
  "Tank_2_Pellets"      = "#08519C",
  "Tank_3_Pellets"      = "#238B45"
)

# Subsets
keep_types_HM  <- c("Tank_1_HM_Contrl_0%", "Tank_2_HM_Poro_1%", "Tank_3_HM_Poro_2_2%")
keep_types_HD  <- c("Tank_1_HD_Contrl_0%", "Tank_2_HD_Poro_1%", "Tank_3_HD_Poro_2_2%")
keep_types_env <- c("Tank_1_Water","Tank_2_Water","Tank_3_Water",
                    "Tank_1_Pellets","Tank_2_Pellets","Tank_3_Pellets")
keep_types_all <- c(keep_types_HM, keep_types_HD, keep_types_env)

# Silencing patterns (genus-level), optional family-level safety net
MALACO_PATTERNS <- c(
  "(?i)\\b(?:candidatus\\s+)?malacoplasma\\b",
  "(?i)\\bmycoplasma\\b"
)
SILENCE_FAMILY   <- TRUE
FAMILY_PATTERNS  <- c("(?i)\\bmycoplasmataceae\\b")

## 3 │ Helpers ---------------------------------------------------------------

load_microeco_object <- function(path){
  if (!file.exists(path)) stop("File not found: ", path)
  ext <- tolower(file_ext(path))
  if (ext == "rds") {
    obj <- readRDS(path)
    if (!inherits(obj, "microtable"))
      stop("The .rds file does not contain a 'microtable'. Got: ", paste(class(obj), collapse=", "))
    return(obj)
  }
  if (ext %in% c("rdata","rda")) {
    env <- new.env(parent = emptyenv())
    loaded <- load(path, envir = env)
    for (nm in loaded) if (inherits(env[[nm]], "microtable")) return(env[[nm]])
    stop("No 'microtable' object found in .RData. Objects: ", paste(loaded, collapse=", "))
  }
  stop("Unsupported file extension: ", ext)
}

rank_groups <- list(
  Family  = c("Family","Families"),
  Genus   = c("Genus","Genera")
)
pick_rank_label <- function(group_key, available_cols){
  opts <- rank_groups[[group_key]]
  hit <- opts[opts %in% available_cols]
  if (length(hit)) return(hit[1])
  NA_character_
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
  if (is.na(genus_col)) {
    message("No Genus/Genera column in tax_table. Skipping silencing.")
    return(mt_new)
  }
  
  g_txt <- normalize_tax(tt[[genus_col]])
  g_hit <- rep(FALSE, length(g_txt))
  for (pat in genus_patterns) g_hit <- g_hit | str_detect(g_txt, regex(pat))
  g_hit <- replace_na(g_hit, FALSE)
  
  f_hit <- rep(FALSE, nrow(tt))
  if (use_family && !is.na(family_col) && length(family_patterns)) {
    f_txt <- normalize_tax(tt[[family_col]])
    for (pat in family_patterns) f_hit <- f_hit | str_detect(f_txt, regex(pat))
    f_hit <- replace_na(f_hit, FALSE)
  }
  
  hit <- g_hit | f_hit
  otu_ids <- intersect(rownames(tt)[hit], rownames(mt_new$otu_table))
  if (!length(otu_ids)) {
    message("• No OTUs matched Malacoplasma/Mycoplasma patterns. Nothing silenced.")
    return(mt_new)
  }
  
  mt_new$otu_table[otu_ids, ] <- 0L
  mt_new$tidy_dataset()
  mt_new
}

scale_legends_consistent <- function(colour_vec, breaks, title = "Sample types", include_shape = FALSE) {
  sc <- list(
    scale_color_manual(values = colour_vec, breaks = breaks, name = title, drop = FALSE),
    scale_fill_manual(values  = colour_vec, breaks = breaks, name = title, drop = FALSE),
    guides(
      color = guide_legend(order = 1),
      fill  = guide_legend(order = 1, override.aes = list(alpha = 0.5))
    )
  )
  if (!include_shape) sc <- c(sc, list(guides(shape = "none")))
  sc
}

save_plot <- function(p, stub, w=8, h=6, dpi=600) {
  ggsave(paste0(stub, ".png"), p, width=w, height=h, dpi=dpi, bg="white")
  ggsave(paste0(stub, ".pdf"), p, width=w, height=h, device = "pdf")
}
save_tsv <- function(df, file) {
  if (!is.null(df)) suppressWarnings(readr::write_tsv(as.data.frame(df), file))
}

## 4 │ Load data
mt <- load_microeco_object(input_rdata)

# ✅ check column exists
st0 <- mt$sample_table
if (!group_var %in% colnames(st0)) {
  stop(glue("❌ '{group_var}' not found in sample_table. Available columns: {paste(colnames(st0), collapse=', ')}"))
}

# Trim + show levels (debug-friendly)
mt$sample_table[[group_var]] <- trimws(as.character(mt$sample_table[[group_var]]))
message("Levels present in Sample_types: ",
        paste(sort(unique(mt$sample_table[[group_var]])), collapse = " | "))

## 5 │ Core runner
run_beta_export <- function(microtable_obj, out_dir, keep_levels, use_shape = use_shapes) {
  message(glue("▶ Running β-diversity export for: {out_dir}"))
  figs_dir <- file.path(out_dir, "figures"); dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)
  tabs_dir <- file.path(out_dir, "tables");  dir.create(tabs_dir,  recursive = TRUE, showWarnings = FALSE)
  
  mtab <- microtable_obj$clone(deep = TRUE)
  
  # subset & lock order
  st <- as.data.frame(mtab$sample_table)
  st[[group_var]] <- trimws(as.character(st[[group_var]]))
  st <- st |> dplyr::filter(.data[[group_var]] %in% keep_levels)
  if (!nrow(st)) stop(glue("❌ No samples after filtering to: {paste(keep_levels, collapse=', ')}"))
  
  # lock factor levels in your requested order
  st[[group_var]] <- factor(st[[group_var]], levels = keep_levels)
  
  # apply subset to microtable
  keep_ids <- rownames(st)
  mtab$sample_table <- st
  mtab$otu_table    <- mtab$otu_table[, keep_ids, drop = FALSE]
  mtab$tidy_dataset()
  
  present <- levels(st[[group_var]])
  colour_vec <- manual_colors[present]
  if (any(is.na(colour_vec))) {
    miss <- present[is.na(colour_vec)]
    warning(glue("⚠️ Missing colors for: {paste(miss, collapse=', ')} — using grey"))
    colour_vec[is.na(colour_vec)] <- "#999999"
  }
  names(colour_vec) <- present
  
  ord_plots <- list()
  
  for (meas in beta_measures) {
    mtab$cal_betadiv(measure = meas)
    tb <- trans_beta$new(dataset = mtab, measure = meas, group = group_var)
    
    # PERMANOVA
    tb$cal_manova(permutations = 999, manova_all = TRUE, na.action = stats::na.omit)
    perm <- tb$res_manova
    save_tsv(perm, file.path(tabs_dir, glue("PERMANOVA_{meas}.tsv")))
    
    annot <- if (!is.null(perm) && nrow(perm) >= 1) {
      glue("PERMANOVA: R²={sprintf('%.3f', perm$R2[1])}  F={sprintf('%.2f', perm$F[1])}  p={format.pval(perm$`Pr(>F)`[1], 2)}")
    } else {
      "PERMANOVA: NA"
    }
    
    make_ordination_plot <- function(method) {
      if (identical(method, "NMDS")) {
        set.seed(SEED_NMDS)
        tb$cal_ordination(method = "NMDS", ncomp = 2, NMDS_matrix = TRUE, trymax = 200)
      } else {
        tb$cal_ordination(method = method, ncomp = 2)
      }
      
      p <- tb$plot_ordination(
        plot_type               = c("point", "chull", "centroid"),
        plot_color              = group_var,
        plot_shape              = if (use_shape) group_var else NULL,
        color_values            = colour_vec,
        NMDS_stress_pos         = NULL,
        NMDS_stress_text_prefix = ""
      ) +
        annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.2,
                 label = annot, size = 4, fontface = "bold") +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 14) +
        theme(plot.margin = margin(5.5, 35, 5.5, 5.5)) +
        labs(title = sprintf("%s (%s)", toupper(method), toupper(meas)))
      
      p + scale_legends_consistent(colour_vec, present, title = "Sample types", include_shape = use_shape)
    }
    
    # Build & save (PCoA then NMDS)
    p_pcoa <- make_ordination_plot("PCoA")
    p_nmds <- make_ordination_plot("NMDS")
    
    save_plot(p_pcoa, file.path(figs_dir, glue("PCoA_{meas}")))
    save_plot(p_nmds, file.path(figs_dir, glue("NMDS_{meas}")))
    
    ord_plots[[paste0("PCoA_", meas)]] <- p_pcoa
    ord_plots[[paste0("NMDS_", meas)]] <- p_nmds
  }
  
  # Composite (A/B/C/D)
  lab_vec <- LETTERS[seq_len(length(ord_plots))]
  comp <- ggpubr::ggarrange(
    plotlist      = ord_plots,
    labels        = lab_vec,
    label.x       = 0.01,
    label.y       = 0.99,
    hjust         = 0,
    vjust         = 1,
    font.label    = list(size = 16, face = "bold"),
    ncol          = 2,
    nrow          = ceiling(length(ord_plots)/2),
    align         = "hv",
    common.legend = TRUE,
    legend        = "bottom"
  )
  save_plot(comp, file.path(figs_dir, "beta_composite"), w = 12, h = 10, dpi = 600)
  
  # Save samples used (nice for supplement)
  samples_used <- data.frame(Sample_ID = rownames(st), Sample_types = st[[group_var]])
  save_tsv(samples_used, file.path(tabs_dir, "Samples_used.tsv"))
}

## 6 │ Branches: With vs Without Malacoplasma
mt_with    <- mt
mt_without <- silence_taxa(mt, MALACO_PATTERNS, use_family = SILENCE_FAMILY, family_patterns = FAMILY_PATTERNS)

branches <- list(
  With_Malacoplasma    = mt_with,
  Without_Malacoplasma = mt_without
)

## 7 │ Run analyses (separate folders per branch + subset)
for (bname in names(branches)) {
  bmt   <- branches[[bname]]
  broot <- file.path(output_root, bname)
  
  run_beta_export(bmt, file.path(broot, "All"), keep_types_all, use_shape = use_shapes)
  run_beta_export(bmt, file.path(broot, "HM"),  keep_types_HM,  use_shape = use_shapes)
  run_beta_export(bmt, file.path(broot, "HD"),  keep_types_HD,  use_shape = use_shapes)
  
  # ✅ NEW: ENV only (Water + Pellets)
  run_beta_export(bmt, file.path(broot, "ENV"), keep_types_env, use_shape = use_shapes)
}

message(glue("✅ β-diversity complete: {normalizePath(output_root)}"))
