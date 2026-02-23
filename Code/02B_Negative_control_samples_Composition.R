#!/usr/bin/env Rscript

# =============================================================================
# Script: Step_4_Negative_Controls_Composition_Plots_AllFormats.R
# Author: Jayalal K Jayanthan (extended outputs; FIXED htmlwidgets/plotly export)
# Date  : 2025-08-22
# Goal  : 1) Load microeco object
#         2) Subset negative controls (Lab1, Lab2, Both_labs)
#         3) Produce nested bar plots (Genera_level_OTUs by Classes)
#         4) Save in PNG, PDF, SVG, EPS, TIFF, JPEG + interactive HTML
# Notes : Fixes htmlwidgets::saveWidget errors on Windows by:
#         - using a lab-specific libdir inside out_dir when pandoc not available
#         - falling back to selfcontained=FALSE when needed
# =============================================================================

set.seed(123)

# ---- 0) Packages ------------------------------------------------------------
suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(cowplot)
  library(svglite)
  library(plotly)
  library(htmlwidgets)
})

# ---- 1) I/O -----------------------------------------------------------------
rds_in   <- "results/Step_02_remove_spikein_taxa/microeco_spikein_removed.rds"
out_base <- "results/Step_03B_Negative_Controls_Composition_Plots"

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(rds_in)) stop("❌ Input RDS not found: ", rds_in)
mt <- readRDS(rds_in)
if (!inherits(mt, "microtable")) stop("❌ Input object is not a microeco::microtable")
message("✅ Loaded microeco object: ", rds_in)

# ---- 2) Column checks + rank resolver --------------------------------------
tt_cols <- colnames(mt$tax_table)
st_cols <- colnames(mt$sample_table)

# Allow either singular/plural column naming (your earlier scripts used both)
rank_candidates <- list(
  Genus   = c("Genera_level_OTUs", "Genera_level_OTU", "Genus_level_OTUs", "Genus_level_OTU"),
  Class   = c("Classes", "Class")
)

pick_col <- function(cands, available) {
  hit <- cands[cands %in% available]
  if (length(hit)) hit[1] else NA_character_
}

tax_genus_col <- pick_col(rank_candidates$Genus, tt_cols)
tax_class_col <- pick_col(rank_candidates$Class, tt_cols)

if (is.na(tax_genus_col) || is.na(tax_class_col)) {
  stop(
    "❌ Missing required taxonomy columns.\n",
    "Needed: one of {", paste(rank_candidates$Genus, collapse = ", "), "} AND one of {",
    paste(rank_candidates$Class, collapse = ", "), "}\n",
    "Found tax_table cols: ", paste(tt_cols, collapse = ", ")
  )
}

if (!("Groups" %in% st_cols)) {
  stop("❌ sample_table must contain column 'Groups' to select blanks. Found: ",
       paste(st_cols, collapse = ", "))
}

# Optional facet vars (only those present will be used)
facet_vars_default <- c("Labs", "DNA_Extraction_Batch_ID", "Groups")
facet_vars <- facet_vars_default[facet_vars_default %in% st_cols]
if (!length(facet_vars)) {
  warning("⚠️ None of facet vars found: ", paste(facet_vars_default, collapse = ", "),
          " — plots will be generated without facetting.")
  facet_vars <- NULL
}

message("✅ Using taxrank (genus): ", tax_genus_col)
message("✅ Using high_level (class): ", tax_class_col)
message("✅ Facet vars: ", if (is.null(facet_vars)) "NONE" else paste(facet_vars, collapse = ", "))

# ---- 3) Define negative-control sets ---------------------------------------
labs_neg <- list(
  Lab1      = "Blank1",
  Lab2      = c("Blank2","Blank3","Blank4","Blank5","Blank6","Blank7",
                "Blank8","Blank9","Blank10","Blank11","Blank12",
                "Blank13","Blank14","Blank15","Blank16","Blank17",
                "Blank18","Blank19"),
  Both_labs = c("Blank1","Blank2","Blank3","Blank4","Blank5","Blank6","Blank7",
                "Blank8","Blank9","Blank10","Blank11","Blank12","Blank13",
                "Blank14","Blank15","Blank16","Blank17","Blank18","Blank19")
)

# ---- 4) Plot builder: nested barplot ---------------------------------------
make_nested_plot <- function(trans_obj, facets = NULL) {
  
  # NOTE: microeco uses ggnested if ggnested package is installed.
  # If it is not installed, plot_bar(ggnested=TRUE) can error.
  # We handle that by trying ggnested=TRUE and falling back to ggnested=FALSE.
  
  p <- try(
    trans_obj$plot_bar(
      ggnested = TRUE,
      facet    = facets
    ),
    silent = TRUE
  )
  
  if (inherits(p, "try-error")) {
    warning("⚠️ ggnested=TRUE failed (ggnested package missing?). Falling back to ggnested=FALSE.")
    p <- trans_obj$plot_bar(
      ggnested = FALSE,
      facet    = facets
    )
  }
  
  p +
    theme_classic() +
    theme(
      axis.title      = element_text(size = 16, face = "bold", family = "Times New Roman"),
      axis.text.x     = element_text(size = 14, angle = 65, hjust = 1,
                                     face = "bold", family = "Times New Roman"),
      axis.text.y     = element_text(size = 16, family = "Times New Roman"),
      strip.text      = element_text(size = 16, face = "bold", family = "Times New Roman"),
      legend.title    = element_text(size = 13, face = "bold", family = "Times New Roman"),
      legend.text     = element_text(size = 12, face = "bold", family = "Times New Roman"),
      legend.key.size = unit(0.4, "cm"),
      legend.spacing  = unit(0.05, "cm"),
      plot.title      = element_text(size = 14, face = "bold", hjust = 0.5,
                                     family = "Times New Roman"),
      plot.margin     = margin(3, 3, 3, 3)
    ) +
    labs(
      x = "Samples",
      y = "Relative abundance (%)",
      fill = tax_genus_col
    ) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE), colour = "none")
}

# ---- 5) Export helper: PNG/PDF/SVG/EPS/TIFF/JPEG + HTML --------------------
save_all_formats <- function(p, out_dir, lab_name) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Common width/height
  w <- 26; h <- 14; dpi <- 1200
  
  formats <- list(
    png  = list(ext = "png",  dev = "png"),
    pdf  = list(ext = "pdf",  dev = cairo_pdf),
    svg  = list(ext = "svg",  dev = svglite::svglite),
    eps  = list(ext = "eps",  dev = cairo_ps),
    tiff = list(ext = "tiff", dev = "tiff"),
    jpeg = list(ext = "jpeg", dev = "jpeg")
  )
  
  for (fmt in names(formats)) {
    f <- formats[[fmt]]
    out_file <- file.path(out_dir, paste0("nested_barplot_", lab_name, ".", f$ext))
    ggsave(out_file, plot = p, width = w, height = h, dpi = dpi, device = f$dev)
  }
  
  # Interactive HTML (FIXED)
  html_file <- file.path(out_dir, paste0("nested_barplot_", lab_name, ".html"))
  widget <- plotly::ggplotly(p)
  
  has_pandoc <- requireNamespace("rmarkdown", quietly = TRUE) &&
    isTRUE(rmarkdown::pandoc_available())
  
  if (has_pandoc) {
    # One-file HTML
    htmlwidgets::saveWidget(widget, file = html_file, selfcontained = TRUE)
  } else {
    # No pandoc => write dependencies to a controlled libdir in out_dir
    libdir <- file.path(out_dir, paste0("nested_barplot_", lab_name, "_files"))
    dir.create(libdir, recursive = TRUE, showWarnings = FALSE)
    
    # selfcontained=FALSE avoids normalizePath mustWork crashes on Windows
    htmlwidgets::saveWidget(widget, file = html_file, selfcontained = FALSE, libdir = libdir)
  }
  
  message("✅ Saved ", lab_name, " plots → ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
}

# ---- 6) Main loop over Lab1/Lab2/Both --------------------------------------
for (lab_name in names(labs_neg)) {
  
  negs <- labs_neg[[lab_name]]
  message("\n==============================")
  message("▶ Processing: ", lab_name)
  message("==============================")
  
  # Deep-clone microtable (safe copy for R6 object)
  mt_sub <- mt$clone(deep = TRUE)
  
  # Subset to blanks by Groups
  mt_sub$sample_table <- mt_sub$sample_table %>%
    dplyr::filter(.data$Groups %in% negs)
  
  # If nothing left, skip safely
  if (nrow(mt_sub$sample_table) == 0) {
    warning("⚠️ No samples found for ", lab_name, " (Groups in: ", paste(negs, collapse = ", "), "). Skipping.")
    next
  }
  
  # Tidy to sync otu_table/sample_table
  mt_sub$tidy_dataset()
  
  # Build transformed abundance object
  trans_obj <- trans_abund$new(
    dataset                = mt_sub,
    taxrank                = tax_genus_col,
    ntaxa                  = 130,
    show                   = 0,
    high_level             = tax_class_col,
    delete_taxonomy_prefix = FALSE
  )
  
  # Plot
  p <- make_nested_plot(trans_obj, facet_vars) +
    ggtitle(paste0("Non-template DNA-extraction blanks (", lab_name, ")"))
  
  # Save
  save_all_formats(p, file.path(out_base, lab_name), lab_name)
}

message("\n🎉 Done! All Lab1, Lab2, and Both_labs plots exported (plus HTML).")
message("Output base: ", normalizePath(out_base, winslash = "/", mustWork = FALSE))
