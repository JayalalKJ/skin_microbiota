#!/usr/bin/env Rscript
# =============================================================================
# Runs the full pipeline (00 → 05) in order, with logs + fail-fast behavior.
#
# Usage:
#   Rscript run_all_pipeline.R
#
# Assumes scripts are in: scripts/
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- Config ----------------------------------------------------------------
scripts_dir <- "scripts_2"
log_dir     <- file.path("results", "_pipeline_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

pipeline <- c(
  #"00_project_setup.R",
  #"01_import_mothur_to_microeco.R",
  #"02_remove_spikein_taxa.R",
  #"03A_Blanks_visual_proofpack.R",
  #"03B_Negative_control_samples_Composition.R",
  #"03C_contaminant_removal_from_blanks.R",
  #"03D_batch_diagnostics_before_conqur_twosets_Script_1.R",
  #"03D_batch_correction_conqur_set6groups_Script_2.R",
  #"03E_convert_to_microeco_obj.R",
  #"04_rarefaction_HMHD_for_alpha.R",
  #"05_alpha_diversity_HMHD.R",
  #"06_Beta_diversity_HMHD.R",
  #"07_Beta_diversity_HMHD_ENV_feed_water.R",
  #"08_composition_analysis_HMHD_feed_water.R",
  "09_Beta_diversity_HM_HD_Fish_ID.R",
  "09A_Beta_diversity_Fish_HD_HM_feed_water.R"
  #"12_network_SparCC_hindgut.R",
  #"12B_LEFSE_LDA_GeneraLevelOTUs_HMHD_6TankGroups.R",
  #"beta_diversity_subset_export_Aitchison_locked63.R"
  
)

# ---- Helpers ---------------------------------------------------------------
stop_missing <- function(path) {
  if (!file.exists(path)) stop("❌ Missing file: ", path)
}

stamp <- function() format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

run_one <- function(script_name) {
  script_path <- file.path(scripts_dir, script_name)
  stop_missing(script_path)

  tag <- paste0(gsub("\\.R$", "", script_name), "__", stamp())
  out_log <- file.path(log_dir, paste0(tag, ".out.log"))
  err_log <- file.path(log_dir, paste0(tag, ".err.log"))

  message("\n──────────────────────────────────────────────────────────────")
  message("▶ Running: ", script_name)
  message("  out: ", out_log)
  message("  err: ", err_log)

  # Use the same R that is running this driver
  rbin <- file.path(R.home("bin"), "Rscript")
  if (!file.exists(rbin)) rbin <- "Rscript"

  exit_code <- system2(
    command = rbin,
    args    = c("--vanilla", script_path),
    stdout  = out_log,
    stderr  = err_log
  )

  if (!identical(exit_code, 0L)) {
    message("❌ Failed: ", script_name)
    message("   Check logs:\n   - ", out_log, "\n   - ", err_log)
    stop("Pipeline aborted (exit code = ", exit_code, ").")
  }

  message("✅ Done: ", script_name)
  invisible(list(script = script_name, out_log = out_log, err_log = err_log))
}

# ---- Run pipeline ----------------------------------------------------------
results <- lapply(pipeline, run_one)

# ---- Manifest --------------------------------------------------------------
manifest <- data.frame(
  step     = seq_along(pipeline),
  script   = vapply(results, `[[`, character(1), "script"),
  out_log  = vapply(results, `[[`, character(1), "out_log"),
  err_log  = vapply(results, `[[`, character(1), "err_log"),
  stringsAsFactors = FALSE
)

manifest_file <- file.path(log_dir, paste0("pipeline_manifest__", stamp(), ".tsv"))
write.table(manifest, file = manifest_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("\n🎉 Pipeline complete.")
message("• Manifest: ", normalizePath(manifest_file))
message("• Logs dir: ", normalizePath(log_dir))
