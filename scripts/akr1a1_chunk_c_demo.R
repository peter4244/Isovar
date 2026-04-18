#!/usr/bin/env Rscript
# Chunk C review artifact — run isoscope against the HAEC-185 basal
# and NMD lungcells sources for AKR1A1, classify isoforms by which
# implicated splice sites they use, and surface per-source provenance.

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(tidyr); library(readr); library(tibble); library(cli); library(jsonlite)
})

isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

cli::cli_h1("AKR1A1 — isoscope runs across long-read sources")

cfgs <- list(
  haec185 = file.path(isovar_dir, "inst/isoscope_configs/haec185_basal.R"),
  nmd     = file.path(isovar_dir, "inst/isoscope_configs/nmd_lungcells.R")
)

# Sites to check: the donors that splaire implicated for the two indels.
# + strand, so donor = junction start position.
implicated_sites <- c(
  donor_ref_intron_5ss   = 45551155,  # rs61467610 INS — canonical phenotype donor
  donor_shifted_INS_5ss  = 45551161,  # rs61467610 INS — +6 bp shifted
  donor_upstream_ref     = 45547554,  # rs66922050 DEL — cryptic donor
  donor_upstream_shifted = 45547544   # rs66922050 DEL — -10 bp shifted
)

results <- list()
for (src in names(cfgs)) {
  cli::cli_h2("source: {src}")
  iso <- runIsoscopeGene("AKR1A1", cfgs[[src]],
                         output_dir = file.path(tempdir(),
                                                paste0("isovar_iso_", src)))
  iso_classified <- classifyIsoformsBySiteUsage(
    iso, site_positions = unname(implicated_sites),
    kind = "junction_start"
  )

  n_gen <- sum(iso_classified$source == "GENCODE", na.rm = TRUE)
  n_sqa <- sum(iso_classified$source == "SQANTI",  na.rm = TRUE)
  cli::cli_alert_info("isoforms: {nrow(iso_classified)} total \\
                      (GENCODE={n_gen}, SQANTI={n_sqa})")

  usage <- tibble(
    site_label = names(implicated_sites),
    pos        = unname(implicated_sites)
  )
  usage$n_total <- vapply(usage$pos, function(p) {
    col <- paste0("uses_", p, "_jstart")
    if (col %in% names(iso_classified)) sum(iso_classified[[col]], na.rm = TRUE) else NA_integer_
  }, integer(1))
  usage$n_sqanti <- vapply(usage$pos, function(p) {
    col <- paste0("uses_", p, "_jstart")
    if (col %in% names(iso_classified))
      sum(iso_classified[[col]] & iso_classified$source == "SQANTI", na.rm = TRUE) else NA_integer_
  }, integer(1))
  cat("\nSite usage:\n")
  print(as.data.frame(usage), row.names = FALSE)

  results[[src]] <- iso_classified
}

cli::cli_h1("Side-by-side site usage (implicated splice sites)")
side_by_side <- tibble(
  site_label = names(implicated_sites),
  pos        = unname(implicated_sites),
  haec185_n_total  = vapply(unname(implicated_sites), function(p) {
    col <- paste0("uses_", p, "_jstart")
    if (col %in% names(results$haec185)) sum(results$haec185[[col]], na.rm = TRUE) else NA_integer_
  }, integer(1)),
  nmd_n_total = vapply(unname(implicated_sites), function(p) {
    col <- paste0("uses_", p, "_jstart")
    if (col %in% names(results$nmd)) sum(results$nmd[[col]], na.rm = TRUE) else NA_integer_
  }, integer(1))
)
print(as.data.frame(side_by_side), row.names = FALSE)

cli::cli_h1("Per-source provenance")
for (src in names(results)) {
  m <- getIsovarMeta(results[[src]])
  cat(sprintf("\n[%s]\n", src))
  cat(sprintf("  cohort     : %s\n", m$sources[[1]]$cohort))
  cat(sprintf("  tissue     : %s\n", m$sources[[1]]$tissue))
  cat(sprintf("  n_samples  : %s\n", m$sources[[1]]$n_samples))
  cat(sprintf("  build      : %s   gtf: %s\n",
              m$genome_build, m$gtf_version))
  cat(sprintf("  conditions : %s\n",
              if (is.null(m$sources[[1]]$conditions)) "—"
              else paste(m$sources[[1]]$conditions, collapse = ", ")))
}
