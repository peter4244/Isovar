#!/usr/bin/env Rscript
# One-time fixture build for Chunk B haplotype tests.
#
# Runs a live LDlinkR::LDmatrix call for all AKR1A1 credible-set rsIDs
# and saves the returned data.frame under tests/testthat/fixtures/.
# Tests then read from the fixture so CI doesn't hit LDlink.
#
# Requires: env var LDLINK_TOKEN set to your LDlink API token.
# Register free at https://ldlink.nci.nih.gov/?tab=apiaccess
#
# Run:
#   cd ~/claude_projects/isovar
#   LDLINK_TOKEN=<your_token> Rscript scripts/build_haplotype_fixture.R

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(tidyr); library(readr); library(tibble); library(cli); library(jsonlite)
  library(LDlinkR)
})

isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

token <- Sys.getenv("LDLINK_TOKEN")
if (!nzchar(token)) {
  cli::cli_abort(c(
    "LDLINK_TOKEN not set.",
    "i" = "Register free at https://ldlink.nci.nih.gov/?tab=apiaccess",
    "i" = "Then re-run: LDLINK_TOKEN=<your_token> Rscript scripts/build_haplotype_fixture.R"
  ))
}

# Build the AKR1A1 variant table (splaire + gnomAD rsIDs; GWAS not needed here).
splaire_fixture <- file.path(isovar_dir, "inst/extdata/akr1a1_splaire_test.tsv.gz")
ranked <- rankSplaireVariants(splaire_fixture, gene = "AKR1A1")
ann <- annotateGnomad(ranked$variant_id, cache_dir = "/tmp/isovar_cache")

rsids <- unique(stats::na.omit(ann$rsid))
cli::cli_alert_info("Querying LDlink for {length(rsids)} AKR1A1 rsIDs...")

ld_df <- LDlinkR::LDmatrix(
  snps         = rsids,
  pop          = "EUR",
  r2d          = "r2",
  token        = token,
  genome_build = "grch38"
)

out_dir <- file.path(isovar_dir, "tests/testthat/fixtures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, "ldlink_akr1a1_r2.rds")
saveRDS(ld_df, out_path)
cli::cli_alert_success("Wrote {.path {out_path}}  ({nrow(ld_df)} x {ncol(ld_df)})")
