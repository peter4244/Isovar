#!/usr/bin/env Rscript
# Chunk A / A.3 review artifact — AKR1A1 credible set merged with
# splaire splice-prediction deltas, gnomAD rsID / AF (default fields),
# and COPD GWAS effect estimates (rsID + position match) as a
# nested credible-set tibble with propagated provenance metadata.
#
# Usage:
#   cd ~/claude_projects/isovar
#   Rscript scripts/akr1a1_chunk_a_demo.R
#
# Writes:
#   runs/akr1a1/credible_sets.tsv (+ .meta.json)
#   runs/akr1a1/variants.tsv       (+ .meta.json)

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(tidyr); library(readr); library(tibble); library(cli); library(jsonlite)
})

isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

splaire_fixture <- file.path(isovar_dir, "inst/extdata/akr1a1_splaire_test.tsv.gz")
gwas_fixture    <- file.path(isovar_dir, "inst/extdata/gwas_akr1a1_test.tsv.gz")

cli::cli_h1("AKR1A1 — Chunk A merged credible-set table (nested)")

nested <- buildAnnotatedCredibleSet(
  sm_predictions   = splaire_fixture,
  gwas             = gwas_fixture,
  gene             = "AKR1A1",
  gnomad_cache_dir = "/tmp/isovar_cache"
)

cli::cli_h2("Nested shape")
cli::cli_alert_info("credible-set rows: {nrow(nested)}   \\
                    top-level cols: {ncol(nested)}")
cli::cli_alert_info("per-variant tibble nrow: {nested$n_variants[1]}, \\
                    ncol: {ncol(nested$variants[[1]])}")
cli::cli_alert_info("per-variant columns:")
cat("  ", paste(names(nested$variants[[1]]), collapse = ", "), "\n")

cli::cli_h2("Provenance metadata (attached to result)")
meta <- getIsovarMeta(nested)
cat("genome_build : ", meta$genome_build,  "\n")
cat("gtf_version  : ", meta$gtf_version,   "\n")
cat("generated_at : ", meta$generated_at,  "\n")
cat("sources      :\n")
for (s in meta$sources) {
  cat(sprintf("  - [%s] %s  (%s)\n",
              s$kind, s$id,
              s$genome_build %||% s$version %||% ""))
}

# Write TSV pair + sidecars
out_dir <- file.path(isovar_dir, "runs/akr1a1")
paths <- write_tsv_pair(nested, out_dir)
cli::cli_alert_success("Wrote:")
cat(sprintf("  credible_sets : %s\n", paths$credible_sets))
cat(sprintf("  variants       : %s\n", paths$variants))
cat(sprintf("  sidecar metas  : .meta.json next to each\n"))

# Flat preview
flat <- as_flat(nested)
cli::cli_h2("Top 10 variants by |splaire delta|")
display <- flat %>%
  arrange(desc(top_abs_delta)) %>%
  mutate(kind = ifelse(nchar(ref) == 1L & nchar(alt) == 1L, "SNP", "indel"),
         splaire = sprintf("%.3f", top_abs_delta),
         pip_   = sprintf("%.3f", pip),
         beta   = ifelse(is.na(gwas_beta), "     —", sprintf("%+.4f", gwas_beta)),
         p      = ifelse(is.na(gwas_p),    "      —", format(gwas_p, digits = 2, scientific = TRUE))) %>%
  select(variant_id, rsid, kind, pip = pip_, splaire, head = top_head, dir = top_direction,
         af, af_nfe, af_eas, grpmax, beta, p, mb = gwas_matched_by)
print(as.data.frame(head(display, 10)), row.names = FALSE)
