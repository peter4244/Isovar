#!/usr/bin/env Rscript
# Chunk A review artifact — AKR1A1 credible set merged with splaire
# splice-prediction deltas, gnomAD rsID / AF, and COPD GWAS effect estimates
# in a single table via buildAnnotatedCredibleSet().
#
# Usage:
#   cd ~/claude_projects/isovar
#   Rscript scripts/akr1a1_chunk_a_demo.R
#
# Writes:
#   runs/akr1a1_chunk_a.tsv   — full merged table
#
# Uses local fixtures for the test data (no 190 MB GWAS load); swap
# fixture paths for the real ones to run against the full dataset.

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(readr); library(tibble); library(cli)
})

isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

splaire_fixture <- file.path(isovar_dir, "inst/extdata/akr1a1_splaire_test.tsv.gz")
gwas_fixture    <- file.path(isovar_dir, "inst/extdata/gwas_akr1a1_test.tsv.gz")

cli::cli_h1("AKR1A1 Chunk A — merged credible-set table")

tbl <- buildAnnotatedCredibleSet(
  splaire_path     = splaire_fixture,
  gwas             = gwas_fixture,
  gene             = "AKR1A1",
  gnomad_cache_dir = "/tmp/isovar_cache"
)

cli::cli_h2("Shape")
cli::cli_alert_info("rows: {nrow(tbl)}  cols: {ncol(tbl)}")
cli::cli_alert_info("variants with rsID       : {sum(!is.na(tbl$rsid))} / {nrow(tbl)}")
cli::cli_alert_info("matched in GWAS (fixture): {sum(!is.na(tbl$gwas_beta))} / {nrow(tbl)}")
cli::cli_alert_info("GWAS matched_by          : \\
                    rsid={sum(tbl$gwas_matched_by == 'rsid', na.rm=TRUE)}, \\
                    position={sum(tbl$gwas_matched_by == 'position', na.rm=TRUE)}, \\
                    none={sum(is.na(tbl$gwas_matched_by))}")

out_dir <- file.path(isovar_dir, "runs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, "akr1a1_chunk_a.tsv")
write_tsv(tbl, out_path)
cli::cli_alert_success("Wrote {.path {out_path}}")

# Human-readable preview (top by |splaire delta|)
cli::cli_h2("Top 10 by |splaire delta|")
display <- tbl %>%
  arrange(desc(top_abs_delta)) %>%
  mutate(kind = ifelse(nchar(ref) == 1L & nchar(alt) == 1L, "SNP", "indel"),
         splaire = sprintf("%.3f", top_abs_delta),
         pip     = sprintf("%.3f", posterior_inclusion_probability),
         af      = sprintf("%.3f", af),
         beta    = ifelse(is.na(gwas_beta), "     —", sprintf("%+.4f", gwas_beta)),
         p       = ifelse(is.na(gwas_p),    "      —", format(gwas_p, digits = 2, scientific = TRUE))) %>%
  select(variant_id, rsid, kind, pip, splaire, head = top_head, dir = top_direction,
         af, beta, p, matched_by = gwas_matched_by) %>%
  head(10)
print(as.data.frame(display), row.names = FALSE)

cli::cli_h2("GWAS-significant tag SNPs (p < 1e-6)")
sig <- tbl %>% filter(!is.na(gwas_p) & gwas_p < 1e-6) %>%
  arrange(gwas_p) %>%
  mutate(p = format(gwas_p, digits = 2, scientific = TRUE),
         beta = sprintf("%+.4f", gwas_beta),
         splaire = sprintf("%.3f", top_abs_delta)) %>%
  select(variant_id, rsid, splaire, beta, p, eaf = gwas_eaf,
         matched_by = gwas_matched_by)
print(as.data.frame(sig), row.names = FALSE)
