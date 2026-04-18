#!/usr/bin/env Rscript
# Chunk B review artifact — AKR1A1 haplotype grouping via LDlinkR.
# Uses the pre-built fixture at tests/testthat/fixtures/ldlink_akr1a1_r2.rds
# (run scripts/build_haplotype_fixture.R once to regenerate).

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(tidyr); library(readr); library(tibble); library(cli); library(jsonlite)
})

isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

flat <- as_flat(buildAnnotatedCredibleSet(
  sm_predictions   = file.path(isovar_dir, "inst/extdata/akr1a1_splaire_test.tsv.gz"),
  gwas             = file.path(isovar_dir, "inst/extdata/gwas_akr1a1_test.tsv.gz"),
  gene             = "AKR1A1",
  gnomad_cache_dir = "/tmp/isovar_cache"
))

fx <- file.path(isovar_dir, "tests/testthat/fixtures/ldlink_akr1a1_r2.rds")

cli::cli_h1("AKR1A1 haplotype grouping at multiple r² thresholds")

for (thresh in c(0.8, 0.9, 0.95, 0.99)) {
  out <- groupHaplotypes(flat, r2_threshold = thresh, fixture = fx)
  hap_counts <- sort(table(out$haplotype_id), decreasing = TRUE)
  cli::cli_h2("r² ≥ {thresh}")
  cli::cli_alert_info("haplotypes: {length(hap_counts)}   \\
                      sizes: {paste(hap_counts, collapse = ', ')}")
  cc <- out[out$causal_candidate, ] %>%
    arrange(haplotype_id)
  show <- cc %>%
    mutate(splaire = sprintf("%.3f", top_abs_delta),
           beta    = ifelse(is.na(gwas_beta), "    —", sprintf("%+.4f", gwas_beta)),
           p       = ifelse(is.na(gwas_p),    "     —", format(gwas_p, digits = 2, scientific = TRUE))) %>%
    select(haplotype_id, rsid, splaire, beta, p, assigned = haplotype_assigned_by)
  print(as.data.frame(show), row.names = FALSE)
}

# Default run (r² >= 0.8) with detailed view of the two indels
cli::cli_h1("Detailed view at r² ≥ 0.8 — both indels + their haplotype tags")
out <- groupHaplotypes(flat, r2_threshold = 0.8, fixture = fx)
focus_rs <- c("rs66922050", "rs61467610", "rs11444006", "rs9147", "rs4660861",
              "rs11211129", "rs2993263")
show <- out %>%
  filter(rsid %in% focus_rs) %>%
  arrange(pos) %>%
  mutate(splaire = sprintf("%.3f", top_abs_delta),
         r2_to_causal = ifelse(is.na(ld_r2_to_causal), "   —",
                               sprintf("%.3f", ld_r2_to_causal)),
         beta = ifelse(is.na(gwas_beta), "    —", sprintf("%+.4f", gwas_beta)),
         p    = ifelse(is.na(gwas_p),    "     —", format(gwas_p, digits = 2, scientific = TRUE)),
         causal = ifelse(causal_candidate, "yes", "")) %>%
  select(rsid, pos, splaire, haplotype_id, causal, assigned = haplotype_assigned_by,
         r2_to_causal, beta, p)
print(as.data.frame(show), row.names = FALSE)

cli::cli_h1("Resolution selection — GWAS coherence per r² level")

out_multi <- groupHaplotypes(flat, r2_thresholds = c(0.8, 0.9, 0.95), fixture = fx)
hap_sum   <- summarizeHaplotypesByResolution(out_multi)
rollup    <- rollupHaplotypeResolutions(hap_sum)

cat("Per-haplotype summary (β aligned for swapped orientation):\n\n")
disp <- hap_sum %>%
  arrange(resolution, desc(n_variants)) %>%
  mutate(
    causal = sprintf("%s (|Δ|=%.3f)", causal_rsid, causal_abs_delta),
    mean_b = ifelse(is.na(mean_beta), "     —", sprintf("%+.4f", mean_beta)),
    sd_b   = ifelse(is.na(sd_beta), "    —", sprintf("%.4f", sd_beta)),
    sc     = ifelse(is.na(sign_concordance), " —", sprintf("%.2f", sign_concordance)),
    min_p  = ifelse(is.na(min_p), "      —", format(min_p, digits = 2, scientific = TRUE))
  ) %>%
  select(res = resolution, hap = haplotype_id, n = n_variants, n_gw = n_gwas_variants,
         causal, mean_b, sd_b, sc, min_p, n_gws)
print(as.data.frame(disp), row.names = FALSE)

cat("\nAggregate per-resolution coherence:\n\n")
print(as.data.frame(rollup %>% mutate(
  mean_within_hap_beta_sd = sprintf("%.4f", mean_within_hap_beta_sd),
  mean_sign_concordance   = sprintf("%.3f", mean_sign_concordance),
  frac_significant_haps   = sprintf("%.2f", frac_significant_haps)
)), row.names = FALSE)

cli::cli_h1("Provenance (sources attached to output)")
meta <- getIsovarMeta(out)
for (s in meta$sources)
  cat(sprintf("  [%s] %s  (%s)\n", s$kind, s$id,
              s$genome_build %||% s$version %||% s$population %||% ""))
