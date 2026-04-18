#!/usr/bin/env Rscript
# AKR1A1 short-read Leafcutter sQTL tag-SNP effects on the phenotype
# intron (clu_33775: 45551155 → 45552436 on chr1). Indels are absent
# from this sQTL call set, so indel effects must be inferred via their
# tag SNPs + LD (handled by the haplotype layer in Chunk B).

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(tidyr); library(readr); library(tibble); library(cli); library(jsonlite); library(arrow)
})
isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

splaire <- "/Users/petecastaldi/claude_projects/HAEC185/HAEC185_sQTL_credible_sets_splaire_scores.tsv.gz"
gwas    <- "/Users/petecastaldi/claude_projects/copd/gwas/icgcUkb-20251213.dn8.gz"
sqtl_cfg <- file.path(isovar_dir, "inst/sqtl_configs/haec185_basal_leafcutter.R")

cli::cli_h1("AKR1A1 — short-read Leafcutter sQTL on phenotype intron")

# Build the merged credible set (splaire + gnomAD + GWAS)
flat <- as_flat(buildAnnotatedCredibleSet(
  sm_predictions   = splaire,
  gwas             = gwas,
  gene             = "AKR1A1",
  gnomad_cache_dir = "/tmp/isovar_cache"
))

# Pull sQTL rows for all credible-set variant_ids × the phenotype intron
sqtl_rows <- loadSqtlResults(
  sqtl_cfg, chr = "chr1",
  variant_ids   = flat$variant_id,
  phenotype_ids = "1:45551155:45552436:clu_33775"
)
cli::cli_alert_info("sQTL rows loaded: {nrow(sqtl_rows)} \\
                    (out of {length(flat$variant_id)} variants tested against phenotype)")

# Join onto the flat table
annotated <- annotateSqtl(flat, sqtl_rows)

# Collapse to one row per variant (sqtl_rows is already per-variant because
# we filtered to one phenotype_id, but be safe).
display <- annotated %>%
  arrange(desc(top_abs_delta)) %>%
  mutate(
    kind    = ifelse(nchar(ref) == 1L & nchar(alt) == 1L, "SNP", "indel"),
    splaire = sprintf("%.3f", top_abs_delta),
    beta    = ifelse(is.na(sqtl_slope), "    —", sprintf("%+.4f", sqtl_slope)),
    se      = ifelse(is.na(sqtl_slope_se), "   —", sprintf("%.4f", sqtl_slope_se)),
    p       = ifelse(is.na(sqtl_p),    "      —", format(sqtl_p, digits = 2, scientific = TRUE)),
    af      = ifelse(is.na(sqtl_af), "  —", sprintf("%.3f", sqtl_af))
  ) %>%
  select(variant_id, rsid, kind, pip,
         splaire, sqtl_slope = beta, sqtl_se = se, sqtl_p = p, sqtl_af = af,
         gwas_p)

cli::cli_h2("Top 12 by |splaire delta|, with sQTL on the phenotype intron")
print(as.data.frame(head(display, 12)), row.names = FALSE)

n_with_sqtl    <- sum(!is.na(annotated$sqtl_slope))
n_sqtl_sig     <- sum(!is.na(annotated$sqtl_p) & annotated$sqtl_p < 1e-6)
n_indels_total <- sum(nchar(annotated$ref) > 1 | nchar(annotated$alt) > 1)
n_indels_miss  <- sum((nchar(annotated$ref) > 1 | nchar(annotated$alt) > 1) &
                       is.na(annotated$sqtl_slope))

cli::cli_h2("Coverage summary")
cat(sprintf("  variants total            : %d\n", nrow(annotated)))
cat(sprintf("  with sQTL on phenotype    : %d\n", n_with_sqtl))
cat(sprintf("  sQTL significant (p<1e-6) : %d\n", n_sqtl_sig))
cat(sprintf("  indels (all)              : %d\n", n_indels_total))
cat(sprintf("  indels missing in sQTL    : %d   (expected — indel calls excluded)\n",
            n_indels_miss))

cli::cli_h2("Provenance")
meta <- getIsovarMeta(annotated)
for (s in meta$sources) {
  cat(sprintf("  [%s] %s  (cohort=%s, tissue=%s)\n",
              s$kind, s$id,
              s$cohort %||% "—",
              s$tissue %||% "—"))
}
