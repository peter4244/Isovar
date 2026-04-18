# isovar sQTL source config — HAEC-185 basal Leafcutter sQTL (short-read).
#
# Consumed by loadSqtlResults() to attach per-intron effect sizes to
# credible-set variants.

ISOVAR_SOURCE_META <- list(
  source_id        = "haec185_basal_leafcutter_sqtl",
  source_label     = "HAEC-185 basal Leafcutter sQTL (short-read, 185-donor cohort)",
  cohort           = "HAEC-185",
  tissue           = "basal airway epithelium",
  n_samples        = 185,
  genome_build     = "GRCh38",
  phenotype        = "leafcutter_intron_psi",
  tool             = "tensorqtl cis nominal",
  parquet_format   = "cis_qtl_pairs",
  includes_indels  = TRUE,   # verified 2026-04-18: chr1 parquet has per-intron β for the 3 AKR1A1 indels
  notes            = "Indels included. For cohorts where indels are absent, the haplotype layer should propagate tag-SNP β to the causal candidate."
)

# Per-chromosome parquet files
SQTL_PARQUET <- list(
  chr1 = "/Users/petecastaldi/claude_projects/HAEC185/shortread/sQTL/sQTL_chr1.cis_qtl_pairs.1.parquet"
  # chr2, chr3, ... added as they arrive
)

# Column schema is tensorqtl's cis_qtl_pairs:
#   phenotype_id, variant_id, tss_distance, af, ma_samples, ma_count,
#   pval_nominal, slope, slope_se
