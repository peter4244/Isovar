test_that(".strip_phenotype_strand removes trailing _+/_-", {
  expect_equal(.strip_phenotype_strand("1:100:200:clu_5_-"), "1:100:200:clu_5")
  expect_equal(.strip_phenotype_strand("1:100:200:clu_5_+"), "1:100:200:clu_5")
  expect_equal(.strip_phenotype_strand("1:100:200:clu_5"),   "1:100:200:clu_5")
})

test_that("loadSqtlResults errors on missing parquet", {
  expect_error(loadSqtlResults("/tmp/__no_such_parquet__"),
               "sQTL parquet not found")
})

test_that("loadSqtlResults errors on invalid config", {
  tmp <- tempfile(fileext = ".R")
  writeLines("# no SQTL_PARQUET here", tmp)
  expect_error(loadSqtlResults(tmp, chr = "chr1"),
               "missing.*SQTL_PARQUET")
})

test_that("loadSqtlResults and annotateSqtl work against the AKR1A1 region", {
  skip_on_cran()
  cfg <- file.path(.isovar_pkgroot(), "inst/sqtl_configs/haec185_basal_leafcutter.R")
  if (!file.exists(cfg)) skip("sQTL config not present")

  parquet_path <- "/Users/petecastaldi/claude_projects/HAEC185/shortread/sQTL/sQTL_chr1.cis_qtl_pairs.1.parquet"
  if (!file.exists(parquet_path)) skip("sQTL parquet not on this laptop")

  # Pull sQTL data for a handful of AKR1A1 credible-set variant IDs
  variant_ids <- c(
    "chr1:45551015:T:C",  # rs9147
    "chr1:45480964:G:T",  # rs4660861
    "chr1:45517408:C:T"   # rs11211129
  )
  sq <- loadSqtlResults(cfg, chr = "chr1", variant_ids = variant_ids,
                        phenotype_ids = "1:45551155:45552436:clu_33775")
  expect_s3_class(sq, "tbl_df")
  # Should have at most 3 rows (one per variant × phenotype match)
  expect_true(nrow(sq) >= 1L)
  expect_true(all(!is.na(sq$slope)))
  # Metadata carries source + tissue
  meta <- getIsovarMeta(sq)
  src <- meta$sources[[1]]
  expect_equal(src$kind,   "short_read_sqtl")
  expect_equal(src$tissue, "basal airway epithelium")

  # annotateSqtl joins onto a variant tibble
  variants <- tibble::tibble(variant_id = variant_ids)
  out <- annotateSqtl(variants, sq)
  expect_true("sqtl_slope" %in% names(out))
  expect_equal(length(unique(out$variant_id)), length(variant_ids))
})
