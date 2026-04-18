akr1a1_gwas_fixture <- function() {
  file.path(.isovar_pkgroot(), "inst/extdata/gwas_akr1a1_test.tsv.gz")
}

test_that("loadCopdGwas reads the icgcUkb schema", {
  g <- loadCopdGwas(akr1a1_gwas_fixture())
  expect_s3_class(g, "tbl_df")
  expect_true(all(c("rsid","effectallele","otherallele","beta","se","p","eaf")
                  %in% names(g)))
  expect_gt(nrow(g), 0L)
  expect_type(g$beta, "double")
  expect_type(g$pos,  "integer")
})

test_that("loadCopdGwas errors on missing file", {
  expect_error(loadCopdGwas("/tmp/__no_such_gwas__.tsv.gz"),
               "GWAS file not found")
})

test_that("annotateGwas joins by rsID and preserves input order", {
  g <- loadCopdGwas(akr1a1_gwas_fixture())
  variants <- tibble::tibble(
    variant_id = c("chr1:45551015:T:C",
                   "chr1:45480964:G:T",
                   "chr1:45549863:A:ATATCTG",
                   "chr1:00000000:A:T"),
    rsid = c("rs9147", "rs4660861", "rs61467610", NA_character_),
    chr  = c("chr1", "chr1", "chr1", "chr1"),
    pos  = c(45551015L, 45480964L, 45549863L, 0L),
    ref  = c("T", "G", "A", "A"),
    alt  = c("C", "T", "ATATCTG", "T")
  )
  # With no liftover chain supplied and different builds, position
  # fallback is skipped — matches come from rsID only.
  ann <- annotateGwas(variants, g)
  expect_equal(ann$variant_id, variants$variant_id)
  expect_equal(sum(!is.na(ann$gwas_beta)), 2L)
  expect_false(is.na(ann$gwas_beta[1]))
  expect_false(is.na(ann$gwas_beta[2]))
  expect_true(is.na(ann$gwas_beta[3]))
  expect_true(is.na(ann$gwas_beta[4]))
  # gwas_matched_by correctly reports "rsid" for the two matches
  expect_equal(ann$gwas_matched_by, c("rsid", "rsid", NA_character_, NA_character_))
})

test_that("annotateGwas attaches GWAS fields for the AKR1A1 tag SNPs present in icgcUkb", {
  g <- loadCopdGwas(akr1a1_gwas_fixture())
  variants <- tibble::tibble(
    rsid = c("rs4660861", "rs9147"),
    chr  = c("chr1", "chr1"),
    pos  = c(45480964L, 45551015L),
    ref  = c("G", "T"),
    alt  = c("T", "C")
  )
  ann <- annotateGwas(variants, g)
  expect_true(all(!is.na(ann$gwas_beta)))
  expect_lt(ann$gwas_p[ann$rsid == "rs4660861"], 1e-6)
  expect_true(all(ann$gwas_matched_by == "rsid"))
})

test_that("annotateGwas position fallback fills in a row when rsID disagrees", {
  # Synthetic GWAS tibble on the same build as variants (no liftover needed).
  gwas <- tibble::tibble(
    rsid          = c("rs9147", "rs999_different_id"),
    chr           = c("chr1", "chr1"),
    pos           = c(45551015L, 45480964L),
    marker        = c("1:45551015:T:C", "1:45480964:G:T"),
    effectallele  = c("C", "T"),
    otherallele   = c("T", "G"),
    eaf           = c(0.5, 0.4),
    beta          = c(0.05, -0.05),
    se            = c(0.01, 0.01),
    p             = c(1e-6, 5e-8)
  )
  variants <- tibble::tibble(
    variant_id = c("chr1:45551015:T:C", "chr1:45480964:G:T"),
    rsid       = c("rs9147", "rs4660861"),   # rs4660861 won't hit rsid join
    chr        = c("chr1", "chr1"),
    pos        = c(45551015L, 45480964L),
    ref        = c("T", "G"),
    alt        = c("C", "T")
  )
  # Same build on both sides -> position fallback runs without a chain.
  ann <- annotateGwas(variants, gwas,
                      variants_build = "GRCh38",
                      gwas_build     = "GRCh38")
  expect_equal(ann$gwas_matched_by, c("rsid", "position"))
  expect_equal(ann$gwas_p[2], 5e-8)
})

test_that("checkAlleleAlignment classifies aligned / swapped / strand_flip / mismatch", {
  # Use non-palindromic SNPs (avoid A/T and C/G pairs) so strand_flip is
  # unambiguous.
  df <- tibble::tibble(
    ref                = c("A",  "A", "A", "A", "A"),
    alt                = c("G",  "G", "C", "C", "G"),
    gwas_effectallele  = c("G",  "A", "G", "A", NA_character_),
    gwas_otherallele   = c("A",  "G", "T", "G", NA_character_)
  )
  # Row 1: aligned (alt=effect, ref=other)
  # Row 2: swapped (ref=effect, alt=other)
  # Row 3: strand_flip (A->C on + strand = T->G on - strand; gwas has G,T)
  # Row 4: mismatch (A->C, gwas reports A,G — neither aligned, swapped, nor
  #                  a valid strand flip: ref_c=T, alt_c=G; T!=G, T!=A, G!=G? -
  #                  actually G == other, T != effect. Nope on aligned_flip.
  #                  swapped_flip: T==effect(A)? no. so mismatch)
  # Row 5: no GWAS match (NA effect/other)
  a <- checkAlleleAlignment(df)
  expect_equal(a$alignment_status,
               c("aligned", "swapped", "strand_flip", "mismatch", NA_character_))
  expect_equal(a$alignment_ok, c(TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("checkAlleleAlignment handles indels without a strand-flip path", {
  df <- tibble::tibble(
    ref                = c("A",       "TTTGAACTTCG"),
    alt                = c("ATATCTG", "T"),
    gwas_effectallele  = c("ATATCTG", "T"),
    gwas_otherallele   = c("A",       "TTTGAACTTCG")
  )
  a <- checkAlleleAlignment(df)
  expect_equal(a$alignment_status, c("aligned", "aligned"))
  expect_true(all(a$alignment_ok))
})

test_that("end-to-end from isovar rank -> gwas join works on AKR1A1 fixture", {
  r <- rankSplaireVariants(akr1a1_fixture(), gene = "AKR1A1")
  rs_map <- c(
    "chr1:45549863:A:ATATCTG"     = "rs61467610",
    "chr1:45542866:TTTGAACTTCG:T" = "rs66922050",
    "chr1:45551015:T:C"           = "rs9147",
    "chr1:45480964:G:T"           = "rs4660861"
  )
  r$rsid <- unname(rs_map[r$variant_id])
  r <- r[!is.na(r$rsid), ]

  g <- loadCopdGwas(akr1a1_gwas_fixture())
  ann <- annotateGwas(r, g)
  expect_equal(nrow(ann), length(rs_map))
  expect_equal(sum(!is.na(ann$gwas_beta)), 2L)
  expect_true(any(ann$gwas_p < 1e-7, na.rm = TRUE))
})

test_that("buildAnnotatedCredibleSet composes the full pipeline into one table", {
  skip_on_cran()
  skip_if_no_internet()
  tbl <- buildAnnotatedCredibleSet(
    sm_predictions   = akr1a1_fixture(),
    gwas             = akr1a1_gwas_fixture(),
    gene             = "AKR1A1",
    gnomad_cache_dir = "/tmp/isovar_cache"
  )
  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 50L)
  # Columns from all three layers should be present
  expect_true(all(c("top_abs_delta", "posterior_inclusion_probability",  # splaire
                    "rsid", "af", "af_nfe", "grpmax",                    # gnomAD
                    "gwas_beta", "gwas_p", "gwas_matched_by")            # GWAS
                  %in% names(tbl)))
  # Every variant got an rsID (gnomAD coverage is complete for AKR1A1)
  expect_true(all(!is.na(tbl$rsid)))
  # The known GWAS tag SNPs have beta populated
  rs_hits <- tbl[tbl$rsid %in% c("rs4660861", "rs9147"), ]
  expect_equal(nrow(rs_hits), 2L)
  expect_true(all(!is.na(rs_hits$gwas_beta)))
  # The known indels have NA GWAS (absent in icgcUkb)
  indel_misses <- tbl[tbl$rsid %in% c("rs61467610", "rs66922050"), ]
  expect_equal(nrow(indel_misses), 2L)
  expect_true(all(is.na(indel_misses$gwas_beta)))
})
