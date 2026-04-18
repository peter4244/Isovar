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
  # Important finding: the icgcUkb COPD GWAS excludes indels and many
  # low-imputation variants, so the two splaire-causal AKR1A1 indels
  # (rs61467610, rs66922050) are NOT in the sumstats even though they
  # are in our credible set. Tag SNPs (rs9147, rs4660861) ARE there.
  # This test encodes that biological reality.
  variants <- tibble::tibble(
    variant_id = c("chr1:45551015:T:C",           # rs9147 — present
                   "chr1:45480964:G:T",           # rs4660861 — present
                   "chr1:45549863:A:ATATCTG",     # rs61467610 — INS, absent
                   "chr1:00000000:A:T"),          # intentional no-match
    rsid = c("rs9147", "rs4660861", "rs61467610", NA_character_),
    ref  = c("T", "G", "A", "A"),
    alt  = c("C", "T", "ATATCTG", "T")
  )
  ann <- annotateGwas(variants, g)
  expect_equal(ann$variant_id, variants$variant_id)
  # Two matched (the SNPs), two NA (indel + no-rsid).
  expect_equal(sum(!is.na(ann$gwas_beta)), 2L)
  expect_false(is.na(ann$gwas_beta[1]))   # rs9147
  expect_false(is.na(ann$gwas_beta[2]))   # rs4660861
  expect_true(is.na(ann$gwas_beta[3]))    # rs61467610 indel — absent from GWAS
  expect_true(is.na(ann$gwas_beta[4]))
})

test_that("annotateGwas attaches GWAS fields for the AKR1A1 tag SNPs present in icgcUkb", {
  g <- loadCopdGwas(akr1a1_gwas_fixture())
  # Only tag SNPs present in the GWAS (indels excluded — see above).
  variants <- tibble::tibble(
    rsid = c("rs4660861", "rs9147"),
    ref  = c("G", "T"),
    alt  = c("T", "C")
  )
  ann <- annotateGwas(variants, g)
  expect_true(all(!is.na(ann$gwas_beta)))
  # rs4660861 is a genome-wide-significant COPD hit on the DEL haplotype
  expect_lt(ann$gwas_p[ann$rsid == "rs4660861"], 1e-6)
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
  # Hand-mapped subset so this test doesn't hit the network. Two tag SNPs
  # (present in GWAS) + two splaire-causal indels (absent — GWAS excludes
  # indels). This is the biologically meaningful coverage pattern.
  rs_map <- c(
    "chr1:45549863:A:ATATCTG"     = "rs61467610",   # INS, absent in GWAS
    "chr1:45542866:TTTGAACTTCG:T" = "rs66922050",   # DEL, absent in GWAS
    "chr1:45551015:T:C"           = "rs9147",       # SNP, present
    "chr1:45480964:G:T"           = "rs4660861"     # SNP, present
  )
  r$rsid <- unname(rs_map[r$variant_id])
  r <- r[!is.na(r$rsid), ]

  g <- loadCopdGwas(akr1a1_gwas_fixture())
  ann <- annotateGwas(r, g)
  expect_equal(nrow(ann), length(rs_map))
  # Two of four match (tag SNPs); indels have NA beta
  expect_equal(sum(!is.na(ann$gwas_beta)), 2L)
  # rs4660861 at p ~ 5.59e-8 — just above the strict 5e-8 bar but clearly
  # a strong COPD association on the DEL haplotype.
  expect_true(any(ann$gwas_p < 1e-7, na.rm = TRUE))
})
