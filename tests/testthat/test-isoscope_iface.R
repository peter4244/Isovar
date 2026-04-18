akr1a1_isoscope_fixture <- function() {
  file.path(.isovar_pkgroot(), "tests/testthat/fixtures/isoscope_akr1a1_haec185.tsv")
}

skip_if_no_isoscope <- function() {
  script <- path.expand("~/claude_projects/isoscope/gene_isoform_annotation.R")
  if (!file.exists(script))
    testthat::skip(paste0("isoscope checkout not found at ", script))
}

skip_if_no_isoscope_fixture <- function() {
  p <- akr1a1_isoscope_fixture()
  if (!file.exists(p))
    testthat::skip(paste0("No isoscope fixture at ", p,
                          " — run an isoscope run + copy the TSV to the fixtures/ dir."))
}

test_that(".junction_starts / .junction_ends parse isoscope junctions correctly", {
  j <- "45551155_45561789_+;45561878_45566569_+;45566688_45566869_+"
  expect_equal(.junction_starts(j), c(45551155L, 45561878L, 45566688L))
  expect_equal(.junction_ends(j),   c(45561789L, 45566569L, 45566869L))
  expect_equal(.junction_starts(""),   integer(0))
  expect_equal(.junction_starts(NA),   integer(0))
})

test_that("classifyIsoformsBySiteUsage flags isoforms by junction-start usage", {
  iso <- tibble::tibble(
    isoform_id = c("A", "B", "C"),
    junctions  = c(
      "45551155_45561789_+;45561878_45566569_+",   # uses 45551155 as a donor
      "45551161_45561789_+;45561878_45566569_+",   # uses 45551161 (shifted)
      NA_character_                                 # single-exon
    )
  )
  out <- classifyIsoformsBySiteUsage(iso,
                                     site_positions = c(45551155, 45551161))
  expect_true("uses_45551155_jstart" %in% names(out))
  expect_true("uses_45551161_jstart" %in% names(out))
  expect_equal(out$uses_45551155_jstart, c(TRUE,  FALSE, FALSE))
  expect_equal(out$uses_45551161_jstart, c(FALSE, TRUE,  FALSE))
})

test_that("classifyIsoformsBySiteUsage kind='junction_end' uses acceptors", {
  iso <- tibble::tibble(
    junctions = c("45551155_45561789_+;45561878_45566569_+")
  )
  out <- classifyIsoformsBySiteUsage(iso, site_positions = c(45561789, 99999),
                                     kind = "junction_end")
  expect_true("uses_45561789_jend" %in% names(out))
  expect_true(out$uses_45561789_jend)
  expect_false(out$uses_99999_jend)
})

test_that("runIsoscopeGene errors cleanly on missing config or script", {
  expect_error(
    runIsoscopeGene("AKR1A1", "/tmp/__no_such_config__.R"),
    "isoscope config not found"
  )
  expect_error(
    runIsoscopeGene("AKR1A1", akr1a1_isoscope_fixture(),
                    isoscope_dir = "/tmp/__no_isoscope__"),
    "isoscope script not found"
  )
})

test_that("runIsoscopeGene validates the gene identifier", {
  # Valid identifiers include HGNC symbols and Ensembl IDs; reject shell
  # / regex injection.
  skip_if_no_isoscope()
  cfg <- file.path(.isovar_pkgroot(), "inst/isoscope_configs/haec185_basal.R")
  expect_error(runIsoscopeGene("A;B", cfg), "Invalid gene identifier")
})

test_that("runIsoscopeGene on AKR1A1 produces annotation with propagated metadata", {
  skip_on_cran()
  skip_if_no_isoscope()
  cfg <- file.path(.isovar_pkgroot(), "inst/isoscope_configs/haec185_basal.R")
  if (!file.exists(cfg))
    testthat::skip("haec185_basal config not present")
  iso <- runIsoscopeGene("AKR1A1", cfg, output_dir = tempfile("isovar_iso_"))
  expect_s3_class(iso, "tbl_df")
  expect_true(nrow(iso) > 0L)
  expect_true("junctions" %in% names(iso))
  # Metadata attached
  meta <- getIsovarMeta(iso)
  expect_equal(meta$genome_build,   "GRCh38")
  expect_equal(meta$gtf_version,    "GENCODE_v49")
  src <- meta$sources[[1]]
  expect_equal(src$kind,   "long_read_isoform_annotation")
  expect_equal(src$cohort, "HAEC-185")
})
