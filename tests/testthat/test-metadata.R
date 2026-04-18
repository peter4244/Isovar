test_that("setIsovarMeta / getIsovarMeta round-trip", {
  x <- tibble::tibble(a = 1:3)
  x <- setIsovarMeta(x,
                     genome_build = "GRCh38",
                     gtf_version  = "GENCODE_v49",
                     sources      = list(list(id = "test", kind = "test")))
  m <- getIsovarMeta(x)
  expect_equal(m$genome_build, "GRCh38")
  expect_equal(m$gtf_version,  "GENCODE_v49")
  expect_equal(m$sources[[1]]$id, "test")
  expect_false(is.null(m$isovar_version))
  expect_false(is.null(m$generated_at))
})

test_that("mergeIsovarMeta concatenates sources and dedups by id", {
  a <- setIsovarMeta(tibble::tibble(x = 1),
                     genome_build = "GRCh38",
                     sources = list(list(id = "s1", kind = "a"),
                                    list(id = "s2", kind = "b")))
  b <- setIsovarMeta(tibble::tibble(x = 1),
                     genome_build = "GRCh38",
                     sources = list(list(id = "s2", kind = "b"),     # dup
                                    list(id = "s3", kind = "c")))
  m <- mergeIsovarMeta(a, b)
  ids <- vapply(m$sources, function(s) s$id, character(1))
  expect_equal(sort(ids), c("s1", "s2", "s3"))
})

test_that("mergeIsovarMeta flags mixed builds", {
  a <- setIsovarMeta(tibble::tibble(x=1), genome_build = "GRCh38")
  b <- setIsovarMeta(tibble::tibble(x=1), genome_build = "GRCh37")
  m <- mergeIsovarMeta(a, b)
  expect_equal(sort(m$mixed_builds), c("GRCh37", "GRCh38"))
})

test_that("writeMeta / readMeta round-trip via JSON sidecar", {
  x <- setIsovarMeta(tibble::tibble(a = 1:2),
                     genome_build = "GRCh38", gtf_version = "GENCODE_v49")
  tmp <- tempfile(fileext = ".tsv")
  readr::write_tsv(tibble::tibble(a = 1:2), tmp)
  writeMeta(x, tmp)
  expect_true(file.exists(paste0(tmp, ".meta.json")))
  m <- readMeta(tmp)
  expect_equal(m$genome_build, "GRCh38")
  expect_equal(m$gtf_version,  "GENCODE_v49")
})

test_that("getIsovarMeta(require=FALSE) returns NULL when no metadata", {
  x <- tibble::tibble(a = 1:3)
  expect_null(getIsovarMeta(x, require = FALSE))
  expect_error(getIsovarMeta(x, require = TRUE), "no .*isovar_meta")
})

test_that("metadata is attached by rankSplaireVariants, annotateGnomad, and loadCopdGwas", {
  r <- rankSplaireVariants(akr1a1_fixture(), gene = "AKR1A1")
  m <- getIsovarMeta(r)
  expect_equal(m$genome_build, "GRCh38")
  expect_true(any(vapply(m$sources,
                         function(s) s$kind %||% "" == "splicing_model_predictions",
                         logical(1))))

  g <- loadCopdGwas(file.path(.isovar_pkgroot(),
                              "inst/extdata/gwas_akr1a1_test.tsv.gz"))
  mg <- getIsovarMeta(g)
  expect_equal(mg$genome_build, "GRCh37")
  expect_equal(mg$sources[[1]]$kind, "gwas_sumstats")
})
