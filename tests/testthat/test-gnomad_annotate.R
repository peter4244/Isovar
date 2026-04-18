test_that("annotateGnomad returns known rsIDs and AFs for AKR1A1 top hits", {
  skip_on_cran()
  skip_if_no_internet()
  cache <- tempfile("isovar_cache_")
  ann <- annotateGnomad(
    c("chr1:45549863:A:ATATCTG", "chr1:45542866:TTTGAACTTCG:T"),
    cache_dir = cache
  )
  expect_equal(nrow(ann), 2L)
  r1 <- ann[ann$variant_id == "chr1:45549863:A:ATATCTG", ]
  r2 <- ann[ann$variant_id == "chr1:45542866:TTTGAACTTCG:T", ]
  expect_equal(r1$rsid, "rs61467610")
  expect_equal(r2$rsid, "rs66922050")
  expect_equal(r1$filter, "PASS")
  # Common variants — AF between 0.4 and 0.6
  expect_true(r1$af > 0.4 && r1$af < 0.6)
  expect_true(r2$af > 0.4 && r2$af < 0.6)
  # Per-ancestry columns are populated
  expect_false(is.na(r1$af_nfe))
  expect_false(is.na(r1$af_eas))
})

test_that("input order is preserved", {
  skip_on_cran(); skip_if_no_internet()
  ids <- c("chr1:45542866:TTTGAACTTCG:T",
           "chr1:45549863:A:ATATCTG",
           "chr1:45508256:G:A")
  ann <- annotateGnomad(ids, cache_dir = tempfile("isovar_cache_"))
  expect_equal(ann$variant_id, ids)
})

test_that("a nonexistent variant yields NA annotation rather than an error", {
  skip_on_cran(); skip_if_no_internet()
  ids <- c("chr1:45549863:A:ATATCTG", "chr1:45549863:A:AAAAAAA")
  ann <- annotateGnomad(ids, cache_dir = tempfile("isovar_cache_"))
  expect_equal(nrow(ann), 2L)
  expect_false(is.na(ann$rsid[1]))
  expect_true(is.na(ann$rsid[2]))
})

test_that("cache_dir actually caches (second call is fast and offline-safe)", {
  skip_on_cran(); skip_if_no_internet()
  cache <- tempfile("isovar_cache_")
  ids <- c("chr1:45549863:A:ATATCTG", "chr1:45542866:TTTGAACTTCG:T")
  ann1 <- annotateGnomad(ids, cache_dir = cache)
  # Cache should contain one RDS file for the chr1 region.
  cache_files <- list.files(cache, pattern = "\\.rds$")
  expect_length(cache_files, 1L)
  # Second call should hit cache — time it (generous tolerance).
  t0 <- Sys.time()
  ann2 <- annotateGnomad(ids, cache_dir = cache)
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  expect_lt(dt, 3)
  expect_equal(ann1, ann2)
})

test_that("malformed variant_ids are rejected with a clear error", {
  expect_error(annotateGnomad("chr1_45549863_A_ATATCTG"),
               "chr:pos:ref:alt")
})

test_that("chr prefix normalization works (input without 'chr')", {
  skip_on_cran(); skip_if_no_internet()
  ann <- annotateGnomad("1:45549863:A:ATATCTG",
                        cache_dir = tempfile("isovar_cache_"))
  expect_equal(ann$rsid, "rs61467610")
  expect_equal(ann$chr, "chr1")
})
