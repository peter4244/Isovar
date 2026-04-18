test_that("rankSplaireVariants returns one row per input and ranks by |delta|", {
  r <- rankSplaireVariants(akr1a1_fixture(), gene = "AKR1A1")
  expect_s3_class(r, "tbl_df")
  expect_equal(nrow(r), 50L)
  expect_equal(unique(r$gene), "AKR1A1")
  expect_true(all(diff(r$top_abs_delta) <= 0))  # sorted descending
  expect_true(all(r$top_abs_delta >= 0))
  expect_setequal(r$top_direction, c("inc", "dec"))
  expect_setequal(r$top_head, c("don", "acc", "ssu"))
})

test_that("the two known high-impact AKR1A1 variants are flagged", {
  r <- rankSplaireVariants(akr1a1_fixture(), gene = "AKR1A1",
                           magnitude_threshold = 0.5)
  hi <- r[r$high_impact, ]
  expect_true("chr1:45549863:A:ATATCTG"       %in% hi$variant_id)
  expect_true("chr1:45542866:TTTGAACTTCG:T"   %in% hi$variant_id)
  # Both should be splaireVar donor-driven and have |delta| >= 0.97.
  top2 <- r[r$variant_id %in% c("chr1:45549863:A:ATATCTG",
                                "chr1:45542866:TTTGAACTTCG:T"), ]
  expect_true(all(top2$top_abs_delta > 0.97))
  expect_true(all(top2$top_head == "don"))
})

test_that("variant_id parsing into chr/pos/ref/alt is correct", {
  r <- rankSplaireVariants(akr1a1_fixture(), gene = "AKR1A1")
  expect_true(all(r$chr == "chr1"))
  expect_equal(r$pos[r$variant_id == "chr1:45549863:A:ATATCTG"], 45549863L)
  expect_equal(r$ref[r$variant_id == "chr1:45549863:A:ATATCTG"], "A")
  expect_equal(r$alt[r$variant_id == "chr1:45549863:A:ATATCTG"], "ATATCTG")
})

test_that("filtering by gene works and empty gene gives a warning + 0 rows", {
  expect_warning(
    r <- rankSplaireVariants(akr1a1_fixture(), gene = "NOT_A_GENE"),
    "No rows match gene"
  )
  expect_equal(nrow(r), 0L)
})

test_that("models argument restricts which family is considered", {
  r_all  <- rankSplaireVariants(akr1a1_fixture(), gene = "AKR1A1")
  r_ref  <- rankSplaireVariants(akr1a1_fixture(), gene = "AKR1A1",
                                models = "splaire")
  r_var  <- rankSplaireVariants(akr1a1_fixture(), gene = "AKR1A1",
                                models = "splaireVar")
  expect_true(all(r_ref$top_model == "splaire"))
  expect_true(all(r_var$top_model == "splaireVar"))
  # All-models result's winners must come from one of the two.
  expect_true(all(r_all$top_model %in% c("splaire", "splaireVar")))
})

test_that("heads argument restricts which splice-site head is considered", {
  r_don <- rankSplaireVariants(akr1a1_fixture(), gene = "AKR1A1",
                               heads = "don")
  expect_true(all(r_don$top_head == "don"))
})

test_that("errors on malformed input", {
  expect_error(rankSplaireVariants(123), "data.frame or a single file path")
  expect_error(rankSplaireVariants("/tmp/__nonexistent__.tsv.gz"),
               "File not found")
})
