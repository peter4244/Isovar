test_that("motif_shift is detected for the AKR1A1 insertion (rs61467610)", {
  sites <- implicatedSpliceSites(akr1a1_fixture(), gene = "AKR1A1")
  ins <- sites[sites$variant_id == "chr1:45549863:A:ATATCTG", ]
  # Both donor motif shifts at the clu_33775 5' boundary:
  don <- ins[ins$head == "don" & startsWith(ins$interpretation, "motif_shift"), ]
  expect_equal(nrow(don), 2L)
  expect_setequal(don$site_pos, c(45551155L, 45551161L))
  expect_true(all(don$shift_bp == 6L))
  expect_true(all(don$magnitude > 0.9))
  # loss/gain must cross-reference each other
  loss <- don[don$interpretation == "motif_shift_loss", ]
  gain <- don[don$interpretation == "motif_shift_gain", ]
  expect_equal(loss$paired_site_pos, gain$site_pos)
  expect_equal(gain$paired_site_pos, loss$site_pos)

  # The acceptor boundary of the same intron should also be flagged:
  acc <- ins[ins$head == "acc" & startsWith(ins$interpretation, "motif_shift"), ]
  expect_equal(nrow(acc), 2L)
  expect_setequal(acc$site_pos, c(45552436L, 45552442L))
  expect_true(all(acc$shift_bp == 6L))
})

test_that("motif_shift is detected for the AKR1A1 deletion (rs66922050)", {
  sites <- implicatedSpliceSites(akr1a1_fixture(), gene = "AKR1A1")
  del <- sites[sites$variant_id == "chr1:45542866:TTTGAACTTCG:T" &
                 startsWith(sites$interpretation, "motif_shift"), ]
  expect_true(nrow(del) >= 4L)
  expect_true(all(del$shift_bp == 10L))
  # Donor shift at ~chr1:45547544/45547554
  don <- del[del$head == "don", ]
  expect_setequal(don$site_pos, c(45547544L, 45547554L))
})

test_that("weak events emit a row but magnitude is below threshold", {
  sites <- implicatedSpliceSites(akr1a1_fixture(), gene = "AKR1A1",
                                 magnitude_threshold = 0.1)
  weak <- sites[sites$interpretation == "weak", ]
  expect_true(all(weak$magnitude < 0.1 | is.na(weak$magnitude)))
})

test_that("model switch changes classification source", {
  s_ref <- implicatedSpliceSites(akr1a1_fixture(), gene = "AKR1A1",
                                 model = "splaire")
  s_var <- implicatedSpliceSites(akr1a1_fixture(), gene = "AKR1A1",
                                 model = "splaireVar")
  expect_true(all(s_ref$model == "splaire"))
  expect_true(all(s_var$model == "splaireVar"))
  # Shapes should match (same 50 variants × 3 heads, plus motif_shift doubles)
  expect_equal(sort(unique(s_ref$variant_id)),
               sort(unique(s_var$variant_id)))
})

test_that("motif_shift_tolerance=0 still catches the AKR1A1 indels exactly", {
  sites <- implicatedSpliceSites(akr1a1_fixture(), gene = "AKR1A1",
                                 motif_shift_tolerance = 0L)
  ms <- sites[startsWith(sites$interpretation, "motif_shift"), ]
  # Insertion (6bp) and deletion (10bp) both have exact offsets
  expect_true(all(ms$shift_bp %in% c(6L, 10L)))
})

test_that("SNP gets a site_gain classification for the known AKR1A1 hit", {
  sites <- implicatedSpliceSites(akr1a1_fixture(), gene = "AKR1A1",
                                 magnitude_threshold = 0.1)
  snp_hit <- sites[sites$variant_id == "chr1:45508256:G:A" &
                     sites$interpretation != "weak", ]
  # Single-row site_gain expected (donor)
  expect_equal(nrow(snp_hit), 1L)
  expect_equal(snp_hit$interpretation, "site_gain")
  expect_equal(snp_hit$head, "don")
})
