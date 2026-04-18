akr1a1_ld_fixture <- function() {
  file.path(.isovar_pkgroot(), "tests/testthat/fixtures/ldlink_akr1a1_r2.rds")
}

skip_if_no_ld_fixture <- function() {
  p <- akr1a1_ld_fixture()
  if (!file.exists(p))
    testthat::skip(paste0("No LD fixture at ", p,
                          " — run scripts/build_haplotype_fixture.R with ",
                          "LDLINK_TOKEN set to generate it."))
}

test_that("groupHaplotypes emits suffixed multi-resolution columns", {
  skip_if_no_ld_fixture()
  flat <- as_flat(buildAnnotatedCredibleSet(
    sm_predictions   = akr1a1_fixture(),
    gwas             = file.path(.isovar_pkgroot(), "inst/extdata/gwas_akr1a1_test.tsv.gz"),
    gene             = "AKR1A1",
    gnomad_cache_dir = "/tmp/isovar_cache"
  ))
  out <- groupHaplotypes(flat, fixture = akr1a1_ld_fixture())

  # Default r² levels: 0.8, 0.9, 0.95 -> suffixes _r080, _r090, _r095
  for (sfx in c("_r080", "_r090", "_r095")) {
    expect_true(all(c(paste0("haplotype_id", sfx),
                      paste0("causal_candidate", sfx),
                      paste0("ld_r2_to_causal", sfx)) %in% names(out)))
  }
  expect_true("haplotype_assigned_by" %in% names(out))

  # All variants get a haplotype at every resolution
  expect_true(all(!is.na(out$haplotype_id_r080)))
  expect_true(all(!is.na(out$haplotype_id_r090)))
  expect_true(all(!is.na(out$haplotype_id_r095)))

  # Finer resolution -> >= haplotype count
  n08 <- length(unique(out$haplotype_id_r080))
  n09 <- length(unique(out$haplotype_id_r090))
  n95 <- length(unique(out$haplotype_id_r095))
  expect_lte(n08, n09)
  expect_lte(n09, n95)

  # Exactly one causal per haplotype per resolution
  for (sfx in c("_r080", "_r090", "_r095")) {
    per_hap <- table(out[[paste0("haplotype_id", sfx)]][out[[paste0("causal_candidate", sfx)]]])
    expect_true(all(per_hap == 1L))
  }
})

test_that("AKR1A1 indels: DEL LD-clustered, INS proximity-attached, separate at r² >= 0.9", {
  skip_if_no_ld_fixture()
  flat <- as_flat(buildAnnotatedCredibleSet(
    sm_predictions   = akr1a1_fixture(),
    gwas             = file.path(.isovar_pkgroot(), "inst/extdata/gwas_akr1a1_test.tsv.gz"),
    gene             = "AKR1A1",
    gnomad_cache_dir = "/tmp/isovar_cache"
  ))
  out <- groupHaplotypes(flat, fixture = akr1a1_ld_fixture())

  rs6692 <- out[out$rsid == "rs66922050", ]   # 10-bp DEL
  rs6146 <- out[out$rsid == "rs61467610", ]   # 6-bp INS

  # Assignment method is resolution-invariant
  expect_equal(rs6692$haplotype_assigned_by, "ld_cluster")
  expect_equal(rs6146$haplotype_assigned_by, "proximity")

  # At r² >= 0.8 the two indels merge; at r² >= 0.9 they separate.
  expect_equal(rs6692$haplotype_id_r080, rs6146$haplotype_id_r080)
  expect_true(rs6692$causal_candidate_r080)
  expect_false(rs6146$causal_candidate_r080)

  expect_false(identical(rs6692$haplotype_id_r090, rs6146$haplotype_id_r090))
  expect_true(rs6692$causal_candidate_r090)
  expect_true(rs6146$causal_candidate_r090)

  expect_false(identical(rs6692$haplotype_id_r095, rs6146$haplotype_id_r095))
  expect_true(rs6692$causal_candidate_r095)
  expect_true(rs6146$causal_candidate_r095)
})

test_that("groupHaplotypes works on nested credible-set input", {
  skip_if_no_ld_fixture()
  nested <- buildAnnotatedCredibleSet(
    sm_predictions   = akr1a1_fixture(),
    gwas             = file.path(.isovar_pkgroot(), "inst/extdata/gwas_akr1a1_test.tsv.gz"),
    gene             = "AKR1A1",
    gnomad_cache_dir = "/tmp/isovar_cache"
  )
  out <- groupHaplotypes(nested, fixture = akr1a1_ld_fixture())
  expect_true(is.list(out$variants))
  inner <- out$variants[[1]]
  expect_true("haplotype_id_r080" %in% names(inner))
  expect_true("haplotype_id_r090" %in% names(inner))
  expect_true("haplotype_id_r095" %in% names(inner))
  expect_true("causal_candidate_r090" %in% names(inner))
})

test_that("groupHaplotypes accepts a single r² threshold via r2_thresholds", {
  skip_if_no_ld_fixture()
  flat <- as_flat(buildAnnotatedCredibleSet(
    sm_predictions   = akr1a1_fixture(),
    gwas             = file.path(.isovar_pkgroot(), "inst/extdata/gwas_akr1a1_test.tsv.gz"),
    gene             = "AKR1A1",
    gnomad_cache_dir = "/tmp/isovar_cache"
  ))
  out <- groupHaplotypes(flat, r2_thresholds = 0.9,
                         fixture = akr1a1_ld_fixture())
  expect_true("haplotype_id_r090" %in% names(out))
  expect_false("haplotype_id_r080" %in% names(out))
})

test_that("missing token and missing fixture produce a clear error", {
  flat <- tibble::tibble(
    variant_id = c("chr1:100:A:G", "chr1:200:C:T"),
    rsid       = c("rs1", "rs2"),
    chr        = c("chr1", "chr1"),
    pos        = c(100L, 200L),
    ref        = c("A", "C"),
    alt        = c("G", "T"),
    top_abs_delta = c(0.5, 0.3)
  )
  # Point secrets_path at a nonexistent file so the auto-loader can't
  # pick up a real token from ~/.config/isovar/.
  empty_secrets <- tempfile(fileext = ".env")
  withr::with_envvar(c(LDLINK_TOKEN = ""), {
    expect_error(
      groupHaplotypes(flat, r2_thresholds = 0.9,
                      fixture = NULL, secrets_path = empty_secrets),
      "No LDlink token"
    )
  })
})

test_that(".cluster_ld produces 1 haplotype for perfectly linked variants", {
  m <- matrix(1, nrow = 3, ncol = 3,
              dimnames = list(c("a","b","c"), c("a","b","c")))
  haps <- .cluster_ld(m, r2_threshold = 0.8)
  expect_length(haps, 3L)
  expect_equal(length(unique(haps)), 1L)
})

test_that(".cluster_ld produces N haplotypes for independent variants", {
  m <- diag(3)
  rownames(m) <- colnames(m) <- c("a","b","c")
  haps <- .cluster_ld(m, r2_threshold = 0.8)
  expect_equal(length(unique(haps)), 3L)
})
