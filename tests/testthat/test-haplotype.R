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

test_that("groupHaplotypes assigns haplotype_id and causal_candidate", {
  skip_if_no_ld_fixture()
  flat <- as_flat(buildAnnotatedCredibleSet(
    sm_predictions   = akr1a1_fixture(),
    gwas             = file.path(.isovar_pkgroot(), "inst/extdata/gwas_akr1a1_test.tsv.gz"),
    gene             = "AKR1A1",
    gnomad_cache_dir = "/tmp/isovar_cache"
  ))
  flat <- flat[, setdiff(names(flat), c("variants"))]     # drop list-col if present
  out <- groupHaplotypes(flat, fixture = akr1a1_ld_fixture())

  expect_true(all(c("haplotype_id", "causal_candidate",
                    "ld_r2_to_causal", "haplotype_assigned_by") %in% names(out)))
  # All variants get a haplotype (either clustered or proximity-attached)
  expect_true(all(!is.na(out$haplotype_id)))
  # Exactly one causal_candidate per haplotype
  per_hap <- table(out$haplotype_id[out$causal_candidate])
  expect_true(all(per_hap == 1L))
  # Expected number of haplotypes for AKR1A1 under r2 >= 0.8 (pending real data):
  # at least 2 (indel-tagged pair) and not more than a handful.
  n_hap <- length(unique(out$haplotype_id))
  expect_gte(n_hap, 2L)
  expect_lte(n_hap, 8L)
})

test_that("AKR1A1 indels are correctly placed by LD / proximity", {
  skip_if_no_ld_fixture()
  flat <- as_flat(buildAnnotatedCredibleSet(
    sm_predictions   = akr1a1_fixture(),
    gwas             = file.path(.isovar_pkgroot(), "inst/extdata/gwas_akr1a1_test.tsv.gz"),
    gene             = "AKR1A1",
    gnomad_cache_dir = "/tmp/isovar_cache"
  ))
  out <- groupHaplotypes(flat, fixture = akr1a1_ld_fixture())

  rs6692 <- out[out$rsid == "rs66922050", ]
  rs6146 <- out[out$rsid == "rs61467610", ]
  expect_equal(nrow(rs6692), 1L)
  expect_equal(nrow(rs6146), 1L)

  # rs66922050 (10-bp DEL) is in 1000 Genomes → LD-clustered.
  # rs61467610 (6-bp INS) is not → proximity-attached.
  expect_equal(rs6692$haplotype_assigned_by, "ld_cluster")
  expect_equal(rs6146$haplotype_assigned_by, "proximity")

  # At r² >= 0.8 EUR, AKR1A1 credible-set variants mostly collapse into one
  # haplotype block. Both indels land in the same cluster, so rs66922050
  # (larger |delta|) is its causal candidate; rs61467610 is NOT because
  # it's in the same haplotype. causal_candidate is per-haplotype-unique.
  expect_true(rs6692$causal_candidate)
  expect_equal(rs6692$haplotype_id, rs6146$haplotype_id)
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
  expect_true("haplotype_id" %in% names(inner))
  expect_true("causal_candidate" %in% names(inner))
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
      groupHaplotypes(flat, fixture = NULL, secrets_path = empty_secrets),
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
