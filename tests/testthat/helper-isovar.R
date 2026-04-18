## testthat helpers — loaded automatically by testthat::test_dir().
## Source the R/ functions so tests work without installing a package.

.isovar_pkgroot <- function() {
  # testthat runs helpers with cwd at the tests/testthat directory.
  normalizePath(file.path("..", ".."))
}

for (f in list.files(file.path(.isovar_pkgroot(), "R"), full.names = TRUE)) {
  source(f, local = FALSE)
}

akr1a1_fixture <- function() {
  file.path(.isovar_pkgroot(), "inst/extdata/akr1a1_splaire_test.tsv.gz")
}

skip_if_no_internet <- function() {
  testthat::skip_if_offline("gnomad-public-us-east-1.s3.amazonaws.com")
}
