## Run tests with:
##   cd /path/to/isovar && Rscript tests/testthat.R
## or:
##   testthat::test_dir("tests/testthat")

suppressPackageStartupMessages({
  library(testthat)
  library(tibble)
  library(dplyr)
})

# Source all function files so tests don't rely on a package namespace.
isovar_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."))
tryCatch(setwd(isovar_root), error = function(e) invisible())
for (f in list.files("R", full.names = TRUE)) source(f)

# Expose the fixture path as a test-global.
fixture_akr1a1 <- "inst/extdata/akr1a1_splaire_test.tsv.gz"

test_dir("tests/testthat")
