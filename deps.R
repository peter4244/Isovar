## isovar dependencies
##
## Source this file (or run check_deps()) to install anything missing.
## Required versions follow the "known good" set on the author's laptop
## (R 4.5.2, Bioconductor 3.20). Loosen only when we verify older versions work.

.isovar_deps <- list(
  cran = c("httr2", "cli", "cachem", "memoise", "readr", "dplyr",
           "tidyr", "stringr", "tibble", "rlang", "vctrs",
           "LDlinkR", "jsonlite", "withr"),
  bioc = c("Rsamtools", "VariantAnnotation", "GenomicRanges", "GenomeInfoDb",
           "rtracklayer", "Biostrings", "S4Vectors")
)

check_deps <- function(install_missing = FALSE) {
  missing_cran <- .isovar_deps$cran[!vapply(
    .isovar_deps$cran, requireNamespace, logical(1), quietly = TRUE
  )]
  missing_bioc <- .isovar_deps$bioc[!vapply(
    .isovar_deps$bioc, requireNamespace, logical(1), quietly = TRUE
  )]

  if (length(missing_cran) == 0L && length(missing_bioc) == 0L) {
    message("All isovar dependencies satisfied.")
    return(invisible(TRUE))
  }

  if (length(missing_cran)) {
    message("Missing CRAN packages: ", paste(missing_cran, collapse = ", "))
    if (install_missing) install.packages(missing_cran)
  }
  if (length(missing_bioc)) {
    message("Missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    if (install_missing) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
    }
  }
  invisible(length(missing_cran) == 0L && length(missing_bioc) == 0L)
}
