#' Load COPD GWAS summary statistics
#'
#' Reads a tab-separated GWAS sumstats file with the icgcUkb schema:
#' `marker chr pos rsid effectallele otherallele eaf beta se p z imputersq
#' n pheno mode tag`. Gzipped input is handled transparently by `readr`.
#'
#' The canonical COPD GWAS for the isovar workflow is
#' `icgcUkb-20251213.dn8.gz` at
#' `/Users/petecastaldi/claude_projects/copd/gwas/` — **on GRCh37**, which
#' is why the isovar join strategy is rsID-based rather than by
#' coordinate (HAEC-185 + gnomAD v4.1 + splaire are GRCh38).
#'
#' @param path Path to a sumstats file (`.tsv` or `.tsv.gz`).
#' @param cols Optional character vector restricting which columns to load
#'   (all are read by default). Large full-GWAS reads benefit from
#'   trimming.
#' @return A tibble with the GWAS columns, typed sensibly.
#' @export
loadCopdGwas <- function(path, cols = NULL) {
  stopifnot(is.character(path), length(path) == 1L)
  if (!file.exists(path))
    cli::cli_abort("GWAS file not found: {.path {path}}")

  types <- readr::cols(
    marker = readr::col_character(),
    chr = readr::col_character(),
    pos = readr::col_integer(),
    rsid = readr::col_character(),
    effectallele = readr::col_character(),
    otherallele = readr::col_character(),
    eaf = readr::col_double(),
    beta = readr::col_double(),
    se = readr::col_double(),
    p = readr::col_double(),
    z = readr::col_double(),
    imputersq = readr::col_double(),
    n = readr::col_integer(),
    pheno = readr::col_character(),
    mode = readr::col_character(),
    tag = readr::col_character()
  )
  out <- readr::read_tsv(path, col_types = types, progress = FALSE)
  if (!is.null(cols)) {
    missing <- setdiff(cols, names(out))
    if (length(missing))
      cli::cli_abort("Requested columns not present: {.val {missing}}")
    out <- out[, cols, drop = FALSE]
  }
  out
}

#' Annotate variants with COPD GWAS summary statistics (rsID join)
#'
#' Attaches GWAS effect estimates to a variant table via left-join on
#' rsID. rsID is build-independent, sidestepping the GRCh37/38 mismatch
#' between the icgcUkb GWAS and the HAEC-185 / gnomAD / splaire layers.
#' Coordinate-based liftover is intentionally **not** implemented here
#' — add only if rsID coverage proves insufficient.
#'
#' Variants without a matching rsID receive `NA` for the GWAS fields.
#' Allele alignment (effect allele orientation) is checked with
#' [checkAlleleAlignment()] and surfaced as a column — beta is **never
#' silently flipped**.
#'
#' @param variants A tibble with at minimum `rsid`, `ref`, `alt`. (Output
#'   of [annotateGnomad()] + [rankSplaireVariants()] already has these.)
#' @param gwas A tibble returned by [loadCopdGwas()].
#' @return `variants` augmented with:
#'   `gwas_effectallele`, `gwas_otherallele`, `gwas_eaf`,
#'   `gwas_beta`, `gwas_se`, `gwas_p`, `gwas_z`, `gwas_n`,
#'   `gwas_imputersq`, `gwas_marker`, `gwas_alignment_status` (one of
#'   `"aligned"`, `"swapped"`, `"strand_flip"`, `"mismatch"`, `NA`),
#'   `gwas_allele_alignment_ok` (logical).
#' @export
annotateGwas <- function(variants, gwas) {
  required <- c("rsid", "ref", "alt")
  missing <- setdiff(required, names(variants))
  if (length(missing))
    cli::cli_abort("{.arg variants} is missing columns: {.val {missing}}")
  gwas_required <- c("rsid", "effectallele", "otherallele", "eaf",
                     "beta", "se", "p")
  missing_gwas <- setdiff(gwas_required, names(gwas))
  if (length(missing_gwas))
    cli::cli_abort("{.arg gwas} is missing columns: {.val {missing_gwas}}")

  # Reduce GWAS to needed columns, rename to gwas_* to avoid clashes.
  g <- gwas[, intersect(names(gwas),
                        c("rsid", "marker", "effectallele", "otherallele",
                          "eaf", "beta", "se", "p", "z", "imputersq", "n")),
            drop = FALSE]
  names(g)[names(g) != "rsid"] <- paste0("gwas_", names(g)[names(g) != "rsid"])

  joined <- dplyr::left_join(variants, g, by = "rsid")

  align <- checkAlleleAlignment(joined,
                                ref_col = "ref",
                                alt_col = "alt",
                                effect_col = "gwas_effectallele",
                                other_col = "gwas_otherallele")
  joined$gwas_alignment_status <- align$alignment_status
  joined$gwas_allele_alignment_ok <- align$alignment_ok
  joined
}

#' Check allele alignment between a variant table and a GWAS join
#'
#' Classifies each row's ref/alt orientation against the GWAS effect /
#' other alleles. Does **not** flip beta; this is a diagnostic, and
#' downstream code should act on `alignment_status` explicitly.
#'
#' Classification:
#' - `"aligned"` — `alt == effectallele && ref == otherallele` (the
#'   common, expected case; beta refers to the alt allele)
#' - `"swapped"` — `ref == effectallele && alt == otherallele` (beta
#'   is on the reference allele; flip sign if using)
#' - `"strand_flip"` — alleles match after complementing both sides
#'   (either aligned-flipped or swapped-flipped)
#' - `"mismatch"` — none of the above; possibly a multi-allelic
#'   ambiguity or a true allele discrepancy. Do not use.
#' - `NA` — GWAS fields are `NA` (no match).
#'
#' `alignment_ok` is `TRUE` only for `"aligned"`.
#'
#' @param df A tibble containing the columns passed by name below.
#' @param ref_col,alt_col,effect_col,other_col Column names.
#' @return A tibble with `alignment_status` (character) and
#'   `alignment_ok` (logical), one row per input.
#' @export
checkAlleleAlignment <- function(df,
                                 ref_col = "ref",
                                 alt_col = "alt",
                                 effect_col = "gwas_effectallele",
                                 other_col = "gwas_otherallele") {
  ref    <- df[[ref_col]]
  alt    <- df[[alt_col]]
  effect <- df[[effect_col]]
  other  <- df[[other_col]]

  status <- rep(NA_character_, nrow(df))
  has_gwas <- !is.na(effect) & !is.na(other)

  aligned <- has_gwas & ref == other   & alt == effect
  swapped <- has_gwas & ref == effect & alt == other
  status[aligned] <- "aligned"
  status[swapped & is.na(status)] <- "swapped"

  # Strand-flip check applies only to biallelic SNPs (single base each).
  flip_candidates <- has_gwas & is.na(status) &
                       nchar(ref) == 1L & nchar(alt) == 1L &
                       nchar(effect) == 1L & nchar(other) == 1L
  if (any(flip_candidates, na.rm = TRUE)) {
    comp <- c(A = "T", T = "A", C = "G", G = "C")
    ref_c    <- comp[ref]
    alt_c    <- comp[alt]
    aligned_flip <- flip_candidates & !is.na(ref_c) & !is.na(alt_c) &
                      ref_c == other & alt_c == effect
    swapped_flip <- flip_candidates & !is.na(ref_c) & !is.na(alt_c) &
                      ref_c == effect & alt_c == other
    status[aligned_flip | swapped_flip] <- "strand_flip"
  }

  status[has_gwas & is.na(status)] <- "mismatch"

  tibble::tibble(
    alignment_status = status,
    alignment_ok     = !is.na(status) & status == "aligned"
  )
}
