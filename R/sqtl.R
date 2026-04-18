#' Load Leafcutter sQTL nominal-pass results (tensorqtl `cis_qtl_pairs`)
#'
#' Streams a tensorqtl `*.cis_qtl_pairs.*.parquet` with column-predicate
#' pushdown via `arrow`. Filter by `variant_ids`, `phenotype_ids`, or
#' a coordinate range — each pushes down to the parquet reader so
#' multi-GB files don't get materialized unnecessarily.
#'
#' See `inst/sqtl_configs/haec185_basal_leafcutter.R` for the canonical
#' HAEC-185 config. Attach an `ISOVAR_SOURCE_META` list to the returned
#' tibble via [setIsovarMeta()] using the source config's metadata.
#'
#' @param source Path to a parquet file, OR path to an isovar sQTL
#'   config file (`inst/sqtl_configs/*.R`) that defines `SQTL_PARQUET`
#'   (a named list of chr → parquet path) and `ISOVAR_SOURCE_META`.
#' @param chr Chromosome string (e.g. `"chr1"` or `"1"`). Required
#'   when `source` is a config file (to pick the right parquet).
#' @param variant_ids Optional character vector of `chr:pos:ref:alt`
#'   variant IDs. Filters applied as parquet predicate.
#' @param phenotype_ids Optional character vector of Leafcutter intron
#'   phenotype IDs. Strand suffix is stripped for matching.
#' @param columns Which columns to keep. Default returns the standard
#'   9 cis_qtl_pairs fields.
#' @return A tibble of filtered sQTL rows with an `isovar_meta`
#'   attribute populated from the config's `ISOVAR_SOURCE_META`.
#' @export
loadSqtlResults <- function(source,
                            chr = NULL,
                            variant_ids  = NULL,
                            phenotype_ids = NULL,
                            columns = c("phenotype_id", "variant_id",
                                        "tss_distance", "af",
                                        "ma_samples", "ma_count",
                                        "pval_nominal", "slope", "slope_se")) {
  if (!requireNamespace("arrow", quietly = TRUE))
    cli::cli_abort("Package {.pkg arrow} is required to read tensorqtl parquet files.")

  meta <- list()
  if (is.character(source) && length(source) == 1L &&
      grepl("\\.R$", source, ignore.case = TRUE)) {
    env <- new.env(parent = baseenv())
    tryCatch(sys.source(source, envir = env),
             error = function(e)
               cli::cli_abort("Failed to source sQTL config {.path {source}}: {e$message}"))
    if (!exists("SQTL_PARQUET", envir = env, inherits = FALSE))
      cli::cli_abort("Config {.path {source}} missing {.field SQTL_PARQUET}.")
    parquet_map <- get("SQTL_PARQUET", envir = env)
    if (is.null(chr))
      cli::cli_abort("Config-based {.arg source} requires a {.arg chr} argument.")
    chr_norm <- if (startsWith(chr, "chr")) chr else paste0("chr", chr)
    parquet_path <- parquet_map[[chr_norm]] %||% parquet_map[[sub("^chr", "", chr_norm)]]
    if (is.null(parquet_path))
      cli::cli_abort("No parquet mapped for {.val {chr}} in {.path {source}}.")
    if (exists("ISOVAR_SOURCE_META", envir = env, inherits = FALSE))
      meta <- get("ISOVAR_SOURCE_META", envir = env)
  } else {
    parquet_path <- source
  }

  if (!file.exists(parquet_path))
    cli::cli_abort("sQTL parquet not found: {.path {parquet_path}}")

  ds <- arrow::open_dataset(parquet_path, format = "parquet")

  pheno_norm_requested <- if (!is.null(phenotype_ids))
    .strip_phenotype_strand(phenotype_ids) else NULL

  query <- ds
  if (!is.null(variant_ids))
    query <- dplyr::filter(query, variant_id %in% variant_ids)
  if (!is.null(pheno_norm_requested)) {
    # arrow doesn't support regex filter-pushdown; expand to both strand
    # suffixes to catch the leafcutter convention mismatch (splaire uses
    # _+ for a + strand gene, sQTL may have _- for the same intron).
    expanded <- c(paste0(pheno_norm_requested, "_+"),
                  paste0(pheno_norm_requested, "_-"))
    query <- dplyr::filter(query, phenotype_id %in% expanded)
  }

  df <- as.data.frame(
    dplyr::collect(dplyr::select(query, dplyr::all_of(columns)))
  )

  out <- tibble::as_tibble(df)
  out$phenotype_id_norm <- .strip_phenotype_strand(out$phenotype_id)

  # Propagate the source's metadata
  source_record <- list(
    id           = meta$source_id %||% "haec185_leafcutter_sqtl",
    kind         = "short_read_sqtl",
    label        = meta$source_label %||% NA_character_,
    cohort       = meta$cohort %||% NA_character_,
    tissue       = meta$tissue %||% NA_character_,
    n_samples    = meta$n_samples %||% NA_integer_,
    genome_build = meta$genome_build %||% "GRCh38",
    tool         = meta$tool %||% "tensorqtl cis nominal",
    includes_indels = meta$includes_indels %||% NA,
    notes        = meta$notes %||% NA_character_,
    parquet_path = normalizePath(parquet_path, mustWork = FALSE)
  )
  setIsovarMeta(out,
                genome_build = meta$genome_build %||% "GRCh38",
                sources = list(source_record))
}

#' Annotate variants with per-intron sQTL effect sizes
#'
#' Left-joins sQTL nominal-pass results onto a variant table, producing
#' one row per (variant, intron) with effect size, SE, p-value, and
#' the source's n_samples. Strand-suffix on `phenotype_id` is stripped
#' for matching (see the sQTL vs splaire strand-mismatch note in
#' `docs/schemas.md`).
#'
#' @param variants A tibble with at minimum `variant_id`.
#' @param sqtl Output of [loadSqtlResults()] (must carry
#'   `phenotype_id_norm`).
#' @return A long tibble with columns `variant_id`,
#'   `phenotype_id_norm`, `sqtl_slope`, `sqtl_slope_se`, `sqtl_p`,
#'   `sqtl_af`, `sqtl_tss_distance`. Metadata is merged from both inputs.
#' @export
annotateSqtl <- function(variants, sqtl) {
  if (!"variant_id" %in% names(variants))
    cli::cli_abort("{.arg variants} must have a {.field variant_id} column.")
  if (!"phenotype_id_norm" %in% names(sqtl))
    cli::cli_abort("{.arg sqtl} must come from {.fn loadSqtlResults}.")

  sq <- sqtl %>%
    dplyr::select(variant_id,
                  phenotype_id,
                  phenotype_id_norm,
                  sqtl_slope       = slope,
                  sqtl_slope_se    = slope_se,
                  sqtl_p           = pval_nominal,
                  sqtl_af          = af,
                  sqtl_tss_distance = tss_distance)
  joined <- dplyr::left_join(variants, sq, by = "variant_id", relationship = "many-to-many")
  setIsovarMeta(joined, mergeIsovarMeta(variants, sqtl))
}

# -- internals ---------------------------------------------------------

.strip_phenotype_strand <- function(ids) {
  # Leafcutter phenotype_id: "chr:start:end:clu_N_+"  or  "..._-"
  sub("(_[+-])$", "", ids)
}
