#' Rank credible-set variants by predicted splice-site effect
#'
#' Given a splaire / splaireVar score table (see
#' `HAEC185_sQTL_credible_sets_splaire_scores.tsv.gz` for the canonical
#' column schema), return one row per input variant-phenotype with the
#' single largest-magnitude predicted effect across both model families
#' (`splaire`, `splaireVar`) and all three splice-site output heads
#' (donor `don`, acceptor `acc`, splice-site usage `ssu`).
#'
#' Each row of the input has twelve candidate deltas:
#' `{splaire,splaireVar}_{don,acc,ssu}_max_{inc,dec}` (each accompanied
#' by an offset and genomic position). This function reduces those to
#' the single winning delta per variant by absolute value, and returns
#' enough context to reproduce the choice (model, head, direction,
#' signed delta, genomic position, offset from variant).
#'
#' For indel variants, note that paired `max_inc` and `max_dec` values
#' offset by exactly the indel length typically reflect a splice-motif
#' *shift* rather than genuine gain/loss — see
#' [implicatedSpliceSites()] for detection of that pattern.
#'
#' @param sm_predictions Splicing-model variant-effect predictions:
#'   either a path to a splaire-scores TSV (optionally gzipped), or a
#'   pre-loaded data.frame with the splaire column schema. Named
#'   generically because future work will feed this function from
#'   non-splaire models (SpliceAI, Pangolin, splaireVar-retrained) and
#'   from multiple tissues — the argument name should not hard-code the
#'   model. The schema is currently splaire-v1-specific; generalizing
#'   the schema (or routing different schemas to per-model parsers) is
#'   a planned extension.
#' @param gene Optional character vector of HGNC gene symbols to
#'   restrict to. `NULL` (default) keeps all genes.
#' @param magnitude_threshold Numeric in `[0, 1]`. Variants whose
#'   winning |delta| falls below this are still returned, but are
#'   flagged `FALSE` in the returned `high_impact` column. Default 0.1.
#' @param models Subset of `c("splaire", "splaireVar")` to consider.
#'   Default `"both"` (equivalent to both families).
#' @param heads Subset of `c("don", "acc", "ssu")` to consider.
#'   Default all three.
#' @return A tibble with one row per input variant-phenotype, columns:
#'   `gene`, `phenotype`, `ensembl_id`, `variant_id`, `chr`, `pos`,
#'   `ref`, `alt`, `strand`, `credible_set_number`,
#'   `posterior_inclusion_probability`, `top_model`, `top_head`,
#'   `top_direction` (`"inc"`/`"dec"`), `top_delta` (signed),
#'   `top_abs_delta`, `top_pos`, `top_off`, `high_impact` (logical).
#'   Results are sorted by `top_abs_delta` descending.
#' @examples
#' \dontrun{
#' r <- rankSplaireVariants(
#'   "HAEC185_sQTL_credible_sets_splaire_scores.tsv.gz",
#'   gene = "AKR1A1"
#' )
#' head(r, 10)
#' }
#' @export
rankSplaireVariants <- function(sm_predictions,
                                gene = NULL,
                                magnitude_threshold = 0.1,
                                models = c("splaire", "splaireVar"),
                                heads = c("don", "acc", "ssu")) {
  stopifnot(
    is.numeric(magnitude_threshold),
    length(magnitude_threshold) == 1L,
    magnitude_threshold >= 0, magnitude_threshold <= 1
  )
  models <- match.arg(models, several.ok = TRUE)
  heads  <- match.arg(heads,  several.ok = TRUE)

  df <- .load_splaire_table(sm_predictions)
  if (!is.null(gene)) {
    df <- df[df$gene %in% gene, , drop = FALSE]
    if (nrow(df) == 0L)
      cli::cli_warn("No rows match gene(s): {.val {gene}}")
  }

  # Build the 12 (model x head x direction) candidate columns and scan
  # each row for the largest |delta|. Column names follow the splaire
  # schema verbatim: e.g. "splaireVar_don_max_inc".
  combos <- expand.grid(
    model     = models,
    head      = heads,
    direction = c("inc", "dec"),
    stringsAsFactors = FALSE
  )
  value_cols <- sprintf("%s_%s_max_%s", combos$model, combos$head, combos$direction)
  pos_cols   <- sprintf("%s_pos", value_cols)
  off_cols   <- sprintf("%s_off", value_cols)

  missing_cols <- setdiff(c(value_cols, pos_cols, off_cols), names(df))
  if (length(missing_cols))
    cli::cli_abort("Input is missing expected splaire columns: {.val {missing_cols}}")

  vals <- as.matrix(df[, value_cols, drop = FALSE])
  abs_vals <- abs(vals)

  # Index of winning column per row; ties resolved by column order.
  winner_idx <- max.col(abs_vals, ties.method = "first")
  row_seq <- seq_len(nrow(df))

  top_delta     <- vals[cbind(row_seq, winner_idx)]
  top_abs_delta <- abs_vals[cbind(row_seq, winner_idx)]
  top_pos       <- as.matrix(df[, pos_cols])[cbind(row_seq, winner_idx)]
  top_off       <- as.matrix(df[, off_cols])[cbind(row_seq, winner_idx)]

  # variant_id convention: "chrN:pos:ref:alt"
  vid_parts <- .split_variant_id(df$variant_id)

  out <- tibble::tibble(
    gene                             = df$gene,
    phenotype                        = df$phenotype,
    ensembl_id                       = df$ensembl_id,
    variant_id                       = df$variant_id,
    chr                              = vid_parts$chr,
    pos                              = vid_parts$pos,
    ref                              = vid_parts$ref,
    alt                              = vid_parts$alt,
    strand                           = df$strand,
    credible_set_number              = df$credible_set_number,
    posterior_inclusion_probability  = df$posterior_inclusion_probability,
    top_model     = combos$model[winner_idx],
    top_head      = combos$head[winner_idx],
    top_direction = combos$direction[winner_idx],
    top_delta     = top_delta,
    top_abs_delta = top_abs_delta,
    top_pos       = top_pos,
    top_off       = top_off,
    high_impact   = top_abs_delta >= magnitude_threshold
  )
  out <- out[order(-out$top_abs_delta), , drop = FALSE]
  out
}

# -- internals ---------------------------------------------------------

.load_splaire_table <- function(sm_predictions) {
  if (is.data.frame(sm_predictions)) return(sm_predictions)
  if (!is.character(sm_predictions) || length(sm_predictions) != 1L)
    cli::cli_abort("{.arg sm_predictions} must be a data.frame or a single file path.")
  if (!file.exists(sm_predictions))
    cli::cli_abort("File not found: {.path {sm_predictions}}")
  readr::read_tsv(sm_predictions, show_col_types = FALSE, progress = FALSE)
}

.split_variant_id <- function(ids) {
  parts <- strsplit(ids, ":", fixed = TRUE)
  if (any(lengths(parts) != 4L))
    cli::cli_abort(
      "variant_id values must have the form {.val chr:pos:ref:alt}; got e.g. {.val {ids[lengths(parts) != 4L][1]}}"
    )
  m <- do.call(rbind, parts)
  list(
    chr = m[, 1],
    pos = as.integer(m[, 2]),
    ref = m[, 3],
    alt = m[, 4]
  )
}
