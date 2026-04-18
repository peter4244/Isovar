#' Extract structured splice-site objects from splaire rows
#'
#' Takes a splaire-scores table (raw, or the output of
#' [rankSplaireVariants()] joined back to the raw columns) and returns
#' a long-format tibble of implicated splice sites: one row per
#' (variant, head, event), giving the genomic coordinate of each
#' affected splice site and a classification of the event.
#'
#' ### Event classification
#'
#' For each variant and each splaire output head (donor `don`,
#' acceptor `acc`, splice-site usage `ssu`), the function examines the
#' paired `max_inc` / `max_dec` positions from either model family. It
#' classifies the event as:
#'
#' - **`motif_shift_loss`** / **`motif_shift_gain`** — paired large
#'   gain and loss at positions differing by exactly the indel length
#'   (within `motif_shift_tolerance` bp), indicating the splice motif
#'   has been displaced by the indel rather than destroyed or newly
#'   created. Loss and gain are emitted as two linked rows
#'   (`paired_site_pos` cross-references the companion).
#' - **`site_loss`** — large negative delta without a corresponding
#'   paired gain nearby. Suggests genuine splice-site disruption.
#' - **`site_gain`** — large positive delta without a corresponding
#'   paired loss nearby. Suggests cryptic / novel splice-site creation.
#' - **`weak`** — neither direction exceeds `magnitude_threshold`; a
#'   single summary row is still emitted so that filter logic has
#'   something to see (filter downstream with `magnitude >= threshold`).
#'
#' The classification uses the specified `model` family (default:
#' `splaireVar`, which was trained with variants and usually gives the
#' more informative variant-effect signal).
#'
#' @param sm_predictions Splicing-model variant-effect predictions —
#'   either a data.frame with the splaire schema or a file path readable
#'   by [rankSplaireVariants()]. See that function for the rationale on
#'   the generic argument name.
#' @param gene Optional HGNC gene symbol(s) to restrict to.
#' @param model Which model family to use for classification. Default
#'   `"splaireVar"`. `"splaire"` uses the reference-trained family.
#' @param heads Splice-site heads to evaluate. Default all three.
#' @param magnitude_threshold Minimum |delta| to call a site
#'   `site_loss` / `site_gain` / `motif_shift_*`. Default 0.1.
#' @param motif_shift_threshold Minimum |delta| that *both* the paired
#'   `max_inc` and `max_dec` must exceed to qualify as a motif shift.
#'   Default 0.5 — a relatively strict bar, since motif-shift
#'   interpretations drive Step 2 of the workflow.
#' @param motif_shift_tolerance Integer; allowed deviation (in bp)
#'   between `|max_inc_pos - max_dec_pos|` and the indel length for a
#'   motif-shift call. Default 1L. Use 0L to require exact match.
#' @return A tibble with one row per (variant, head) event; motif
#'   shifts contribute two rows (loss + gain) linked by
#'   `paired_site_pos`. Columns: `variant_id`, `chr`, `variant_pos`,
#'   `ref`, `alt`, `posterior_inclusion_probability`, `gene`,
#'   `phenotype`, `model`, `head`, `interpretation`, `site_pos`,
#'   `delta`, `magnitude`, `offset_from_variant`, `paired_site_pos`,
#'   `shift_bp`.
#' @examples
#' \dontrun{
#' sites <- implicatedSpliceSites(
#'   "HAEC185_sQTL_credible_sets_splaire_scores.tsv.gz",
#'   gene = "AKR1A1"
#' )
#' subset(sites, magnitude >= 0.5)
#' }
#' @export
implicatedSpliceSites <- function(sm_predictions,
                                  gene = NULL,
                                  model = c("splaireVar", "splaire"),
                                  heads = c("don", "acc", "ssu"),
                                  magnitude_threshold = 0.1,
                                  motif_shift_threshold = 0.5,
                                  motif_shift_tolerance = 1L) {
  model <- match.arg(model)
  heads <- match.arg(heads, several.ok = TRUE)
  stopifnot(
    is.numeric(magnitude_threshold), magnitude_threshold >= 0,
    is.numeric(motif_shift_threshold), motif_shift_threshold >= 0,
    is.numeric(motif_shift_tolerance), motif_shift_tolerance >= 0
  )
  motif_shift_tolerance <- as.integer(motif_shift_tolerance)

  df <- .load_splaire_table(sm_predictions)
  if (!is.null(gene)) df <- df[df$gene %in% gene, , drop = FALSE]
  if (nrow(df) == 0L) return(.empty_sites_df())

  parts <- .split_variant_id(df$variant_id)
  indel_len <- abs(nchar(parts$alt) - nchar(parts$ref))

  out <- vector("list", nrow(df) * length(heads))
  k <- 0L
  for (i in seq_len(nrow(df))) {
    base <- list(
      variant_id                        = df$variant_id[i],
      chr                               = parts$chr[i],
      variant_pos                       = parts$pos[i],
      ref                               = parts$ref[i],
      alt                               = parts$alt[i],
      posterior_inclusion_probability   = df$posterior_inclusion_probability[i],
      gene                              = df$gene[i],
      phenotype                         = df$phenotype[i],
      model                             = model
    )
    for (h in heads) {
      inc_v <- df[[sprintf("%s_%s_max_inc", model, h)]][i]
      inc_p <- df[[sprintf("%s_%s_max_inc_pos", model, h)]][i]
      inc_o <- df[[sprintf("%s_%s_max_inc_off", model, h)]][i]
      dec_v <- df[[sprintf("%s_%s_max_dec", model, h)]][i]
      dec_p <- df[[sprintf("%s_%s_max_dec_pos", model, h)]][i]
      dec_o <- df[[sprintf("%s_%s_max_dec_off", model, h)]][i]

      events <- .classify_head(
        inc_v = inc_v, inc_p = inc_p, inc_o = inc_o,
        dec_v = dec_v, dec_p = dec_p, dec_o = dec_o,
        indel_len             = indel_len[i],
        magnitude_threshold   = magnitude_threshold,
        motif_shift_threshold = motif_shift_threshold,
        motif_shift_tolerance = motif_shift_tolerance
      )
      for (ev in events) {
        k <- k + 1L
        out[[k]] <- c(base, list(head = h), ev)
      }
    }
  }
  length(out) <- k
  dplyr::bind_rows(lapply(out, tibble::as_tibble_row))
}

# -- internals ---------------------------------------------------------

.classify_head <- function(inc_v, inc_p, inc_o,
                           dec_v, dec_p, dec_o,
                           indel_len,
                           magnitude_threshold,
                           motif_shift_threshold,
                           motif_shift_tolerance) {
  # NaNs propagate through comparisons; treat as weak.
  inc_mag <- if (is.na(inc_v)) 0 else abs(inc_v)
  dec_mag <- if (is.na(dec_v)) 0 else abs(dec_v)

  is_indel <- !is.na(indel_len) && indel_len > 0L
  is_motif_shift <- FALSE
  if (is_indel &&
      inc_mag >= motif_shift_threshold &&
      dec_mag >= motif_shift_threshold &&
      !is.na(inc_p) && !is.na(dec_p)) {
    is_motif_shift <- abs(abs(inc_p - dec_p) - indel_len) <= motif_shift_tolerance
  }

  if (is_motif_shift) {
    return(list(
      list(interpretation = "motif_shift_loss",
           site_pos       = as.integer(dec_p),
           delta          = dec_v,
           magnitude      = dec_mag,
           offset_from_variant = as.integer(dec_o),
           paired_site_pos = as.integer(inc_p),
           shift_bp       = as.integer(indel_len)),
      list(interpretation = "motif_shift_gain",
           site_pos       = as.integer(inc_p),
           delta          = inc_v,
           magnitude      = inc_mag,
           offset_from_variant = as.integer(inc_o),
           paired_site_pos = as.integer(dec_p),
           shift_bp       = as.integer(indel_len))
    ))
  }

  # Not a motif shift: pick the stronger direction, emit one event.
  if (max(inc_mag, dec_mag) < magnitude_threshold) {
    # Weak — still record the stronger side for completeness.
    stronger <- if (inc_mag >= dec_mag) "inc" else "dec"
    return(list(.weak_row(inc_v, inc_p, inc_o, dec_v, dec_p, dec_o, stronger)))
  }

  if (inc_mag >= dec_mag) {
    return(list(list(
      interpretation        = "site_gain",
      site_pos              = as.integer(inc_p),
      delta                 = inc_v,
      magnitude             = inc_mag,
      offset_from_variant   = as.integer(inc_o),
      paired_site_pos       = NA_integer_,
      shift_bp              = NA_integer_
    )))
  }
  list(list(
    interpretation        = "site_loss",
    site_pos              = as.integer(dec_p),
    delta                 = dec_v,
    magnitude             = dec_mag,
    offset_from_variant   = as.integer(dec_o),
    paired_site_pos       = NA_integer_,
    shift_bp              = NA_integer_
  ))
}

.weak_row <- function(inc_v, inc_p, inc_o, dec_v, dec_p, dec_o, which) {
  if (which == "inc") {
    list(interpretation = "weak", site_pos = as.integer(inc_p),
         delta = inc_v, magnitude = if (is.na(inc_v)) 0 else abs(inc_v),
         offset_from_variant = as.integer(inc_o),
         paired_site_pos = NA_integer_, shift_bp = NA_integer_)
  } else {
    list(interpretation = "weak", site_pos = as.integer(dec_p),
         delta = dec_v, magnitude = if (is.na(dec_v)) 0 else abs(dec_v),
         offset_from_variant = as.integer(dec_o),
         paired_site_pos = NA_integer_, shift_bp = NA_integer_)
  }
}

.empty_sites_df <- function() {
  tibble::tibble(
    variant_id = character(), chr = character(), variant_pos = integer(),
    ref = character(), alt = character(),
    posterior_inclusion_probability = numeric(),
    gene = character(), phenotype = character(),
    model = character(), head = character(),
    interpretation = character(), site_pos = integer(),
    delta = numeric(), magnitude = numeric(),
    offset_from_variant = integer(), paired_site_pos = integer(),
    shift_bp = integer()
  )
}
