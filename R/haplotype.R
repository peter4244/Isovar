#' Group variants into haplotypes by LD structure
#'
#' Uses the LDlink REST API (via the `LDlinkR` CRAN wrapper) to compute
#' pairwise r² across the variants' rsIDs, then hierarchically clusters
#' the r² matrix to assign each variant to a haplotype. Within each
#' haplotype, the variant with the largest splaire |delta| is flagged
#' as the `causal_candidate`.
#'
#' ### Design notes
#'
#' LDlink is backed by 1000 Genomes Phase 3, so variants absent from
#' the 1KG panel — notably most indels — will not have LD data. For
#' those, the function falls back to **proximity assignment**: the
#' orphan variant inherits the haplotype of its nearest-position
#' credible-set neighbour. This is deliberately simple and should be
#' replaced with a local-cohort LD path when HAEC-185 genotypes are
#' available; the `haplotype_assigned_by` column surfaces which rows
#' were LD-clustered vs. proximity-attached.
#'
#' `LDlinkR` requires an API token (free from
#' <https://ldlink.nci.nih.gov>). Set it in the environment variable
#' `LDLINK_TOKEN`; the function errors cleanly if absent and no
#' `fixture` is provided.
#'
#' For CI and offline use, supply a pre-computed r² matrix via
#' `fixture` (an RDS path or the loaded data.frame). Build one with
#' `scripts/build_haplotype_fixture.R`.
#'
#' @param variants Either a flat per-variant tibble (from
#'   [buildAnnotatedCredibleSet()] → [as_flat()]) or a nested
#'   credible-set tibble. Must contain `rsid`, `pos`, and
#'   `top_abs_delta` columns. For nested input, haplotype grouping is
#'   performed independently within each credible set.
#' @param population LDlink population code (default `"EUR"`). Use
#'   [LDlinkR::list_pop()] to see options. Multi-population queries
#'   like `c("CEU","YRI","CHB")` are supported by LDlink.
#' @param r2_thresholds Numeric vector of r² bars for clustering.
#'   Default `c(0.8, 0.9, 0.95)` — a three-level resolution view from
#'   broad LD blocks (0.8) to fine structure (0.95). For a single
#'   threshold, pass a length-1 vector (e.g. `0.9`).
#' @param token LDlink API token. Defaults to `Sys.getenv("LDLINK_TOKEN")`.
#'   If empty, [loadIsovarSecrets()] is invoked on `secrets_path` and
#'   the env var is re-checked. Ignored when `fixture` is supplied.
#' @param fixture Optional: an RDS path, or a pre-loaded LDlinkR
#'   `LDmatrix` data.frame. When supplied, skips the live LDlink call.
#' @param secrets_path Path to the isovar secrets file for auto-load
#'   of `LDLINK_TOKEN`. Default `~/.config/isovar/secrets.env`. Pass a
#'   nonexistent path to disable auto-load (useful in tests).
#' @return The input augmented with, for each r² level in
#'   `r2_thresholds`, four columns suffixed by `_r{int(100*r²)}` (e.g.
#'   `_r080`, `_r090`, `_r095`):
#'   - `haplotype_id_r080` (character, `"hap_1"`, `"hap_2"`, …)
#'   - `causal_candidate_r080` (logical; TRUE for the variant with the
#'     largest `top_abs_delta` within each haplotype)
#'   - `ld_r2_to_causal_r080` (numeric; r² of this variant to its
#'     haplotype's causal candidate; `NA` for proximity-attached rows)
#'   And one resolution-invariant column:
#'   - `haplotype_assigned_by` (`"ld_cluster"` or `"proximity"`) —
#'     the same regardless of threshold, since it tracks whether the
#'     variant is in the LD matrix at all.
#' @export
groupHaplotypes <- function(variants,
                            population = "EUR",
                            r2_thresholds = c(0.8, 0.9, 0.95),
                            token = Sys.getenv("LDLINK_TOKEN"),
                            fixture = NULL,
                            secrets_path = "~/.config/isovar/secrets.env") {
  stopifnot(
    is.numeric(r2_thresholds),
    length(r2_thresholds) >= 1L,
    all(r2_thresholds > 0 & r2_thresholds <= 1)
  )
  is_nested <- "variants" %in% names(variants) && is.list(variants$variants)

  if (is_nested) {
    variants$variants <- lapply(
      variants$variants,
      function(v) .group_haplotypes_multi(v, population, r2_thresholds, token,
                                          fixture, secrets_path)
    )
    return(variants)
  }
  .group_haplotypes_multi(variants, population, r2_thresholds, token, fixture,
                          secrets_path)
}

.group_haplotypes_multi <- function(variants,
                                    population,
                                    r2_thresholds,
                                    token,
                                    fixture,
                                    secrets_path) {
  # Call the single-threshold function once per r² level, merging the
  # per-threshold columns. The LD matrix is fetched / read once (via
  # .get_ld_matrix caching) and the clustering is cheap, so multi-res
  # adds negligible cost.
  for (thresh in r2_thresholds) {
    one <- .group_haplotypes_flat(variants, population, thresh, token,
                                  fixture, secrets_path)
    suffix <- sprintf("_r%03d", as.integer(round(100 * thresh)))
    for (col in c("haplotype_id", "causal_candidate", "ld_r2_to_causal")) {
      variants[[paste0(col, suffix)]] <- one[[col]]
    }
    # haplotype_assigned_by is threshold-invariant; attach once
    if (!"haplotype_assigned_by" %in% names(variants))
      variants$haplotype_assigned_by <- one$haplotype_assigned_by
  }
  # Carry the metadata from the last call (any of them has the same sources)
  meta <- getIsovarMeta(one, require = FALSE)
  if (!is.null(meta)) {
    meta$sources[[length(meta$sources)]]$r2_thresholds <- r2_thresholds
    setIsovarMeta(variants, meta)
  } else {
    variants
  }
}

# -- internals ---------------------------------------------------------

.group_haplotypes_flat <- function(variants,
                                   population,
                                   r2_threshold,
                                   token,
                                   fixture,
                                   secrets_path) {
  required <- c("rsid", "pos", "top_abs_delta")
  missing <- setdiff(required, names(variants))
  if (length(missing))
    cli::cli_abort("{.arg variants} is missing columns: {.val {missing}}")

  variants$haplotype_id         <- NA_character_
  variants$causal_candidate     <- FALSE
  variants$ld_r2_to_causal      <- NA_real_
  variants$haplotype_assigned_by <- NA_character_

  rsids <- unique(variants$rsid[!is.na(variants$rsid)])
  if (length(rsids) < 2L) {
    cli::cli_warn("Fewer than 2 variants with rsIDs; haplotype grouping skipped.")
    return(variants)
  }

  ld_df  <- .get_ld_matrix(rsids, population, token, fixture, secrets_path)
  ld_mat <- .ld_df_to_matrix(ld_df)

  # Keep only variants that are actually in the LD matrix
  rsids_in_ld <- intersect(rownames(ld_mat), rsids)
  ld_sub <- ld_mat[rsids_in_ld, rsids_in_ld, drop = FALSE]

  clusters <- .cluster_ld(ld_sub, r2_threshold)   # named vec rsid -> "hap_1" etc.

  # Assign LD-cluster rows
  idx_ld <- which(variants$rsid %in% names(clusters))
  variants$haplotype_id[idx_ld]          <- clusters[variants$rsid[idx_ld]]
  variants$haplotype_assigned_by[idx_ld] <- "ld_cluster"

  # Proximity fallback for orphan rsIDs and variants without an rsID
  idx_orphan <- which(is.na(variants$haplotype_id))
  if (length(idx_orphan) > 0L && any(!is.na(variants$haplotype_id))) {
    known_rows <- which(!is.na(variants$haplotype_id))
    for (i in idx_orphan) {
      pos_i <- variants$pos[i]
      if (is.na(pos_i)) next
      dists <- abs(variants$pos[known_rows] - pos_i)
      nearest <- known_rows[which.min(dists)]
      variants$haplotype_id[i]           <- variants$haplotype_id[nearest]
      variants$haplotype_assigned_by[i]  <- "proximity"
    }
  }

  # Flag causal_candidate within each haplotype
  for (h in stats::na.omit(unique(variants$haplotype_id))) {
    rows <- which(variants$haplotype_id == h)
    if (length(rows) == 0L) next
    deltas <- variants$top_abs_delta[rows]
    if (all(is.na(deltas))) next
    winner <- rows[which.max(deltas)]
    variants$causal_candidate[winner] <- TRUE
  }

  # ld_r2_to_causal for LD-clustered rows
  for (h in stats::na.omit(unique(variants$haplotype_id))) {
    rows_h <- which(variants$haplotype_id == h)
    causal_row <- rows_h[variants$causal_candidate[rows_h]][1]
    if (is.na(causal_row)) next
    causal_rsid <- variants$rsid[causal_row]
    if (is.na(causal_rsid) || !causal_rsid %in% rownames(ld_mat)) next
    for (r in rows_h) {
      v_rsid <- variants$rsid[r]
      if (is.na(v_rsid) || !v_rsid %in% rownames(ld_mat)) next
      variants$ld_r2_to_causal[r] <- ld_mat[v_rsid, causal_rsid]
    }
  }

  existing <- getIsovarMeta(variants, require = FALSE)
  ld_meta <- list(
    sources = list(list(
      id           = "ldlink_1000g",
      kind         = "haplotype_ld",
      source       = "LDlinkR / LDlink REST API (1000 Genomes Phase 3)",
      population   = population,
      r2_threshold = r2_threshold
    ))
  )
  merged <- if (is.null(existing)) ld_meta else
    mergeIsovarMeta(existing, ld_meta)
  setIsovarMeta(variants, merged)
}

.get_ld_matrix <- function(rsids, population, token, fixture,
                           secrets_path = "~/.config/isovar/secrets.env") {
  if (!is.null(fixture)) {
    if (is.character(fixture) && length(fixture) == 1L) {
      if (!file.exists(fixture))
        cli::cli_abort("LD fixture not found: {.path {fixture}}")
      return(readRDS(fixture))
    }
    if (is.data.frame(fixture)) return(fixture)
    cli::cli_abort("{.arg fixture} must be an RDS path or a data.frame.")
  }

  # Auto-load from secrets file if env var isn't set.
  if (!nzchar(token)) {
    loadIsovarSecrets(path = secrets_path, quiet = TRUE)
    token <- Sys.getenv("LDLINK_TOKEN")
  }
  if (!nzchar(token))
    cli::cli_abort(c(
      "No LDlink token available.",
      "i" = "Set env var {.envvar LDLINK_TOKEN}, add {.val LDLINK_TOKEN=...} to {.path ~/.config/isovar/secrets.env}, or pass a precomputed {.arg fixture}.",
      "i" = "Register free at {.url https://ldlink.nci.nih.gov/?tab=apiaccess}."
    ))
  if (!requireNamespace("LDlinkR", quietly = TRUE))
    cli::cli_abort("Package {.pkg LDlinkR} is required.")

  cli::cli_inform(c("i" = "Calling LDlink LDmatrix: {length(rsids)} rsIDs, pop={.val {population}}"))
  LDlinkR::LDmatrix(
    snps     = rsids,
    pop      = population,
    r2d      = "r2",
    token    = token,
    genome_build = "grch38"
  )
}

.ld_df_to_matrix <- function(ld_df) {
  if (!is.data.frame(ld_df))
    cli::cli_abort("LD matrix must be a data.frame.")
  rs_col <- intersect(c("RS_number", "RSnumber", "RS", "rs"), names(ld_df))[1]
  if (is.na(rs_col))
    cli::cli_abort("Could not find an rsID column in the LD matrix.")
  rn <- as.character(ld_df[[rs_col]])
  m  <- as.matrix(ld_df[, setdiff(names(ld_df), rs_col), drop = FALSE])
  rownames(m) <- rn
  colnames(m) <- sub("^X", "", colnames(m))  # LDlinkR sometimes prefixes colnames
  storage.mode(m) <- "numeric"
  diag(m) <- 1
  keep <- intersect(rn, colnames(m))
  m[keep, keep, drop = FALSE]
}

#' Summarize haplotypes across resolutions with GWAS coherence metrics
#'
#' For each r² level in a multi-resolution [groupHaplotypes()] output,
#' summarize each haplotype cluster with GWAS consistency statistics.
#' The goal is to surface **which resolution produces the most
#' biologically coherent haplotypes** — ones whose variants share a
#' GWAS effect direction and magnitude.
#'
#' A haplotype is "cohesive" when its GWAS-matched variants agree in
#' β sign and cluster tightly in β magnitude. Splitting a real
#' haplotype too finely loses coherence (more clusters each with fewer
#' variants); merging too coarsely mixes variants with opposite effect
#' directions (high within-cluster β SD).
#'
#' The returned per-haplotype tibble is the long-format view. Use
#' [rollupHaplotypeResolutions()] to reduce it to one row per
#' resolution for threshold selection.
#'
#' @param x A flat per-variant tibble augmented by [groupHaplotypes()]
#'   (must contain the `_r0XX` suffix columns). For nested
#'   credible-set input, unnest with [as_flat()] first.
#' @param gwas_beta_col,gwas_p_col Column names for GWAS β and p.
#'   Defaults match [buildAnnotatedCredibleSet()] output.
#' @return A tibble with one row per (resolution, haplotype) and
#'   columns: `resolution` (e.g. `"r080"`), `r2_threshold`,
#'   `haplotype_id`, `n_variants`, `n_gwas_variants`,
#'   `causal_rsid`, `causal_abs_delta`, `mean_beta`, `sd_beta`,
#'   `sign_concordance` (fraction of GWAS variants whose β sign
#'   matches the causal / mode), `min_p`, `n_gws` (p < 5e-8).
#' @export
summarizeHaplotypesByResolution <- function(x,
                                            gwas_beta_col = "gwas_beta",
                                            gwas_p_col    = "gwas_p") {
  res_cols <- grep("^haplotype_id_r", names(x), value = TRUE)
  if (length(res_cols) == 0L)
    cli::cli_abort("Input has no {.field haplotype_id_r*} columns — did you call {.fn groupHaplotypes}?")

  parts <- vector("list", length(res_cols))
  for (i in seq_along(res_cols)) {
    hid_col    <- res_cols[i]
    sfx        <- sub("^haplotype_id", "", hid_col)   # e.g. "_r080"
    causal_col <- paste0("causal_candidate", sfx)
    resolution <- sub("^_", "", sfx)
    r2_thresh  <- as.numeric(sub("^r", "", resolution)) / 100

    df <- x
    df$.hap    <- df[[hid_col]]
    df$.causal <- df[[causal_col]]
    df$.beta   <- if (gwas_beta_col %in% names(df)) df[[gwas_beta_col]] else NA_real_
    df$.p      <- if (gwas_p_col    %in% names(df)) df[[gwas_p_col]]    else NA_real_

    # Adjust β for allele alignment (flip on 'swapped' rows if column present)
    if ("gwas_alignment_status" %in% names(df)) {
      flip <- !is.na(df$gwas_alignment_status) & df$gwas_alignment_status == "swapped"
      df$.beta[flip] <- -df$.beta[flip]
    }

    agg <- df %>%
      dplyr::group_by(haplotype_id = .data$.hap) %>%
      dplyr::summarise(
        n_variants       = dplyr::n(),
        n_gwas_variants  = sum(!is.na(.data$.beta)),
        causal_rsid      = .data$rsid[.data$.causal][1],
        causal_abs_delta = max(.data$top_abs_delta, na.rm = TRUE),
        mean_beta        = mean(.data$.beta, na.rm = TRUE),
        sd_beta          = stats::sd(.data$.beta, na.rm = TRUE),
        sign_concordance = {
          vals <- .data$.beta[!is.na(.data$.beta)]
          if (length(vals) == 0L) NA_real_
          else {
            mode_sign <- sign(sum(sign(vals)))
            if (mode_sign == 0) 0.5 else
              mean(sign(vals) == mode_sign)
          }
        },
        min_p            = suppressWarnings(min(.data$.p, na.rm = TRUE)),
        n_gws            = sum(!is.na(.data$.p) & .data$.p < 5e-8),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        resolution   = resolution,
        r2_threshold = r2_thresh,
        min_p        = ifelse(is.infinite(.data$min_p), NA_real_, .data$min_p)
      ) %>%
      dplyr::relocate(resolution, r2_threshold)
    parts[[i]] <- agg
  }
  out <- dplyr::bind_rows(parts)
  meta <- getIsovarMeta(x, require = FALSE)
  if (!is.null(meta)) attr(out, "isovar_meta") <- meta
  out
}

#' Roll haplotype summaries up to one row per resolution
#'
#' Takes the long per-haplotype summary from
#' [summarizeHaplotypesByResolution()] and reduces to one row per r²
#' level with aggregate cohesion metrics. Intended for **threshold
#' selection**: pick the resolution where haplotypes are most
#' internally consistent in GWAS signal without over-splitting.
#'
#' Metrics (weighted by haplotype size):
#' - `mean_within_hap_beta_sd` — lower is better (cohesive clusters)
#' - `mean_sign_concordance` — higher is better (consistent β signs
#'   within clusters). Bounded `[0.5, 1.0]`; 1.0 = perfect.
#' - `frac_significant_haps` — fraction of haplotypes with at least
#'   one GWAS-significant (p < 5e-8) variant.
#' - `n_haplotypes` — total haplotype count at that resolution.
#'
#' @param hap_summary Output of [summarizeHaplotypesByResolution()].
#' @return A tibble with one row per resolution.
#' @export
rollupHaplotypeResolutions <- function(hap_summary) {
  hap_summary %>%
    dplyr::group_by(resolution, r2_threshold) %>%
    dplyr::summarise(
      n_haplotypes           = dplyr::n(),
      n_gwas_significant_hap = sum(.data$n_gws > 0, na.rm = TRUE),
      frac_significant_haps  = sum(.data$n_gws > 0, na.rm = TRUE) / dplyr::n(),
      mean_within_hap_beta_sd = stats::weighted.mean(
        .data$sd_beta, w = .data$n_gwas_variants, na.rm = TRUE
      ),
      mean_sign_concordance  = stats::weighted.mean(
        .data$sign_concordance, w = .data$n_gwas_variants, na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(r2_threshold)
}

.cluster_ld <- function(ld_mat, r2_threshold) {
  if (nrow(ld_mat) == 0L) return(character())
  if (nrow(ld_mat) == 1L) {
    out <- "hap_1"
    names(out) <- rownames(ld_mat)
    return(out)
  }
  dist_mat <- as.dist(1 - ld_mat)
  hc <- stats::hclust(dist_mat, method = "average")
  clusters <- stats::cutree(hc, h = 1 - r2_threshold)
  hap_ids <- paste0("hap_", clusters)
  names(hap_ids) <- rownames(ld_mat)
  hap_ids
}
