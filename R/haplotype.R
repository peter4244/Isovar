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
#' @param r2_threshold r² bar for grouping into the same haplotype
#'   (default 0.8). Variants with `r² >= threshold` cluster together.
#' @param token LDlink API token. Defaults to `Sys.getenv("LDLINK_TOKEN")`.
#'   Ignored when `fixture` is supplied.
#' @param fixture Optional: an RDS path, or a pre-loaded LDlinkR
#'   `LDmatrix` data.frame. When supplied, skips the live LDlink call.
#' @return The input augmented with:
#'   - `haplotype_id` (character, `"hap_1"`, `"hap_2"`, …)
#'   - `causal_candidate` (logical; TRUE for the variant with the
#'     largest `top_abs_delta` within each haplotype)
#'   - `ld_r2_to_causal` (numeric; r² of this variant to its
#'     haplotype's causal candidate; `NA` for proximity-attached rows)
#'   - `haplotype_assigned_by` (`"ld_cluster"` or `"proximity"`)
#' @export
groupHaplotypes <- function(variants,
                            population = "EUR",
                            r2_threshold = 0.8,
                            token = Sys.getenv("LDLINK_TOKEN"),
                            fixture = NULL) {
  is_nested <- "variants" %in% names(variants) && is.list(variants$variants)

  if (is_nested) {
    variants$variants <- lapply(
      variants$variants,
      function(v) .group_haplotypes_flat(v, population, r2_threshold, token, fixture)
    )
    return(variants)
  }
  .group_haplotypes_flat(variants, population, r2_threshold, token, fixture)
}

# -- internals ---------------------------------------------------------

.group_haplotypes_flat <- function(variants,
                                   population,
                                   r2_threshold,
                                   token,
                                   fixture) {
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

  ld_df  <- .get_ld_matrix(rsids, population, token, fixture)
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

.get_ld_matrix <- function(rsids, population, token, fixture) {
  if (!is.null(fixture)) {
    if (is.character(fixture) && length(fixture) == 1L) {
      if (!file.exists(fixture))
        cli::cli_abort("LD fixture not found: {.path {fixture}}")
      return(readRDS(fixture))
    }
    if (is.data.frame(fixture)) return(fixture)
    cli::cli_abort("{.arg fixture} must be an RDS path or a data.frame.")
  }

  if (!nzchar(token))
    cli::cli_abort(c(
      "No LDlink token available.",
      "i" = "Set env var {.envvar LDLINK_TOKEN} or pass a precomputed {.arg fixture}.",
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
