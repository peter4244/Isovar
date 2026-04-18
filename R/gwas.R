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
  setIsovarMeta(
    out,
    genome_build = "GRCh37",
    sources = list(list(
      id           = "copd_gwas",
      kind         = "gwas_sumstats",
      path         = normalizePath(path, mustWork = FALSE),
      phenotype    = if ("pheno" %in% names(out) && nrow(out) > 0L)
                       as.character(out$pheno[1]) else NA_character_,
      study_tag    = if ("tag" %in% names(out) && nrow(out) > 0L)
                       as.character(out$tag[1]) else NA_character_,
      genome_build = "GRCh37"
    ))
  )
}

#' Annotate variants with COPD GWAS summary statistics
#'
#' Attaches GWAS effect estimates to a variant table via a two-pass
#' match: first by rsID (build-independent), then — for any rows that
#' missed the rsID join — by genomic position after optional liftover
#' to the GWAS genome build. This catches variants whose rsIDs drift
#' across annotation releases or whose representation differs between
#' databases, without silently relying on a single key.
#'
#' rsID is the primary join key because it is build-independent and
#' robust to the GRCh37/GRCh38 mismatch between the icgcUkb GWAS and
#' HAEC-185 / gnomAD v4.1 / splaire. Position matching is a secondary
#' check — it requires a liftover chain when the builds differ.
#'
#' Allele alignment (effect allele orientation) is checked with
#' [checkAlleleAlignment()] and surfaced as a column — beta is **never
#' silently flipped**.
#'
#' @param variants A tibble with at minimum `rsid`, `chr`, `pos`,
#'   `ref`, `alt`. (Output of [rankSplaireVariants()] +
#'   [annotateGnomad()] has all of these.)
#' @param gwas A tibble returned by [loadCopdGwas()]. Must carry `chr`,
#'   `pos`, `rsid`, `effectallele`, `otherallele`.
#' @param liftover_chain Optional; either a `Chain` object from
#'   [rtracklayer::import.chain()] or a path to one, mapping
#'   `variants`' build to `gwas`' build (e.g. hg38 → hg19 for the
#'   icgcUkb GWAS). If `NULL`, position-match is done only when both
#'   sides are already on the same build. Use
#'   [getLiftoverChain()] to fetch and cache a UCSC chain.
#' @param variants_build,gwas_build Character labels, only used to
#'   decide whether to invoke `liftover_chain`. Default `"GRCh38"` and
#'   `"GRCh37"` matches our canonical setup.
#' @return `variants` augmented with:
#'   `gwas_effectallele`, `gwas_otherallele`, `gwas_eaf`,
#'   `gwas_beta`, `gwas_se`, `gwas_p`, `gwas_z`, `gwas_n`,
#'   `gwas_imputersq`, `gwas_marker`,
#'   `gwas_matched_by` (one of `"rsid"`, `"position"`, `NA`),
#'   `gwas_alignment_status`, `gwas_allele_alignment_ok`.
#' @export
annotateGwas <- function(variants, gwas,
                         liftover_chain = NULL,
                         variants_build = "GRCh38",
                         gwas_build = "GRCh37") {
  required <- c("rsid", "chr", "pos", "ref", "alt")
  missing <- setdiff(required, names(variants))
  if (length(missing))
    cli::cli_abort("{.arg variants} is missing columns: {.val {missing}}")
  gwas_required <- c("rsid", "chr", "pos", "effectallele", "otherallele",
                     "eaf", "beta", "se", "p")
  missing_gwas <- setdiff(gwas_required, names(gwas))
  if (length(missing_gwas))
    cli::cli_abort("{.arg gwas} is missing columns: {.val {missing_gwas}}")

  gwas_cols <- intersect(
    names(gwas),
    c("rsid", "chr", "pos", "marker", "effectallele", "otherallele",
      "eaf", "beta", "se", "p", "z", "imputersq", "n")
  )
  g <- gwas[, gwas_cols, drop = FALSE]
  g_named <- g
  names(g_named)[names(g_named) != "rsid"] <-
    paste0("gwas_", names(g_named)[names(g_named) != "rsid"])

  # --- pass 1: rsID join --------------------------------------------
  joined <- dplyr::left_join(
    variants,
    g_named[, c("rsid", setdiff(names(g_named), "rsid")), drop = FALSE],
    by = "rsid"
  )
  joined$gwas_matched_by <- ifelse(!is.na(joined$gwas_beta), "rsid", NA_character_)

  # --- pass 2: position-based fallback ------------------------------
  unmatched_rows <- which(is.na(joined$gwas_matched_by))
  if (length(unmatched_rows) > 0L) {
    lifted <- .lift_positions(
      chr   = joined$chr[unmatched_rows],
      pos   = joined$pos[unmatched_rows],
      chain = liftover_chain,
      from  = variants_build,
      to    = gwas_build
    )
    norm_chr_v <- sub("^chr", "", lifted$chr)
    norm_chr_g <- sub("^chr", "", g$chr)

    for (k in seq_along(unmatched_rows)) {
      i <- unmatched_rows[k]
      v_chr <- norm_chr_v[k]
      v_pos <- lifted$pos[k]
      if (is.na(v_pos) || is.na(v_chr)) next

      hit <- which(norm_chr_g == v_chr & g$pos == v_pos)
      if (length(hit) == 0L) next

      aln_at_pos <- checkAlleleAlignment(
        data.frame(
          ref = joined$ref[i],
          alt = joined$alt[i],
          gwas_effectallele = g$effectallele[hit],
          gwas_otherallele  = g$otherallele[hit],
          stringsAsFactors = FALSE
        )
      )
      ok <- !is.na(aln_at_pos$alignment_status) &
              aln_at_pos$alignment_status %in% c("aligned", "swapped", "strand_flip")
      if (!any(ok)) next

      h <- hit[which(ok)[1]]
      for (nm in setdiff(names(g_named), "rsid"))
        joined[[nm]][i] <- g_named[[nm]][h]
      joined$gwas_matched_by[i] <- "position"
    }
  }

  align <- checkAlleleAlignment(joined,
                                ref_col = "ref",
                                alt_col = "alt",
                                effect_col = "gwas_effectallele",
                                other_col = "gwas_otherallele")
  joined$gwas_alignment_status <- align$alignment_status
  joined$gwas_allele_alignment_ok <- align$alignment_ok
  setIsovarMeta(joined, mergeIsovarMeta(variants, gwas))
}

# Lift GRCh38 -> GRCh37 positions using a chain. Returns a list with
# chr + pos (same length as input; NA where unmappable).
.lift_positions <- function(chr, pos, chain, from, to) {
  n <- length(chr)
  if (identical(from, to) || is.null(chain)) {
    return(list(chr = chr, pos = pos))
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE))
    cli::cli_abort("rtracklayer + GenomicRanges are required for liftover.")

  if (is.character(chain)) {
    if (!file.exists(chain))
      cli::cli_abort("Chain file not found: {.path {chain}}")
    chain <- rtracklayer::import.chain(chain)
  }

  # Add "chr" prefix if missing (UCSC-style chains expect it)
  chr_q <- ifelse(startsWith(chr, "chr"), chr, paste0("chr", chr))
  gr <- GenomicRanges::GRanges(chr_q, IRanges::IRanges(start = pos, width = 1))
  lifted <- rtracklayer::liftOver(gr, chain)

  out_chr <- rep(NA_character_, n)
  out_pos <- rep(NA_integer_, n)
  lens <- lengths(lifted)
  ok <- lens > 0L
  if (any(ok)) {
    # Use the first lifted range per input
    first <- unlist(lifted[ok])[cumsum(lens[ok]) - (lens[ok] - 1L)]
    out_chr[ok] <- as.character(GenomicRanges::seqnames(first))
    out_pos[ok] <- as.integer(GenomicRanges::start(first))
  }
  list(chr = out_chr, pos = out_pos)
}

#' Fetch (and cache) a UCSC liftOver chain file
#'
#' Downloads a UCSC chain file on first call and caches it under
#' `cache_dir`. Subsequent calls read from cache. Returns an absolute
#' path to the uncompressed chain file suitable for
#' [rtracklayer::import.chain()] or for passing as `liftover_chain` to
#' [annotateGwas()].
#'
#' @param from,to Assembly labels as they appear in the UCSC chain
#'   filename, e.g. `from = "hg38"`, `to = "hg19"`.
#' @param cache_dir Directory for cached chain files. Default is a
#'   user-level cache path.
#' @return Absolute path to the uncompressed chain file.
#' @export
getLiftoverChain <- function(from = "hg38", to = "hg19",
                             cache_dir = file.path("~", ".cache", "isovar", "liftover")) {
  cache_dir <- path.expand(cache_dir)
  if (!dir.exists(cache_dir))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  fname_base <- sprintf("%sTo%s.over.chain",
                        from, paste0(toupper(substr(to, 1, 1)),
                                     substr(to, 2, nchar(to))))
  local <- file.path(cache_dir, fname_base)
  gzlocal <- paste0(local, ".gz")
  if (file.exists(local)) return(local)

  if (!file.exists(gzlocal)) {
    url <- sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/%s/liftOver/%s.gz",
                   from, fname_base)
    cli::cli_inform(c("i" = "Downloading chain: {.url {url}}"))
    utils::download.file(url, gzlocal, mode = "wb", quiet = TRUE)
  }

  # Decompress; prefer R.utils, fall back to system gunzip.
  if (requireNamespace("R.utils", quietly = TRUE)) {
    R.utils::gunzip(gzlocal, destname = local, remove = FALSE)
  } else {
    exit <- system2("gunzip", c("-c", gzlocal), stdout = local)
    if (!identical(exit, 0L))
      cli::cli_abort("Failed to decompress {.path {gzlocal}} (install R.utils or gunzip).")
  }
  local
}

#' Build a merged annotated credible-set table (Chunk A composed output)
#'
#' Composes the Step-1 data products into a single per-variant tibble:
#' fine-mapped sQTL credible-set rows (splaire TSV) + splice-prediction
#' deltas (splaire / splaireVar) + gnomAD rsID + allele frequencies +
#' COPD GWAS effect estimates (rsID then position fallback).
#'
#' This is the canonical Chunk A deliverable — one function, one
#' merged table per gene.
#'
#' @param sm_predictions Splicing-model variant-effect predictions
#'   (path or data.frame). Currently expects the splaire v1 schema;
#'   multi-tissue / multi-model composition is a planned extension —
#'   see "Design notes" in `docs/methods.md`.
#' @param gwas A GWAS tibble from [loadCopdGwas()], OR a path to a
#'   sumstats file (loaded transparently).
#' @param gene Optional HGNC gene symbol(s) to restrict to.
#' @param gnomad_cache_dir Cache directory for [annotateGnomad()].
#' @param liftover_chain Optional chain path or object; passed to
#'   [annotateGwas()]. If `NULL` and `variants_build != gwas_build`,
#'   position-matching is skipped and only rsID matches are used.
#' @param variants_build,gwas_build Build labels. Defaults assume
#'   splaire=GRCh38 and icgcUkb=GRCh37.
#' @param fields `"default"` (returns only canonical columns from each
#'   upstream layer per `docs/schemas.md`) or `"all"` (carries every
#'   upstream column forward, equivalent to prior behavior).
#' @return A **nested tibble**: one row per credible set, with a
#'   list-column `variants` holding the per-variant tibble. Carries a
#'   merged `isovar_meta` attribute with provenance from all three
#'   upstream sources. Use [as_flat()] to unnest into the flat form,
#'   and [write_tsv_pair()] to write credible_sets + variants TSVs.
#' @export
buildAnnotatedCredibleSet <- function(sm_predictions,
                                      gwas,
                                      gene = NULL,
                                      gnomad_cache_dir = NULL,
                                      liftover_chain = NULL,
                                      variants_build = "GRCh38",
                                      gwas_build = "GRCh37",
                                      fields = c("default", "all")) {
  fields <- match.arg(fields)

  ranked <- rankSplaireVariants(sm_predictions, gene = gene)
  if (nrow(ranked) == 0L) return(nestCredibleSets(ranked))

  gnom <- annotateGnomad(ranked$variant_id, cache_dir = gnomad_cache_dir,
                         fields = fields)
  gnom_keep <- gnom[, setdiff(names(gnom),
                              c("chr", "pos", "ref", "alt"))]
  variants <- dplyr::left_join(ranked, gnom_keep, by = "variant_id")
  variants <- setIsovarMeta(variants, mergeIsovarMeta(ranked, gnom))

  if (is.character(gwas) && length(gwas) == 1L)
    gwas <- loadCopdGwas(gwas)

  flat <- annotateGwas(variants, gwas,
                       liftover_chain = liftover_chain,
                       variants_build = variants_build,
                       gwas_build     = gwas_build)

  # Drop GWAS "extended" columns when fields == "default".
  if (fields == "default") {
    drop_gwas <- c("gwas_z", "gwas_n", "gwas_imputersq", "gwas_marker",
                   "gwas_chr", "gwas_pos")
    flat <- flat[, setdiff(names(flat), drop_gwas), drop = FALSE]
  }

  nested <- nestCredibleSets(flat)
  setIsovarMeta(nested, mergeIsovarMeta(flat, list(sources = list())))
}

#' Re-shape a flat per-variant tibble into a nested credible-set tibble
#'
#' Groups the input by `(gene, phenotype, ensembl_id, strand,
#' credible_set_number)` and packs each group's variant-level fields
#' into a list-column `variants`. Produces the canonical isovar
#' credible-set output shape.
#'
#' @param flat A tibble from the Step-1 pipeline (output of
#'   `rankSplaireVariants()` augmented with gnomAD + GWAS columns).
#' @return A nested tibble with one row per credible set and columns:
#'   `credible_set_id`, `gene`, `phenotype`, `ensembl_id`, `strand`,
#'   `credible_set_number`, `n_variants`, `variants` (list-column).
#'   The input's `isovar_meta` attribute is preserved.
#' @export
nestCredibleSets <- function(flat) {
  key_cols <- c("gene", "phenotype", "ensembl_id", "strand",
                "credible_set_number")
  missing <- setdiff(key_cols, names(flat))
  if (length(missing))
    cli::cli_abort("Input missing required credible-set columns: {.val {missing}}")

  # Pull 'pip' alias into the variant-level tibble for cleaner schema.
  flat2 <- flat
  if ("posterior_inclusion_probability" %in% names(flat2)) {
    flat2$pip <- flat2$posterior_inclusion_probability
    flat2$posterior_inclusion_probability <- NULL
  }

  variant_cols <- setdiff(names(flat2), key_cols)

  nested <- flat2 %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
    tidyr::nest(variants = dplyr::all_of(variant_cols)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      credible_set_id = paste(gene, phenotype, sep = "__"),
      n_variants      = vapply(variants, nrow, integer(1))
    ) %>%
    dplyr::relocate(credible_set_id) %>%
    dplyr::relocate(n_variants, .after = "credible_set_number")

  meta <- getIsovarMeta(flat, require = FALSE)
  if (!is.null(meta)) attr(nested, "isovar_meta") <- meta
  nested
}

#' Unnest a credible-set tibble back into the flat per-variant shape
#'
#' Inverse of [nestCredibleSets()]. Credible-set-level columns are
#' duplicated across each CS's variants.
#'
#' @param x A nested credible-set tibble.
#' @return A flat tibble with one row per variant.
#' @export
as_flat <- function(x) {
  if (!"variants" %in% names(x))
    cli::cli_abort("Input does not have a {.field variants} list-column.")
  flat <- tidyr::unnest(x, cols = "variants")
  meta <- getIsovarMeta(x, require = FALSE)
  if (!is.null(meta)) attr(flat, "isovar_meta") <- meta
  flat
}

#' Write a credible-set tibble as a pair of TSVs with metadata sidecars
#'
#' Produces `<dir>/credible_sets.tsv` (one row per CS, no list-column)
#' and `<dir>/variants.tsv` (one row per variant, with
#' `credible_set_id` joining back). Metadata sidecars `.meta.json` are
#' written next to each TSV.
#'
#' @param x A nested credible-set tibble from [buildAnnotatedCredibleSet()].
#' @param dir Output directory; created if absent.
#' @return Invisibly, a named list of paths written.
#' @export
write_tsv_pair <- function(x, dir) {
  if (!"variants" %in% names(x))
    cli::cli_abort("Input does not have a {.field variants} list-column.")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  # Credible-set-level: strip variants, strip its meta so it writes plain
  cs_tbl <- x
  cs_tbl$variants <- NULL
  cs_path <- file.path(dir, "credible_sets.tsv")
  readr::write_tsv(cs_tbl, cs_path)

  # Per-variant: unnest + add credible_set_id first
  var_tbl <- as_flat(x)
  var_path <- file.path(dir, "variants.tsv")
  readr::write_tsv(var_tbl, var_path)

  # Sidecar metadata
  meta <- getIsovarMeta(x, require = FALSE)
  if (!is.null(meta)) {
    writeMeta(meta, cs_path)
    writeMeta(meta, var_path)
  }

  invisible(list(credible_sets = cs_path, variants = var_path,
                 credible_sets_meta = paste0(cs_path, ".meta.json"),
                 variants_meta = paste0(var_path, ".meta.json")))
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
