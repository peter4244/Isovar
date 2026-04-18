#' Provenance metadata for isovar result objects
#'
#' Every isovar-produced data object carries a structured metadata
#' attribute (`attr(x, "isovar_meta")`) that describes its provenance:
#' isovar version, genome build, GTF version (when applicable), upstream
#' sources, and generation timestamp. Metadata propagates through every
#' transformation and can be written to a sidecar `.meta.json` next to
#' saved results.
#'
#' See `docs/schemas.md` ("Provenance metadata") for the canonical
#' schema.
#'
#' @name isovar-metadata
NULL

#' Attach or update provenance metadata on an isovar object
#'
#' Overwrites any existing `isovar_meta` attribute. Use
#' [mergeIsovarMeta()] when combining with pre-existing metadata.
#'
#' Accepts either a single metadata list (from [mergeIsovarMeta()] or
#' [getIsovarMeta()]) or named field arguments:
#'
#' ```
#' setIsovarMeta(x, merged_meta_list)               # full replacement
#' setIsovarMeta(x, genome_build = "GRCh38", ...)   # named fields
#' ```
#'
#' @param x Any R object (typically a tibble).
#' @param meta Optional metadata list (e.g. from [mergeIsovarMeta()]).
#'   When supplied, takes precedence over `...`.
#' @param ... Named fields to include in the metadata list. Common:
#'   `isovar_version`, `generated_at`, `genome_build`, `gtf_version`,
#'   `sources` (list of source descriptors).
#' @return `x`, with `isovar_meta` attribute set.
#' @export
setIsovarMeta <- function(x, meta = NULL, ...) {
  if (is.null(meta)) meta <- list(...)
  attr(x, "isovar_meta") <- .normalize_meta(meta)
  x
}

#' Retrieve the provenance metadata from an isovar object
#'
#' @param x An isovar-produced object (e.g., a tibble returned by one
#'   of the annotation functions).
#' @param require If `TRUE` (default), error when no metadata is
#'   attached; if `FALSE`, return `NULL`.
#' @return The metadata list, or `NULL` if absent and `require=FALSE`.
#' @export
getIsovarMeta <- function(x, require = TRUE) {
  meta <- attr(x, "isovar_meta")
  if (is.null(meta) && require)
    cli::cli_abort("Object has no {.field isovar_meta} attribute.")
  meta
}

#' Merge provenance metadata from multiple isovar objects
#'
#' Preserves all `sources` from every input; concatenates into a single
#' list, dedup'd by `id`. `genome_build` and `gtf_version` are carried
#' from the first argument; disagreements across inputs are flagged in
#' a `mixed_builds` / `mixed_gtf` field rather than silently resolved.
#'
#' @param ... isovar-produced objects or bare metadata lists.
#' @return A merged metadata list suitable for passing to [setIsovarMeta()].
#' @export
mergeIsovarMeta <- function(...) {
  xs <- list(...)
  metas <- lapply(xs, function(x) {
    if (is.list(x) && !is.data.frame(x)) x else getIsovarMeta(x, require = FALSE)
  })
  metas <- metas[!vapply(metas, is.null, logical(1))]
  if (length(metas) == 0L) return(.default_meta())

  first <- metas[[1]]
  builds <- unique(stats::na.omit(vapply(metas,
                                         function(m) m$genome_build %||% NA_character_,
                                         character(1))))
  gtfs <- unique(stats::na.omit(vapply(metas,
                                       function(m) m$gtf_version %||% NA_character_,
                                       character(1))))

  all_sources <- unlist(lapply(metas, function(m) m$sources), recursive = FALSE)
  if (!is.null(all_sources)) {
    ids <- vapply(all_sources, function(s) s$id %||% "", character(1))
    all_sources <- all_sources[!duplicated(ids)]
  }

  out <- list(
    isovar_version = first$isovar_version %||% .isovar_version(),
    generated_at   = first$generated_at   %||% .now_iso8601(),
    genome_build   = if (length(builds)) builds[1] else NA_character_,
    gtf_version    = if (length(gtfs))   gtfs[1]   else NA_character_,
    sources        = all_sources
  )
  if (length(builds) > 1L) out$mixed_builds <- builds
  if (length(gtfs)   > 1L) out$mixed_gtf    <- gtfs
  .normalize_meta(out)
}

#' Write provenance metadata as a sidecar JSON next to a data file
#'
#' Produces `<path>.meta.json` with the metadata object.
#'
#' @param x An isovar object carrying metadata (or a metadata list).
#' @param data_path Path to the associated data file. The sidecar is
#'   written next to it as `<data_path>.meta.json`.
#' @return Invisibly, the sidecar path.
#' @export
writeMeta <- function(x, data_path) {
  meta <- if (is.list(x) && !is.data.frame(x)) x else getIsovarMeta(x)
  out_path <- paste0(data_path, ".meta.json")
  if (!requireNamespace("jsonlite", quietly = TRUE))
    cli::cli_abort("jsonlite is required to write sidecar metadata.")
  jsonlite::write_json(meta, out_path, pretty = TRUE, auto_unbox = TRUE, null = "null")
  invisible(out_path)
}

#' Read a provenance metadata sidecar JSON
#'
#' @param data_path Path to the data file. Reads `<data_path>.meta.json`.
#' @return The metadata list.
#' @export
readMeta <- function(data_path) {
  side <- paste0(data_path, ".meta.json")
  if (!file.exists(side))
    cli::cli_abort("No sidecar found at {.path {side}}.")
  if (!requireNamespace("jsonlite", quietly = TRUE))
    cli::cli_abort("jsonlite is required to read sidecar metadata.")
  jsonlite::read_json(side, simplifyVector = FALSE)
}

# -- internals ---------------------------------------------------------

.isovar_version <- function() {
  tryCatch(as.character(utils::packageVersion("isovar")),
           error = function(e) "0.0.2+dev")
}

.now_iso8601 <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

.default_meta <- function() {
  list(
    isovar_version = .isovar_version(),
    generated_at   = .now_iso8601(),
    genome_build   = NA_character_,
    gtf_version    = NA_character_,
    sources        = list()
  )
}

.normalize_meta <- function(meta) {
  base <- .default_meta()
  for (nm in names(meta)) base[[nm]] <- meta[[nm]]
  if (is.null(base$sources)) base$sources <- list()
  base
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || (is.atomic(a) && is.na(a[1]))) b else a
