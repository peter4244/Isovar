#' Run isoscope's gene_isoform_annotation against a source config
#'
#' Thin wrapper around the canonical
#' `gene_isoform_annotation.R` in a local isoscope checkout. Invokes
#' it as a subprocess with `--config <path>` so each long-read source
#' can have its own config without mutating the isoscope clone. Parses
#' the resulting TSV and attaches an `isovar_meta` attribute derived
#' from the config's `ISOVAR_SOURCE_META` list.
#'
#' See `inst/isoscope_configs/` for example configs and `docs/methods.md`
#' for the multi-source long-read design.
#'
#' @param gene HGNC symbol or versioned Ensembl ID (e.g. `"AKR1A1"`,
#'   `"ENSG00000117448.15"`).
#' @param isoscope_config_path Absolute path to a `config.R` file
#'   (preferably one of `inst/isoscope_configs/*.R`). Must set the
#'   isoscope paths and ideally an `ISOVAR_SOURCE_META` list.
#' @param output_dir Where isoscope will write its per-gene TSV. If
#'   `NULL`, a temporary directory is used. Re-runs with the same
#'   `output_dir` are idempotent — isoscope will simply overwrite.
#' @param isoscope_dir Path to the canonical isoscope checkout (holds
#'   `gene_isoform_annotation.R` and `gencode_v49_genes.gtf`). Default
#'   `~/claude_projects/isoscope`.
#' @param no_expr If `TRUE` (default), passes `--no-expr` to isoscope
#'   (skips expression computation; Chunk D handles per-condition
#'   counts externally).
#' @param extra_args Optional character vector of extra CLI args to
#'   pass through (e.g. `c("--gtf", "--fasta")`).
#' @return A tibble of the parsed isoform annotation (see
#'   `docs/schemas.md` — isoscope output columns), with an
#'   `isovar_meta` attribute carrying the source's provenance.
#' @export
runIsoscopeGene <- function(gene,
                            isoscope_config_path,
                            output_dir = NULL,
                            isoscope_dir = "~/claude_projects/isoscope",
                            no_expr = TRUE,
                            extra_args = NULL) {
  stopifnot(is.character(gene), length(gene) == 1L, nzchar(gene))
  if (!grepl("^[A-Za-z0-9._-]+$", gene))
    cli::cli_abort("Invalid gene identifier: {.val {gene}}")

  isoscope_dir  <- path.expand(isoscope_dir)
  script        <- file.path(isoscope_dir, "gene_isoform_annotation.R")
  if (!file.exists(script))
    cli::cli_abort("isoscope script not found at {.path {script}}. Clone from github.com/peter4244/isoscope.")

  isoscope_config_path <- path.expand(isoscope_config_path)
  if (!file.exists(isoscope_config_path))
    cli::cli_abort("isoscope config not found: {.path {isoscope_config_path}}")

  if (is.null(output_dir)) output_dir <- tempfile(paste0("isovar_isoscope_", gene, "_"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  src_meta <- .read_isovar_source_meta(isoscope_config_path)

  args <- c(gene,
            "--config", isoscope_config_path,
            "--output-dir", output_dir)
  if (isTRUE(no_expr)) args <- c(args, "--no-expr")
  if (!is.null(extra_args)) args <- c(args, extra_args)

  cli::cli_inform(c("i" = "Running isoscope for {.val {gene}} via {.path {basename(isoscope_config_path)}}..."))
  t0 <- Sys.time()
  # isoscope resolves ./gencode_v49_genes.gtf relative to CWD; run from
  # its own directory so the lookup path resolves.
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(isoscope_dir)
  status <- system2("Rscript", c(shQuote(script), args),
                    stdout = FALSE, stderr = FALSE)
  setwd(oldwd)
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  if (!identical(status, 0L))
    cli::cli_abort("isoscope exited with status {status}. Inspect output_dir: {.path {output_dir}}")
  cli::cli_alert_success("isoscope completed in {sprintf('%.1f', dt)}s")

  # isoscope writes isoform_annotation_<GENE_NAME>.tsv. Gene-symbol inputs
  # produce a file named after the symbol; Ensembl-id inputs name it
  # after the resolved HGNC symbol — search the dir either way.
  tsv_candidates <- list.files(output_dir,
                               pattern = "^isoform_annotation_.*\\.tsv$",
                               full.names = TRUE)
  if (length(tsv_candidates) == 0L)
    cli::cli_abort("No isoform_annotation_*.tsv produced in {.path {output_dir}}.")
  tsv_path <- tsv_candidates[1]

  iso <- readr::read_tsv(tsv_path, show_col_types = FALSE, progress = FALSE)
  iso <- setIsovarMeta(
    iso,
    genome_build = src_meta$genome_build %||% NA_character_,
    gtf_version  = src_meta$gencode_version %||% NA_character_,
    sources      = list(list(
      id           = src_meta$source_id %||% "unknown_source",
      kind         = "long_read_isoform_annotation",
      label        = src_meta$source_label %||% NA_character_,
      cohort       = src_meta$cohort %||% NA_character_,
      tissue       = src_meta$tissue %||% NA_character_,
      n_samples    = src_meta$n_samples %||% NA_integer_,
      genome_build = src_meta$genome_build %||% NA_character_,
      gencode_version = src_meta$gencode_version %||% NA_character_,
      sqanti_run_date = src_meta$sqanti_run_date %||% NA_character_,
      conditions   = src_meta$conditions,
      isoscope_tsv = normalizePath(tsv_path, mustWork = FALSE)
    ))
  )
  iso
}

#' Classify isoforms by whether they use specified splice-site positions
#'
#' For each genomic position in `site_positions`, adds a logical column
#' `uses_<pos>_<kind>` indicating whether that isoform's `junctions`
#' column contains a tuple anchored at that position. Default
#' `kind = "junction_start"` treats a site as "used" when the isoform
#' has a junction whose **start** equals the site (i.e. donor on +
#' strand / acceptor on − strand). `kind = "junction_end"` uses the
#' junction's end (acceptor on + strand).
#'
#' The isoscope `junctions` column schema is
#' `start_end_strand;start_end_strand;...` — this function parses it.
#'
#' @param iso_df A tibble produced by [runIsoscopeGene()] (must have
#'   a `junctions` column in the isoscope format).
#' @param site_positions Integer vector of genomic positions.
#' @param kind `"junction_start"` or `"junction_end"`.
#' @return `iso_df` with one new logical column per site.
#' @export
classifyIsoformsBySiteUsage <- function(iso_df, site_positions,
                                        kind = c("junction_start", "junction_end")) {
  kind <- match.arg(kind)
  if (!"junctions" %in% names(iso_df))
    cli::cli_abort("Input has no {.field junctions} column.")
  if (!is.numeric(site_positions) || length(site_positions) == 0L)
    cli::cli_abort("{.arg site_positions} must be a non-empty numeric vector.")
  site_positions <- as.integer(site_positions)

  extract <- if (kind == "junction_start") .junction_starts else .junction_ends

  junc_anchors <- lapply(iso_df$junctions, extract)
  for (p in site_positions) {
    col <- paste0("uses_", p, "_",
                  if (kind == "junction_start") "jstart" else "jend")
    iso_df[[col]] <- vapply(junc_anchors,
                            function(v) !is.null(v) && p %in% v,
                            logical(1))
  }
  iso_df
}

# -- internals ---------------------------------------------------------

.read_isovar_source_meta <- function(config_path) {
  # Source the isovar config into an isolated env (parented on baseenv
  # so base operators like `<-` resolve) so it doesn't clobber the
  # caller's workspace; isoscope sources it independently inside the
  # subprocess.
  env <- new.env(parent = baseenv())
  tryCatch(
    sys.source(config_path, envir = env),
    error = function(e)
      cli::cli_abort("Failed to source isoscope config {.path {config_path}}: {e$message}")
  )
  if (!exists("ISOVAR_SOURCE_META", envir = env, inherits = FALSE)) {
    cli::cli_warn("Config {.path {config_path}} has no {.field ISOVAR_SOURCE_META}; metadata will be minimal.")
    return(list())
  }
  get("ISOVAR_SOURCE_META", envir = env)
}

.junction_starts <- function(junctions) {
  if (is.na(junctions) || !nzchar(junctions)) return(integer(0))
  parts <- strsplit(junctions, ";", fixed = TRUE)[[1]]
  suppressWarnings(as.integer(sub("_.*", "", parts)))
}

.junction_ends <- function(junctions) {
  if (is.na(junctions) || !nzchar(junctions)) return(integer(0))
  parts <- strsplit(junctions, ";", fixed = TRUE)[[1]]
  m <- regmatches(parts, regexec("^[^_]+_([^_]+)_", parts))
  ends <- vapply(m, function(mm) if (length(mm) >= 2L) mm[2] else NA_character_, character(1))
  suppressWarnings(as.integer(ends))
}
