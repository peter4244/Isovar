#' Parse isoform structures across GENCODE + long-read GTFs
#'
#' Calls [Isopair::parseIsoformStructures()] once per input GTF, filters
#' to the requested `isoform_ids`, and binds the results into a single
#' structures tibble. A new `source_gtf` column tags each row with the
#' GTF it came from (e.g. `"gencode"`, `"nmd_sqanti"`), so downstream
#' analyses can distinguish reference vs. long-read origin.
#'
#' This is the bridge between Chunk C/D's annotation layer and
#' Isopair's structural analyses (buildUnionExons, buildProfiles,
#' traceReferenceAtg, etc.) — Isopair needs a single `structures`
#' tibble spanning every isoform in the pair set.
#'
#' @param gtfs Named list of GTF paths. Names become the `source_gtf`
#'   tag (e.g. `list(gencode = ..., nmd_sqanti = ..., nmd_isocall = ...)`).
#' @param isoform_ids Character vector of isoform IDs to parse. IDs that
#'   don't appear in a given GTF are silently skipped for that GTF.
#' @return A tibble with the Isopair structures schema plus a
#'   `source_gtf` column. Rows are unique on `isoform_id` — if the same
#'   ID appears in multiple GTFs (an ENST resolved both in GENCODE and
#'   in the long-read run), the first GTF in `gtfs` wins.
#' @export
parseStructuresMultiGtf <- function(gtfs, isoform_ids) {
  stopifnot(is.list(gtfs), length(gtfs) > 0L, !is.null(names(gtfs)))
  if (!requireNamespace("Isopair", quietly = TRUE))
    cli::cli_abort("Package {.pkg Isopair} is required.")

  parts <- list()
  seen  <- character(0)
  for (tag in names(gtfs)) {
    gtf_path <- gtfs[[tag]]
    if (!file.exists(gtf_path))
      cli::cli_abort("GTF not found: {.path {gtf_path}}")
    todo <- setdiff(isoform_ids, seen)
    if (length(todo) == 0L) next

    cli::cli_inform(c("i" = "Parsing {length(todo)} isoform(s) from {.val {tag}}..."))
    structs <- Isopair::parseIsoformStructures(gtf_path, todo, verbose = FALSE)
    if (!is.null(structs) && nrow(structs) > 0L) {
      structs$source_gtf <- tag
      parts[[tag]] <- structs
      seen <- c(seen, structs$isoform_id)
    }
  }

  missing <- setdiff(isoform_ids, seen)
  if (length(missing) > 0L)
    cli::cli_warn("Could not parse structures for {length(missing)} isoform(s); \\
                  first few: {.val {head(missing, 5)}}.")

  dplyr::bind_rows(parts)
}

#' Extract CDS annotations across GENCODE + long-read GTFs
#'
#' Companion to [parseStructuresMultiGtf()]. Calls
#' [Isopair::extractCdsAnnotations()] per GTF, preserves a `source_gtf`
#' tag, and binds. Used for [Isopair::computePtcStatus()] and
#' [Isopair::traceReferenceAtg()] inputs.
#'
#' @param gtfs Named list of GTF paths. Same shape as
#'   `parseStructuresMultiGtf()`.
#' @param isoform_ids Character vector of isoform IDs.
#' @return A tibble with Isopair's CDS schema plus `source_gtf`.
#' @export
extractCdsMultiGtf <- function(gtfs, isoform_ids) {
  stopifnot(is.list(gtfs), length(gtfs) > 0L)
  if (!requireNamespace("Isopair", quietly = TRUE))
    cli::cli_abort("Package {.pkg Isopair} is required.")

  parts <- list()
  seen  <- character(0)
  for (tag in names(gtfs)) {
    gtf_path <- gtfs[[tag]]
    todo <- setdiff(isoform_ids, seen)
    if (length(todo) == 0L) next
    cds <- Isopair::extractCdsAnnotations(gtf_path, todo, verbose = FALSE)
    if (!is.null(cds) && nrow(cds) > 0L) {
      cds$source_gtf <- tag
      parts[[tag]] <- cds
      seen <- c(seen, cds$isoform_id)
    }
  }
  dplyr::bind_rows(parts)
}

#' Load transcript sequences into the named-character-vector shape
#' `Isopair::traceReferenceAtg()` expects
#'
#' For each FASTA path in `fastas`, reads with
#' `Biostrings::readDNAStringSet`, filters to the requested
#' `isoform_ids`, and coerces to a named character vector with
#' **uppercase** DNA sequence values.
#'
#' FASTA header formats across GENCODE / SQANTI / isocall differ.
#' GENCODE uses pipe-delimited headers with the ENST as the first
#' token; SQANTI / isocall typically use the bare transcript ID. This
#' function extracts the first whitespace-delimited token as the key,
#' then strips any trailing `|…` to match either style.
#'
#' @param fastas Character vector of FASTA paths (gzipped OK). Earlier
#'   entries take precedence on collisions.
#' @param isoform_ids Character vector of IDs to retain.
#' @return A named character vector (length = number of matched IDs,
#'   values uppercase DNA).
#' @export
loadTranscriptSequences <- function(fastas, isoform_ids) {
  if (!requireNamespace("Biostrings", quietly = TRUE))
    cli::cli_abort("Package {.pkg Biostrings} is required.")
  if (length(fastas) == 0L) return(character(0))

  combined <- character(0)
  seen     <- character(0)
  for (fa in fastas) {
    if (!file.exists(fa))
      cli::cli_abort("FASTA not found: {.path {fa}}")
    dss <- Biostrings::readDNAStringSet(fa)
    # Normalize names: first whitespace-delimited token, then strip |-tail
    raw <- names(dss)
    keys <- sub("\\|.*$", "", sub("\\s.*$", "", raw))
    names(dss) <- keys
    want <- setdiff(intersect(keys, isoform_ids), seen)
    if (length(want) == 0L) next
    subset <- dss[want]
    seq_chr <- toupper(as.character(subset))
    names(seq_chr) <- want
    combined <- c(combined, seq_chr)
    seen <- c(seen, want)
  }

  missing <- setdiff(isoform_ids, seen)
  if (length(missing) > 0L)
    cli::cli_warn("{length(missing)} isoform(s) had no FASTA sequence; \\
                  first few: {.val {head(missing, 5)}}.")
  combined
}

#' Classify hidden PTCs via `Isopair::traceReferenceAtg`
#'
#' Thin wrapper. For each (reference, comparator) pair, traces the
#' reference isoform's ATG through the comparator's coordinates and
#' calls a hidden-PTC category (`effectively_ptc`, `truncated_no_ejc`,
#' `ref_atg_lost`, `no_downstream_ejc`, `no_ref_cds`,
#' `mapping_failed`). Essential for comparators whose own CDS
#' prediction failed — i.e. most novel IR isoforms in SQANTI.
#'
#' @param pairs A tibble with `reference_isoform_id` and
#'   `comparator_isoform_id` columns (and typically `gene_id`).
#' @param structures Output of [parseStructuresMultiGtf()].
#' @param cds CDS annotations (output of [extractCdsMultiGtf()]).
#' @param sequences Named char vector from [loadTranscriptSequences()].
#' @param ejc_threshold EJC distance cutoff (default 50, the 50-nt rule).
#' @param resolve_alt_start Logical; when TRUE (isovar default), comparators
#'   whose reference ATG isn't exonic are re-evaluated from the first viable
#'   alternative ATG in the comparator (Isopair's `resolve_alt_start` path).
#'   Defaults to TRUE here because isovar's workflow always wants a resolved
#'   NMD verdict — `ref_atg_lost` as a terminal label silently drops
#'   translated-but-alt-start isoforms from downstream NMD analysis.
#' @param min_alt_orf_nt Minimum alt-start ORF length (nt) to count as
#'   viable. Default 30 (10 aa). Passed through to Isopair.
#' @return Isopair's per-pair hidden-PTC tibble, unchanged in schema; when
#'   `resolve_alt_start = TRUE` the tibble also contains `alt_start_tx_pos`,
#'   `alt_start_orf_length`, and may emit `alt_start_effectively_ptc`,
#'   `alt_start_no_downstream_ejc`, `ref_atg_lost_no_viable_start` categories.
#' @export
hiddenPtcStatus <- function(pairs, structures, cds, sequences,
                            ejc_threshold = 50L,
                            resolve_alt_start = TRUE,
                            min_alt_orf_nt = 30L) {
  if (!requireNamespace("Isopair", quietly = TRUE))
    cli::cli_abort("Package {.pkg Isopair} is required.")
  Isopair::traceReferenceAtg(pairs, structures, cds, sequences,
                             ejc_threshold     = ejc_threshold,
                             resolve_alt_start = resolve_alt_start,
                             min_alt_orf_nt    = min_alt_orf_nt)
}

#' Enumerate and classify every viable ORF in a set of transcripts
#'
#' Thin wrapper around \code{Isopair::enumerateOrfs()} — returns one row
#' per (isoform, viable ORF), with per-ORF NMD classification. This is
#' the per-ORF companion to [hiddenPtcStatus()]'s per-(ref, comp)
#' output: while `hiddenPtcStatus` evaluates one ORF per pair (ref's
#' ATG traced through comp, with optional alt-start fallback),
#' `enumerateComparatorOrfs` enumerates *every* viable ORF in each
#' comparator transcript and classifies each independently.
#'
#' See the `feedback_per_orf_classifications` memory note for the
#' design principle: transcript-level PTC / NMD labels are rollups over
#' per-ORF verdicts, not intrinsic properties of transcripts.
#'
#' @param structures Output of [parseStructuresMultiGtf()] — or any
#'   Isopair structures tibble.
#' @param cds Output of [extractCdsMultiGtf()].
#' @param sequences Output of [loadTranscriptSequences()].
#' @param min_orf_nt Minimum ORF length (nt); default 30.
#' @param ejc_threshold Downstream-EJC distance cutoff; default 50.
#' @param include_no_stop Emit rows for ATGs whose frame runs off the
#'   transcript (category `no_stop_in_frame`). Default `TRUE`.
#' @param kozak_filter Logical; when `TRUE` (isovar default), only ATGs
#'   whose Kozak context scores at/above `kozak_threshold` are kept
#'   before ORF tracing. Enumerating every ATG regardless of Kozak
#'   produces mostly noise — plausible translation starts need
#'   initiation-competent context.
#' @param kozak_threshold Numeric log-odds threshold. Default 0 (above
#'   random). For a data-driven threshold, pass the output of
#'   [Isopair::empiricalKozakThreshold()] on an annotated-CDS training
#'   set (e.g. GENCODE coding isoforms).
#' @return A long-format tibble with one row per plausibly-translated
#'   (isoform, ORF): `isoform_id`, `atg_tx_pos`, `kozak_score`,
#'   `stop_tx_pos`, `orf_length`, `n_downstream_ejc`, `is_annotated_cds`,
#'   `category`.
#' @export
enumerateComparatorOrfs <- function(structures, cds, sequences,
                                    min_orf_nt = 30L,
                                    ejc_threshold = 50L,
                                    include_no_stop = TRUE,
                                    kozak_filter = TRUE,
                                    kozak_threshold = 0) {
  if (!requireNamespace("Isopair", quietly = TRUE))
    cli::cli_abort("Package {.pkg Isopair} is required.")
  Isopair::enumerateOrfs(
    structures     = structures,
    cds_metadata   = cds,
    sequences      = sequences,
    min_orf_nt     = min_orf_nt,
    ejc_threshold  = ejc_threshold,
    include_no_stop = include_no_stop,
    kozak_filter   = kozak_filter,
    kozak_threshold = kozak_threshold
  )
}

#' Roll per-ORF classifications up to a per-transcript NMD verdict
#'
#' Convenience aggregator over the output of
#' [enumerateComparatorOrfs()]. An isoform is flagged as an NMD
#' substrate when *any* of its ORFs is `effectively_ptc`. Use this
#' when a report needs a transcript-level verdict while keeping the
#' per-ORF detail accessible.
#'
#' @param orfs Output of [enumerateComparatorOrfs()] /
#'   [Isopair::enumerateOrfs()].
#' @return A tibble with one row per isoform: `isoform_id`, `n_orfs`,
#'   `n_ptc_orfs`, `any_ptc` (logical), `has_annotated_cds_ptc`.
#' @export
summarizeOrfsToTranscript <- function(orfs) {
  orfs %>%
    dplyr::group_by(isoform_id) %>%
    dplyr::summarise(
      n_orfs               = dplyr::n(),
      n_ptc_orfs           = sum(category == "effectively_ptc", na.rm = TRUE),
      n_no_stop            = sum(category == "no_stop_in_frame", na.rm = TRUE),
      any_ptc              = any(category == "effectively_ptc", na.rm = TRUE),
      has_annotated_cds_ptc = any(is_annotated_cds &
                                   category == "effectively_ptc", na.rm = TRUE),
      .groups = "drop"
    )
}

#' Identify dominant isoforms across long-read sources via Isopair
#'
#' Wraps `Isopair::identifyDominantIsoforms()` with the multi-source
#' count-matrix + gene-map shape produced by [loadLongreadEvidence()].
#' Returns one dominant isoform per (source, gene) at the requested
#' dominance threshold (default 0.5 = must account for > 50% of gene
#' expression summed across samples).
#'
#' @param lr_evidence Output of [loadLongreadEvidence()].
#' @param source_id Name of the source to use (default `"nmd_sqanti"`).
#' @param threshold Dominance threshold passed to Isopair (default 0.5).
#' @param gene_id Optional gene identifier to tag every isoform with in
#'   the gene_map. isoscope's per-gene TSV doesn't carry a `gene_id`
#'   column, but every row corresponds to a single queried gene — pass
#'   that gene's Ensembl ID or HGNC symbol here. Defaults to the value
#'   of the first non-NA `gene_id` / `gene_name` column if one exists,
#'   otherwise `"UNKNOWN_GENE"`.
#' @return Tibble with one row per gene: `gene_id`, `dominant_isoform_id`,
#'   `dominant_cpm_fraction`.
#' @export
computeDominantIsoform <- function(lr_evidence,
                                   source_id = "nmd_sqanti",
                                   threshold = 0.5,
                                   gene_id = NULL) {
  if (!requireNamespace("Isopair", quietly = TRUE))
    cli::cli_abort("Package {.pkg Isopair} is required.")
  idx <- which(lr_evidence$source_id == source_id)
  if (length(idx) == 0L)
    cli::cli_abort("source_id {.val {source_id}} not in {.arg lr_evidence}.")

  ann <- lr_evidence$annotation[[idx]]
  cnt <- lr_evidence$counts[[idx]]
  smp <- lr_evidence$samples[[idx]]
  if (is.null(cnt) || nrow(cnt) == 0L)
    cli::cli_abort("Counts missing for source {.val {source_id}}.")

  if (is.null(gene_id)) {
    if ("gene_id" %in% names(ann) && any(!is.na(ann$gene_id))) {
      gene_id <- ann$gene_id[which(!is.na(ann$gene_id))[1]]
    } else if ("gene_name" %in% names(ann) && any(!is.na(ann$gene_name))) {
      gene_id <- ann$gene_name[which(!is.na(ann$gene_name))[1]]
    } else {
      gene_id <- "UNKNOWN_GENE"
    }
  }

  common <- intersect(ann$isoform_id, rownames(cnt))
  gene_map <- tibble::tibble(
    isoform_id = common,
    gene_id    = rep(gene_id, length(common))
  )

  Isopair::identifyDominantIsoforms(
    expression_matrix = cnt[common, , drop = FALSE],
    gene_map          = gene_map,
    samples           = smp$column_name,
    threshold         = threshold
  )
}
