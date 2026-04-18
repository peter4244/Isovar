#' Catalogue the long-read sources available to isovar
#'
#' Returns a list of source descriptors, each containing the isoscope
#' config path, the count-matrix path, and the regex used to parse the
#' count-matrix column names into per-sample metadata (donor, celltype,
#' treatment). The defaults cover the four canonical long-read views
#' we currently support:
#'
#' * `haec185_sqanti`  — SQANTI-corrected, 9 basal donors
#' * `haec185_isocall` — pre-SQANTI isocall, same 9 donors
#' * `nmd_sqanti`      — SQANTI-corrected, 38 samples (6 cell types × Smg1i/DMSO)
#' * `nmd_isocall`     — pre-SQANTI isocall, same 38 samples
#'
#' Each source also specifies a `sample_schema` (`"donor_only"` for
#' HAEC-185 or `"sample_celltype_donor_treatment"` for NMD) that tells
#' [parseSampleColumns()] how to decode its count-matrix header.
#'
#' @param isovar_dir Path to the isovar repo (used to resolve
#'   `inst/isoscope_configs/*.R` paths). Default the repo root.
#' @return A named list of source descriptors.
#' @export
longreadSources <- function(isovar_dir = "/Users/petecastaldi/claude_projects/isovar") {
  cfg <- function(f) file.path(isovar_dir, "inst/isoscope_configs", f)

  list(
    haec185_sqanti = list(
      label             = "HAEC-185 basal (SQANTI)",
      cohort            = "HAEC-185",
      tissue            = "basal airway epithelium",
      source_kind       = "long_read_isoform_annotation",
      calling_stage     = "sqanti",
      isoscope_config   = cfg("haec185_basal.R"),
      count_matrix_path = "/Users/petecastaldi/claude_projects/HAEC185/longread/2024_9samples/isocall/sqanti/20260304/LR_RNA_Basal_filtered.count_matrix.txt",
      count_delim       = ",",
      sample_schema     = "donor_only"
    ),
    haec185_isocall = list(
      label             = "HAEC-185 basal (isocall, pre-SQANTI)",
      cohort            = "HAEC-185",
      tissue            = "basal airway epithelium",
      source_kind       = "long_read_isoform_annotation",
      calling_stage     = "isocall_uncorrected",
      isoscope_config   = cfg("haec185_basal_isocall.R"),
      count_matrix_path = "/Users/petecastaldi/claude_projects/HAEC185/longread/2024_9samples/isocall/20260304/call/LR_RNA_Basal_isocall.count_matrix.txt",
      count_delim       = ",",
      sample_schema     = "donor_only"
    ),
    nmd_sqanti = list(
      label             = "NMD lungcells (SQANTI)",
      cohort            = "NMD lung cell lines",
      tissue            = "lung cell lines (DD, DD_ALI, AT, DO, FB, MV)",
      source_kind       = "long_read_isoform_annotation",
      calling_stage     = "sqanti",
      isoscope_config   = cfg("nmd_lungcells.R"),
      count_matrix_path = "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_filtered.count_matrix.txt",
      count_delim       = ",",
      sample_schema     = "sample_celltype_donor_treatment"
    ),
    nmd_isocall = list(
      label             = "NMD lungcells (isocall, pre-SQANTI)",
      cohort            = "NMD lung cell lines",
      tissue            = "lung cell lines (DD, DD_ALI, AT, DO, FB, MV)",
      source_kind       = "long_read_isoform_annotation",
      calling_stage     = "isocall_uncorrected",
      isoscope_config   = cfg("nmd_lungcells_isocall.R"),
      count_matrix_path = "/Users/petecastaldi/claude_projects/nmd/isocall/nmd_lungcells/results/call/nmd_isocall.count_matrix.txt",
      count_delim       = ",",
      sample_schema     = "sample_celltype_donor_treatment"
    )
  )
}

#' Parse count-matrix column names into per-sample metadata
#'
#' Decodes each sample column in a long-read count matrix into
#' `donor`, `celltype`, `treatment`, and `sample_num` fields according
#' to the source's `sample_schema`. Supported schemas:
#'
#' * **`donor_only`** — each column is a bare donor ID (HAEC-185
#'   style: `DD071Q`, `DD032OP2`, …). `celltype` / `treatment` are
#'   `NA`; `condition` is recorded as `"baseline"`.
#' * **`sample_celltype_donor_treatment`** — NMD style:
#'   `^Sample(\d+)_(CELLTYPE)_(DONOR)_(TREATMENT)$`, with `CELLTYPE`
#'   matching `DD_ALI|DD|AT|DO|FB|MV` (DD_ALI must precede DD in the
#'   alternation so the longer match wins) and `TREATMENT` matching
#'   `DMSO|Smg1i`.
#'
#' Columns whose names are exactly `"id"` are dropped (isoform ID
#' column, not a sample). Other unparseable columns emit a warning
#' and land with all fields `NA`.
#'
#' @param column_names Character vector of count-matrix column names.
#' @param sample_schema One of `"donor_only"` or
#'   `"sample_celltype_donor_treatment"`.
#' @return A tibble with one row per input column: `column_name`,
#'   `sample_num`, `donor`, `celltype`, `treatment`, `condition`.
#' @export
parseSampleColumns <- function(column_names, sample_schema) {
  sample_schema <- match.arg(sample_schema,
                             c("donor_only", "sample_celltype_donor_treatment"))
  keep <- column_names != "id"
  col  <- column_names[keep]

  out <- tibble::tibble(
    column_name = col,
    sample_num  = NA_integer_,
    donor       = NA_character_,
    celltype    = NA_character_,
    treatment   = NA_character_,
    condition   = NA_character_
  )

  if (sample_schema == "donor_only") {
    out$donor     <- col
    out$condition <- "baseline"
  } else if (sample_schema == "sample_celltype_donor_treatment") {
    # DD_ALI alternation must come first so the longer token wins.
    pat <- "^Sample(\\d+)_(DD_ALI|DD|AT|DO|FB|MV)_([^_]+)_(DMSO|Smg1i)$"
    m <- regmatches(col, regexec(pat, col))
    for (i in seq_along(m)) {
      grps <- m[[i]]
      if (length(grps) == 5L) {
        out$sample_num[i] <- as.integer(grps[2])
        out$celltype[i]   <- grps[3]
        out$donor[i]      <- grps[4]
        out$treatment[i]  <- grps[5]
      }
    }
    out$condition <- dplyr::coalesce(
      paste(out$celltype, out$treatment, sep = "_"),
      NA_character_
    )
    n_unparsed <- sum(is.na(out$celltype))
    if (n_unparsed > 0L)
      cli::cli_warn("{n_unparsed}/{nrow(out)} NMD-style columns did not parse.")
  }
  out
}

#' Load per-source long-read evidence for a gene
#'
#' For each source, runs [runIsoscopeGene()] to get the per-gene
#' isoform annotation, then loads the corresponding count matrix and
#' parses its sample columns. Returns a nested tibble with one row per
#' source and list-columns for annotation, counts (wide, isoform ×
#' sample matrix), and samples (per-column metadata).
#'
#' `celltypes` and `treatments` filters apply to the NMD source's
#' count matrix at load time — useful for restricting to DD / DD_ALI
#' under Smg1i / DMSO per the user's focus. HAEC-185 sources ignore
#' these filters (they have only "baseline" condition).
#'
#' @param gene HGNC symbol or Ensembl ID.
#' @param sources Named list of source descriptors; default
#'   [longreadSources()].
#' @param celltypes Optional character vector of cell types to retain
#'   (NMD only).
#' @param treatments Optional character vector of treatments to retain
#'   (NMD only).
#' @return A tibble with one row per source and columns:
#'   `source_id`, `cohort`, `tissue`, `calling_stage`, `n_isoforms`,
#'   `n_samples`, `annotation` (list of tibbles),
#'   `counts` (list of matrices), `samples` (list of tibbles).
#'   Carries a merged `isovar_meta` attribute.
#' @export
loadLongreadEvidence <- function(gene,
                                 sources = longreadSources(),
                                 celltypes = NULL,
                                 treatments = NULL) {
  rows <- vector("list", length(sources))
  names(rows) <- names(sources)
  metas <- list()

  for (src_id in names(sources)) {
    s <- sources[[src_id]]
    cli::cli_inform(c("i" = "Loading long-read evidence from {.val {src_id}}"))

    # 1. isoscope annotation
    ann <- runIsoscopeGene(gene, s$isoscope_config,
                           output_dir = file.path(tempdir(),
                                                  paste0("isovar_lr_", src_id, "_", gene)))
    meta <- getIsovarMeta(ann, require = FALSE)
    if (!is.null(meta)) metas[[src_id]] <- meta

    # 2. count matrix
    if (!file.exists(s$count_matrix_path)) {
      cli::cli_warn("Count matrix not found for {.val {src_id}}: {.path {s$count_matrix_path}}")
      rows[[src_id]] <- tibble::tibble(
        source_id = src_id, cohort = s$cohort, tissue = s$tissue,
        calling_stage = s$calling_stage,
        n_isoforms = nrow(ann), n_samples = 0L,
        annotation = list(ann), counts = list(NULL), samples = list(NULL)
      )
      next
    }
    cm <- readr::read_delim(s$count_matrix_path, delim = s$count_delim,
                            show_col_types = FALSE, progress = FALSE)

    # Filter count matrix to this gene's isoforms (annotation has source column)
    iso_ids <- ann$isoform_id
    cm_gene <- cm[cm$id %in% iso_ids, , drop = FALSE]

    # 3. sample metadata from column names
    samples <- parseSampleColumns(names(cm_gene), s$sample_schema)
    if (s$sample_schema == "sample_celltype_donor_treatment") {
      if (!is.null(celltypes))
        samples <- samples[samples$celltype %in% celltypes, , drop = FALSE]
      if (!is.null(treatments))
        samples <- samples[samples$treatment %in% treatments, , drop = FALSE]
    }

    # 4. counts matrix (isoform × filtered sample)
    sample_cols <- intersect(samples$column_name, names(cm_gene))
    counts_mat  <- as.matrix(cm_gene[, sample_cols, drop = FALSE])
    rownames(counts_mat) <- cm_gene$id

    rows[[src_id]] <- tibble::tibble(
      source_id     = src_id,
      cohort        = s$cohort,
      tissue        = s$tissue,
      calling_stage = s$calling_stage,
      n_isoforms    = nrow(ann),
      n_samples     = length(sample_cols),
      annotation    = list(ann),
      counts        = list(counts_mat),
      samples       = list(samples[samples$column_name %in% sample_cols, , drop = FALSE])
    )
  }

  out <- dplyr::bind_rows(rows)
  if (length(metas)) {
    merged <- Reduce(mergeIsovarMeta, metas)
    setIsovarMeta(out, merged)
  } else {
    out
  }
}
