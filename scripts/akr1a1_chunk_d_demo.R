#!/usr/bin/env Rscript
# Chunk D review artifact — AKR1A1 long-read evidence across all 4
# sources (HAEC-185 basal × {SQANTI, isocall}, NMD lungcells × {SQANTI,
# isocall}), filtered to DD + DD_ALI under Smg1i + DMSO for the NMD
# cohort. Shows per-source isoform counts, phenotype-intron usage, and
# per-condition sample totals.

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(tidyr); library(readr); library(tibble); library(cli); library(jsonlite)
})
isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

cli::cli_h1("AKR1A1 long-read evidence — 4 sources, per-condition counts")

lr <- loadLongreadEvidence("AKR1A1",
                           celltypes  = c("DD", "DD_ALI"),
                           treatments = c("DMSO", "Smg1i"))

# Compact top-line summary
cli::cli_h2("Per-source shape")
print(as.data.frame(lr %>%
  select(source_id, cohort, tissue, calling_stage,
         n_isoforms, n_samples)), row.names = FALSE)

# For each source, compute total reads on the phenotype intron and on
# the shifted sites, aggregated across samples (filtered to the
# celltype / treatment subset where relevant).
cli::cli_h2("Per-source, per-condition phenotype-intron read totals")
rows_out <- list()
for (i in seq_len(nrow(lr))) {
  ann <- lr$annotation[[i]]
  cnt <- lr$counts[[i]]
  smp <- lr$samples[[i]]
  if (is.null(cnt) || nrow(cnt) == 0L) next

  uses_phenotype <- sapply(ann$junctions, function(j) {
    if (is.na(j) || !nzchar(j)) return(FALSE)
    any(grepl("^45551155_45552436_[+-]",
              strsplit(j, ";", fixed = TRUE)[[1]]))
  })
  iso_ids_phen <- ann$isoform_id[uses_phenotype]
  rows_in_cnt  <- intersect(iso_ids_phen, rownames(cnt))
  if (length(rows_in_cnt) == 0L) next
  sub_cnt <- cnt[rows_in_cnt, , drop = FALSE]

  # Aggregate by (celltype, treatment) — use parsed sample metadata
  # directly rather than re-parsing a concatenated condition string,
  # which breaks on underscores inside celltype names like DD_ALI.
  if (all(smp$condition == "baseline" | is.na(smp$condition))) {
    cols <- intersect(smp$column_name, colnames(sub_cnt))
    rows_out[[length(rows_out) + 1L]] <- tibble(
      source_id       = lr$source_id[i],
      calling_stage   = lr$calling_stage[i],
      celltype        = NA_character_,
      treatment       = NA_character_,
      n_samples       = length(cols),
      n_iso_phenotype = length(rows_in_cnt),
      total_reads     = as.integer(sum(sub_cnt[, cols, drop = FALSE]))
    )
  } else {
    strata <- dplyr::distinct(smp, celltype, treatment)
    for (j in seq_len(nrow(strata))) {
      cols <- smp$column_name[smp$celltype  == strata$celltype[j] &
                              smp$treatment == strata$treatment[j]]
      cols <- intersect(cols, colnames(sub_cnt))
      if (length(cols) == 0L) next
      rows_out[[length(rows_out) + 1L]] <- tibble(
        source_id       = lr$source_id[i],
        calling_stage   = lr$calling_stage[i],
        celltype        = strata$celltype[j],
        treatment       = strata$treatment[j],
        n_samples       = length(cols),
        n_iso_phenotype = length(rows_in_cnt),
        total_reads     = as.integer(sum(sub_cnt[, cols, drop = FALSE]))
      )
    }
  }
}
dd <- dplyr::bind_rows(rows_out)
print(as.data.frame(dd %>% arrange(source_id, celltype, treatment)), row.names = FALSE)

# NMD Smg1i vs DMSO contrast on phenotype-intron reads (DD + DD_ALI).
cli::cli_h2("NMD Smg1i vs DMSO contrast on phenotype-intron reads (DD + DD_ALI)")
nmd <- dd %>% filter(grepl("nmd", source_id))
if (nrow(nmd)) {
  tidy <- nmd %>%
    group_by(source_id, calling_stage, celltype, treatment) %>%
    summarise(total_reads = sum(total_reads),
              n_samples   = sum(n_samples), .groups = "drop") %>%
    arrange(source_id, celltype, treatment)
  print(as.data.frame(tidy), row.names = FALSE)
}

cli::cli_h2("Provenance")
meta <- getIsovarMeta(lr)
for (s in meta$sources)
  cat(sprintf("  [%s] %s  (cohort=%s, tissue=%s, stage=%s)\n",
              s$kind %||% "—", s$id %||% "—",
              s$cohort %||% "—", s$tissue %||% "—",
              s$calling_stage %||% "—"))
