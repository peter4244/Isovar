#!/usr/bin/env Rscript
# Chunk E.2 review artifact — per-ORF classifications for the AKR1A1
# isoforms most rescued by Smg1i in DD_ALI. Compares per-ORF verdicts
# (Isopair::enumerateOrfs via enumerateComparatorOrfs) with the
# per-transcript verdict from hiddenPtcStatus.

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(readr); library(tibble); library(cli); library(jsonlite)
  library(Isopair); library(Biostrings)
})
isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

gencode_gtf <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.sorted.gtf.gz"
gencode_fa  <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.transcripts.fa.gz"
nmd_gtf     <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.sorted.gtf.gz"
nmd_fa      <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"

lr <- loadLongreadEvidence("AKR1A1", celltypes = c("DD","DD_ALI"),
                           treatments = c("DMSO","Smg1i"))
idx <- which(lr$source_id == "nmd_sqanti")
ann <- lr$annotation[[idx]]; cnt <- lr$counts[[idx]]; smp <- lr$samples[[idx]]
all_ids <- ann$isoform_id

gtfs <- list(gencode = gencode_gtf, nmd_sqanti = nmd_gtf)
structs <- parseStructuresMultiGtf(gtfs, all_ids)
cds     <- extractCdsMultiGtf(gtfs, all_ids)
seqs    <- loadTranscriptSequences(c(gencode_fa, nmd_fa), all_ids)

# Top 10 by DD_ALI Smg1i/DMSO ratio — pick them straight from the count matrix
dmso <- smp$column_name[smp$celltype == "DD_ALI" & smp$treatment == "DMSO"]
smg  <- smp$column_name[smp$celltype == "DD_ALI" & smp$treatment == "Smg1i"]
iso_in_cnt <- intersect(all_ids, rownames(cnt))
per_iso <- tibble(
  isoform_id = iso_in_cnt,
  dmso_mean  = rowMeans(cnt[iso_in_cnt, dmso, drop = FALSE]),
  smg_mean   = rowMeans(cnt[iso_in_cnt, smg,  drop = FALSE])
) %>%
  mutate(ratio = ifelse(dmso_mean > 0, smg_mean / dmso_mean, NA)) %>%
  filter(is.finite(ratio)) %>%
  arrange(desc(ratio)) %>%
  head(10)

cli::cli_h1("AKR1A1 — per-ORF enumeration on the top 10 Smg1i-rescued isoforms")
top_structs <- structs[structs$isoform_id %in% per_iso$isoform_id, ]
top_cds     <- cds[cds$isoform_id %in% per_iso$isoform_id, ]
top_seqs    <- seqs[names(seqs) %in% per_iso$isoform_id]

orfs <- enumerateComparatorOrfs(top_structs, top_cds, top_seqs)

cat(sprintf("\ntotal ORFs across %d isoforms: %d\n",
            length(unique(orfs$isoform_id)), nrow(orfs)))
cat("category counts:\n")
print(as.data.frame(as.table(table(orfs$category))))

cat("\n--- per-isoform ORF counts + roll-up verdict ---\n")
summary <- summarizeOrfsToTranscript(orfs) %>%
  left_join(per_iso, by = "isoform_id") %>%
  arrange(desc(ratio)) %>%
  mutate(
    dmso_mean = sprintf("%.1f", dmso_mean),
    smg_mean  = sprintf("%.1f", smg_mean),
    ratio     = sprintf("%.1fx", ratio)
  ) %>%
  select(isoform_id, n_orfs, n_ptc_orfs, n_no_stop,
         any_ptc, dmso_mean, smg_mean, ratio)
print(as.data.frame(summary), row.names = FALSE)

# Write artifacts
out_dir <- file.path(isovar_dir, "runs/akr1a1/chunk_e2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_tsv(orfs,    file.path(out_dir, "per_orf_classifications_top10.tsv"))
write_tsv(summary, file.path(out_dir, "orf_rollup_vs_rescue_top10.tsv"))
cli::cli_alert_success("Artifacts written to {.path {out_dir}}")
