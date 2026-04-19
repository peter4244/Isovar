#!/usr/bin/env Rscript
# Chunk E review artifact — end-to-end AKR1A1 run through Isopair.
#   - Pull the 72 AKR1A1 isoforms (GENCODE + NMD-SQANTI novels).
#   - Pick the dominant isoform in NMD (identifyDominantIsoforms).
#   - Pair dominant vs each coding-tractable comparator.
#   - Parse structures + CDS across both GTFs; load transcript sequences
#     from both FASTAs.
#   - Compute standard PTC status (computePtcStatus).
#   - Compute hidden PTC categories (traceReferenceAtg) for comparators
#     whose own CDS was missing or PTC-negative.
#   - Cross-reference with the Chunk D DD_ALI Smg1i rescue signal.

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(tidyr); library(readr); library(tibble); library(cli); library(jsonlite)
  library(Isopair); library(Biostrings)
})

isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

out_dir <- file.path(isovar_dir, "runs/akr1a1/chunk_e")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gencode_gtf <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.sorted.gtf.gz"
gencode_fa  <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.transcripts.fa.gz"
nmd_gtf     <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.sorted.gtf.gz"
nmd_fa      <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"

cli::cli_h1("AKR1A1 Isopair end-to-end on NMD-SQANTI")

# 1. Long-read evidence (Chunk D)
lr <- loadLongreadEvidence("AKR1A1",
                           celltypes = c("DD","DD_ALI"),
                           treatments = c("DMSO","Smg1i"))
nmd_idx <- which(lr$source_id == "nmd_sqanti")
ann <- lr$annotation[[nmd_idx]]
cnt <- lr$counts[[nmd_idx]]
smp <- lr$samples[[nmd_idx]]
all_ids <- ann$isoform_id

cli::cli_h2("Dominant isoform (Isopair::identifyDominantIsoforms)")
dominant <- computeDominantIsoform(lr, source_id = "nmd_sqanti", threshold = 0.3)
print(as.data.frame(dominant), row.names = FALSE)
ref_id <- dominant$dominant_isoform_id[1]
cli::cli_alert_success("Reference isoform: {.val {ref_id}}")

# 2. Structures + CDS across GENCODE + SQANTI
cli::cli_h2("Parse structures + CDS across GENCODE + SQANTI GTFs")
gtfs  <- list(gencode = gencode_gtf, nmd_sqanti = nmd_gtf)
structs <- parseStructuresMultiGtf(gtfs, all_ids)
cds     <- extractCdsMultiGtf(gtfs, all_ids)
cli::cli_alert_info("structures: {nrow(structs)}   cds: {nrow(cds)}")
cli::cli_alert_info("  coding: {sum(cds$coding_status == 'coding')}   \\
                    non-coding: {sum(cds$coding_status != 'coding')}")

# 3. Sequences — GENCODE fasta for known ENSTs, SQANTI fasta for novels
cli::cli_h2("Load transcript sequences")
seqs <- loadTranscriptSequences(c(gencode_fa, nmd_fa), all_ids)
cli::cli_alert_info("sequences loaded: {length(seqs)} / {length(all_ids)}")

# 4. Standard PTC status (50-nt rule on each isoform's own CDS)
cli::cli_h2("computePtcStatus (standard 50-nt rule)")
ptc_std <- Isopair::computePtcStatus(structs, cds)
cli::cli_alert_info("PTC+ isoforms (standard): {sum(ptc_std$has_ptc, na.rm = TRUE)} / {nrow(ptc_std)}")

# 5. Build pairs: dominant vs every other isoform
comp_ids <- setdiff(all_ids, ref_id)
pairs <- tibble(
  gene_id               = "ENSG00000117448.15",
  reference_isoform_id  = ref_id,
  comparator_isoform_id = comp_ids
)
# Filter to comparators we have structures for
pairs <- pairs[pairs$comparator_isoform_id %in% structs$isoform_id, ]
cli::cli_alert_info("pairs built: {nrow(pairs)} (dominant vs each comparator)")

# 6. Hidden PTC via reference ATG trace
cli::cli_h2("hiddenPtcStatus (traceReferenceAtg on all pairs)")
hidden <- hiddenPtcStatus(pairs, structs, cds, seqs)
cat("category counts:\n")
print(table(hidden$category, useNA = "always"))

# 7. Cross-check with Smg1i rescue signal from Chunk D
cli::cli_h2("Cross-check: hidden_ptc vs DD_ALI Smg1i rescue")
# Compute per-isoform Smg1i/DMSO ratio in DD_ALI
dmso_cols <- smp$column_name[smp$celltype == "DD_ALI" & smp$treatment == "DMSO"]
smg1i_cols <- smp$column_name[smp$celltype == "DD_ALI" & smp$treatment == "Smg1i"]
iso_in_cnt <- intersect(all_ids, rownames(cnt))
per_iso <- tibble(
  isoform_id = iso_in_cnt,
  dmso_mean  = rowMeans(cnt[iso_in_cnt, dmso_cols,  drop = FALSE]),
  smg1i_mean = rowMeans(cnt[iso_in_cnt, smg1i_cols, drop = FALSE])
) %>%
  mutate(smg1i_over_dmso = ifelse(dmso_mean > 0, smg1i_mean / dmso_mean, NA_real_))

# Join hidden-PTC category
merged <- hidden %>%
  select(comparator_isoform_id, category, comp_orf_length, n_downstream_ejc) %>%
  left_join(per_iso, by = c("comparator_isoform_id" = "isoform_id"))

# Show top isoforms ranked by Smg1i rescue magnitude
cli::cli_h2("Top 10 isoforms by DD_ALI Smg1i/DMSO rescue")
show <- merged %>%
  filter(!is.na(smg1i_over_dmso) & is.finite(smg1i_over_dmso)) %>%
  arrange(desc(smg1i_over_dmso)) %>%
  head(10) %>%
  mutate(dmso  = sprintf("%.1f", dmso_mean),
         smg1i = sprintf("%.1f", smg1i_mean),
         ratio = sprintf("%.2fx", smg1i_over_dmso)) %>%
  select(isoform = comparator_isoform_id, category,
         orf_len = comp_orf_length, n_ejc_down = n_downstream_ejc,
         dmso, smg1i, ratio)
print(as.data.frame(show), row.names = FALSE)

# 8. Write outputs for review
write_tsv(hidden, file.path(out_dir, "hidden_ptc_status.tsv"))
write_tsv(merged, file.path(out_dir, "hidden_ptc_vs_smg1i_rescue.tsv"))
write_tsv(ptc_std, file.path(out_dir, "standard_ptc_status.tsv"))
cli::cli_alert_success("Artifacts written to {.path {out_dir}}")
