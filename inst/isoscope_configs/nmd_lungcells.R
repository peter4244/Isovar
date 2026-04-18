# isovar -> isoscope source config: NMD lungcells long-read run.
#
# Single merged SQANTI3 run across 38 samples × 6 cell types
# (DD, DD_ALI, AT, DO, FB, MV) × 2 treatments (DMSO, Smg1i).
# Per-cell-type / per-condition splits happen at the count-matrix
# column level in Chunk D, not at the isoscope level.

ISOVAR_SOURCE_META <- list(
  source_id       = "nmd_lungcells",
  source_label    = "NMD lungcells (38-sample LR SQANTI3, multi-cell-type, Smg1i/DMSO)",
  cohort          = "NMD lung cell lines (Randell_Lung_Cells_2025)",
  tissue          = "lung cell lines (DD, DD_ALI, AT, DO, FB, MV)",
  n_samples       = 38,
  genome_build    = "GRCh38",
  gencode_version = "GENCODE_v49",
  sqanti_run_date = NA_character_,
  conditions      = c("DMSO", "Smg1i")
)

GENCODE_GTF_INDEXED   <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.sorted.gtf.gz"
GENCODE_FASTA         <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.transcripts.fa.gz"

SQANTI_CLASSIFICATION <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_classification.txt"
SQANTI_GTF_INDEXED    <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.sorted.gtf.gz"
SQANTI_FASTA          <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.fasta"
SQANTI_PROTEIN_FASTA  <- "/Users/petecastaldi/claude_projects/nmd/sqanti/nmd_lungcells/results/nmd_lungcells_corrected.faa"

DGELIST_RDS           <- NULL
TRANSCRIPT_ID_COLUMN  <- "transcript_id"
STRATIFY_BY           <- NULL
LEVEL_LABELS          <- NULL
