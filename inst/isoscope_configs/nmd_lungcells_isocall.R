# isovar -> isoscope source config: NMD lungcells ISOCALL (pre-SQANTI).

SOURCE_KIND          <- "isocall"

ISOVAR_SOURCE_META <- list(
  source_id       = "nmd_lungcells_isocall",
  source_label    = "NMD lungcells (38-sample LR isocall, pre-SQANTI)",
  cohort          = "NMD lung cell lines (Randell_Lung_Cells_2025)",
  tissue          = "lung cell lines (DD, DD_ALI, AT, DO, FB, MV)",
  n_samples       = 38,
  genome_build    = "GRCh38",
  gencode_version = "GENCODE_v49",
  sqanti_run_date = NA_character_,
  conditions      = c("DMSO", "Smg1i"),
  calling_stage   = "isocall_uncorrected"
)

GENCODE_GTF_INDEXED   <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.sorted.gtf.gz"
GENCODE_FASTA         <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.transcripts.fa.gz"

SQANTI_CLASSIFICATION <- NULL
SQANTI_GTF_INDEXED    <- "/Users/petecastaldi/claude_projects/nmd/isocall/nmd_lungcells/results/call/nmd_isocall.isoforms.sorted.gtf.gz"
SQANTI_FASTA          <- NULL
SQANTI_PROTEIN_FASTA  <- NULL

DGELIST_RDS           <- NULL
TRANSCRIPT_ID_COLUMN  <- "transcript_id"
STRATIFY_BY           <- NULL
LEVEL_LABELS          <- NULL
