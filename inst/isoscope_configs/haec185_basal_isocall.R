# isovar -> isoscope source config: HAEC-185 basal ISOCALL (pre-SQANTI).
#
# Uses isoscope's isocall mode (SOURCE_KIND = "isocall") — enumerates
# isoforms directly from the long-read GTF without a SQANTI classification
# file. Useful for confirming whether SQANTI correction has altered
# junctions relative to the raw isocall calls.

SOURCE_KIND          <- "isocall"

ISOVAR_SOURCE_META <- list(
  source_id       = "haec185_basal_2024_9samples_isocall",
  source_label    = "HAEC-185 basal (2024 9-sample LR isocall, pre-SQANTI)",
  cohort          = "HAEC-185",
  tissue          = "basal airway epithelium",
  n_samples       = 9,
  genome_build    = "GRCh38",
  gencode_version = "GENCODE_v49",
  sqanti_run_date = NA_character_,
  conditions      = NULL,
  calling_stage   = "isocall_uncorrected"
)

GENCODE_GTF_INDEXED   <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.sorted.gtf.gz"
GENCODE_FASTA         <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.transcripts.fa.gz"

# No SQANTI_CLASSIFICATION in isocall mode.
SQANTI_CLASSIFICATION <- NULL
SQANTI_GTF_INDEXED    <- "/Users/petecastaldi/claude_projects/HAEC185/longread/2024_9samples/isocall/20260304/call/LR_RNA_Basal_isocall.isoforms.sorted.gtf.gz"
SQANTI_FASTA          <- NULL
SQANTI_PROTEIN_FASTA  <- NULL

DGELIST_RDS           <- NULL
TRANSCRIPT_ID_COLUMN  <- "transcript_id"
STRATIFY_BY           <- NULL
LEVEL_LABELS          <- NULL
