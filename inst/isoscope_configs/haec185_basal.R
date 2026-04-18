# isovar -> isoscope source config: HAEC-185 basal long-read run.
#
# Consumed by gene_isoform_annotation.R (canonical isoscope) via:
#   Rscript gene_isoform_annotation.R <GENE> --config <this_file>
# and by isovar via runIsoscopeGene(gene, "inst/isoscope_configs/haec185_basal.R").

ISOVAR_SOURCE_META <- list(
  source_id       = "haec185_basal_2024_9samples",
  source_label    = "HAEC-185 basal (2024 9-sample LR SQANTI3, 20260304)",
  cohort          = "HAEC-185",
  tissue          = "basal airway epithelium",
  n_samples       = 9,
  genome_build    = "GRCh38",
  gencode_version = "GENCODE_v49",
  sqanti_run_date = "2026-03-04",
  conditions      = NULL
)

GENCODE_GTF_INDEXED   <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.sorted.gtf.gz"
GENCODE_FASTA         <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.transcripts.fa.gz"

SQANTI_CLASSIFICATION <- "/Users/petecastaldi/claude_projects/HAEC185/longread/2024_9samples/isocall/sqanti/20260304/LR_RNA_Basal_classification.txt"
SQANTI_GTF_INDEXED    <- "/Users/petecastaldi/claude_projects/HAEC185/longread/2024_9samples/isocall/sqanti/20260304/LR_RNA_Basal_corrected.sorted.gtf.gz"
SQANTI_FASTA          <- "/Users/petecastaldi/claude_projects/HAEC185/longread/2024_9samples/isocall/sqanti/20260304/LR_RNA_Basal_corrected.fasta"
SQANTI_PROTEIN_FASTA  <- "/Users/petecastaldi/claude_projects/HAEC185/longread/2024_9samples/isocall/sqanti/20260304/LR_RNA_Basal_corrected.faa"

# Expression stratification handled in Chunk D via count-matrix masks
# — isoscope itself runs here with --no-expr.
DGELIST_RDS           <- NULL
TRANSCRIPT_ID_COLUMN  <- "transcript_id"
STRATIFY_BY           <- NULL
LEVEL_LABELS          <- NULL
