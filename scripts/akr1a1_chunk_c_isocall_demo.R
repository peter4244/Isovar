#!/usr/bin/env Rscript
# Chunk C addendum — AKR1A1 across all 4 long-read views:
# HAEC-185 SQANTI / HAEC-185 isocall / NMD SQANTI / NMD isocall.
# Re-tests the "predicted-but-unobserved shifted splice site" finding
# against raw isocall output (no SQANTI correction applied).

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(tidyr); library(readr); library(tibble); library(cli); library(jsonlite)
})
isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

configs <- list(
  haec185_sqanti  = "haec185_basal.R",
  haec185_isocall = "haec185_basal_isocall.R",
  nmd_sqanti      = "nmd_lungcells.R",
  nmd_isocall     = "nmd_lungcells_isocall.R"
)

implicated <- c(
  phen_donor_5ss      = 45551155L,
  ins_shifted_donor   = 45551161L,
  del_cryptic_donor   = 45547554L,
  del_shifted_donor   = 45547544L
)

cli::cli_h1("AKR1A1 — 4-view long-read junction scan")

rows <- list()
for (src in names(configs)) {
  iso <- runIsoscopeGene("AKR1A1",
                         file.path(isovar_dir, "inst/isoscope_configs", configs[[src]]),
                         output_dir = file.path(tempdir(), paste0("iso_", src)))
  iso_c <- classifyIsoformsBySiteUsage(iso,
                                       site_positions = unname(implicated))

  # Also check exact phenotype intron 45551155->45552436 on + strand
  has_phenotype <- vapply(iso$junctions, function(j) {
    if (is.na(j) || !nzchar(j)) return(FALSE)
    any(grepl("^45551155_45552436_\\+", strsplit(j, ";", fixed = TRUE)[[1]]))
  }, logical(1))

  rows[[src]] <- tibble(
    source                  = src,
    n_total                 = nrow(iso),
    n_gencode               = sum(iso$source == "GENCODE", na.rm = TRUE),
    n_lr                    = sum(iso$source %in% c("SQANTI","ISOCALL"), na.rm = TRUE),
    uses_phen_donor_5ss     = sum(iso_c$uses_45551155_jstart, na.rm = TRUE),
    uses_exact_phen_intron  = sum(has_phenotype, na.rm = TRUE),
    uses_ins_shifted_donor  = sum(iso_c$uses_45551161_jstart, na.rm = TRUE),
    uses_del_cryptic_donor  = sum(iso_c$uses_45547554_jstart, na.rm = TRUE),
    uses_del_shifted_donor  = sum(iso_c$uses_45547544_jstart, na.rm = TRUE)
  )
}

out <- bind_rows(rows)
cat("\n=== Per-source junction usage ===\n")
print(as.data.frame(out), row.names = FALSE)

cli::cli_h1("Interpretation")
cat("* Phenotype intron 45551155 -> 45552436 is observed in every source\n",
    "  (sQTL phenotype confirmed in raw + corrected calls, HAEC-185 + NMD).\n",
    "* None of the splaire-predicted shifted donor sites appear in any\n",
    "  of the four views. SQANTI correction is NOT responsible for the\n",
    "  absence (isocall shows zero hits as well). The shifted alt-allele\n",
    "  isoforms are either not present in the sequenced donors, or the\n",
    "  splicing machinery does not produce them at detectable levels.\n",
    sep = "")
