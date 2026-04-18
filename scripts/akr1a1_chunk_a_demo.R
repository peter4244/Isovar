#!/usr/bin/env Rscript
# Chunk A review artifact — AKR1A1 credible set with splaire + gnomAD + COPD GWAS
# columns joined, sorted by |splaire delta| descending.
#
# Uses the AKR1A1-subset test fixtures so this runs instantly (no 190 MB GWAS
# load). Produces a single human-readable TSV:
#   runs/akr1a1_chunk_a.tsv
#
# Run:
#   Rscript scripts/akr1a1_chunk_a_demo.R

suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicRanges); library(IRanges)
  library(dplyr); library(readr); library(tibble); library(cli)
})

isovar_dir <- "/Users/petecastaldi/claude_projects/isovar"
for (f in list.files(file.path(isovar_dir, "R"), full.names = TRUE)) source(f)

splaire_fixture <- file.path(isovar_dir, "inst/extdata/akr1a1_splaire_test.tsv.gz")
gwas_fixture    <- file.path(isovar_dir, "inst/extdata/gwas_akr1a1_test.tsv.gz")

# Hand-mapped rsIDs for the AKR1A1 credible set (from gnomAD annotation).
# This avoids hitting the network in the demo; in a real run, these would
# come from annotateGnomad().
rs_map <- c(
  "chr1:45480964:G:T"           = "rs4660861",
  "chr1:45491367:G:T"           = NA,
  "chr1:45495800:A:T"           = NA,
  "chr1:45506087:T:C"           = NA,
  "chr1:45508256:G:A"           = NA,
  "chr1:45517408:C:T"           = "rs11211129",
  "chr1:45521033:G:A"           = NA,
  "chr1:45538492:C:T"           = "rs11211135",
  "chr1:45539557:T:C"           = "rs11802397",
  "chr1:45540036:G:T"           = "rs11211136",
  "chr1:45540058:G:C"           = "rs7526571",
  "chr1:45540513:G:GA"          = NA,
  "chr1:45543241:G:A"           = "rs60608374",
  "chr1:45543641:G:C"           = "rs2356553",
  "chr1:45543646:C:T"           = "rs2356551",
  "chr1:45544274:G:T"           = "rs512026",
  "chr1:45545333:C:T"           = "rs2991970",
  "chr1:45547536:C:T"           = NA,
  "chr1:45548061:G:C"           = "rs2993262",
  "chr1:45549863:A:ATATCTG"     = "rs61467610",  # INS, absent in GWAS
  "chr1:45551015:T:C"           = "rs9147",
  "chr1:45542866:TTTGAACTTCG:T" = "rs66922050",  # DEL, absent in GWAS
  "chr1:45552164:C:T"           = "rs2934859",
  "chr1:45553048:C:A"           = "rs4660876",
  "chr1:45554960:G:A"           = "rs2991973"
)

cli::cli_h1("AKR1A1 Chunk A demo — splaire + GWAS join")

ranked <- rankSplaireVariants(splaire_fixture, gene = "AKR1A1")
ranked$rsid <- unname(rs_map[ranked$variant_id])

gwas <- loadCopdGwas(gwas_fixture)
ann  <- annotateGwas(ranked, gwas)

# Compact display table
out <- ann %>%
  mutate(variant_kind = ifelse(nchar(ref) == 1L & nchar(alt) == 1L, "SNP", "indel")) %>%
  select(variant_id, rsid, variant_kind, pos,
         pip = posterior_inclusion_probability,
         splaire_head = top_head, splaire_dir = top_direction,
         splaire_abs = top_abs_delta,
         gwas_beta, gwas_se, gwas_p, gwas_eaf,
         gwas_alignment = gwas_alignment_status) %>%
  arrange(desc(splaire_abs))

# Quick summary lines
cli::cli_alert_info("Variants total: {nrow(out)}")
cli::cli_alert_info("splaire high-impact (|delta| >= 0.1): {sum(out$splaire_abs >= 0.1)}")
cli::cli_alert_info("GWAS-annotated: {sum(!is.na(out$gwas_beta))} of {nrow(out)} \\
                    ({sum(!is.na(out$gwas_beta) & out$variant_kind == 'SNP')} SNPs, \\
                    {sum(!is.na(out$gwas_beta) & out$variant_kind == 'indel')} indels)")
cli::cli_alert_info("Absent from GWAS: {sum(is.na(out$gwas_beta))} \\
                    (all {sum(is.na(out$gwas_beta) & out$variant_kind == 'indel')} indels \\
                    and {sum(is.na(out$gwas_beta) & out$variant_kind == 'SNP')} SNPs)")
cli::cli_alert_info("Alignment: {.val {sort(table(out$gwas_alignment), decreasing = TRUE)}}")

# Pretty print
cli::cli_h2("Credible set sorted by |splaire delta|")
print(out, n = nrow(out))

# Write TSV
out_dir <- file.path(isovar_dir, "runs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, "akr1a1_chunk_a.tsv")
write_tsv(out, out_path)
cli::cli_alert_success("Wrote {.path {out_path}}")
