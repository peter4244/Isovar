akr_gencode <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.sorted.gtf.gz"
akr_gencode_fa <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.transcripts.fa.gz"

test_that("parseStructuresMultiGtf handles a single GTF with tagging", {
  skip_on_cran()
  skip_if(!file.exists(akr_gencode), "GENCODE GTF not present")
  skip_if(!requireNamespace("Isopair", quietly = TRUE), "Isopair not installed")

  structs <- parseStructuresMultiGtf(
    gtfs        = list(gencode = akr_gencode),
    isoform_ids = c("ENST00000372070.7", "ENST00000621846.4")
  )
  expect_equal(nrow(structs), 2L)
  expect_true("source_gtf" %in% names(structs))
  expect_true(all(structs$source_gtf == "gencode"))
  expect_setequal(structs$isoform_id,
                  c("ENST00000372070.7", "ENST00000621846.4"))
})

test_that("parseStructuresMultiGtf warns on unparseable IDs", {
  skip_on_cran()
  skip_if(!file.exists(akr_gencode), "GENCODE GTF not present")
  skip_if(!requireNamespace("Isopair", quietly = TRUE), "Isopair not installed")

  expect_warning(
    out <- parseStructuresMultiGtf(
      gtfs        = list(gencode = akr_gencode),
      isoform_ids = c("ENST00000372070.7", "FAKE_ID_FOR_TEST")
    ),
    "Could not parse"
  )
  expect_equal(nrow(out), 1L)
})

test_that("loadTranscriptSequences returns a named character vector", {
  skip_on_cran()
  skip_if(!file.exists(akr_gencode_fa), "GENCODE FASTA not present")
  skip_if(!requireNamespace("Biostrings", quietly = TRUE))

  seqs <- loadTranscriptSequences(
    fastas      = akr_gencode_fa,
    isoform_ids = c("ENST00000372070.7", "ENST00000621846.4")
  )
  expect_type(seqs, "character")
  expect_true(length(seqs) >= 1L)
  expect_true(all(grepl("^[ACGTN]+$", seqs[1])))  # upper-case DNA
  expect_true(all(names(seqs) %in% c("ENST00000372070.7", "ENST00000621846.4")))
})

test_that("extractCdsMultiGtf returns the Isopair CDS schema + tag", {
  skip_on_cran()
  skip_if(!file.exists(akr_gencode), "GENCODE GTF not present")
  skip_if(!requireNamespace("Isopair", quietly = TRUE), "Isopair not installed")

  cds <- extractCdsMultiGtf(
    gtfs        = list(gencode = akr_gencode),
    isoform_ids = c("ENST00000372070.7", "ENST00000621846.4")
  )
  expect_true(all(c("isoform_id", "coding_status", "cds_start", "cds_stop",
                    "orf_length", "source_gtf") %in% names(cds)))
})

test_that("enumerateComparatorOrfs returns per-(isoform, ORF) rows", {
  skip_on_cran()
  skip_if(!file.exists(akr_gencode), "GENCODE GTF not present")
  skip_if(!file.exists(akr_gencode_fa), "GENCODE FASTA not present")
  skip_if(!requireNamespace("Isopair", quietly = TRUE), "Isopair not installed")

  ids <- c("ENST00000372070.7", "ENST00000621846.4")
  structs <- parseStructuresMultiGtf(list(gencode = akr_gencode), ids)
  cds     <- extractCdsMultiGtf(list(gencode = akr_gencode), ids)
  seqs    <- loadTranscriptSequences(akr_gencode_fa, ids)
  orfs    <- enumerateComparatorOrfs(structs, cds, seqs)

  expect_true(all(c("isoform_id", "atg_tx_pos", "orf_length",
                    "n_downstream_ejc", "is_annotated_cds", "category")
                  %in% names(orfs)))
  expect_true(all(orfs$isoform_id %in% ids))
  # Each isoform contributes multiple ORFs
  expect_gt(nrow(orfs), length(ids))
  # At least one annotated-CDS ORF should be present
  expect_true(any(orfs$is_annotated_cds))
})

test_that("summarizeOrfsToTranscript collapses to per-isoform NMD verdict", {
  orfs <- tibble::tibble(
    isoform_id       = c("A","A","A","B","B"),
    atg_tx_pos       = c(10L, 50L, 90L, 20L, 80L),
    stop_tx_pos      = c(40L, 80L, NA_integer_, 50L, 110L),
    orf_length       = c(30L, 30L, NA_integer_, 30L, 30L),
    n_downstream_ejc = c(2L,  0L,  NA_integer_, 1L,  0L),
    is_annotated_cds = c(FALSE, TRUE, FALSE, FALSE, FALSE),
    category         = c("effectively_ptc", "no_downstream_ejc",
                         "no_stop_in_frame",
                         "effectively_ptc", "no_downstream_ejc")
  )
  summ <- summarizeOrfsToTranscript(orfs)
  expect_equal(nrow(summ), 2L)
  expect_equal(summ$n_orfs,     c(3L, 2L))
  expect_equal(summ$n_ptc_orfs, c(1L, 1L))
  expect_equal(summ$any_ptc,    c(TRUE, TRUE))
  expect_equal(summ$n_no_stop,  c(1L, 0L))
})

test_that("hiddenPtcStatus forwards to Isopair::traceReferenceAtg", {
  skip_on_cran()
  skip_if(!file.exists(akr_gencode), "GENCODE GTF not present")
  skip_if(!file.exists(akr_gencode_fa), "GENCODE FASTA not present")
  skip_if(!requireNamespace("Isopair", quietly = TRUE), "Isopair not installed")

  ids <- c("ENST00000372070.7", "ENST00000621846.4")
  structs <- parseStructuresMultiGtf(list(gencode = akr_gencode), ids)
  cds     <- extractCdsMultiGtf(list(gencode = akr_gencode), ids)
  seqs    <- loadTranscriptSequences(akr_gencode_fa, ids)

  pairs <- tibble::tibble(
    gene_id               = "ENSG00000117448.15",
    reference_isoform_id  = "ENST00000372070.7",
    comparator_isoform_id = "ENST00000621846.4"
  )
  out <- hiddenPtcStatus(pairs, structs, cds, seqs)
  expect_true("category" %in% names(out))
  expect_equal(nrow(out), 1L)
})
