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
