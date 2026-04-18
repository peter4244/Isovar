test_that("longreadSources returns 4 canonical sources", {
  s <- longreadSources()
  expect_setequal(names(s), c("haec185_sqanti", "haec185_isocall",
                              "nmd_sqanti", "nmd_isocall"))
  for (nm in names(s)) {
    cfg <- s[[nm]]$isoscope_config
    expect_true(file.exists(cfg), info = paste("missing config:", cfg))
    expect_true(cfg != "")
    expect_true(s[[nm]]$sample_schema %in%
                  c("donor_only", "sample_celltype_donor_treatment"))
  }
})

test_that("parseSampleColumns donor_only decodes HAEC-185 style", {
  cols <- c("id", "DD071Q", "DD032OP2", "DD066J", "NJ45", "T36")
  parsed <- parseSampleColumns(cols, "donor_only")
  expect_equal(nrow(parsed), 5L)             # id column dropped
  expect_equal(parsed$donor, c("DD071Q", "DD032OP2", "DD066J", "NJ45", "T36"))
  expect_true(all(is.na(parsed$celltype)))
  expect_true(all(is.na(parsed$treatment)))
  expect_true(all(parsed$condition == "baseline"))
})

test_that("parseSampleColumns handles every real NMD column name", {
  nmd_columns <- c(
    "id",
    "Sample1_DD_ALI_001V_Smg1i",   "Sample2_DD_010R_Smg1i",
    "Sample3_DD_ALI_027U_Smg1i",   "Sample4_DD_ALI_029T_Smg1i",
    "Sample5_AT_001V_DMSO",        "Sample6_AT_001V_Smg1i",
    "Sample7_AT_027U_DMSO",        "Sample8_AT_027U_Smg1i",
    "Sample9_AT_029T_DMSO",        "Sample10_AT_029T_Smg1i",
    "Sample11_DD_ALI_001V_DMSO",   "Sample12_DD_010R_DMSO",
    "Sample13_DD_017Q_DMSO",       "Sample14_DD_017Q_Smg1i",
    "Sample15_DD_ALI_027U_DMSO",   "Sample16_DD_ALI_029T_DMSO",
    "Sample17_DD_047N_DMSO",       "Sample18_DD_047N_Smg1i",
    "Sample19_DD_066Q_DMSO",       "Sample20_DD_066Q_Smg1i",
    "Sample21_DO_001V_DMSO",       "Sample22_DO_001V_Smg1i",
    "Sample23_DO_027U_DMSO",       "Sample24_DO_027U_Smg1i",
    "Sample25_DO_029T_DMSO",       "Sample26_DO_029T_Smg1i",
    "Sample27_FB_001V_DMSO",       "Sample28_FB_001V_Smg1i",
    "Sample29_FB_027U_DMSO",       "Sample30_FB_027U_Smg1i",
    "Sample31_FB_029T_DMSO",       "Sample32_FB_029T_Smg1i",
    "Sample33_MV_001V_DMSO",       "Sample34_MV_001V_Smg1i",
    "Sample35_MV_027U_DMSO",       "Sample36_MV_027U_Smg1i",
    "Sample37_MV_029T_DMSO",       "Sample38_MV_029T_Smg1i"
  )
  parsed <- parseSampleColumns(nmd_columns, "sample_celltype_donor_treatment")
  expect_equal(nrow(parsed), 38L)
  expect_true(all(!is.na(parsed$celltype)))
  expect_true(all(!is.na(parsed$donor)))
  expect_true(all(!is.na(parsed$treatment)))
  # DD_ALI must not collapse into DD. Counts verified by hand from the 38
  # NMD sample columns: 6 DD_ALI, 8 DD, 6 AT, 6 DO, 6 FB, 6 MV.
  expect_equal(sum(parsed$celltype == "DD_ALI"), 6L)
  expect_equal(sum(parsed$celltype == "DD"),     8L)
  # All 6 cell types represented
  expect_setequal(unique(parsed$celltype),
                  c("DD_ALI", "DD", "AT", "DO", "FB", "MV"))
  # Both treatments
  expect_setequal(unique(parsed$treatment), c("DMSO", "Smg1i"))
})

test_that("loadLongreadEvidence returns per-source structure on AKR1A1", {
  skip_on_cran()
  isoscope_script <- path.expand("~/claude_projects/isoscope/gene_isoform_annotation.R")
  if (!file.exists(isoscope_script)) skip("isoscope checkout not available")
  # This hits all 4 isoscope runs; skip on network/load-limited environments.
  skip_on_ci()

  lr <- loadLongreadEvidence("AKR1A1")
  expect_s3_class(lr, "tbl_df")
  expect_equal(nrow(lr), 4L)
  expect_setequal(lr$source_id,
                  c("haec185_sqanti", "haec185_isocall", "nmd_sqanti", "nmd_isocall"))
  # All sources should have annotation rows
  expect_true(all(sapply(lr$annotation, nrow) > 0))
  # NMD source's samples should include both treatments
  nmd_idx <- which(lr$source_id == "nmd_sqanti")
  nmd_samples <- lr$samples[[nmd_idx]]
  expect_setequal(unique(nmd_samples$treatment), c("DMSO", "Smg1i"))
})

test_that("celltypes filter drops non-matching NMD samples", {
  skip_on_cran()
  isoscope_script <- path.expand("~/claude_projects/isoscope/gene_isoform_annotation.R")
  if (!file.exists(isoscope_script)) skip("isoscope checkout not available")
  skip_on_ci()

  lr <- loadLongreadEvidence("AKR1A1",
                             celltypes = c("DD", "DD_ALI"),
                             treatments = c("DMSO", "Smg1i"))
  nmd_idx <- which(lr$source_id == "nmd_sqanti")
  nmd_samples <- lr$samples[[nmd_idx]]
  expect_true(all(nmd_samples$celltype %in% c("DD", "DD_ALI")))
  # DD (8) + DD_ALI (6) across DMSO + Smg1i in the 38-column matrix = 14 samples
  expect_equal(nrow(nmd_samples), 14L)
})
