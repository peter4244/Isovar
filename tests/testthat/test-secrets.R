test_that("loadIsovarSecrets reads KEY=VALUE lines and sets env vars", {
  tmp <- tempfile(fileext = ".env")
  writeLines(c(
    "# comment line",
    "",
    "ISOVAR_TEST_KEY_A=alpha",
    "ISOVAR_TEST_KEY_B = beta value",
    "  ISOVAR_TEST_KEY_C=\"gamma in quotes\"  "
  ), tmp)
  withr::with_envvar(c(
    ISOVAR_TEST_KEY_A = "",
    ISOVAR_TEST_KEY_B = "",
    ISOVAR_TEST_KEY_C = ""
  ), {
    keys <- loadIsovarSecrets(tmp, quiet = TRUE)
    expect_setequal(keys, c("ISOVAR_TEST_KEY_A",
                            "ISOVAR_TEST_KEY_B",
                            "ISOVAR_TEST_KEY_C"))
    expect_equal(Sys.getenv("ISOVAR_TEST_KEY_A"), "alpha")
    expect_equal(Sys.getenv("ISOVAR_TEST_KEY_B"), "beta value")
    expect_equal(Sys.getenv("ISOVAR_TEST_KEY_C"), "gamma in quotes")
  })
})

test_that("loadIsovarSecrets does not overwrite existing env vars by default", {
  tmp <- tempfile(fileext = ".env")
  writeLines("ISOVAR_TEST_PRESET=from_file", tmp)
  withr::with_envvar(c(ISOVAR_TEST_PRESET = "from_shell"), {
    keys <- loadIsovarSecrets(tmp, quiet = TRUE)
    expect_equal(Sys.getenv("ISOVAR_TEST_PRESET"), "from_shell")
    expect_false("ISOVAR_TEST_PRESET" %in% keys)
  })
})

test_that("overwrite=TRUE forces replacement", {
  tmp <- tempfile(fileext = ".env")
  writeLines("ISOVAR_TEST_PRESET2=from_file", tmp)
  withr::with_envvar(c(ISOVAR_TEST_PRESET2 = "from_shell"), {
    keys <- loadIsovarSecrets(tmp, overwrite = TRUE, quiet = TRUE)
    expect_equal(Sys.getenv("ISOVAR_TEST_PRESET2"), "from_file")
    expect_true("ISOVAR_TEST_PRESET2" %in% keys)
  })
})

test_that("loadIsovarSecrets handles a missing file gracefully", {
  expect_silent(
    keys <- loadIsovarSecrets("/tmp/__no_such_secrets_file__.env", quiet = TRUE)
  )
  expect_equal(length(keys), 0L)
})

test_that("malformed lines are skipped, not raised", {
  tmp <- tempfile(fileext = ".env")
  writeLines(c("VALID=1", "not_a_kv_line", "=novalue", "123BAD=start"), tmp)
  keys <- loadIsovarSecrets(tmp, overwrite = TRUE, quiet = TRUE)
  expect_equal(keys, "VALID")
})
