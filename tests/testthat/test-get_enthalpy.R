# tests/testthat/test-get_enthalpy.R
# Tests for get_enthalpy() — multi-source retrieval without merging.

test_that("get_enthalpy returns a data.frame with expected columns", {
  res <- get_enthalpy("CH4", phase = "g")
  expect_s3_class(res, "data.frame")
  expected_cols <- c(
    "source", "name", "formula", "CAS", "phase", "compound_nature",
    "delta_fH", "delta_fG", "S", "Cp",
    "temperature", "uncertainty", "unit",
    "source_id", "reference"
  )
  expect_true(all(expected_cols %in% names(res)))
})

test_that("get_enthalpy finds methane in at least one source", {
  res <- get_enthalpy("CH4", phase = "g")
  expect_gt(nrow(res), 0)
  expect_true(any(!is.na(res$delta_fH)))
})

test_that("get_enthalpy is case-insensitive (exact = TRUE)", {
  res_lower <- get_enthalpy("methane", phase = "g", exact = TRUE)
  res_upper <- get_enthalpy("Methane", phase = "g", exact = TRUE)
  expect_equal(nrow(res_lower), nrow(res_upper))
})

test_that("get_enthalpy with source argument restricts the query", {
  res_atct <- get_enthalpy("CH4", source = "ATcT", phase = "g")
  expect_true(all(res_atct$source == "ATcT"))
})

test_that("get_enthalpy returns one row per source when species exists", {
  res <- get_enthalpy("CH4", phase = "g")
  # Au moins 2 sources parmi les 3 doivent avoir CH4
  expect_gte(length(unique(res$source)), 2)
})

test_that("get_enthalpy returns empty result for unknown compound", {
  res <- suppressMessages(get_enthalpy("XxXxXx_not_a_compound"))
  expect_true(is.null(res) || NROW(res) == 0)
})

test_that("get_enthalpy preserves NA for unavailable properties", {
  res <- get_enthalpy("CH4", source = "ATcT", phase = "g")
  # ATcT ne fournit pas delta_fG, S, Cp
  expect_true(all(is.na(res$delta_fG)))
  expect_true(all(is.na(res$S)))
  expect_true(all(is.na(res$Cp)))
})

test_that("get_enthalpy filters by phase", {
  res_g  <- get_enthalpy("water", phase = "g")
  res_l  <- get_enthalpy("water", phase = "l")
  res_aq <- get_enthalpy("water", phase = "aq")
  if (nrow(res_g)  > 0) expect_true(all(res_g$phase  == "g"))
  if (nrow(res_l)  > 0) expect_true(all(res_l$phase  == "l"))
  if (nrow(res_aq) > 0) expect_true(all(res_aq$phase == "aq"))
})