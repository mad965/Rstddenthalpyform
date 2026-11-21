# tests/testthat/test-prediction-model.R
# Tests for the linear regression model on elemental counts.

test_that("count_elements returns a data.frame sorted by frequency", {
  df <- data.frame(
    formula = c("CH4", "C2H6", "H2O", "CO2"),
    delta_fH = c(-74.5, -84.0, -285.8, -393.5),
    stringsAsFactors = FALSE
  )
  res <- count_elements(df)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("element", "n_composes") %in% names(res)))
  expect_true("C" %in% res$element)
  expect_true("H" %in% res$element)
  expect_true("O" %in% res$element)
  # Trié par fréquence décroissante
  expect_true(all(diff(res$n_composes) <= 0))
})

test_that("build_design_matrix produces the right dimensions", {
  df <- data.frame(
    formula = c("CH4", "C2H6", "H2O"),
    delta_fH = c(-74.5, -84.0, -285.8),
    phase = c("g", "g", "l"),
    stringsAsFactors = FALSE
  )
  elements <- c("C", "H", "O")
  dm <- build_design_matrix(df, elements)
  expect_equal(nrow(dm$X), 3)
  expect_equal(ncol(dm$X), 3)
  expect_equal(colnames(dm$X), elements)
  # CH4 : C=1, H=4, O=0
  expect_equal(dm$X[1, "C"], 1)
  expect_equal(dm$X[1, "H"], 4)
  expect_equal(dm$X[1, "O"], 0)
})

test_that("fit_enthalpy_model fits a linear model without intercept", {
  df <- data.frame(
    formula = c("CH4", "C2H6", "C3H8", "C4H10"),
    delta_fH = c(-74.5, -84.0, -104.0, -125.0),
    phase = rep("g", 4),
    stringsAsFactors = FALSE
  )
  m <- fit_enthalpy_model(df, elements = c("C", "H"), intercept = FALSE)
  expect_s3_class(m, "lm")
  expect_equal(length(coef(m)), 2)
  expect_equal(attr(m, "elements"), c("C", "H"))
})

test_that("fit_enthalpy_model errors with too few compounds", {
  df <- data.frame(
    formula = "CH4",
    delta_fH = -74.5,
    phase = "g",
    stringsAsFactors = FALSE
  )
  expect_error(
    fit_enthalpy_model(df, elements = c("C", "H", "O")),
    "Pas assez de composés"
  )
})

test_that("predict_enthalpy returns a numeric value", {
  df <- data.frame(
    formula = c("CH4", "C2H6", "C3H8", "C4H10"),
    delta_fH = c(-74.5, -84.0, -104.0, -125.0),
    phase = rep("g", 4),
    stringsAsFactors = FALSE
  )
  m <- fit_enthalpy_model(df, elements = c("C", "H"))
  p <- predict_enthalpy(m, "CH4")
  expect_type(p, "double")
  expect_false(is.na(p))
  # La prédiction pour CH4 (dans le train) doit être proche de -74.5
  expect_lt(abs(p - (-74.5)), 5)
})

test_that("fit_all_models returns 3 models named ATcT, RMG, CHNOSZ", {
  # Utilise les vraies bases (chargées depuis data/)
  skip_if_not(file.exists("data/ATcT.rda"), "Data files not available")
  load("data/ATcT.rda"); load("data/RMG.rda"); load("data/CHNOSZ.rda")
  
  models <- fit_all_models(atct_harmonized, rmg_harmonized, chnosz_harmonized)
  expect_type(models, "list")
  expect_named(models, c("ATcT", "RMG", "CHNOSZ"))
  for (m in models) expect_s3_class(m, "lm")
})

