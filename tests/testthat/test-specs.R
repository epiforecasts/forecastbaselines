# Tests for model, interval method and transformation specs

test_that("constructors describe the Julia call without starting Julia", {
  expect_identical(
    juliabridge::as_julia(LSDModel(s = 12, window_width = 2)),
    "ForecastBaselines.LSDModel(; s = 12, w = 2)"
  )
  expect_identical(
    juliabridge::as_julia(ARMAModel(p = 1, s = 52, k = 2)),
    "arma_model(; p = 1, q = 0, s = 52, k = 2, trend = false)"
  )
})

test_that("arguments reach Julia under their Julia names", {
  expect_identical(LSDModel(s = 12, window_width = 2)$kwargs$w, 2L)
  expect_identical(OLSModel(degree = 2)$kwargs$p, 2L)
  expect_true(ARMAModel(p = 1, include_drift = TRUE)$kwargs$trend)
  expect_identical(ETSModel(error_type = "M")$kwargs$error, "N")
})

test_that("roles keep models and interval methods apart", {
  expect_error(
    forecast(NULL, interval_method = ConstantModel()),
    "must be an interval spec"
  )
})
