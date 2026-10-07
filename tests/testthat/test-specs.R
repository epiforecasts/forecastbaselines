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
  # A constructor passed without its brackets is refused in R too.
  expect_error(fit_baseline(1:10, ConstantModel), "must be a model spec")
  expect_error(
    forecast(NULL, interval_method = EmpiricalInterval),
    "must be an interval spec"
  )
})

test_that("models and interval methods already in Julia are accepted", {
  skip_if_no_julia()

  fitted <- fit_baseline(
    as.numeric(1:20), as_julia_object(ConstantModel())
  )
  fc <- forecast(fitted, interval_method = as_julia_object(NoInterval()))
  expect_s3_class(fc, "ForecastBaselines_Forecast")
  method <- as_julia_object(EmpiricalInterval(seed = 1))
  expect_length(interval_forecast(fitted, method)$point, 1)
})
