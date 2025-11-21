# Tests for Forecasting Functions

# fit_baseline Tests ----------------------------------------------------

test_that("fit_baseline works with ConstantModel", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(50, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  expect_true(!is.null(fitted))
})

test_that("fit_baseline works with ARMAModel", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(50, mean = 100, sd = 10)
  model <- ARMAModel(p = 1)
  fitted <- fit_baseline(data, model)

  expect_true(!is.null(fitted))
})

test_that("fit_baseline works with temporal_info", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(50, mean = 100, sd = 10)
  model <- ConstantModel()
  temp_info <- TemporalInfo(start = 1, resolution = 1)
  fitted <- fit_baseline(data, model, temporal_info = temp_info)

  expect_true(!is.null(fitted))
})

# point_forecast Tests --------------------------------------------------

test_that("point_forecast generates forecasts for single horizon", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(50, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  forecasts <- point_forecast(fitted, horizon = 1)

  expect_type(forecasts, "double")
  expect_length(forecasts, 1)
})

test_that("point_forecast generates forecasts for multiple horizons", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(50, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  forecasts <- point_forecast(fitted, horizon = 1:12)

  expect_type(forecasts, "double")
  expect_length(forecasts, 12)
})

# forecast Tests --------------------------------------------------------

test_that("forecast creates Forecast object", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(50, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:5)

  expect_s3_class(fc, "ForecastBaselines_Forecast")
  expect_true(!is.null(fc$horizon))
  expect_true(!is.null(fc$mean))
})

test_that("forecast works with single horizon", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1L)

  expect_s3_class(fc, "ForecastBaselines_Forecast")
  expect_length(fc$mean, 1)
})

test_that("forecast works with multiple horizons", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:10)

  expect_s3_class(fc, "ForecastBaselines_Forecast")
  expect_length(fc$mean, 10)
})

test_that("forecast works with truth values", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(50, mean = 100, sd = 10)
  truth_vals <- rnorm(5, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:5,
    truth = truth_vals
  )

  expect_s3_class(fc, "ForecastBaselines_Forecast")
  expect_true(!is.null(fc$truth))
  expect_length(fc$truth, 5)
})

test_that("forecast works with model_name", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:5,
    model_name = "TestModel"
  )

  expect_s3_class(fc, "ForecastBaselines_Forecast")
  expect_equal(fc$model_name, "TestModel")
})

# TemporalInfo Tests ----------------------------------------------------

test_that("TemporalInfo creates object with integer start", {
  skip_if_no_julia()

  temp_info <- TemporalInfo(start = 1, resolution = 1)
  expect_true(!is.null(temp_info))
})

test_that("TemporalInfo creates object with Date start", {
  skip_if_no_julia()

  start_date <- as.Date("2024-01-01")
  temp_info <- TemporalInfo(start = start_date, resolution = 1)
  expect_true(!is.null(temp_info))
})

# Print method Tests ----------------------------------------------------

test_that("print.ForecastBaselines_Forecast works", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:5,
    model_name = "TestModel"
  )

  expect_output(print(fc), "ForecastBaselines Forecast Object")
  expect_output(print(fc), "TestModel")
  expect_output(print(fc), "Horizon")
})
