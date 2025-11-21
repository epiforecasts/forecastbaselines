# Tests for Utility Functions

# Query Functions -------------------------------------------------------

test_that("has_horizon checks forecast horizon", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:5)

  result <- has_horizon(fc)
  expect_type(result, "logical")
  expect_length(result, 1)
  expect_true(result)
})

test_that("has_mean checks forecast mean", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:5)

  result <- has_mean(fc)
  expect_type(result, "logical")
  expect_length(result, 1)
  expect_true(result)
})

test_that("has_median checks forecast median", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:5,
    include_median = TRUE
  )

  result <- has_median(fc)
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("has_intervals checks forecast intervals", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:5)

  result <- has_intervals(fc)
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("has_truth checks forecast truth values", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  truth_vals <- rnorm(5, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:5,
    truth = truth_vals
  )

  result <- has_truth(fc)
  expect_type(result, "logical")
  expect_length(result, 1)
  expect_true(result)
})

test_that("has_trajectories checks forecast trajectories", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:5)

  result <- has_trajectories(fc)
  expect_type(result, "logical")
  expect_length(result, 1)
})

# Modification Functions ------------------------------------------------

test_that("add_truth adds truth values to forecast", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  truth_vals <- rnorm(5, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:5)

  fc_with_truth <- add_truth(fc, truth_vals)

  expect_s3_class(fc_with_truth, "ForecastBaselines_Forecast")
  expect_true(has_truth(fc_with_truth))
  expect_equal(fc_with_truth$truth, truth_vals, tolerance = 1e-10)
})

test_that("add_median adds median to forecast", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  median_vals <- rnorm(5, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:5,
    include_median = FALSE
  )

  fc_with_median <- add_median(fc, median_vals)

  expect_s3_class(fc_with_median, "ForecastBaselines_Forecast")
  expect_true(has_median(fc_with_median))
})

# Filter/Subset Functions -----------------------------------------------

test_that("truncate_horizon truncates forecast", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:10)

  fc_truncated <- truncate_horizon(fc, max_h = 5)

  expect_s3_class(fc_truncated, "ForecastBaselines_Forecast")
  expect_true(max(fc_truncated$horizon) <= 5)
})

test_that("filter_horizons filters specific horizons", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:10)

  fc_filtered <- filter_horizons(fc, horizons = c(1, 3, 5, 7))

  expect_s3_class(fc_filtered, "ForecastBaselines_Forecast")
  expect_length(fc_filtered$horizon, 4)
})

# Integration Tests -----------------------------------------------------

test_that("utility functions work together", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  truth_vals <- rnorm(10, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  # Create forecast
  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:10)

  # Check components
  expect_true(has_horizon(fc))
  expect_true(has_mean(fc))
  expect_false(has_truth(fc))

  # Add truth
  fc <- add_truth(fc, truth_vals)
  expect_true(has_truth(fc))

  # Truncate
  fc <- truncate_horizon(fc, max_h = 5)
  expect_true(max(fc$horizon) <= 5)

  # Filter
  fc <- filter_horizons(fc, horizons = c(1, 2, 3))
  expect_length(fc$horizon, 3)
})

test_that("has_* functions handle missing components correctly", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:5,
    include_median = FALSE
  )

  # Should have these
  expect_true(has_horizon(fc))
  expect_true(has_mean(fc))

  # Should not have these
  expect_false(has_truth(fc))
  expect_false(has_trajectories(fc))
})

test_that("truncate_horizon preserves forecast structure", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  truth_vals <- rnorm(10, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:10,
    truth = truth_vals,
    model_name = "TestModel"
  )

  fc_truncated <- truncate_horizon(fc, max_h = 5)

  # Check that structure is preserved
  expect_s3_class(fc_truncated, "ForecastBaselines_Forecast")
  expect_equal(fc_truncated$model_name, "TestModel")
  expect_true(has_truth(fc_truncated))
  expect_true(has_mean(fc_truncated))
})

test_that("filter_horizons preserves forecast properties", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(30, mean = 100, sd = 10)
  truth_vals <- rnorm(10, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:10,
    truth = truth_vals,
    model_name = "TestModel"
  )

  fc_filtered <- filter_horizons(fc, horizons = c(2, 4, 6, 8))

  # Check that properties are preserved
  expect_s3_class(fc_filtered, "ForecastBaselines_Forecast")
  expect_equal(fc_filtered$model_name, "TestModel")
  expect_true(has_truth(fc_filtered))
  expect_length(fc_filtered$mean, 4)
  expect_length(fc_filtered$truth, 4)
})
