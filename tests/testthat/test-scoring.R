# Tests for Scoring Rules and Evaluation using scoringutils

# Conversion Tests ------------------------------------------------------

test_that("as_forecast_point S3 method converts forecast correctly", {
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

  # Use scoringutils S3 method
  fc_point <- as_forecast_point(fc)

  expect_s3_class(fc_point, "forecast_point")
  expect_true("observed" %in% names(fc_point))
  expect_true("predicted" %in% names(fc_point))
  expect_true("horizon" %in% names(fc_point))
  expect_true("model" %in% names(fc_point))
  expect_equal(nrow(fc_point), 5)
  expect_equal(fc_point$observed, truth_vals)
})

test_that("as_forecast_point requires truth values", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(50, mean = 100, sd = 10)

  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  # Forecast without truth
  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:5
  )

  expect_error(
    as_forecast_point(fc),
    "truth values"
  )
})

# Score Function Tests --------------------------------------------------

test_that("score works with converted forecast", {
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

  # Convert and score using scoringutils
  fc_point <- as_forecast_point(fc)
  scores <- score(fc_point)
  scores_summary <- summarise_scores(scores, by = "model")

  expect_s3_class(scores_summary, "data.frame")
  expect_true("model" %in% names(scores_summary))
  # Point forecast should have ae_point, se_point, ape
  expect_true(any(c("ae_point", "se_point", "ape") %in% names(scores_summary)))
})

test_that("score returns individual scores without summarising", {
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

  # Convert and score without summarising
  fc_point <- as_forecast_point(fc)
  scores <- score(fc_point)

  expect_s3_class(scores, "data.frame")
  expect_equal(nrow(scores), 5) # One row per horizon
  expect_true("horizon" %in% names(scores))
})

# Metric Access Tests --------------------------------------------------

test_that("score returns correct metrics that can be accessed", {
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

  # Convert and score
  fc_point <- as_forecast_point(fc)
  scores <- score(fc_point)
  scores_summary <- summarise_scores(scores, by = "model")

  # Check that key metrics are accessible
  expect_true("ae_point" %in% names(scores))
  expect_true("se_point" %in% names(scores))
  expect_true("ape" %in% names(scores))

  # Check values are numeric and non-negative
  expect_type(scores$ae_point, "double")
  expect_true(all(scores$ae_point >= 0))
  expect_true(all(scores$se_point >= 0))
  expect_true(all(scores$ape >= 0))

  # Check RMSE calculation
  rmse <- sqrt(scores$se_point)
  expect_type(rmse, "double")
  expect_true(all(rmse >= 0))
})

# Comparative Tests -----------------------------------------------------

test_that("MAE and RMSE have expected relationship", {
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

  # Convert and score
  fc_point <- as_forecast_point(fc)
  scores <- score(fc_point)
  scores_summary <- summarise_scores(scores, by = "model")

  mae_score <- scores_summary$ae_point
  mse_score <- scores_summary$se_point
  rmse_score <- sqrt(scores_summary$se_point)

  # RMSE should equal sqrt(MSE)
  expect_equal(rmse_score, sqrt(mse_score), tolerance = 1e-10)

  # By Cauchy-Schwarz, RMSE >= MAE (with equality only when all errors are equal)
  expect_true(rmse_score >= mae_score - 1e-10)
})

# Multiple Models Comparison --------------------------------------------

test_that("score works with multiple horizons and produces correct output", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(456)
  data <- rnorm(50, mean = 50, sd = 5)
  truth_vals <- rnorm(10, mean = 50, sd = 5)

  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted,
    interval_method = NoInterval(),
    horizon = 1:10,
    truth = truth_vals,
    model_name = "ConstantModel"
  )

  # Convert and score
  fc_point <- as_forecast_point(fc)
  scores <- score(fc_point)

  # Check individual scores
  expect_equal(nrow(scores), 10) # One row per horizon

  # Get summarised scores
  scores_sum <- summarise_scores(scores, by = "model")
  expect_equal(nrow(scores_sum), 1) # One row for the model
})
