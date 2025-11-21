# Tests for Model Constructors

# ConstantModel Tests ---------------------------------------------------

test_that("ConstantModel creates valid model", {
  skip_if_no_julia()

  model <- ConstantModel()
  expect_true(!is.null(model))
})

# MarginalModel Tests ---------------------------------------------------

test_that("MarginalModel creates valid model with default p", {
  skip_if_no_julia()

  model <- MarginalModel()
  expect_true(!is.null(model))
})

test_that("MarginalModel creates valid model with specified p", {
  skip_if_no_julia()

  model <- MarginalModel(p = 10)
  expect_true(!is.null(model))
})

# KDEModel Tests --------------------------------------------------------

test_that("KDEModel creates valid model with defaults", {
  skip_if_no_julia()

  model <- KDEModel()
  expect_true(!is.null(model))
})

test_that("KDEModel creates valid model with bandwidth", {
  skip_if_no_julia()

  model <- KDEModel(bandwidth = 0.5)
  expect_true(!is.null(model))
})

# LSDModel Tests --------------------------------------------------------

test_that("LSDModel creates valid model with seasonal period", {
  skip_if_no_julia()

  model <- LSDModel(s = 7)
  expect_true(!is.null(model))
})

test_that("LSDModel creates valid model with window_width", {
  skip_if_no_julia()

  model <- LSDModel(s = 12, window_width = 2)
  expect_true(!is.null(model))
})

test_that("LSDModel creates valid model with trend_correction", {
  skip_if_no_julia()

  model <- LSDModel(s = 7, trend_correction = TRUE)
  expect_true(!is.null(model))
})

# OLSModel Tests --------------------------------------------------------

test_that("OLSModel creates valid model with linear trend", {
  skip_if_no_julia()

  model <- OLSModel(degree = 1)
  expect_true(!is.null(model))
})

test_that("OLSModel creates valid model with quadratic trend", {
  skip_if_no_julia()

  model <- OLSModel(degree = 2)
  expect_true(!is.null(model))
})

test_that("OLSModel creates valid model with differencing", {
  skip_if_no_julia()

  # Need degree >= 2 for differencing = 1 (sufficient temporal lag)
  model <- OLSModel(degree = 2, differencing = 1)
  expect_true(!is.null(model))
})

# IDSModel Tests --------------------------------------------------------

test_that("IDSModel creates valid model with defaults", {
  skip_if_no_julia()

  model <- IDSModel()
  expect_true(!is.null(model))
})

test_that("IDSModel creates valid model with custom parameters", {
  skip_if_no_julia()

  model <- IDSModel(threshold = 0.1, window_size = 5)
  expect_true(!is.null(model))
})

# STLModel Tests --------------------------------------------------------

test_that("STLModel creates valid model with seasonal period", {
  skip_if_no_julia()

  model <- STLModel(s = 12)
  expect_true(!is.null(model))
})

test_that("STLModel creates valid model with trend", {
  skip_if_no_julia()

  model <- STLModel(s = 12, trend = TRUE)
  expect_true(!is.null(model))
})

test_that("STLModel creates valid model with robust option", {
  skip_if_no_julia()

  model <- STLModel(s = 12, robust = TRUE)
  expect_true(!is.null(model))
})

# ARMAModel Tests -------------------------------------------------------

test_that("ARMAModel creates AR(1) model", {
  skip_if_no_julia()

  model <- ARMAModel(p = 1)
  expect_true(!is.null(model))
})

test_that("ARMAModel creates MA(1) model", {
  skip_if_no_julia()

  model <- ARMAModel(q = 1)
  expect_true(!is.null(model))
})

test_that("ARMAModel creates ARMA(2,1) model", {
  skip_if_no_julia()

  model <- ARMAModel(p = 2, q = 1)
  expect_true(!is.null(model))
})

test_that("ARMAModel creates model with seasonality", {
  skip_if_no_julia()

  model <- ARMAModel(p = 1, q = 1, s = 12)
  expect_true(!is.null(model))
})

test_that("ARMAModel creates model with mean and drift", {
  skip_if_no_julia()

  model <- ARMAModel(p = 1, include_mean = TRUE, include_drift = TRUE)
  expect_true(!is.null(model))
})

# INARCHModel Tests -----------------------------------------------------

test_that("INARCHModel creates valid model", {
  skip_if_no_julia()

  model <- INARCHModel(p = 1)
  expect_true(!is.null(model))
})

test_that("INARCHModel creates valid model with higher order", {
  skip_if_no_julia()

  model <- INARCHModel(p = 2)
  expect_true(!is.null(model))
})

# ETSModel Tests --------------------------------------------------------

test_that("ETSModel creates simple exponential smoothing (A,N,N)", {
  skip_if_no_julia()

  model <- ETSModel(error_type = "A", trend_type = "N", season_type = "N")
  expect_true(!is.null(model))
})

test_that("ETSModel creates Holt's linear trend (A,A,N)", {
  skip_if_no_julia()

  model <- ETSModel(error_type = "A", trend_type = "A", season_type = "N")
  expect_true(!is.null(model))
})

test_that("ETSModel creates Holt-Winters additive (A,A,A)", {
  skip_if_no_julia()

  model <- ETSModel(error_type = "A", trend_type = "A", season_type = "A", s = 12)
  expect_true(!is.null(model))
})

test_that("ETSModel creates Holt-Winters multiplicative (M,M,M)", {
  skip_if_no_julia()

  model <- ETSModel(error_type = "M", trend_type = "M", season_type = "M", s = 12)
  expect_true(!is.null(model))
})

test_that("ETSModel validates error_type", {
  skip_if_no_julia()

  expect_error(
    ETSModel(error_type = "X", trend_type = "N", season_type = "N"),
    "error_type must be one of"
  )
})

test_that("ETSModel validates trend_type", {
  skip_if_no_julia()

  expect_error(
    ETSModel(error_type = "A", trend_type = "X", season_type = "N"),
    "trend_type must be one of"
  )
})

test_that("ETSModel validates season_type", {
  skip_if_no_julia()

  expect_error(
    ETSModel(error_type = "A", trend_type = "N", season_type = "X"),
    "season_type must be one of"
  )
})

test_that("ETSModel requires s for seasonal models", {
  skip_if_no_julia()

  expect_error(
    ETSModel(error_type = "A", trend_type = "N", season_type = "A"),
    "Seasonal period 's' must be provided"
  )
})
