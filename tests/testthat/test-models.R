# Tests for Model Functions

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

test_that("KDEModel takes no arguments", {
  expect_error(KDEModel(bandwidth = 0.5), "unused argument")
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

test_that("LSDModel rejects unsupported trend_correction", {
  expect_error(LSDModel(s = 7, trend_correction = TRUE), "unused argument")
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

test_that("IDSModel creates valid model with custom window size", {
  skip_if_no_julia()

  model <- IDSModel(window_size = 5)
  expect_true(!is.null(model))
})

test_that("IDSModel rejects unsupported threshold", {
  expect_error(IDSModel(threshold = 0.1), "unused argument")
})

# STLModel Tests --------------------------------------------------------

test_that("STLModel creates valid model with seasonal period", {
  skip_if_no_julia()

  model <- STLModel(s = 12)
  expect_true(!is.null(model))
})

test_that("STLModel rejects unsupported trend and robust options", {
  expect_error(STLModel(s = 12, trend = TRUE), "unused argument")
  expect_error(STLModel(s = 12, robust = TRUE), "unused argument")
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

test_that("ARMAModel creates model with drift", {
  skip_if_no_julia()

  model <- ARMAModel(p = 1, include_drift = TRUE)
  expect_true(!is.null(model))
})

test_that("ARMAModel drift reaches the Julia model", {
  skip_if_no_julia()

  # A drift term adds one parameter to the mean function; were it dropped,
  # the two models would be identical.
  expect_equal(mean_dim(ARMAModel(p = 1)), 1L)
  expect_equal(mean_dim(ARMAModel(p = 1, include_drift = TRUE)), 2L)
})

test_that("ARMAModel harmonics reach the Julia model", {
  skip_if_no_julia()

  # Intercept plus 2k harmonic coefficients, plus one more for the drift.
  expect_equal(mean_dim(ARMAModel(p = 1, s = 52)), 3L)
  expect_equal(mean_dim(ARMAModel(p = 1, s = 52, k = 4)), 9L)
  expect_equal(
    mean_dim(ARMAModel(p = 1, s = 52, k = 4, include_drift = TRUE)), 10L
  )
})

test_that("ARMAModel rejects unsupported include_mean", {
  expect_error(ARMAModel(p = 1, include_mean = FALSE), "unused argument")
})

test_that("ARMAModel rejects harmonics without a seasonal period", {
  expect_error(ARMAModel(p = 1, k = 4), "require a seasonal period")
})

test_that("ARMAModel rejects a negative seasonal period", {
  # A negative period would reach Julia as a non-seasonal model, silently
  # dropping both the period and the harmonics.
  expect_error(ARMAModel(p = 1, s = -12, k = 3), "non-negative seasonal period")
})

test_that("ARMAModel rejects fewer than one harmonic wave", {
  # k = 0 would otherwise drop the seasonal period on the Julia side.
  expect_error(ARMAModel(p = 1, s = 12, k = 0), "at least 1 harmonic wave")
})

test_that("ARMAModel rejects a non-logical include_drift", {
  expect_error(ARMAModel(p = 1, include_drift = "yes"), "TRUE or FALSE")
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

test_that("INARCHModel warns about harmonics beyond the first", {
  skip_if_no_julia()

  expect_warning(
    INARCHModel(p = 1, s = 52, k = 2),
    "not carried into the forecasts"
  )
})

test_that("INARCHModel passes seasonality and negative binomial through", {
  skip_if_no_julia()

  model <- suppressWarnings(INARCHModel(p = 1, s = 52, k = 4, nb = TRUE))
  expect_equal(as.integer(model$s), 52L)
  expect_equal(as.integer(model$k), 4L)
  expect_true(as.logical(model$nb))
})

test_that("INARCHModel rejects harmonics without a seasonal period", {
  expect_error(INARCHModel(p = 1, k = 4), "require a seasonal period")
})

test_that("INARCHModel rejects fewer than one harmonic wave", {
  expect_error(INARCHModel(p = 1, s = 52, k = 0), "at least 1 harmonic wave")
})

test_that("INARCHModel rejects a non-logical nb", {
  expect_error(INARCHModel(p = 1, nb = 1), "TRUE or FALSE")
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

  model <- ETSModel(
    error_type = "A", trend_type = "A", season_type = "A", s = 12
  )
  expect_true(!is.null(model))
})

test_that("ETSModel creates Holt-Winters multiplicative (M,M,M)", {
  skip_if_no_julia()

  model <- ETSModel(
    error_type = "M", trend_type = "M", season_type = "M", s = 12
  )
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

test_that("ETSModel rejects unsupported damped argument", {
  expect_error(
    ETSModel(error_type = "A", trend_type = "A", damped = TRUE),
    "unused argument"
  )
})

test_that("ETSModel requires s for seasonal models", {
  skip_if_no_julia()

  expect_error(
    ETSModel(error_type = "A", trend_type = "N", season_type = "A"),
    "Seasonal period 's' must be provided"
  )
})
