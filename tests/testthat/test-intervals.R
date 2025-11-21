# Tests for Interval Methods

# NoInterval Tests ------------------------------------------------------

test_that("NoInterval creates valid interval method", {
  skip_if_no_julia()

  method <- NoInterval()
  expect_true(!is.null(method))
})

# EmpiricalInterval Tests -----------------------------------------------

test_that("EmpiricalInterval creates valid method with defaults", {
  skip_if_no_julia()

  method <- EmpiricalInterval()
  expect_true(!is.null(method))
})

test_that("EmpiricalInterval accepts custom n_trajectories", {
  skip_if_no_julia()

  method <- EmpiricalInterval(n_trajectories = 2000)
  expect_true(!is.null(method))
})

test_that("EmpiricalInterval accepts seed parameter", {
  skip_if_no_julia()

  method <- EmpiricalInterval(n_trajectories = 1000, seed = 123)
  expect_true(!is.null(method))
})

test_that("EmpiricalInterval accepts positivity_correction options", {
  skip_if_no_julia()

  method1 <- EmpiricalInterval(positivity_correction = "none")
  expect_true(!is.null(method1))

  method2 <- EmpiricalInterval(positivity_correction = "post_clip")
  expect_true(!is.null(method2))

  method3 <- EmpiricalInterval(positivity_correction = "truncate")
  expect_true(!is.null(method3))

  method4 <- EmpiricalInterval(positivity_correction = "zero_floor")
  expect_true(!is.null(method4))
})

test_that("EmpiricalInterval validates positivity_correction", {
  skip_if_no_julia()

  expect_error(
    EmpiricalInterval(positivity_correction = "invalid"),
    "positivity_correction must be one of"
  )
})

test_that("EmpiricalInterval accepts symmetry_correction", {
  skip_if_no_julia()

  method <- EmpiricalInterval(symmetry_correction = TRUE)
  expect_true(!is.null(method))
})

test_that("EmpiricalInterval accepts stepwise parameter", {
  skip_if_no_julia()

  method <- EmpiricalInterval(stepwise = TRUE)
  expect_true(!is.null(method))
})

test_that("EmpiricalInterval accepts return_trajectories", {
  skip_if_no_julia()

  method <- EmpiricalInterval(return_trajectories = TRUE)
  expect_true(!is.null(method))
})

test_that("EmpiricalInterval works with multiple options", {
  skip_if_no_julia()

  method <- EmpiricalInterval(
    n_trajectories = 500,
    seed = 456,
    positivity_correction = "post_clip",
    symmetry_correction = TRUE,
    stepwise = FALSE,
    return_trajectories = TRUE
  )
  expect_true(!is.null(method))
})

# ParametricInterval Tests ----------------------------------------------

test_that("ParametricInterval creates valid method with defaults", {
  skip_if_no_julia()

  method <- ParametricInterval()
  expect_true(!is.null(method))
})

test_that("ParametricInterval accepts positivity_correction", {
  skip_if_no_julia()

  method1 <- ParametricInterval(positivity_correction = "none")
  expect_true(!is.null(method1))

  method2 <- ParametricInterval(positivity_correction = "post_clip")
  expect_true(!is.null(method2))
})

test_that("ParametricInterval validates positivity_correction", {
  skip_if_no_julia()

  expect_error(
    ParametricInterval(positivity_correction = "invalid"),
    "positivity_correction must be one of"
  )
})

# ModelTrajectoryInterval Tests -----------------------------------------

test_that("ModelTrajectoryInterval creates valid method with defaults", {
  skip_if_no_julia()

  method <- ModelTrajectoryInterval()
  expect_true(!is.null(method))
})

test_that("ModelTrajectoryInterval accepts custom n_trajectories", {
  skip_if_no_julia()

  method <- ModelTrajectoryInterval(n_trajectories = 2000)
  expect_true(!is.null(method))
})

test_that("ModelTrajectoryInterval accepts seed parameter", {
  skip_if_no_julia()

  method <- ModelTrajectoryInterval(n_trajectories = 1000, seed = 789)
  expect_true(!is.null(method))
})

test_that("ModelTrajectoryInterval accepts positivity_correction options", {
  skip_if_no_julia()

  method1 <- ModelTrajectoryInterval(positivity_correction = "none")
  expect_true(!is.null(method1))

  method2 <- ModelTrajectoryInterval(positivity_correction = "post_clip")
  expect_true(!is.null(method2))

  method3 <- ModelTrajectoryInterval(positivity_correction = "truncate")
  expect_true(!is.null(method3))

  method4 <- ModelTrajectoryInterval(positivity_correction = "zero_floor")
  expect_true(!is.null(method4))
})

test_that("ModelTrajectoryInterval validates positivity_correction", {
  skip_if_no_julia()

  expect_error(
    ModelTrajectoryInterval(positivity_correction = "invalid"),
    "positivity_correction must be one of"
  )
})

test_that("ModelTrajectoryInterval accepts return_trajectories", {
  skip_if_no_julia()

  method <- ModelTrajectoryInterval(return_trajectories = TRUE)
  expect_true(!is.null(method))
})

test_that("ModelTrajectoryInterval works with multiple options", {
  skip_if_no_julia()

  method <- ModelTrajectoryInterval(
    n_trajectories = 500,
    seed = 321,
    positivity_correction = "truncate",
    return_trajectories = TRUE
  )
  expect_true(!is.null(method))
})

# Integration Tests -----------------------------------------------------

test_that("NoInterval produces forecast without intervals", {
  skip_if_no_julia()
  skip_on_cran()

  set.seed(123)
  data <- rnorm(50, mean = 100, sd = 10)
  model <- ConstantModel()
  fitted <- fit_baseline(data, model)

  fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:5)

  expect_s3_class(fc, "ForecastBaselines_Forecast")
  expect_true(!is.null(fc$mean))
})
