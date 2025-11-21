# Tests for Data Transformations

# NoTransform Tests -----------------------------------------------------

test_that("NoTransform creates valid transformation", {
  skip_if_no_julia()

  trans <- NoTransform()
  expect_true(!is.null(trans))
})

# LogTransform Tests ----------------------------------------------------

test_that("LogTransform creates valid transformation", {
  skip_if_no_julia()

  trans <- LogTransform()
  expect_true(!is.null(trans))
})

# LogPlusOneTransform Tests ---------------------------------------------

test_that("LogPlusOneTransform creates valid transformation with default c", {
  skip_if_no_julia()

  trans <- LogPlusOneTransform()
  expect_true(!is.null(trans))
})

test_that("LogPlusOneTransform creates valid transformation with custom c", {
  skip_if_no_julia()

  trans <- LogPlusOneTransform(c = 0.5)
  expect_true(!is.null(trans))
})

# SquareRootTransform Tests ---------------------------------------------

test_that("SquareRootTransform creates valid transformation", {
  skip_if_no_julia()

  trans <- SquareRootTransform()
  expect_true(!is.null(trans))
})

# PowerTransform Tests --------------------------------------------------

test_that("PowerTransform creates valid transformation", {
  skip_if_no_julia()

  trans <- PowerTransform(lambda = 0.5)
  expect_true(!is.null(trans))
})

test_that("PowerTransform works with different lambda values", {
  skip_if_no_julia()

  trans1 <- PowerTransform(lambda = 0.3)
  expect_true(!is.null(trans1))

  trans2 <- PowerTransform(lambda = 1 / 3)
  expect_true(!is.null(trans2))

  trans3 <- PowerTransform(lambda = 2)
  expect_true(!is.null(trans3))
})

# PowerPlusOneTransform Tests -------------------------------------------

test_that("PowerPlusOneTransform creates valid transformation", {
  skip_if_no_julia()

  trans <- PowerPlusOneTransform(lambda = 0.3)
  expect_true(!is.null(trans))
})

test_that("PowerPlusOneTransform accepts custom constant", {
  skip_if_no_julia()

  trans <- PowerPlusOneTransform(lambda = 0.5, constant = 0.5)
  expect_true(!is.null(trans))
})

# transform_data Tests --------------------------------------------------

test_that("transform_data applies NoTransform correctly", {
  skip_if_no_julia()

  data <- c(1, 2, 3, 4, 5)
  trans <- NoTransform()
  transformed <- transform_data(data, trans)

  expect_type(transformed, "double")
  expect_length(transformed, 5)
  expect_equal(transformed, data, tolerance = 1e-10)
})

test_that("transform_data applies LogTransform correctly", {
  skip_if_no_julia()

  data <- c(1, 2, 3, 4, 5)
  trans <- LogTransform()
  transformed <- transform_data(data, trans)

  expect_type(transformed, "double")
  expect_length(transformed, 5)
  expect_equal(transformed, log(data), tolerance = 1e-10)
})

test_that("transform_data applies LogPlusOneTransform correctly", {
  skip_if_no_julia()

  data <- c(0, 1, 2, 3, 4)
  trans <- LogPlusOneTransform(c = 1)
  transformed <- transform_data(data, trans)

  expect_type(transformed, "double")
  expect_length(transformed, 5)
  expect_equal(transformed, log(data + 1), tolerance = 1e-10)
})

test_that("transform_data applies PowerTransform correctly", {
  skip_if_no_julia()

  data <- c(1, 4, 9, 16, 25)
  trans <- PowerTransform(lambda = 0.5)
  transformed <- transform_data(data, trans)

  expect_type(transformed, "double")
  expect_length(transformed, 5)
  expect_equal(transformed, sqrt(data), tolerance = 1e-10)
})

# inverse_transform_data Tests ------------------------------------------

test_that("inverse_transform_data inverts NoTransform correctly", {
  skip_if_no_julia()

  data <- c(1, 2, 3, 4, 5)
  trans <- NoTransform()
  transformed <- transform_data(data, trans)
  original <- inverse_transform_data(transformed, trans)

  expect_equal(original, data, tolerance = 1e-10)
})

test_that("inverse_transform_data inverts LogTransform correctly", {
  skip_if_no_julia()

  data <- c(1, 2, 3, 4, 5)
  trans <- LogTransform()
  transformed <- transform_data(data, trans)
  original <- inverse_transform_data(transformed, trans)

  expect_equal(original, data, tolerance = 1e-10)
})

test_that("inverse_transform_data inverts LogPlusOneTransform correctly", {
  skip_if_no_julia()

  data <- c(0, 1, 2, 3, 4)
  trans <- LogPlusOneTransform(c = 1)
  transformed <- transform_data(data, trans)
  original <- inverse_transform_data(transformed, trans)

  expect_equal(original, data, tolerance = 1e-10)
})

test_that("inverse_transform_data inverts PowerTransform correctly", {
  skip_if_no_julia()

  data <- c(1, 4, 9, 16, 25)
  trans <- PowerTransform(lambda = 0.5)
  transformed <- transform_data(data, trans)
  original <- inverse_transform_data(transformed, trans)

  expect_equal(original, data, tolerance = 1e-10)
})
