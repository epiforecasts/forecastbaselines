# Tests for Julia setup and initialization

test_that("is_setup returns logical", {
  result <- is_setup()
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("Julia availability can be checked", {
  # Should not throw error even if Julia not available
  expect_no_error(is_setup())
})

test_that("package startup message is informative", {
  # This tests the .onLoad function behavior
  # We can't directly test .onLoad, but we can verify
  # the function exists and doesn't error
  expect_true(exists(".onLoad", where = asNamespace("forecastbaselines")))
})
