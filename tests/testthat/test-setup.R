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

test_that("setup_ForecastBaselines exists and is exported", {
  # Setup is now manual / lazy via check_setup() rather than .onLoad,
  # to avoid crashing R during package attach when other compiled
  # backends (e.g. Stan) are present.
  expect_true(exists("setup_ForecastBaselines",
                    where = asNamespace("forecastbaselines"),
                    mode = "function"))
})
