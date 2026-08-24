# Test setup and helpers
# This file is automatically loaded before all tests

# Try to setup Julia and ForecastBaselines automatically for tests
tryCatch(
  {
    if (!is_setup()) {
      setup_ForecastBaselines(verbose = FALSE)
    }

    # Load helper file for tests (since package may not be installed)
    # Try multiple possible locations
    possible_paths <- c(
      file.path("inst", "julia", "forecast_helpers.jl"),
      file.path("..", "..", "inst", "julia", "forecast_helpers.jl"),
      system.file("julia", "forecast_helpers.jl", package = "forecastbaselines")
    )

    helper_loaded <- FALSE
    for (helper_file in possible_paths) {
      if (file.exists(helper_file) && !helper_loaded) {
        juliaready::command_julia(sprintf('include("%s")', normalizePath(helper_file)))
        helper_loaded <- TRUE
        break
      }
    }

    if (!helper_loaded) {
      warning("Could not load forecast_helpers.jl")
    }
  },
  error = function(e) {
    message("Could not automatically setup ForecastBaselines for tests: ", e$message)
  }
)

# Skip helper for Julia/ForecastBaselines availability
skip_if_no_julia <- function() {
  testthat::skip_if_not(is_setup(), "Julia/ForecastBaselines not available")
}

# Number of mean-function parameters of a Julia ARMAModel. The field is named
# "\u03bcDim" (muDim) in Julia; it is written as an escape to keep this file
# ASCII.
mean_dim <- function(model) {
  as.integer(model[["\u03bcDim"]])
}
