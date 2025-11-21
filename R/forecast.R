# Core Forecasting Functions

#' Fit a Baseline Model
#'
#' Fits a baseline forecasting model to observed data.
#'
#' @param x Numeric vector of time series data
#' @param model A model object created by one of the model constructors
#' @param temporal_info Optional TemporalInfo object with start date and
#'   resolution
#'
#' @return A fitted model object that can be used for forecasting
#' @export
#'
#' @examples
#' \dontrun{
#' data <- c(1.2, 2.3, 3.1, 2.8, 3.5, 4.2, 3.9)
#' model <- ARMAModel(p = 1, q = 1)
#' fitted <- fit_baseline(data, model)
#' }
fit_baseline <- function(x, model, temporal_info = NULL) {
  check_setup()

  # Assign data to Julia
  JuliaCall::julia_assign("x_data", as.numeric(x))

  # Assign model to Julia
  JuliaCall::julia_assign("model_obj", model)

  # Fit the model
  if (is.null(temporal_info)) {
    JuliaCall::julia_eval(
      "fitted = ForecastBaselines.fit_baseline(x_data, model_obj)"
    )
  } else {
    JuliaCall::julia_assign("temp_info", temporal_info)
    JuliaCall::julia_eval(
      paste0(
        "fitted = ForecastBaselines.fit_baseline(",
        "x_data, model_obj, temporal_info=temp_info)"
      )
    )
  }

  # Return the fitted model
  JuliaCall::julia_eval("fitted")
}

#' Generate Point Forecasts
#'
#' Generates point forecasts from a fitted model.
#'
#' @param fitted A fitted model object from fit_baseline()
#' @param horizon Integer or vector of integers specifying forecast horizons
#'
#' @return Numeric vector of point forecasts
#' @export
#'
#' @examples
#' \dontrun{
#' forecasts <- point_forecast(fitted, horizon = 1:12)
#' }
point_forecast <- function(fitted, horizon = 1L) {
  check_setup()

  JuliaCall::julia_assign("fitted_obj", fitted)

  if (length(horizon) == 1) {
    JuliaCall::julia_assign("h", as.integer(horizon))
  } else {
    JuliaCall::julia_assign("h", as.integer(horizon))
  }

  result <- JuliaCall::julia_eval(
    "ForecastBaselines.point_forecast(fitted_obj, h)"
  )
  as.numeric(result)
}

#' Generate Complete Forecast with Intervals
#'
#' Generates a complete forecast including point forecasts, prediction
#' intervals, and optionally trajectories.
#'
#' @param fitted A fitted model object from fit_baseline()
#' @param interval_method Interval method object (NoInterval,
#'   EmpiricalInterval, etc.)
#' @param horizon Integer or vector of integers specifying forecast horizons
#' @param levels Numeric vector of confidence levels (default: 0.95)
#' @param include_median Whether to include median forecast (default: TRUE)
#' @param truth Optional vector of true values for evaluation
#' @param model_name Optional name for the model
#'
#' @return A Forecast object (list) containing forecasts and metadata
#' @export
#'
#' @examples
#' \dontrun{
#' # Point forecast only
#' fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:12)
#'
#' # With prediction intervals
#' fc <- forecast(fitted,
#'   interval_method = EmpiricalInterval(n_trajectories = 1000),
#'   horizon = 1:12,
#'   levels = c(0.80, 0.95)
#' )
#'
#' # With truth for evaluation
#' fc <- forecast(fitted,
#'   interval_method = EmpiricalInterval(),
#'   horizon = 1:12,
#'   truth = c(3.6, 3.8, 4.1, ...)
#' )
#' }
forecast <- function(fitted,
                     interval_method = NoInterval(),
                     horizon = 1L,
                     levels = 0.95,
                     include_median = TRUE,
                     truth = NULL,
                     model_name = "") {
  check_setup()

  # Assign objects to Julia
  JuliaCall::julia_assign("fitted_obj", fitted)
  JuliaCall::julia_assign("interval_method_obj", interval_method)
  JuliaCall::julia_assign("h", as.integer(horizon))
  # Ensure levels is always a vector in Julia
  if (length(levels) == 1) {
    # For single values, create Julia vector directly
    JuliaCall::julia_command(sprintf("lvls = [%f]", as.numeric(levels)))
  } else {
    # For multiple values, assign and convert
    JuliaCall::julia_assign("lvls_temp", as.numeric(levels))
    JuliaCall::julia_command("lvls = vec(lvls_temp)")
  }
  JuliaCall::julia_assign("inc_median", include_median)
  JuliaCall::julia_assign("mdl_name", as.character(model_name))

  # Build the forecast call and convert using Julia helper
  # Note: Julia code must be on a single line to avoid parse errors
  # Generate a unique ID for this forecast object
  fc_id <- paste0("fc_", sample.int(.Machine$integer.max, 1))

  if (is.null(truth)) {
    julia_cmd <- paste0(
      fc_id, " = ForecastBaselines.forecast(fitted_obj; ",
      "interval_method = interval_method_obj, horizon = h, levels = lvls, ",
      "include_median = inc_median, model_name = mdl_name)"
    )
    JuliaCall::julia_command(julia_cmd)
    forecast_result <- JuliaCall::julia_eval(
      sprintf("forecast_to_r_dict(%s)", fc_id)
    )
  } else {
    JuliaCall::julia_assign("truth_vals", as.numeric(truth))
    julia_cmd <- paste0(
      fc_id, " = ForecastBaselines.forecast(fitted_obj; ",
      "interval_method = interval_method_obj, horizon = h, levels = lvls, ",
      "include_median = inc_median, truth = truth_vals, model_name = mdl_name)"
    )
    JuliaCall::julia_command(julia_cmd)
    forecast_result <- JuliaCall::julia_eval(
      sprintf("forecast_to_r_dict(%s)", fc_id)
    )
  }

  # Store the Julia variable name for later use (e.g., in score())
  attr(forecast_result, "julia_ref") <- fc_id

  # Add S3 class
  class(forecast_result) <- c("ForecastBaselines_Forecast", "list")
  forecast_result
}

#' Generate Interval Forecasts
#'
#' Generates prediction intervals from a fitted model.
#'
#' @param fitted A fitted model object from fit_baseline()
#' @param method Interval method object
#' @param horizon Integer or vector of integers specifying forecast horizons
#' @param levels Numeric vector of confidence levels
#'
#' @return List containing point forecasts, median, intervals, and trajectories
#' @export
#'
#' @examples
#' \dontrun{
#' intervals <- interval_forecast(fitted,
#'   method = EmpiricalInterval(),
#'   horizon = 1:12,
#'   levels = c(0.80, 0.95)
#' )
#' }
interval_forecast <- function(fitted, method, horizon = 1L, levels = 0.95) {
  check_setup()

  JuliaCall::julia_assign("fitted_obj", fitted)
  JuliaCall::julia_assign("method_obj", method)
  JuliaCall::julia_assign("h", as.integer(horizon))
  # Ensure levels is always a vector in Julia
  if (length(levels) == 1) {
    # For single values, create Julia vector directly
    JuliaCall::julia_command(sprintf("lvls = [%f]", as.numeric(levels)))
  } else {
    # For multiple values, assign and convert
    JuliaCall::julia_assign("lvls_temp", as.numeric(levels))
    JuliaCall::julia_command("lvls = vec(lvls_temp)")
  }

  result <- JuliaCall::julia_eval("
    res = ForecastBaselines.interval_forecast(fitted_obj, method_obj, h, lvls)
    interval_result_to_r_dict(res)
  ")

  # Already converted by Julia helper
  result
}

#' Create Temporal Information Object
#'
#' Creates a TemporalInfo object to track dates and time resolution.
#'
#' @param start Start date (Date, POSIXct, or integer)
#' @param resolution Time resolution (integer for generic, or period object)
#'
#' @return A TemporalInfo object
#' @export
#'
#' @examples
#' \dontrun{
#' # Simple integer indexing
#' temp_info <- TemporalInfo(start = 1, resolution = 1)
#'
#' # Date-based
#' temp_info <- TemporalInfo(start = as.Date("2024-01-01"), resolution = 1)
#' }
TemporalInfo <- function(start = 1, resolution = 1) {
  check_setup()

  # Convert R date to Julia if needed
  if (inherits(start, "Date")) {
    start_str <- format(start, "%Y-%m-%d")
    JuliaCall::julia_assign("start_val", start_str)
    # Import Dates module and create Date
    JuliaCall::julia_eval("using Dates")
    JuliaCall::julia_eval("start_date = Date(start_val)")
    # Convert integer resolution to Day period
    JuliaCall::julia_assign("res_val", as.integer(resolution))
    JuliaCall::julia_eval(
      "ForecastBaselines.TemporalInfo(start_date, Day(res_val))"
    )
  } else {
    JuliaCall::julia_call(
      "ForecastBaselines.TemporalInfo",
      as.integer(start),
      as.integer(resolution)
    )
  }
}

# Internal helper function - no longer needed, conversion done inline
# Kept for backwards compatibility
convert_forecast_to_r <- function(jl_forecast) {
  # If already converted (from new inline approach), just return
  if (is.list(jl_forecast)) {
    class(jl_forecast) <- c("ForecastBaselines_Forecast", "list")
    return(jl_forecast)
  }

  # Fallback: shouldn't reach here with new approach
  warning("Using deprecated conversion path")
  class(jl_forecast) <- c("ForecastBaselines_Forecast", "list")
  jl_forecast
}

#' Print method for Forecast objects
#'
#' @param x A Forecast object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.ForecastBaselines_Forecast <- function(x, ...) {
  cat("ForecastBaselines Forecast Object\n")
  cat("==================================\n\n")

  if (!is.null(x$model_name) && nchar(x$model_name) > 0) {
    cat("Model:", x$model_name, "\n")
  }

  if (!is.null(x$horizon)) {
    cat("Horizon:", min(x$horizon), "to", max(x$horizon), "\n")
  }

  cat("\nComponents:\n")
  if (!is.null(x$mean)) cat("  - Mean forecasts\n")
  if (!is.null(x$median)) cat("  - Median forecasts\n")
  if (!is.null(x$intervals)) cat("  - Prediction intervals\n")
  if (!is.null(x$trajectories)) cat("  - Forecast trajectories\n")
  if (!is.null(x$truth)) cat("  - Truth values\n")

  invisible(x)
}
