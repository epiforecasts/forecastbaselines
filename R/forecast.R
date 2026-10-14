# Core Forecasting Functions

#' Fit a Baseline Model
#'
#' Fits a baseline forecasting model to observed data.
#'
#' @param x Numeric vector of time series data
#' @param model A model created using ConstantModel(), ARMAModel(), ETSModel(),
#'   or other model functions (see ?ConstantModel for available models)
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
  check_role(model, "model")
  check_setup()
  args <- list(
    "ForecastBaselines.fit_baseline", julia_vector(as.numeric(x)),
    as_julia_object(model)
  )
  if (!is.null(temporal_info)) {
    args$temporal_info <- as_julia_object(temporal_info)
  }
  do.call(juliabridge::call_julia, args)
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
  as.numeric(juliabridge::call_julia(
    "ForecastBaselines.point_forecast", fitted, as.integer(horizon)
  ))
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
  check_role(interval_method, "interval")
  check_setup()
  args <- list(
    "ForecastBaselines.forecast", fitted,
    interval_method = as_julia_object(interval_method),
    horizon = as.integer(horizon),
    levels = julia_vector(as.numeric(levels)),
    include_median = include_median,
    model_name = as.character(model_name)
  )
  if (!is.null(truth)) args$truth <- julia_vector(as.numeric(truth))
  new_forecast(do.call(juliabridge::call_julia, args))
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
  check_role(method, "interval")
  check_setup()
  result <- juliabridge::call_julia(
    "ForecastBaselines.interval_forecast",
    fitted, as_julia_object(method),
    as.integer(horizon), julia_vector(as.numeric(levels))
  )
  juliabridge::get_julia(
    juliabridge::call_julia("interval_result_to_r", result)
  )
}

#' Create Temporal Information Object
#'
#' Creates a TemporalInfo object to track dates and time resolution.
#'
#' @param start Start date (Date, POSIXct, or integer)
#' @param resolution Time resolution (integer for generic, or period object)
#'
#' @return A spec for a Julia `TemporalInfo`, built when it is used
#'   (see [juliabridge::julia_spec()])
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
  if (inherits(start, "Date")) {
    start <- juliabridge::julia_spec(
      "ForecastBaselines.Dates.Date", format(start, "%Y-%m-%d")
    )
    resolution <- juliabridge::julia_spec(
      "ForecastBaselines.Dates.Day", as.integer(resolution)
    )
  } else {
    start <- as.integer(start)
    resolution <- as.integer(resolution)
  }
  juliabridge::julia_spec(
    "ForecastBaselines.TemporalInfo", start, resolution,
    .role = "temporal_info"
  )
}

# Internal: translate a Julia Forecast into the R list users work with,
# keeping the Julia object alongside for the functions that act on it.
new_forecast <- function(jl_forecast) {
  result <- juliabridge::get_julia(
    juliabridge::call_julia("forecast_to_r", jl_forecast)
  )
  class(result) <- c("ForecastBaselines_Forecast", "list")
  attr(result, "julia_ref") <- jl_forecast
  result
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
