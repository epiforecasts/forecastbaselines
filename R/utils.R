# Forecast Utility Functions

# Internal: the Julia Forecast behind a forecast from forecast(). It does not
# survive saving, nor a restart of Julia.
julia_forecast <- function(forecast) {
  ref <- attr(forecast, "julia_ref")
  if (!inherits(forecast, "ForecastBaselines_Forecast") || is.null(ref)) {
    stop(
      "`forecast` must be a forecast created by forecast() in this session.",
      call. = FALSE
    )
  }
  ref
}

# Internal: ask a question of the Julia Forecast.
query_forecast <- function(fn, forecast) {
  check_setup()
  as.logical(juliabridge::call_julia(
    paste0("ForecastBaselines.", fn), julia_forecast(forecast)
  ))
}

# Internal: apply a ForecastBaselines.jl function that returns a new Forecast.
update_forecast <- function(fn, forecast, ...) {
  check_setup()
  new_forecast(juliabridge::call_julia(
    paste0("ForecastBaselines.", fn), julia_forecast(forecast), ...
  ))
}

# Query Functions

#' Check if Forecast has Horizon
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_horizon <- function(forecast) {
  query_forecast("has_horizon", forecast)
}

#' Check if Forecast has Mean
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_mean <- function(forecast) {
  query_forecast("has_mean", forecast)
}

#' Check if Forecast has Median
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_median <- function(forecast) {
  query_forecast("has_median", forecast)
}

#' Check if Forecast has Intervals
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_intervals <- function(forecast) {
  query_forecast("has_intervals", forecast)
}

#' Check if Forecast has Truth Values
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_truth <- function(forecast) {
  query_forecast("has_truth", forecast)
}

#' Check if Forecast has Trajectories
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_trajectories <- function(forecast) {
  query_forecast("has_trajectories", forecast)
}

# Modification Functions

#' Add Truth Values to Forecast
#'
#' @param forecast A Forecast object
#' @param truth Numeric vector of observed values
#' @return Updated Forecast object
#' @export
add_truth <- function(forecast, truth) {
  update_forecast("add_truth", forecast, julia_vector(as.numeric(truth)))
}

#' Add Median to Forecast
#'
#' @param forecast A Forecast object
#' @param median Numeric vector of median forecasts
#' @return Updated Forecast object
#' @export
add_median <- function(forecast, median) {
  update_forecast("add_median", forecast, julia_vector(as.numeric(median)))
}

#' Add Intervals to Forecast
#'
#' @param forecast A Forecast object
#' @param intervals Interval data structure
#' @return Updated Forecast object
#' @export
add_intervals <- function(forecast, intervals) {
  update_forecast("add_intervals", forecast, intervals)
}

#' Add Trajectories to Forecast
#'
#' @param forecast A Forecast object
#' @param trajectories Matrix of trajectories
#' @return Updated Forecast object
#' @export
add_trajectories <- function(forecast, trajectories) {
  update_forecast("add_trajectories", forecast, as.matrix(trajectories))
}

#' Add Temporal Information to Forecast
#'
#' Sets the reference date on a forecast; ForecastBaselines.jl derives the
#' target dates from it, the forecast's horizon and `resolution`.
#'
#' @param forecast A Forecast object
#' @param reference_date Reference date (a `Date`, or an integer for
#'   position-based indexing)
#' @param resolution Time resolution: a number of days when `reference_date`
#'   is a `Date`, otherwise an integer step (default: 1)
#' @return Updated Forecast object
#' @export
add_temporal_info <- function(forecast, reference_date, resolution = 1) {
  if (inherits(reference_date, "Date")) {
    reference_date <- as_julia_object(juliabridge::julia_spec(
      "ForecastBaselines.Dates.Date", format(reference_date, "%Y-%m-%d")
    ))
    resolution <- as_julia_object(juliabridge::julia_spec(
      "ForecastBaselines.Dates.Day", as.integer(resolution)
    ))
  } else {
    reference_date <- as.integer(reference_date)
    resolution <- as.integer(resolution)
  }
  update_forecast(
    "add_temporal_info", forecast, reference_date, resolution = resolution
  )
}

# Filter Functions

#' Truncate Forecast Horizon
#'
#' @param forecast A Forecast object
#' @param max_h Maximum horizon to keep
#' @return Truncated Forecast object
#' @export
truncate_horizon <- function(forecast, max_h) {
  update_forecast("truncate_horizon", forecast, as.integer(max_h))
}

#' Filter Forecast to Specific Horizons
#'
#' @param forecast A Forecast object
#' @param horizons Vector of horizons to keep
#' @return Filtered Forecast object
#' @export
filter_horizons <- function(forecast, horizons) {
  update_forecast(
    "filter_horizons", forecast, julia_vector(as.integer(horizons))
  )
}

#' Filter Forecast to Specific Confidence Levels
#'
#' @param forecast A Forecast object
#' @param levels Vector of confidence levels to keep
#' @return Filtered Forecast object
#' @export
filter_levels <- function(forecast, levels) {
  update_forecast(
    "filter_levels", forecast, julia_vector(as.numeric(levels))
  )
}
