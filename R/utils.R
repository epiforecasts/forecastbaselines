# Forecast Utility Functions

# Internal helper to get Julia forecast reference
get_fc_ref <- function(forecast) {
  if (inherits(forecast, "ForecastBaselines_Forecast") &&
    !is.null(attr(forecast, "julia_ref"))) {
    return(attr(forecast, "julia_ref"))
  }
  # Fallback: assign to a temp variable
  juliaready::assign_julia("fc_temp", forecast)
  return("fc_temp")
}

# Query Functions

#' Check if Forecast has Horizon
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_horizon <- function(forecast) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  as.logical(juliaready::eval_julia(
    sprintf("ForecastBaselines.has_horizon(%s)", fc_var)
  ))
}

#' Check if Forecast has Mean
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_mean <- function(forecast) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  as.logical(juliaready::eval_julia(
    sprintf("ForecastBaselines.has_mean(%s)", fc_var)
  ))
}

#' Check if Forecast has Median
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_median <- function(forecast) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  as.logical(juliaready::eval_julia(
    sprintf("ForecastBaselines.has_median(%s)", fc_var)
  ))
}

#' Check if Forecast has Intervals
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_intervals <- function(forecast) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  as.logical(juliaready::eval_julia(
    sprintf("ForecastBaselines.has_intervals(%s)", fc_var)
  ))
}

#' Check if Forecast has Truth Values
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_truth <- function(forecast) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  as.logical(juliaready::eval_julia(
    sprintf("ForecastBaselines.has_truth(%s)", fc_var)
  ))
}

#' Check if Forecast has Trajectories
#'
#' @param forecast A Forecast object
#' @return Logical TRUE/FALSE
#' @export
has_trajectories <- function(forecast) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  as.logical(juliaready::eval_julia(
    sprintf("ForecastBaselines.has_trajectories(%s)", fc_var)
  ))
}

# Modification Functions

#' Add Truth Values to Forecast
#'
#' @param forecast A Forecast object
#' @param truth Numeric vector of observed values
#' @return Updated Forecast object
#' @export
add_truth <- function(forecast, truth) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  juliaready::assign_julia("truth_vals", as.numeric(truth))
  # Store result in Julia and get reference
  fc_id <- paste0("fc_", sample.int(.Machine$integer.max, 1))
  juliaready::command_julia(
    sprintf(
      "%s = ForecastBaselines.add_truth(%s, truth_vals)",
      fc_id, fc_var
    )
  )
  result <- eval_julia_value(sprintf("forecast_to_r(%s)", fc_id))
  class(result) <- c("ForecastBaselines_Forecast", "list")
  attr(result, "julia_ref") <- fc_id
  result
}

#' Add Median to Forecast
#'
#' @param forecast A Forecast object
#' @param median Numeric vector of median forecasts
#' @return Updated Forecast object
#' @export
add_median <- function(forecast, median) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  juliaready::assign_julia("median_vals", as.numeric(median))
  # Store result in Julia and get reference
  fc_id <- paste0("fc_", sample.int(.Machine$integer.max, 1))
  juliaready::command_julia(
    sprintf(
      "%s = ForecastBaselines.add_median(%s, median_vals)",
      fc_id, fc_var
    )
  )
  result <- eval_julia_value(sprintf("forecast_to_r(%s)", fc_id))
  class(result) <- c("ForecastBaselines_Forecast", "list")
  attr(result, "julia_ref") <- fc_id
  result
}

#' Add Intervals to Forecast
#'
#' @param forecast A Forecast object
#' @param intervals Interval data structure
#' @return Updated Forecast object
#' @export
add_intervals <- function(forecast, intervals) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  juliaready::assign_julia("intervals_obj", intervals)
  result <- eval_julia_value(
    sprintf(
      "forecast_to_r(ForecastBaselines.add_intervals(%s, intervals_obj))",
      fc_var
    )
  )
  class(result) <- c("ForecastBaselines_Forecast", "list")
  result
}

#' Add Trajectories to Forecast
#'
#' @param forecast A Forecast object
#' @param trajectories Matrix of trajectories
#' @return Updated Forecast object
#' @export
add_trajectories <- function(forecast, trajectories) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  juliaready::assign_julia("traj_mat", as.matrix(trajectories))
  result <- eval_julia_value(
    sprintf(
      "forecast_to_r(ForecastBaselines.add_trajectories(%s, traj_mat))",
      fc_var
    )
  )
  class(result) <- c("ForecastBaselines_Forecast", "list")
  result
}

#' Add Temporal Information to Forecast
#'
#' @param forecast A Forecast object
#' @param reference_date Reference date
#' @param target_date Target date
#' @param resolution Time resolution
#' @return Updated Forecast object
#' @export
add_temporal_info <- function(forecast, reference_date, target_date,
                              resolution) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  juliaready::assign_julia("ref_date", reference_date)
  juliaready::assign_julia("tgt_date", target_date)
  juliaready::assign_julia("res_val", resolution)
  result <- eval_julia_value(
    sprintf(
      paste0(
        "forecast_to_r(ForecastBaselines.add_temporal_info(",
        "%s, ref_date, tgt_date, res_val))"
      ),
      fc_var
    )
  )
  class(result) <- c("ForecastBaselines_Forecast", "list")
  result
}

# Filter Functions

#' Truncate Forecast Horizon
#'
#' @param forecast A Forecast object
#' @param max_h Maximum horizon to keep
#' @return Truncated Forecast object
#' @export
truncate_horizon <- function(forecast, max_h) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  juliaready::assign_julia("max_horizon", as.integer(max_h))
  # Store result in Julia and get reference
  fc_id <- paste0("fc_", sample.int(.Machine$integer.max, 1))
  juliaready::command_julia(
    sprintf(
      "%s = ForecastBaselines.truncate_horizon(%s, max_horizon)",
      fc_id, fc_var
    )
  )
  result <- eval_julia_value(sprintf("forecast_to_r(%s)", fc_id))
  class(result) <- c("ForecastBaselines_Forecast", "list")
  attr(result, "julia_ref") <- fc_id
  result
}

#' Filter Forecast to Specific Horizons
#'
#' @param forecast A Forecast object
#' @param horizons Vector of horizons to keep
#' @return Filtered Forecast object
#' @export
filter_horizons <- function(forecast, horizons) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  juliaready::assign_julia("h_vec", as.integer(horizons))
  # Store result in Julia and get reference
  fc_id <- paste0("fc_", sample.int(.Machine$integer.max, 1))
  juliaready::command_julia(
    sprintf(
      "%s = ForecastBaselines.filter_horizons(%s, h_vec)",
      fc_id, fc_var
    )
  )
  result <- eval_julia_value(sprintf("forecast_to_r(%s)", fc_id))
  class(result) <- c("ForecastBaselines_Forecast", "list")
  attr(result, "julia_ref") <- fc_id
  result
}

#' Filter Forecast to Specific Confidence Levels
#'
#' @param forecast A Forecast object
#' @param levels Vector of confidence levels to keep
#' @return Filtered Forecast object
#' @export
filter_levels <- function(forecast, levels) {
  check_setup()
  fc_var <- get_fc_ref(forecast)
  juliaready::assign_julia("lvls_vec", as.numeric(levels))
  result <- eval_julia_value(
    sprintf(
      "forecast_to_r(ForecastBaselines.filter_levels(%s, lvls_vec))",
      fc_var
    )
  )
  class(result) <- c("ForecastBaselines_Forecast", "list")
  result
}
