# Interval Methods for Uncertainty Quantification

#' No Interval Method
#'
#' Creates an interval method that produces only point forecasts without
#' prediction intervals.
#'
#' @return A spec for a Julia `NoInterval`, built when it is used
#'   (see [juliabridge::julia_spec()])
#' @export
#'
#' @examples
#' \dontrun{
#' method <- NoInterval()
#' fc <- forecast(fitted, interval_method = method, horizon = 1:12)
#' }
NoInterval <- function() {
  baseline_spec("NoInterval", "interval")
}

#' Empirical Interval Method
#'
#' Creates prediction intervals by bootstrapping from historical forecast
#' errors. This is a non-parametric method that doesn't assume any particular
#' distribution.
#'
#' @param n_trajectories Number of bootstrap samples to generate
#'   (default: 1000)
#' @param min_observation Minimum number of observations required (default: 1)
#' @param bootstrap_distribution Optional distribution to sample from
#'   (default: NULL)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param positivity_correction Method to ensure positive forecasts: "none",
#'   "post_clip", "truncate", or "zero_floor" (default: "none")
#' @param symmetry_correction Whether to apply symmetry correction
#'   (default: FALSE)
#' @param stepwise Whether to use stepwise intervals (default: FALSE)
#' @param return_trajectories Whether to return forecast trajectories
#'   (default: FALSE)
#'
#' @return A spec for a Julia `EmpiricalInterval`, built when it is used
#'   (see [juliabridge::julia_spec()])
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic empirical intervals
#' method <- EmpiricalInterval()
#'
#' # With more trajectories and seed
#' method <- EmpiricalInterval(n_trajectories = 2000, seed = 123)
#'
#' # With positivity correction for count data
#' method <- EmpiricalInterval(
#'   n_trajectories = 1000,
#'   positivity_correction = "post_clip"
#' )
#'
#' # Return trajectories for visualization
#' method <- EmpiricalInterval(return_trajectories = TRUE)
#' }
EmpiricalInterval <- function(n_trajectories = 1000L,
                              min_observation = 1L,
                              bootstrap_distribution = NULL,
                              seed = NULL,
                              positivity_correction = "none",
                              symmetry_correction = FALSE,
                              stepwise = FALSE,
                              return_trajectories = FALSE) {
  check_positivity(
    positivity_correction, c("none", "post_clip", "truncate", "zero_floor")
  )
  baseline_spec(
    "EmpiricalInterval", "interval",
    n_trajectories = as.integer(n_trajectories),
    min_observation = as.integer(min_observation),
    seed = if (!is.null(seed)) as.integer(seed),
    positivity_correction = julia_symbol(positivity_correction),
    symmetry_correction = symmetry_correction,
    stepwise = stepwise,
    return_trajectories = return_trajectories
  )
}

#' Parametric Interval Method
#'
#' Creates prediction intervals using parametric assumptions based on the
#' model's distribution (e.g., assuming normality for ARMA models).
#'
#' @param positivity_correction Method to ensure positive forecasts: "none"
#'   or "post_clip" (default: "none")
#'
#' @return A spec for a Julia `ParametricInterval`, built when it is used
#'   (see [juliabridge::julia_spec()])
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard parametric intervals
#' method <- ParametricInterval()
#'
#' # With positivity correction
#' method <- ParametricInterval(positivity_correction = "post_clip")
#' }
ParametricInterval <- function(positivity_correction = "none") {
  check_positivity(positivity_correction, c("none", "post_clip"))
  baseline_spec(
    "ParametricInterval", "interval",
    positivity_correction = julia_symbol(positivity_correction)
  )
}

#' Model Trajectory Interval Method
#'
#' Creates prediction intervals by simulating trajectories from the fitted
#' model. This method uses the model's own simulation mechanism.
#'
#' @param n_trajectories Number of trajectories to simulate (default: 1000)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param positivity_correction Method to ensure positive forecasts: "none",
#'   "post_clip", "truncate", or "zero_floor" (default: "none")
#' @param return_trajectories Whether to return forecast trajectories
#'   (default: FALSE)
#'
#' @return A spec for a Julia `ModelTrajectoryInterval`, built when it is used
#'   (see [juliabridge::julia_spec()])
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic model-based intervals
#' method <- ModelTrajectoryInterval()
#'
#' # With more trajectories
#' method <- ModelTrajectoryInterval(n_trajectories = 2000, seed = 456)
#'
#' # Return trajectories for analysis
#' method <- ModelTrajectoryInterval(
#'   n_trajectories = 1000,
#'   return_trajectories = TRUE
#' )
#' }
ModelTrajectoryInterval <- function(n_trajectories = 1000L,
                                    seed = NULL,
                                    positivity_correction = "none",
                                    return_trajectories = FALSE) {
  check_positivity(
    positivity_correction, c("none", "post_clip", "truncate", "zero_floor")
  )
  baseline_spec(
    "ModelTrajectoryInterval", "interval",
    n_trajectories = as.integer(n_trajectories),
    seed = if (!is.null(seed)) as.integer(seed),
    positivity_correction = julia_symbol(positivity_correction),
    return_trajectories = return_trajectories
  )
}


# Internal: the positivity correction is checked against the known options
# before it is written into Julia source as a symbol.
check_positivity <- function(x, valid) {
  if (!is.character(x) || length(x) != 1L || !x %in% valid) {
    stop(
      "positivity_correction must be one of: ", paste(valid, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

julia_symbol <- function(x) {
  juliabridge::julia(paste0(":", x))
}
