# Interval Methods for Uncertainty Quantification

#' No Interval Method
#'
#' Creates an interval method that produces only point forecasts without
#' prediction intervals.
#'
#' @return A NoInterval object
#' @export
#'
#' @examples
#' \dontrun{
#' method <- NoInterval()
#' fc <- forecast(fitted, interval_method = method, horizon = 1:12)
#' }
NoInterval <- function() {
  check_setup()
  JuliaCall::julia_eval("ForecastBaselines.NoInterval()")
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
#' @return An EmpiricalInterval object
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
  check_setup()

  # Validate positivity_correction
  valid_pos_corr <- c("none", "post_clip", "truncate", "zero_floor")
  if (!positivity_correction %in% valid_pos_corr) {
    stop(
      "positivity_correction must be one of: ",
      paste(valid_pos_corr, collapse = ", ")
    )
  }

  # Build Julia code string for the constructor
  args <- list()
  args$n_trajectories <- as.integer(n_trajectories)
  args$min_observation <- as.integer(min_observation)

  if (!is.null(seed)) {
    args$seed <- as.integer(seed)
  }

  args$positivity_correction <- sprintf(":%s", positivity_correction)
  args$symmetry_correction <- symmetry_correction
  args$stepwise <- stepwise
  args$return_trajectories <- return_trajectories

  # Construct Julia call
  args_list <- c(
    sprintf("n_trajectories=%d", args$n_trajectories),
    sprintf("min_observation=%d", args$min_observation),
    if (!is.null(seed)) sprintf("seed=%d", args$seed),
    sprintf("positivity_correction=%s", args$positivity_correction),
    sprintf(
      "symmetry_correction=%s",
      tolower(as.character(args$symmetry_correction))
    ),
    sprintf("stepwise=%s", tolower(as.character(args$stepwise))),
    sprintf(
      "return_trajectories=%s",
      tolower(as.character(args$return_trajectories))
    )
  )
  args_str <- paste(args_list, collapse = ", ")

  julia_code <- sprintf("ForecastBaselines.EmpiricalInterval(%s)", args_str)
  JuliaCall::julia_eval(julia_code)
}

#' Parametric Interval Method
#'
#' Creates prediction intervals using parametric assumptions based on the
#' model's distribution (e.g., assuming normality for ARMA models).
#'
#' @param positivity_correction Method to ensure positive forecasts: "none"
#'   or "post_clip" (default: "none")
#'
#' @return A ParametricInterval object
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
  check_setup()

  valid_pos_corr <- c("none", "post_clip")
  if (!positivity_correction %in% valid_pos_corr) {
    stop(
      "positivity_correction must be one of: ",
      paste(valid_pos_corr, collapse = ", ")
    )
  }

  julia_code <- sprintf(
    "ForecastBaselines.ParametricInterval(positivity_correction=:%s)",
    positivity_correction
  )
  JuliaCall::julia_eval(julia_code)
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
#' @return A ModelTrajectoryInterval object
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
  check_setup()

  valid_pos_corr <- c("none", "post_clip", "truncate", "zero_floor")
  if (!positivity_correction %in% valid_pos_corr) {
    stop(
      "positivity_correction must be one of: ",
      paste(valid_pos_corr, collapse = ", ")
    )
  }

  args <- list()
  args$n_trajectories <- as.integer(n_trajectories)

  if (!is.null(seed)) {
    args$seed <- as.integer(seed)
  }

  args$positivity_correction <- sprintf(":%s", positivity_correction)
  args$return_trajectories <- return_trajectories

  # Construct Julia call
  args_list <- c(
    sprintf("n_trajectories=%d", args$n_trajectories),
    if (!is.null(seed)) sprintf("seed=%d", args$seed),
    sprintf("positivity_correction=%s", args$positivity_correction),
    sprintf(
      "return_trajectories=%s",
      tolower(as.character(args$return_trajectories))
    )
  )
  args_str <- paste(args_list, collapse = ", ")

  julia_code <- sprintf(
    "ForecastBaselines.ModelTrajectoryInterval(%s)",
    args_str
  )
  JuliaCall::julia_eval(julia_code)
}
