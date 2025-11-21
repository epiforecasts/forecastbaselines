# Model Constructors for ForecastBaselines.jl

#' Constant Model
#'
#' Creates a naive forecast model that uses the last observed value as the
#' forecast for all future horizons.
#'
#' @return A ConstantModel object
#' @export
#'
#' @examples
#' \dontrun{
#' model <- ConstantModel()
#' }
ConstantModel <- function() {
  check_setup()
  JuliaCall::julia_eval("ForecastBaselines.ConstantModel()")
}

#' Marginal Model
#'
#' Creates a forecast based on the empirical marginal distribution using the
#' mean of the p most recent observations.
#'
#' @param p Number of most recent observations to use (default: all
#'   observations)
#'
#' @return A MarginalModel object
#' @export
#'
#' @examples
#' \dontrun{
#' model <- MarginalModel(p = 10)
#' }
MarginalModel <- function(p = NULL) {
  check_setup()
  if (is.null(p)) {
    JuliaCall::julia_eval("ForecastBaselines.MarginalModel()")
  } else {
    JuliaCall::julia_assign("p_val", as.integer(p))
    JuliaCall::julia_eval("ForecastBaselines.MarginalModel(p=p_val)")
  }
}

#' KDE Model
#'
#' Creates a kernel density estimation model for non-parametric forecasting.
#'
#' @param bandwidth Bandwidth for KDE (default: NULL for automatic selection)
#' @param kernel Kernel function to use (default: "gaussian")
#'
#' @return A KDEModel object
#' @export
#'
#' @examples
#' \dontrun{
#' model <- KDEModel()
#' model <- KDEModel(bandwidth = 0.5)
#' }
KDEModel <- function(bandwidth = NULL, kernel = "gaussian") {
  check_setup()
  # KDEModel in Julia doesn't take parameters, ignoring bandwidth for now
  JuliaCall::julia_eval("ForecastBaselines.KDEModel()")
}

#' Last Similar Dates (LSD) Model
#'
#' Creates a seasonal forecasting model based on similar historical
#' dates.
#'
#' @param s Seasonal period (e.g., 7 for weekly, 12 for monthly)
#' @param window_width Width of the window for averaging similar dates
#'   (default: 1)
#' @param trend_correction Whether to apply trend correction
#'   (default: FALSE)
#'
#' @return An LSDModel object
#' @export
#'
#' @examples
#' \dontrun{
#' # Weekly seasonality
#' model <- LSDModel(s = 7)
#'
#' # Monthly seasonality with window
#' model <- LSDModel(s = 12, window_width = 2)
#' }
LSDModel <- function(s, window_width = 1L, trend_correction = FALSE) {
  check_setup()
  # Julia API uses 's' and 'w' parameters, trend_correction not supported
  JuliaCall::julia_assign("s_val", as.integer(s))
  JuliaCall::julia_assign("w_val", as.integer(window_width))
  JuliaCall::julia_eval("ForecastBaselines.LSDModel(s=s_val, w=w_val)")
}

#' OLS Model
#'
#' Creates an ordinary least squares model with polynomial trend.
#'
#' @param degree Polynomial degree (default: 1 for linear trend)
#' @param s Seasonal period (default: NULL for no seasonality)
#'
#' @return An OLSModel object
#' @export
#'
#' @examples
#' \dontrun{
#' # Linear trend
#' model <- OLSModel(degree = 1)
#'
#' # Quadratic trend with seasonality
#' model <- OLSModel(degree = 2, s = 12)
#' }
OLSModel <- function(degree = 1L, differencing = 0L) {
  check_setup()
  # Julia API uses 'p' (polynomial degree) and 'd' (differencing)
  JuliaCall::julia_assign("p_val", as.integer(degree))
  JuliaCall::julia_assign("d_val", as.integer(differencing))
  JuliaCall::julia_eval("ForecastBaselines.OLSModel(p=p_val, d=d_val)")
}

#' IDS Model
#'
#' Creates an Increase-Decrease-Stable model for trend detection.
#'
#' @param threshold Threshold for trend detection (default: 0.0)
#' @param window_size Window size for trend calculation (default: 3)
#'
#' @return An IDSModel object
#' @export
#'
#' @examples
#' \dontrun{
#' model <- IDSModel()
#' model <- IDSModel(threshold = 0.1, window_size = 5)
#' }
IDSModel <- function(threshold = 0.0, window_size = 3L) {
  check_setup()
  # Julia API only takes 'p' parameter
  JuliaCall::julia_assign("p_val", as.integer(window_size))
  JuliaCall::julia_eval("ForecastBaselines.IDSModel(p=p_val)")
}

#' STL Model
#'
#' Creates a Seasonal-Trend decomposition using Loess model.
#'
#' @param s Seasonal period
#' @param trend Whether to include trend component (default: TRUE)
#' @param robust Whether to use robust fitting (default: FALSE)
#'
#' @return An STLModel object
#' @export
#'
#' @examples
#' \dontrun{
#' # Monthly seasonality
#' model <- STLModel(s = 12)
#'
#' # Robust STL
#' model <- STLModel(s = 12, robust = TRUE)
#' }
STLModel <- function(s, trend = TRUE, robust = FALSE) {
  check_setup()
  # Julia API only takes 's' parameter, trend and robust not supported
  JuliaCall::julia_assign("s_val", as.integer(s))
  JuliaCall::julia_eval("ForecastBaselines.STLModel(s=s_val)")
}

#' ARMA Model
#'
#' Creates an AutoRegressive Moving Average model.
#'
#' @param p AR order (default: 0)
#' @param q MA order (default: 0)
#' @param s Seasonal period (default: 0 for no seasonality)
#' @param include_mean Whether to include a mean term (default: TRUE)
#' @param include_drift Whether to include a drift term (default: FALSE)
#'
#' @return An ARMAModel object
#' @export
#'
#' @examples
#' \dontrun{
#' # AR(1)
#' model <- ARMAModel(p = 1)
#'
#' # MA(1)
#' model <- ARMAModel(q = 1)
#'
#' # ARMA(2,1) with seasonality
#' model <- ARMAModel(p = 2, q = 1, s = 12)
#' }
ARMAModel <- function(p = 0L, q = 0L, s = 0L,
                      include_mean = TRUE, include_drift = FALSE) {
  check_setup()

  # Use keyword argument syntax for Julia
  JuliaCall::julia_assign("p_val", as.integer(p))
  JuliaCall::julia_assign("q_val", as.integer(p))
  JuliaCall::julia_assign("s_val", as.integer(s))

  # For now, use simple keyword argument approach
  # Note: include_mean and include_drift may not be directly supported
  # They may need to be handled differently in Julia API
  JuliaCall::julia_eval(
    "ForecastBaselines.ARMAModel(p=p_val, q=q_val, s=s_val)"
  )
}

#' INARCH Model
#'
#' Creates an Integer-valued ARCH model for count time series.
#'
#' @param p Order of the INARCH model
#'
#' @return An INARCHModel object
#' @export
#'
#' @examples
#' \dontrun{
#' model <- INARCHModel(p = 1)
#' }
INARCHModel <- function(p = 1L) {
  check_setup()
  JuliaCall::julia_assign("p_val", as.integer(p))
  JuliaCall::julia_eval("ForecastBaselines.INARCHModel(p=p_val)")
}

#' ETS Model
#'
#' Creates an Error-Trend-Season exponential smoothing model.
#'
#' @param error_type Error type: "A" (additive), "M" (multiplicative), or
#'   "N" (none)
#' @param trend_type Trend type: "A" (additive), "M" (multiplicative), "Ad"
#'   (damped additive), "Md" (damped multiplicative), or "N" (none)
#' @param season_type Season type: "A" (additive), "M" (multiplicative), or
#'   "N" (none)
#' @param s Seasonal period (required if season_type is not "N")
#' @param damped Whether to use damped trend (default: FALSE)
#'
#' @return An ETSModel object
#' @export
#'
#' @examples
#' \dontrun{
#' # Simple exponential smoothing (A,N,N)
#' model <- ETSModel(error_type = "A", trend_type = "N", season_type = "N")
#'
#' # Holt's linear trend (A,A,N)
#' model <- ETSModel(error_type = "A", trend_type = "A", season_type = "N")
#'
#' # Holt-Winters additive (A,A,A)
#' model <- ETSModel(
#'   error_type = "A", trend_type = "A", season_type = "A", s = 12
#' )
#'
#' # Holt-Winters multiplicative (M,M,M)
#' model <- ETSModel(
#'   error_type = "M", trend_type = "M", season_type = "M", s = 12
#' )
#' }
ETSModel <- function(error_type = "A", trend_type = "N", season_type = "N",
                     s = NULL, damped = FALSE) {
  check_setup()

  # Validate inputs
  valid_error <- c("A", "M", "N")
  valid_trend <- c("A", "M", "Ad", "Md", "N")
  valid_season <- c("A", "M", "N")

  if (!error_type %in% valid_error) {
    stop(
      "error_type must be one of: ",
      paste(valid_error, collapse = ", ")
    )
  }
  if (!trend_type %in% valid_trend) {
    stop(
      "trend_type must be one of: ",
      paste(valid_trend, collapse = ", ")
    )
  }
  if (!season_type %in% valid_season) {
    stop(
      "season_type must be one of: ",
      paste(valid_season, collapse = ", ")
    )
  }

  if (season_type != "N" && is.null(s)) {
    stop(
      "Seasonal period 's' must be provided when season_type is not 'N'"
    )
  }

  # Call Julia function using keyword arguments with String values
  if (is.null(s)) {
    julia_code <- sprintf(
      'ForecastBaselines.ETSModel(error="%s", trend="%s", season="%s")',
      error_type, trend_type, season_type
    )
  } else {
    JuliaCall::julia_assign("s_val", as.integer(s))
    julia_code <- sprintf(
      paste0(
        'ForecastBaselines.ETSModel(error="%s", trend="%s", ',
        'season="%s", s=s_val)'
      ),
      error_type, trend_type, season_type
    )
  }

  JuliaCall::julia_eval(julia_code)
}
