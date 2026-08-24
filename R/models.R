# Model Functions for ForecastBaselines.jl

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
  juliaready::eval_julia("ForecastBaselines.ConstantModel()")
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
    juliaready::eval_julia("ForecastBaselines.MarginalModel()")
  } else {
    juliaready::assign_julia("p_val", as.integer(p))
    juliaready::eval_julia("ForecastBaselines.MarginalModel(p=p_val)")
  }
}

#' KDE Model
#'
#' Creates a kernel density estimation model for non-parametric forecasting.
#'
#' The kernel and the bandwidth selector belong to the estimation settings in
#' ForecastBaselines.jl, which `fit_baseline()` does not expose, so the model
#' itself has no arguments: the kernel is Gaussian and the bandwidth is
#' selected automatically.
#'
#' @return A KDEModel object
#' @export
#'
#' @examples
#' \dontrun{
#' model <- KDEModel()
#' }
KDEModel <- function() {
  check_setup()
  juliaready::eval_julia("ForecastBaselines.KDEModel()")
}

#' Last Similar Dates (LSD) Model
#'
#' Creates a seasonal forecasting model based on similar historical
#' dates.
#'
#' @param s Seasonal period (e.g., 7 for weekly, 12 for monthly)
#' @param window_width Width of the window for averaging similar dates
#'   (default: 1)
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
LSDModel <- function(s, window_width = 1L) {
  check_setup()
  juliaready::assign_julia("s_val", as.integer(s))
  juliaready::assign_julia("w_val", as.integer(window_width))
  juliaready::eval_julia("ForecastBaselines.LSDModel(s=s_val, w=w_val)")
}

#' OLS Model
#'
#' Creates an ordinary least squares model with polynomial trend.
#'
#' @param degree Polynomial degree (default: 1 for linear trend)
#' @param differencing Order of differencing (default: 0 for no differencing)
#'
#' @return An OLSModel object
#' @export
#'
#' @examples
#' \dontrun{
#' # Linear trend
#' model <- OLSModel(degree = 1)
#'
#' # Quadratic trend with differencing
#' model <- OLSModel(degree = 2, differencing = 1)
#' }
OLSModel <- function(degree = 1L, differencing = 0L) {
  check_setup()
  # Julia API uses 'p' (polynomial degree) and 'd' (differencing)
  juliaready::assign_julia("p_val", as.integer(degree))
  juliaready::assign_julia("d_val", as.integer(differencing))
  juliaready::eval_julia("ForecastBaselines.OLSModel(p=p_val, d=d_val)")
}

#' IDS Model
#'
#' Creates an Increase-Decrease-Stable model for trend detection.
#'
#' @param window_size Window size for trend calculation (default: 3)
#'
#' @return An IDSModel object
#' @export
#'
#' @examples
#' \dontrun{
#' model <- IDSModel()
#' model <- IDSModel(window_size = 5)
#' }
IDSModel <- function(window_size = 3L) {
  check_setup()
  juliaready::assign_julia("p_val", as.integer(window_size))
  juliaready::eval_julia("ForecastBaselines.IDSModel(p=p_val)")
}

#' STL Model
#'
#' Creates a Seasonal-Trend decomposition using Loess model.
#'
#' The decomposition always includes a trend component, and the loess
#' smoothing settings, robust iterations among them, belong to the estimation
#' settings in ForecastBaselines.jl, which `fit_baseline()` does not expose.
#' The seasonal period is therefore the only argument.
#'
#' @param s Seasonal period
#'
#' @return An STLModel object
#' @export
#'
#' @examples
#' \dontrun{
#' # Monthly seasonality
#' model <- STLModel(s = 12)
#' }
STLModel <- function(s) {
  check_setup()
  juliaready::assign_julia("s_val", as.integer(s))
  juliaready::eval_julia("ForecastBaselines.STLModel(s=s_val)")
}

#' ARMA Model
#'
#' Creates an AutoRegressive Moving Average model. The series is modelled
#' around a deterministic mean
#' \deqn{\mu_t = \theta_1 +
#'   \sum_{j=1}^k [\theta_{2j} \sin(2\pi j t / s) +
#'                 \theta_{2j+1} \cos(2\pi j t / s)] + \theta_{2k+2} t}
#' where the harmonic terms are present when `s > 0` and the linear term when
#' `include_drift = TRUE`. An intercept is always estimated.
#'
#' @param p AR order (default: 0)
#' @param q MA order (default: 0)
#' @param s Seasonal period (default: 0 for no seasonality)
#' @param k Number of harmonic waves used to represent seasonality
#'   (default: 1). Higher values allow sharper, asymmetric seasonal shapes.
#'   Only meaningful when `s > 0`. Each wave adds two parameters to the mean
#'   function, which are estimated by a derivative-free search started from
#'   zero, so fits become unreliable as `k` grows: on a series with a mean of
#'   100 and two genuine harmonics, `k = 1` recovers the seasonal shape while
#'   `k = 2` and above collapse towards a zero intercept. Check the fitted
#'   mean against the data before trusting a model with several waves.
#' @param include_drift Whether to include a linear trend in the mean
#'   (default: FALSE)
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
#'
#' # Annual seasonality with sharper peaks, plus a linear trend
#' model <- ARMAModel(p = 1, s = 52, k = 2, include_drift = TRUE)
#' }
ARMAModel <- function(p = 0L, q = 0L, s = 0L, k = 1L,
                      include_drift = FALSE) {
  check_harmonics(s, k)
  trend <- julia_bool(include_drift)
  check_setup()

  juliaready::assign_julia("p_val", as.integer(p))
  juliaready::assign_julia("q_val", as.integer(q))
  juliaready::assign_julia("s_val", as.integer(s))
  juliaready::assign_julia("k_val", as.integer(k))

  # arma_model() is defined in inst/julia/forecast_helpers.jl; it builds the
  # mean function for arbitrary k on top of the ARMAModel inner constructor.
  juliaready::eval_julia(sprintf(
    "arma_model(p_val, q_val, s_val, k_val, %s)", trend
  ))
}

#' INARCH Model
#'
#' Creates an Integer-valued ARCH model for count time series. Seasonality
#' enters the conditional mean through `k` harmonic waves of period `s`, and
#' counts follow a Poisson distribution unless `nb = TRUE`.
#'
#' @param p Order of the INARCH model (default: 1)
#' @param s Seasonal period (default: 0 for no seasonality)
#' @param k Number of harmonic waves used to represent seasonality
#'   (default: 1). Only meaningful when `s > 0`. Values above 1 are fitted
#'   correctly but not forecast correctly: `point_forecast()` and
#'   `interval_forecast()` in ForecastBaselines.jl evaluate every harmonic
#'   coefficient at the fundamental frequency, so forecasts do not reflect
#'   the seasonality that was estimated. Constructing such a model warns.
#' @param nb Whether to use a negative binomial conditional distribution
#'   instead of Poisson, allowing for overdispersion (default: FALSE)
#'
#' @return An INARCHModel object
#' @export
#'
#' @examples
#' \dontrun{
#' model <- INARCHModel(p = 1)
#'
#' # Overdispersed counts with weekly data and annual seasonality
#' model <- INARCHModel(p = 1, s = 52, nb = TRUE)
#' }
INARCHModel <- function(p = 1L, s = 0L, k = 1L, nb = FALSE) {
  check_harmonics(s, k)
  if (as.integer(k) > 1L) {
    warning(
      "INARCH forecasts evaluate every harmonic at the fundamental ",
      "frequency, so seasonality estimated with k > 1 is not carried into ",
      "the forecasts.",
      call. = FALSE
    )
  }
  negbin <- julia_bool(nb)
  check_setup()

  juliaready::assign_julia("p_val", as.integer(p))
  juliaready::assign_julia("s_val", as.integer(s))
  juliaready::assign_julia("k_val", as.integer(k))
  juliaready::eval_julia(sprintf(
    "ForecastBaselines.INARCHModel(p=p_val, s=s_val, k=k_val, nb=%s)", negbin
  ))
}

#' ETS Model
#'
#' Creates an Error-Trend-Season exponential smoothing model.
#'
#' @param error_type Error type: "A" (additive), "M" (multiplicative), or
#'   "N" (none)
#' @param trend_type Trend type: "A" (additive), "M" (multiplicative), "Ad"
#'   (damped additive), "Md" (damped multiplicative), or "N" (none). Damping
#'   is requested through "Ad" or "Md".
#' @param season_type Season type: "A" (additive), "M" (multiplicative), or
#'   "N" (none)
#' @param s Seasonal period (required if season_type is not "N")
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
                     s = NULL) {
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
    juliaready::assign_julia("s_val", as.integer(s))
    julia_code <- sprintf(
      paste0(
        'ForecastBaselines.ETSModel(error="%s", trend="%s", ',
        'season="%s", s=s_val)'
      ),
      error_type, trend_type, season_type
    )
  }

  juliaready::eval_julia(julia_code)
}


# Internal: harmonic seasonality only has meaning alongside a seasonal period.
check_harmonics <- function(s, k) {
  if (as.integer(k) > 1L && as.integer(s) == 0L) {
    stop(
      "'k' harmonic waves require a seasonal period: set 's' to a positive ",
      "value, or leave 'k' at 1.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Internal: render an R logical as a Julia boolean literal, rejecting values
# that cannot be honoured rather than silently treating them as FALSE.
julia_bool <- function(x, arg = deparse(substitute(x))) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", arg, "' must be either TRUE or FALSE", call. = FALSE)
  }
  if (x) "true" else "false"
}
