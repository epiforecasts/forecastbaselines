# Scoring Rules and Evaluation Functions using scoringutils

#' @importFrom scoringutils as_forecast_point as_forecast_quantile score
#'   summarise_scores get_metrics
NULL

#' Convert ForecastBaselines Forecast to scoringutils point forecast
#'
#' S3 method to convert a ForecastBaselines_Forecast object to a scoringutils
#' forecast_point object. This allows seamless integration with scoringutils.
#'
#' @param data A ForecastBaselines_Forecast object
#' @param ... Additional arguments (not used)
#'
#' @return A forecast_point object from scoringutils
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert and validate as point forecast
#' fc_point <- scoringutils::as_forecast_point(forecast)
#' }
as_forecast_point.ForecastBaselines_Forecast <- function(data, ...) {
  if (is.null(data$truth) || all(is.na(data$truth))) {
    stop(
      "Forecast must contain truth values for scoring. ",
      "Use add_truth() to add them."
    )
  }

  model_name <- if (!is.null(data$model_name)) data$model_name else "model"

  df <- data.frame(
    observed = data$truth,
    predicted = data$mean,
    horizon = data$horizon,
    model = model_name,
    stringsAsFactors = FALSE
  )

  as_forecast_point(df)
}

#' Convert ForecastBaselines Forecast to scoringutils quantile forecast
#'
#' S3 method to convert a ForecastBaselines_Forecast object to a scoringutils
#' forecast_quantile object. This allows seamless integration with scoringutils.
#'
#' @param data A ForecastBaselines_Forecast object with quantiles
#' @param ... Additional arguments (not used)
#'
#' @return A forecast_quantile object from scoringutils
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert and validate as quantile forecast
#' fc_quantile <- scoringutils::as_forecast_quantile(forecast)
#' }
as_forecast_quantile.ForecastBaselines_Forecast <- function(data, ...) {
  if (is.null(data$truth) || all(is.na(data$truth))) {
    stop(
      "Forecast must contain truth values for scoring. ",
      "Use add_truth() to add them."
    )
  }

  # Check if we have quantiles
  if (is.null(data$quantiles) || length(data$quantiles) == 0) {
    stop(
      "Forecast does not have quantile data. ",
      "Cannot convert to quantile forecast."
    )
  }

  model_name <- if (!is.null(data$model_name)) data$model_name else "model"

  # Assume quantiles is a matrix: rows = horizons, cols = quantile levels
  n_horizons <- length(data$horizon)
  n_quantiles <- ncol(data$quantiles)

  df <- data.frame(
    observed = rep(data$truth, each = n_quantiles),
    predicted = as.vector(t(data$quantiles)),
    quantile_level = rep(data$quantile_levels, times = n_horizons),
    horizon = rep(data$horizon, each = n_quantiles),
    model = model_name,
    stringsAsFactors = FALSE
  )

  as_forecast_quantile(df)
}
