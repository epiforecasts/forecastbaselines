# Scoring Rules and Evaluation Functions using scoringutils

#' @importFrom scoringutils as_forecast_point as_forecast_quantile score
#'   summarise_scores get_metrics
NULL

#' Convert ForecastBaselines Forecast to hubverse format
#'
#' Converts a ForecastBaselines_Forecast object to the hubverse format used by
#' hubVis for visualization. Returns a list with model_output and target_data
#' data frames.
#'
#' @param data A ForecastBaselines_Forecast object with quantiles, or a
#'   scoringutils forecast_quantile object
#' @param origin_date Optional origin date for converting horizons to target_dates.
#'   If NULL (default), uses horizon numbers as target_date values.
#'   Can be a Date object or character string in "YYYY-MM-DD" format.
#' @param horizon_unit Character string specifying the time unit for horizons.
#'   One of "day" (default), "week", or "month". Used when origin_date is provided.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{model_output}{Data frame with columns: target_date, value, model_id,
#'       output_type, output_type_id}
#'     \item{target_data}{Data frame with columns: date, observation}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic conversion (horizons as dates)
#' hubverse_data <- as_hubverse(forecast)
#'
#' # With origin date for proper dates
#' hubverse_data <- as_hubverse(forecast,
#'   origin_date = "2024-01-01",
#'   horizon_unit = "week"
#' )
#'
#' # Visualize with hubVis
#' if (requireNamespace("hubVis", quietly = TRUE)) {
#'   hubVis::plot_step_ahead_model_output(
#'     hubUtils::as_model_out_tbl(hubverse_data$model_output),
#'     hubverse_data$target_data
#'   )
#' }
#' }
as_hubverse <- function(data,
                        origin_date = NULL,
                        horizon_unit = c("day", "week", "month")) {
  UseMethod("as_hubverse")
}

#' @export
as_hubverse.ForecastBaselines_Forecast <- function(data,
                                                    origin_date = NULL,
                                                    horizon_unit = c("day", "week", "month")) {
  # First convert to scoringutils quantile format
  if (is.null(data$quantiles) || length(data$quantiles) == 0) {
    stop(
      "Forecast does not have quantile data. ",
      "Cannot convert to hubverse format."
    )
  }

  if (is.null(data$truth) || all(is.na(data$truth))) {
    stop(
      "Forecast must contain truth values. ",
      "Use add_truth() to add them."
    )
  }

  fc_quantile <- as_forecast_quantile.ForecastBaselines_Forecast(data)
  as_hubverse.default(fc_quantile, origin_date, horizon_unit)
}

#' @export
as_hubverse.default <- function(data,
                                 origin_date = NULL,
                                 horizon_unit = c("day", "week", "month")) {
  horizon_unit <- match.arg(horizon_unit)

  # Check for required columns
  required_cols <- c("observed", "predicted", "quantile_level", "horizon", "model")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Calculate target dates
  if (is.null(origin_date)) {
    # Use horizon as-is
    target_dates <- data$horizon
  } else {
    # Convert origin_date if character
    if (is.character(origin_date)) {
      origin_date <- as.Date(origin_date)
    }

    # Calculate dates based on horizon and unit
    multiplier <- switch(horizon_unit,
      "day" = 1,
      "week" = 7,
      "month" = 30 # Approximate
    )

    target_dates <- origin_date + (data$horizon * multiplier)
  }

  # Create model output table
  model_output <- data.frame(
    target_date = target_dates,
    value = data$predicted,
    model_id = data$model,
    output_type = "quantile",
    output_type_id = data$quantile_level,
    stringsAsFactors = FALSE
  )

  # Create target data table (unique observed values per date)
  target_data <- unique(data.frame(
    date = target_dates,
    observation = data$observed,
    stringsAsFactors = FALSE
  ))

  list(
    model_output = model_output,
    target_data = target_data
  )
}

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
