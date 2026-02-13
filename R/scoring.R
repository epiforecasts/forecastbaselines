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
#' @param start_date Optional start date for the time series.
#'   If NULL (default), uses horizon numbers as target_date values.
#'   Can be a Date object or character string in "YYYY-MM-DD" format.
#' @param horizon_unit Character string specifying the time unit for horizons.
#'   One of "day" (default), "week", or "month". Used when start_date is provided.
#' @param observed_data Optional numeric vector of all observed values (training + test).
#'   If provided, will be used to populate target_data with complete historical context.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{model_output}{Data frame (model_out_tbl) with columns: target_date, value,
#'       model_id, output_type, output_type_id}
#'     \item{target_data}{Data frame with columns: date, observation}
#'   }
#'
#' @details
#' The hubVis package is not on CRAN and must be installed from GitHub:
#' \code{remotes::install_github("hubverse-org/hubVis")}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # With observed data for complete visualization context
#' hubverse_data <- as_hubverse(forecast,
#'   start_date = "2024-01-01",
#'   horizon_unit = "week",
#'   observed_data = full_data
#' )
#'
#' # Install hubVis from GitHub for visualization
#' # remotes::install_github("hubverse-org/hubVis")
#'
#' # Visualize with hubVis (if installed)
#' # Note: as_hubverse() automatically converts model_output to model_out_tbl
#' if (requireNamespace("hubVis", quietly = TRUE)) {
#'   hubVis::plot_step_ahead_model_output(
#'     hubverse_data$model_output,
#'     hubverse_data$target_data
#'   )
#' }
#' }
as_hubverse <- function(data,
                        start_date = NULL,
                        horizon_unit = c("day", "week", "month"),
                        observed_data = NULL) {
  UseMethod("as_hubverse")
}

#' @export
as_hubverse.ForecastBaselines_Forecast <- function(data,
                                                    start_date = NULL,
                                                    horizon_unit = c("day", "week", "month"),
                                                    observed_data = NULL) {
  # Check for quantiles
  if (is.null(data$quantiles) || length(data$quantiles) == 0) {
    stop(
      "Forecast does not have quantile data. ",
      "Cannot convert to hubverse format."
    )
  }

  # If observed_data provided, we can work without truth in the forecast
  if (!is.null(observed_data)) {
    # Extract just the test/forecast period from observed_data to use as "truth"
    n_horizons <- length(data$horizon)
    n_obs <- length(observed_data)

    if (n_obs < n_horizons) {
      stop("observed_data has fewer observations than forecast horizons")
    }

    # Use the last n_horizons observations as truth for conversion
    temp_truth <- observed_data[(n_obs - n_horizons + 1):n_obs]

    # Temporarily add truth for conversion
    data_with_truth <- data
    data_with_truth$truth <- temp_truth

    fc_quantile <- as_forecast_quantile.ForecastBaselines_Forecast(data_with_truth)
  } else {
    # Fall back to requiring truth in the forecast
    if (is.null(data$truth) || all(is.na(data$truth))) {
      stop(
        "Forecast must contain truth values or observed_data must be provided. ",
        "Use add_truth() or provide observed_data parameter."
      )
    }

    fc_quantile <- as_forecast_quantile.ForecastBaselines_Forecast(data)
  }

  as_hubverse.default(fc_quantile, start_date, horizon_unit, observed_data)
}

#' @export
as_hubverse.default <- function(data,
                                 start_date = NULL,
                                 horizon_unit = c("day", "week", "month"),
                                 observed_data = NULL) {
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

  # Calculate target dates for forecasts
  if (is.null(start_date)) {
    # Use horizon as-is
    target_dates <- data$horizon
  } else {
    # Convert start_date if character
    if (is.character(start_date)) {
      start_date <- as.Date(start_date)
    }

    # Calculate dates based on horizon and unit
    multiplier <- switch(horizon_unit,
      "day" = 1,
      "week" = 7,
      "month" = 30 # Approximate
    )

    # If observed_data provided, calculate where forecast starts in the full series
    if (!is.null(observed_data)) {
      # Forecast starts after training data
      forecast_start_index <- length(observed_data) - length(unique(data$horizon)) + 1
      target_dates <- start_date + ((forecast_start_index - 1 + data$horizon - 1) * multiplier)
    } else {
      # Legacy behavior: horizons relative to start_date
      target_dates <- start_date + ((data$horizon - 1) * multiplier)
    }
  }

  # Create model output table
  model_output <- data.frame(
    target_date = target_dates,
    value = data$predicted,
    model_id = data$model,
    output_type = "quantile",
    output_type_id = round(data$quantile_level, 3),
    stringsAsFactors = FALSE
  )

  # Create target data table
  if (!is.null(observed_data)) {
    # Use provided observed_data for complete time series
    if (is.null(start_date)) {
      stop("start_date must be provided when using observed_data")
    }

    multiplier <- switch(horizon_unit,
      "day" = 1,
      "week" = 7,
      "month" = 30
    )

    all_dates <- start_date + (0:(length(observed_data) - 1)) * multiplier
    target_data <- data.frame(
      date = all_dates,
      observation = observed_data,
      stringsAsFactors = FALSE
    )
  } else {
    # Fall back to using observed values from forecast (truth only)
    target_data <- unique(data.frame(
      date = target_dates,
      observation = data$observed,
      stringsAsFactors = FALSE
    ))
  }

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
