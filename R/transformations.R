# Data Transformation Functions

#' No Transformation
#'
#' Creates a transformation that applies no changes to the data.
#'
#' @return A NoTransform object
#' @export
#'
#' @examples
#' \dontrun{
#' trans <- NoTransform()
#' }
NoTransform <- function() {
  check_setup()
  JuliaCall::julia_eval("ForecastBaselines.NoTransform()")
}

#' Log Transformation
#'
#' Creates a natural logarithm transformation. Data must be positive.
#'
#' @return A LogTransform object
#' @export
#'
#' @examples
#' \dontrun{
#' trans <- LogTransform()
#' }
LogTransform <- function() {
  check_setup()
  JuliaCall::julia_eval("ForecastBaselines.LogTransform()")
}

#' Log Plus One Transformation
#'
#' Creates a log(x + c) transformation. Useful for data with zeros.
#'
#' @param c Constant to add before taking log (default: 1)
#'
#' @return A LogPlusOneTransform object
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard log(x + 1)
#' trans <- LogPlusOneTransform()
#'
#' # Custom constant
#' trans <- LogPlusOneTransform(c = 0.5)
#' }
LogPlusOneTransform <- function(c = 1.0) {
  check_setup()
  JuliaCall::julia_call("ForecastBaselines.LogPlusOneTransform", as.numeric(c))
}

#' Square Root Transformation
#'
#' Creates a square root transformation. Useful for count data and variance
#' stabilization.
#'
#' @return A SquareRootTransform object
#' @export
#'
#' @examples
#' \dontrun{
#' trans <- SquareRootTransform()
#' }
SquareRootTransform <- function() {
  check_setup()
  JuliaCall::julia_eval("ForecastBaselines.SquareRootTransform()")
}

#' Power Transformation
#'
#' Creates a power transformation: x^lambda (Box-Cox family).
#'
#' @param lambda Power parameter (lambda = 0 is log, lambda = 0.5 is sqrt)
#'
#' @return A PowerTransform object
#' @export
#'
#' @examples
#' \dontrun{
#' # Square root (equivalent to SquareRootTransform)
#' trans <- PowerTransform(lambda = 0.5)
#'
#' # Cube root
#' trans <- PowerTransform(lambda = 1 / 3)
#' }
PowerTransform <- function(lambda) {
  check_setup()
  JuliaCall::julia_call("ForecastBaselines.PowerTransform", as.numeric(lambda))
}

#' Power Plus One Transformation
#'
#' Creates a power transformation with a constant: (x + c)^lambda.
#' Useful for Box-Cox transformations with zeros.
#'
#' @param lambda Power parameter
#' @param constant Constant to add before transformation (default: 1.0)
#'
#' @return A PowerPlusOneTransform object
#' @export
#'
#' @examples
#' \dontrun{
#' # Box-Cox transformation
#' trans <- PowerPlusOneTransform(lambda = 0.3)
#'
#' # Custom constant
#' trans <- PowerPlusOneTransform(lambda = 0.5, constant = 0.5)
#' }
PowerPlusOneTransform <- function(lambda, constant = 1.0) {
  check_setup()
  JuliaCall::julia_call("ForecastBaselines.PowerPlusOneTransform",
    as.numeric(lambda),
    constant = as.numeric(constant)
  )
}

#' Apply Data Transformation
#'
#' Applies a transformation to data.
#'
#' @param x Numeric vector of data
#' @param transformation A transformation object
#'
#' @return Transformed numeric vector
#' @export
#'
#' @details
#' **Recommended:** Use R's built-in transformation functions instead for
#' better reliability and flexibility. This function has limitations due to
#' bugs in ForecastBaselines.jl (e.g., SquareRootTransform fails).
#'
#' See `vignette("transformations")` for the recommended R approach.
#'
#' @examples
#' \dontrun{
#' # Julia approach (limited):
#' data <- c(1, 2, 3, 4, 5)
#' trans <- LogTransform()
#' transformed <- transform_data(data, trans)
#'
#' # Recommended R approach:
#' transformed <- log(data)
#' }
transform_data <- function(x, transformation) {
  check_setup()

  # Warn about SquareRootTransform bug
  if (inherits(transformation, "SquareRootTransform")) {
    warning(
      "SquareRootTransform has a bug in ForecastBaselines.jl. ",
      "Use sqrt() in R instead, or PowerTransform(lambda=0.5)"
    )
  }

  # Julia transform expects a matrix (n x 1)
  x_mat <- matrix(as.numeric(x), ncol = 1)
  JuliaCall::julia_assign("x_mat", x_mat)
  JuliaCall::julia_assign("trans_obj", transformation)
  result <- JuliaCall::julia_eval(
    "ForecastBaselines.transform(x_mat, trans_obj)"
  )
  as.numeric(result)
}

#' Apply Inverse Transformation
#'
#' Applies the inverse of a transformation to data (e.g., exp for log).
#'
#' @param y Numeric vector of transformed data
#' @param transformation A transformation object
#'
#' @return Original-scale numeric vector
#' @export
#'
#' @examples
#' \dontrun{
#' transformed <- c(0, 0.693, 1.099, 1.386, 1.609)
#' trans <- LogTransform()
#' original <- inverse_transform_data(transformed, trans)
#' }
inverse_transform_data <- function(y, transformation) {
  check_setup()
  # Julia inverse_transform expects a vector (unlike transform which needs
  # a matrix)
  JuliaCall::julia_assign("y_vec", as.numeric(y))
  JuliaCall::julia_assign("trans_obj", transformation)
  result <- JuliaCall::julia_eval(
    "ForecastBaselines.inverse_transform(y_vec, trans_obj)"
  )
  as.numeric(result)
}

#' Apply Transformation to Model
#'
#' Wraps a model with a data transformation.
#'
#' @param model A model object
#' @param transformation A transformation object
#'
#' @return A transformed model object
#' @export
#'
#' @details
#' **Not Implemented:** This function does not work because
#' ForecastBaselines.jl does not implement `transform()` for model types.
#'
#' **Recommended approach:** Transform your data manually in R before fitting:
#'
#' ```r
#' # Instead of transform_model():
#' log_data <- log(data)
#' fitted <- fit_baseline(log_data, model)
#' fc <- forecast(fitted, ...)
#' fc$mean <- exp(fc$mean)  # Back-transform
#' ```
#'
#' See `vignette("transformations")` for complete examples.
#'
#' @examples
#' \dontrun{
#' # This will fail - not implemented in Julia package
#' model <- ARMAModel(p = 1, q = 1)
#' trans <- LogTransform()
#' transformed_model <- transform_model(model, trans) # Error!
#'
#' # Use manual transformation instead:
#' log_data <- log(data)
#' fitted <- fit_baseline(log_data, model)
#' fc <- forecast(fitted, interval_method = NoInterval(), horizon = 1:5)
#' fc$mean <- exp(fc$mean)
#' }
transform_model <- function(model, transformation) {
  check_setup()
  stop("transform_model() is not implemented in ForecastBaselines.jl.\n",
    "Use manual transformation instead:\n",
    "  1. Transform data: log_data <- log(data)\n",
    "  2. Fit model: fitted <- fit_baseline(log_data, model)\n",
    "  3. Forecast: fc <- forecast(fitted, ...)\n",
    "  4. Back-transform: fc$mean <- exp(fc$mean)\n",
    "See vignette('transformations') for details.",
    call. = FALSE
  )
}
