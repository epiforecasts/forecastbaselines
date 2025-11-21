# Basic Usage Example for forecastbaselines
# ==========================================

library(forecastbaselines)

# Initialize Julia and load ForecastBaselines.jl
setup_ForecastBaselines()

# ============================================================================
# Example 1: Simple ARMA Model
# ============================================================================

cat("\n=== Example 1: Simple ARMA Model ===\n\n")

# Create sample data
set.seed(123)
n <- 100
data <- cumsum(rnorm(n, 0, 1)) + 10

# Plot the data
plot(data, type = "l", main = "Time Series Data",
     xlab = "Time", ylab = "Value", col = "blue")

# Create and fit ARMA(1,1) model
cat("Fitting ARMA(1,1) model...\n")
model <- ARMAModel(p = 1, q = 1)
fitted <- fit_baseline(data, model)

# Generate point forecasts
cat("Generating point forecasts...\n")
horizon <- 1:12
point_fc <- point_forecast(fitted, horizon)

cat("Point forecasts:\n")
print(round(point_fc, 2))

# ============================================================================
# Example 2: Forecast with Prediction Intervals
# ============================================================================

cat("\n=== Example 2: Forecast with Prediction Intervals ===\n\n")

# Generate forecasts with 80% and 95% prediction intervals
cat("Generating forecasts with prediction intervals...\n")
fc <- forecast(
  fitted,
  interval_method = EmpiricalInterval(n_trajectories = 1000, seed = 42),
  horizon = 1:12,
  levels = c(0.80, 0.95),
  model_name = "ARMA(1,1)"
)

# Display forecast summary
print(fc)

# ============================================================================
# Example 3: Model Evaluation
# ============================================================================

cat("\n=== Example 3: Model Evaluation ===\n\n")

# Generate some "true" future values for evaluation
set.seed(456)
truth <- cumsum(rnorm(12, 0, 1)) + tail(data, 1)

cat("True values:\n")
print(round(truth, 2))

# Add truth to forecast
fc_with_truth <- add_truth(fc, truth)

# Convert to scoringutils format and score
fc_point <- scoringutils::as_forecast_point(fc_with_truth)
all_scores <- scoringutils::score(fc_point)
scores_summary <- scoringutils::summarise_scores(all_scores, by = "model")

cat("\nForecast scores:\n")
print(scores_summary)

# Access specific scores from the result
cat("\nMAE (ae_point):", round(scores_summary$ae_point, 3), "\n")

# ============================================================================
# Example 4: Different Models Comparison
# ============================================================================

cat("\n=== Example 4: Different Models Comparison ===\n\n")

# Create a list of models to compare
models <- list(
  "Constant" = ConstantModel(),
  "Marginal(p=10)" = MarginalModel(p = 10),
  "ARMA(1,1)" = ARMAModel(p = 1, q = 1),
  "ARMA(2,1)" = ARMAModel(p = 2, q = 1)
)

# Fit each model and compute MAE
results <- data.frame(
  Model = character(),
  MAE = numeric(),
  RMSE = numeric(),
  CRPS = numeric(),
  stringsAsFactors = FALSE
)

for (model_name in names(models)) {
  cat("Fitting", model_name, "...\n")

  fitted <- fit_baseline(data, models[[model_name]])
  fc <- forecast(
    fitted,
    interval_method = EmpiricalInterval(n_trajectories = 500),
    horizon = 1:12,
    truth = truth,
    model_name = model_name
  )

  # Score and extract metrics
  fc_point <- scoringutils::as_forecast_point(fc)
  fc_scores <- scoringutils::score(fc_point)
  fc_summary <- scoringutils::summarise_scores(fc_scores, by = "model")

  results <- rbind(results, data.frame(
    Model = model_name,
    MAE = fc_summary$ae_point,
    RMSE = sqrt(fc_summary$se_point)
  ))
}

# Display results
cat("\nModel Comparison:\n")
print(results)

# Find best model
best_model <- results$Model[which.min(results$MAE)]
cat("\nBest model (by MAE):", best_model, "\n")

# ============================================================================
# Example 5: Visualization
# ============================================================================

cat("\n=== Example 5: Forecast Visualization ===\n\n")

# Create a forecast with intervals
fitted_best <- fit_baseline(data, ARMAModel(p = 2, q = 1))
fc_vis <- forecast(
  fitted_best,
  interval_method = EmpiricalInterval(n_trajectories = 2000, seed = 42),
  horizon = 1:20,
  levels = c(0.50, 0.80, 0.95),
  model_name = "ARMA(2,1)"
)

# Plot historical data and forecasts
n_hist <- length(data)
plot(1:n_hist, data, type = "l", xlim = c(1, n_hist + 20),
     ylim = range(c(data, fc_vis$mean)),
     main = "Forecast with Prediction Intervals",
     xlab = "Time", ylab = "Value", col = "blue", lwd = 2)

# Add forecast
fc_time <- (n_hist + 1):(n_hist + length(fc_vis$mean))
lines(fc_time, fc_vis$mean, col = "red", lwd = 2)

# Add a vertical line at forecast origin
abline(v = n_hist, lty = 2, col = "gray")

legend("topleft",
       legend = c("Historical", "Forecast"),
       col = c("blue", "red"),
       lwd = 2)

cat("\nPlot created successfully!\n")

cat("\n=== Examples Complete ===\n")
