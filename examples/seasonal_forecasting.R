# Seasonal Forecasting Example
# =============================

library(forecastbaselines)

# Initialize Julia and load ForecastBaselines.jl
setup_ForecastBaselines()

# ============================================================================
# Generate Seasonal Data
# ============================================================================

cat("\n=== Generating Seasonal Time Series ===\n\n")

set.seed(123)
n_years <- 4
months_per_year <- 12
n <- n_years * months_per_year

# Time index
time <- 1:n

# Components
trend <- time * 0.2
seasonal <- 3 * sin(2 * pi * time / months_per_year)
noise <- rnorm(n, 0, 0.5)

# Combined series
data <- 20 + trend + seasonal + noise

cat("Generated", n, "months of data with seasonal pattern\n")

# Plot the data
plot(data, type = "l", main = "Seasonal Time Series",
     xlab = "Month", ylab = "Value", col = "blue", lwd = 2)

# ============================================================================
# Method 1: STL Decomposition
# ============================================================================

cat("\n=== Method 1: STL Model ===\n\n")

# Create and fit STL model
model_stl <- STLModel(s = months_per_year, robust = TRUE)
fitted_stl <- fit_baseline(data, model_stl)

# Forecast next 12 months
fc_stl <- forecast(
  fitted_stl,
  interval_method = EmpiricalInterval(n_trajectories = 2000, seed = 42),
  horizon = 1:12,
  levels = c(0.80, 0.95),
  model_name = "STL(12)"
)

cat("STL Forecasts (next 12 months):\n")
print(round(fc_stl$mean, 2))

# ============================================================================
# Method 2: Last Similar Dates (LSD)
# ============================================================================

cat("\n=== Method 2: Last Similar Dates ===\n\n")

# LSD model with monthly seasonality
model_lsd <- LSDModel(s = months_per_year, window_width = 2)
fitted_lsd <- fit_baseline(data, model_lsd)

fc_lsd <- forecast(
  fitted_lsd,
  interval_method = EmpiricalInterval(n_trajectories = 2000),
  horizon = 1:12,
  levels = c(0.80, 0.95),
  model_name = "LSD(12)"
)

cat("LSD Forecasts (next 12 months):\n")
print(round(fc_lsd$mean, 2))

# ============================================================================
# Method 3: Seasonal ARMA
# ============================================================================

cat("\n=== Method 3: Seasonal ARMA ===\n\n")

# ARMA with seasonal component
model_sarma <- ARMAModel(p = 1, q = 1, s = months_per_year)
fitted_sarma <- fit_baseline(data, model_sarma)

fc_sarma <- forecast(
  fitted_sarma,
  interval_method = EmpiricalInterval(n_trajectories = 2000),
  horizon = 1:12,
  levels = c(0.80, 0.95),
  model_name = "ARMA(1,1,12)"
)

cat("Seasonal ARMA Forecasts (next 12 months):\n")
print(round(fc_sarma$mean, 2))

# ============================================================================
# Method 4: ETS Model (Holt-Winters)
# ============================================================================

cat("\n=== Method 4: ETS (Holt-Winters) ===\n\n")

# ETS model with additive seasonality
model_ets <- ETSModel(
  error_type = "A",
  trend_type = "A",
  season_type = "A",
  s = months_per_year
)
fitted_ets <- fit_baseline(data, model_ets)

fc_ets <- forecast(
  fitted_ets,
  interval_method = ParametricInterval(),
  horizon = 1:12,
  levels = c(0.80, 0.95),
  model_name = "ETS(A,A,A)"
)

cat("ETS Forecasts (next 12 months):\n")
print(round(fc_ets$mean, 2))

# ============================================================================
# Model Comparison
# ============================================================================

cat("\n=== Model Comparison ===\n\n")

# Generate true future values
set.seed(789)
future_time <- (n + 1):(n + 12)
future_trend <- future_time * 0.2
future_seasonal <- 3 * sin(2 * pi * future_time / months_per_year)
future_noise <- rnorm(12, 0, 0.5)
truth <- 20 + future_trend + future_seasonal + future_noise

cat("True future values:\n")
print(round(truth, 2))

# Add truth to forecasts
fc_stl <- add_truth(fc_stl, truth)
fc_lsd <- add_truth(fc_lsd, truth)
fc_sarma <- add_truth(fc_sarma, truth)
fc_ets <- add_truth(fc_ets, truth)

# Compare models
forecasts <- list(
  STL = fc_stl,
  LSD = fc_lsd,
  SARMA = fc_sarma,
  ETS = fc_ets
)

# Score all forecasts
all_scores <- lapply(forecasts, function(fc) {
  fc_point <- scoringutils::as_forecast_point(fc)
  scores <- scoringutils::score(fc_point)
  scoringutils::summarise_scores(scores, by = "model")
})

comparison <- data.frame(
  Model = names(forecasts),
  MAE = vapply(all_scores, function(s) s$ae_point, numeric(1)),
  RMSE = vapply(all_scores, function(s) sqrt(s$se_point), numeric(1))
)

cat("\nModel Performance:\n")
print(comparison)

best_mae <- comparison$Model[which.min(comparison$MAE)]
best_crps <- comparison$Model[which.min(comparison$CRPS)]

cat("\nBest model (MAE):", best_mae, "\n")
cat("Best model (CRPS):", best_crps, "\n")

# ============================================================================
# Visualization
# ============================================================================

cat("\n=== Creating Forecast Visualization ===\n\n")

# Create comprehensive plot
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Function to plot forecast
plot_forecast <- function(data, forecast, truth, title) {
  n_hist <- length(data)
  n_fc <- length(forecast$mean)
  fc_time <- (n_hist + 1):(n_hist + n_fc)

  # Set up plot
  plot(1:n_hist, data, type = "l",
       xlim = c(max(1, n_hist - 24), n_hist + n_fc),
       ylim = range(c(tail(data, 24), forecast$mean, truth)),
       main = title,
       xlab = "Month", ylab = "Value",
       col = "blue", lwd = 2)

  # Add forecast
  lines(fc_time, forecast$mean, col = "red", lwd = 2)

  # Add truth
  lines(fc_time, truth, col = "darkgreen", lwd = 2, lty = 2)

  # Add vertical line
  abline(v = n_hist, lty = 2, col = "gray")

  # Legend
  legend("topleft",
         legend = c("Historical", "Forecast", "Truth"),
         col = c("blue", "red", "darkgreen"),
         lwd = 2, lty = c(1, 1, 2),
         cex = 0.7)
}

# Plot each model
plot_forecast(data, fc_stl, truth, "STL Model")
plot_forecast(data, fc_lsd, truth, "LSD Model")
plot_forecast(data, fc_sarma, truth, "SARMA Model")
plot_forecast(data, fc_ets, truth, "ETS Model")

par(mfrow = c(1, 1))

cat("\n=== Seasonal Forecasting Example Complete ===\n")
