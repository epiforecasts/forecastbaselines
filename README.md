# forecastbaselines

> R interface to [ForecastBaselines.jl](https://github.com/ManuelStapper/ForecastBaselines.jl)

An R package providing access to 10 baseline forecasting models from the Julia ForecastBaselines.jl library, with comprehensive uncertainty quantification and seamless integration with the scoringutils evaluation framework.

## Features

- **10 Forecasting Models**: From simple baselines (Constant, Marginal) to advanced time series models (ARMA, ETS, STL)
- **Probabilistic Forecasting**: Multiple methods for prediction intervals (empirical, parametric, model-based)
- **Comprehensive Scoring**: Compatible with [scoringutils](https://epiforecasts.io/scoringutils/) - access to all standard forecast evaluation metrics
- **Data Transformations**: Log, power, Box-Cox transformations with automatic back-transformation

## Available Models

### Simple Baseline Models
- **ConstantModel**: Naive forecast using last observed value
- **MarginalModel**: Forecasts based on empirical marginal distribution
- **KDEModel**: Kernel density estimation

### Seasonal/Trend Models
- **LSDModel**: Last Similar Dates method (seasonal patterns)
- **OLSModel**: Ordinary least squares with polynomial trends
- **IDSModel**: Increase-Decrease-Stable trend detection
- **STLModel**: Seasonal-Trend decomposition using Loess

### Advanced Time Series Models
- **ARMAModel**: Autoregressive Moving Average
- **INARCHModel**: Integer-valued ARCH for count data
- **ETSModel**: Error-Trend-Season exponential smoothing (all 30 variants)

## Installation

### Prerequisites

1. **Julia** (>= 1.9): Download from [julialang.org](https://julialang.org/downloads/)

2. **R** (>= 3.5.0)

3. **JuliaCall R package**:
   ```r
   install.packages("JuliaCall")
   ```

### Installing forecastbaselines

```r
# Install from local directory
devtools::install_local("/path/to/forecastbaselines")

# Or if using GitHub (when published)
# devtools::install_github("ManuelStapper/ForecastBaselines.jl", subdir = "forecastbaselines")
```

### Setup

After installation, initialize Julia and load ForecastBaselines.jl:

```r
library(forecastbaselines)

# Initialize Julia and install/load ForecastBaselines.jl
setup_ForecastBaselines()
```

This only needs to be done once per R session. The function will:
- Initialize Julia
- Install ForecastBaselines.jl if not already installed
- Load the package

## Quick Start

### Basic Workflow

```r
library(forecastbaselines)
setup_ForecastBaselines()

# Your time series data
data <- c(1.2, 2.3, 3.1, 2.8, 3.5, 4.2, 3.9, 4.5, 4.1, 4.8)

# 1. Create a model
model <- ARMAModel(p = 1, q = 1)

# 2. Fit the model
fitted <- fit_baseline(data, model)

# 3. Generate forecasts with quantiles for probabilistic evaluation
fc <- forecast(
  fitted,
  interval_method = EmpiricalInterval(n_trajectories = 1000),
  horizon = 1:5,
  levels = c(0.50, 0.80, 0.95),
  model_name = "ARMA(1,1)"
)

# 4. Evaluate with true values
truth <- c(5.0, 5.2, 5.4, 5.1, 5.3)
fc_with_truth <- add_truth(fc, truth)

# 5. Score the probabilistic forecast
fc_quantile <- scoringutils::as_forecast_quantile(fc_with_truth)
scores <- scoringutils::score(fc_quantile)
scores_summary <- scoringutils::summarise_scores(scores, by = "model")

# Access probabilistic metrics
cat("WIS:", scores_summary$wis, "\n")
cat("Interval Coverage 50%:", scores_summary$interval_coverage_50, "\n")
cat("Interval Coverage 95%:", scores_summary$interval_coverage_95, "\n")
print(scores_summary)
```

## Examples

### Example 1: Seasonal Data with STL

```r
# Simulated monthly data with seasonality
set.seed(123)
months <- 1:48
seasonal <- sin(2 * pi * months / 12) * 2
trend <- months * 0.1
noise <- rnorm(48, 0, 0.3)
data <- 10 + trend + seasonal + noise

# Fit STL model with monthly seasonality
model <- STLModel(s = 12, robust = TRUE)
fitted <- fit_baseline(data, model)

# Forecast next 12 months
fc <- forecast(
  fitted,
  interval_method = EmpiricalInterval(n_trajectories = 2000, seed = 42),
  horizon = 1:12,
  levels = c(0.50, 0.80, 0.95),
  model_name = "STL(12)"
)

print(fc)
```

### Example 2: Comparing Multiple Models

```r
# Create multiple models
models <- list(
  Naive = ConstantModel(),
  ARMA = ARMAModel(p = 2, q = 1),
  ETS = ETSModel(error_type = "A", trend_type = "A", season_type = "N")
)

# Fit and forecast with each model
results <- lapply(names(models), function(name) {
  fitted <- fit_baseline(data, models[[name]])
  fc <- forecast(
    fitted,
    interval_method = EmpiricalInterval(),
    horizon = 1:12,
    truth = truth_values,
    model_name = name
  )

  fc_point <- scoringutils::as_forecast_point(fc)
  fc_scores <- scoringutils::score(fc_point)
  fc_summary <- scoringutils::summarise_scores(fc_scores, by = "model")

  list(
    model = name,
    mae = fc_summary$ae_point,
    rmse = sqrt(fc_summary$se_point)
  )
})

# Print comparison
do.call(rbind, lapply(results, as.data.frame))
```

### Example 3: Data Transformation

```r
# Fit model to log-transformed data
model <- ARMAModel(p = 1)
trans <- LogTransform()
transformed_model <- transform_model(model, trans)

fitted <- fit_baseline(data, transformed_model)

# Forecasts are automatically back-transformed to original scale
fc <- forecast(
  fitted,
  interval_method = ParametricInterval(),
  horizon = 1:10
)
```

### Example 4: Count Data with INARCH

```r
# Count data (e.g., number of events per day)
count_data <- c(5, 3, 7, 4, 6, 8, 5, 7, 9, 6, 4, 5)

# INARCH model for count time series
model <- INARCHModel(p = 1)
fitted <- fit_baseline(count_data, model)

# Forecast with positivity correction
fc <- forecast(
  fitted,
  interval_method = EmpiricalInterval(
    n_trajectories = 1000,
    positivity_correction = "post_clip",
    return_trajectories = TRUE
  ),
  horizon = 1:7
)
```

### Example 5: Calibration Assessment

```r
# Generate forecast with intervals
fc <- forecast(
  fitted,
  interval_method = EmpiricalInterval(n_trajectories = 2000),
  horizon = 1:20,
  levels = seq(0.1, 0.9, by = 0.1),
  truth = observed_values
)

# Check calibration using PIT
pit_values <- PIT_function(fc)
hist(pit_values, breaks = 20, main = "PIT Histogram",
     xlab = "PIT Value", col = "lightblue")
abline(h = length(pit_values) / 20, col = "red", lty = 2)

# Compute calibration divergence
cvm <- CvM_divergence(fc)
cat("Cramér-von Mises divergence:", cvm, "\n")
cat("(0 = perfect calibration)\n")
```

## Interval Methods

### NoInterval
Point forecasts only:
```r
interval_method = NoInterval()
```

### EmpiricalInterval
Bootstrap from historical forecast errors:
```r
interval_method = EmpiricalInterval(
  n_trajectories = 1000,
  seed = 123,
  positivity_correction = "post_clip",
  return_trajectories = TRUE
)
```

### ParametricInterval
Model-based parametric intervals (assumes normality):
```r
interval_method = ParametricInterval()
```

### ModelTrajectoryInterval
Simulate from the fitted model:
```r
interval_method = ModelTrajectoryInterval(
  n_trajectories = 1000,
  seed = 456,
  return_trajectories = TRUE
)
```

## Scoring Rules

**Scoring powered by [scoringutils](https://epiforecasts.io/scoringutils/)**

Use `score(forecast)` to compute all metrics, then access them from the result:

```r
# Get all scores
scores <- score(forecast)

# Access specific metrics
mae <- scores$ae_point        # Mean Absolute Error
mse <- scores$se_point        # Mean Squared Error
rmse <- sqrt(scores$se_point) # Root Mean Squared Error
ape <- scores$ape             # Absolute Percentage Error

# For quantile forecasts
wis <- scores$wis             # Weighted Interval Score
bias <- scores$bias           # Forecast bias
```

**Additional functions:**
- `score(forecast, summarise = FALSE)` - Get per-horizon scores
- `get_available_metrics("point")` - List available point forecast metrics
- `get_available_metrics("quantile")` - List available quantile forecast metrics

## Package Structure

```
forecastbaselines/
├── DESCRIPTION         # Package metadata
├── NAMESPACE          # Exported functions
├── README.md          # This file
├── R/
│   ├── setup.R        # Package initialization and Julia setup
│   ├── models.R       # Model constructors
│   ├── forecast.R     # Forecasting functions
│   ├── intervals.R    # Interval methods
│   ├── scoring.R      # Scoring rules (via scoringutils)
│   ├── utils.R        # Utility functions
│   └── transformations.R  # Data transformations
└── examples/          # Example scripts
```

## Troubleshooting

### Julia not found
If Julia is not automatically detected:
```r
setup_ForecastBaselines(JULIA_HOME = "/path/to/julia/bin")
```

### Package not loading
Try rebuilding the Julia system image:
```r
setup_ForecastBaselines(rebuild = TRUE)
```

### Manual Julia package installation
If automatic installation fails, install ForecastBaselines.jl manually in Julia:
```julia
using Pkg
Pkg.add(url="https://github.com/ManuelStapper/ForecastBaselines.jl")
```

Then in R:
```r
setup_ForecastBaselines(install_package = FALSE)
```

## Citation

If you use this package in your research, please cite the software and the associated preprint:

**Software:**
```
@software{forecastbaselinesr,
  title = {forecastbaselines: R Interface to ForecastBaselines.jl},
  author = {Stapper, Manuel and Funk, Sebastian},
  year = {2025},
  url = {https://github.com/epiforecasts/forecastbaselines}
}
```

**Preprint:**
```
@article{stapper2025baseline,
  title = {Mind the Baseline: The Hidden Impact of Reference Model Selection on Forecast Assessment},
  author = {Stapper, Manuel and Funk, Sebastian},
  year = {2025},
  doi = {10.1101/2025.08.01.25332807},
  url = {https://doi.org/10.1101/2025.08.01.25332807}
}
```

## License

MIT License - see LICENSE file for details

## Contributing

Contributions are welcome! Please open an issue or pull request on GitHub.

## Support

For bugs and feature requests:
- R package issues: [GitHub Issues](https://github.com/epiforecasts/forecastbaselines/issues)
- Julia package issues: [GitHub Issues](https://github.com/ManuelStapper/ForecastBaselines.jl/issues)
