
<!-- README.md is generated from README.Rmd. Please edit that file -->

# forecastbaselines

> R interface to
> [ForecastBaselines.jl](https://github.com/ManuelStapper/ForecastBaselines.jl)

An R package providing access to 10 baseline forecasting models from the
Julia ForecastBaselines.jl library, with comprehensive uncertainty
quantification and seamless integration with the
[scoringutils](https://epiforecasts.io/scoringutils/) evaluation
framework.

## Features

- **10 Forecasting Models**: From simple baselines (Constant, Marginal)
  to advanced time series models (ARMA, ETS, STL)
- **Probabilistic Forecasting**: Multiple methods for prediction
  intervals (empirical, parametric, model-based)
- **Comprehensive Scoring**: Compatible with scoringutils for all
  standard forecast evaluation metrics
- **Data Transformations**: Log, power, Box-Cox transformations with
  automatic back-transformation

## Installation

### Prerequisites

1.  **Julia** (\>= 1.9): Download from
    [julialang.org](https://julialang.org/downloads/)
2.  **R** (\>= 3.5.0)
3.  **JuliaCall R package**: `install.packages("JuliaCall")`

### Installing forecastbaselines

``` r
# Install from GitHub
remotes::install_github("epiforecasts/forecastbaselines")
```

### Setup

After installation, initialize Julia and load ForecastBaselines.jl:

``` r
library(forecastbaselines)

# Initialize Julia and install/load ForecastBaselines.jl
setup_ForecastBaselines()
```

This only needs to be done once per R session.

## Quick Example

``` r
# Your time series data
data <- c(1.2, 2.3, 3.1, 2.8, 3.5, 4.2, 3.9, 4.5, 4.1, 4.8)

# 1. Create and fit a model
model <- ARMAModel(p = 1, q = 1)
fitted <- fit_baseline(data, model)

# 2. Generate forecasts with prediction intervals
fc <- forecast(
  fitted,
  interval_method = EmpiricalInterval(n_trajectories = 1000),
  horizon = 1:5,
  levels = c(0.50, 0.95),
  model_name = "ARMA(1,1)"
)

# 3. Evaluate with true values
truth <- c(5.0, 5.2, 5.4, 5.1, 5.3)
fc_with_truth <- add_truth(fc, truth)

# 4. Score the forecast
fc_point <- as_forecast_point(fc_with_truth)
scores <- score(fc_point)
scores_summary <- summarise_scores(scores, by = "model")

scores_summary[, c("model", "ae_point", "se_point")]
#>        model  ae_point   se_point
#>       <char>     <num>      <num>
#> 1: ARMA(1,1) 0.2368795 0.07611188
```

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
- **ETSModel**: Error-Trend-Season exponential smoothing (all 30
  variants)

## Documentation

For detailed guides and examples:

- `vignette("forecastbaselines")` - Complete introduction with examples
- `vignette("forecast-models")` - Detailed guide to all 10 models
- `vignette("transformations")` - Working with data transformations

## Citation

If you use this package in your research, please cite the software and
the associated preprint:

**Software:**

    @software{forecastbaselinesr,
      title = {forecastbaselines: R Interface to ForecastBaselines.jl},
      author = {Stapper, Manuel and Funk, Sebastian},
      year = {2025},
      url = {https://github.com/epiforecasts/forecastbaselines}
    }

**Preprint:**

    @article{stapper2025baseline,
      title = {Mind the Baseline: The Hidden Impact of Reference Model Selection on Forecast Assessment},
      author = {Stapper, Manuel and Funk, Sebastian},
      year = {2025},
      doi = {10.1101/2025.08.01.25332807},
      url = {https://doi.org/10.1101/2025.08.01.25332807}
    }

## License

MIT License - see LICENSE file for details

## Support

For bugs and feature requests:

- R package issues: [GitHub
  Issues](https://github.com/epiforecasts/forecastbaselines/issues)
- Julia package issues: [ForecastBaselines.jl
  Issues](https://github.com/ManuelStapper/ForecastBaselines.jl/issues)

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

All contributions to this project are gratefully acknowledged using the
[`allcontributors` package](https://github.com/ropensci/allcontributors)
following the [allcontributors](https://allcontributors.org)
specification. Contributions of any kind are welcome!

<a href="https://github.com/epiforecasts/forecastbaselines/commits?author=sbfnk">sbfnk</a>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->
