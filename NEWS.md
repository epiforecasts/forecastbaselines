# forecastbaselines 0.1.0

## Breaking changes

* Model wrappers no longer accept arguments they cannot pass on to
  ForecastBaselines.jl. `LSDModel()` drops `trend_correction`, `IDSModel()`
  drops `threshold` and `STLModel()` drops `trend`, none of which exists in
  Julia; `KDEModel()` drops `bandwidth` and `kernel` and `STLModel()` drops
  `robust`, which are estimation settings that `fit_baseline()` does not
  expose; `ARMAModel()` drops `include_mean`, as the intercept is always
  estimated; and `ETSModel()` drops `damped`, since damping is requested
  through `trend_type = "Ad"` or `"Md"`. These arguments were accepted and
  then discarded, so a model built with them was identical to one built
  without them (#7).

## New features

* `ARMAModel()` passes `include_drift` through to Julia, which adds a linear
  trend to the mean of the series (#7).
* `ARMAModel()` and `INARCHModel()` gain `k`, the number of harmonic waves
  representing seasonality. Values above 1 allow sharper, asymmetric seasonal
  shapes than a single sinusoid, subject to the limits documented on each
  function: ARMA fits with several waves can settle far from the data, and
  INARCH forecasts do not carry the estimated seasonality, which
  `INARCHModel()` warns about (#7).
* `INARCHModel()` gains `s` (seasonal period) and `nb` (negative binomial
  conditional distribution for overdispersed counts) (#7).
