# forecastbaselines 0.1.0

## Breaking changes

* Model wrappers no longer accept arguments that ForecastBaselines.jl cannot
  honour. `KDEModel()` drops `bandwidth` and `kernel`, `LSDModel()` drops
  `trend_correction`, `IDSModel()` drops `threshold`, `STLModel()` drops
  `trend` and `robust`, `ARMAModel()` drops `include_mean`, and `ETSModel()`
  drops `damped` (damping is requested through `trend_type = "Ad"` or `"Md"`).
  These arguments were accepted and then discarded, so a model built with them
  was identical to one built without them (#7).

## New features

* `ARMAModel()` passes `include_drift` through to Julia, which adds a linear
  trend to the mean of the series (#7).
* `ARMAModel()` and `INARCHModel()` gain `k`, the number of harmonic waves
  representing seasonality. Values above 1 allow sharper, asymmetric seasonal
  shapes than a single sinusoid (#7).
* `INARCHModel()` gains `s` (seasonal period) and `nb` (negative binomial
  conditional distribution for overdispersed counts) (#7).
