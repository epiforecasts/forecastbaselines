# Helper functions for converting ForecastBaselines.jl objects to R-friendly formats
# These handle the problematic 0-dimensional arrays and complex types

"""
Convert a Julia value to an R-compatible format.
Handles Nothing, 0-d arrays, and ensures proper vector types.
"""
function to_r_compatible(x::Nothing)
    return nothing
end

function to_r_compatible(x::AbstractArray{T, 0}) where T
    # 0-dimensional array - extract scalar and wrap in 1-element vector
    return T[x[]]
end

function to_r_compatible(x::AbstractArray{T}) where T
    # Multi-dimensional array - ensure it's a proper Vector
    return Vector{T}(vec(x))
end

function to_r_compatible(x::AbstractString)
    return String(x)
end

function to_r_compatible(x::Number)
    return [x]  # Wrap scalars in vectors for consistency
end

function to_r_compatible(x)
    # Pass through anything else
    return x
end

"""
Convert a ForecastBaselines.Forecast object to a NamedTuple that
JuliaConnectoR's `juliaGet()` translates into a plain named R list.
"""
function forecast_to_r(fc)
    quantiles = nothing
    quantile_levels_sorted = nothing

    # Extract quantiles from intervals if present
    if fc.intervals !== nothing && !isempty(fc.intervals)
        # intervals is a Vector{ForecastInterval}, one per horizon
        # Each ForecastInterval has .lower, .upper, .levels (all Vector{Float64})
        # .lower[i] and .upper[i] correspond to .levels[i]

        n_horizons = length(fc.intervals)

        # Get the levels from the first interval (should be same for all)
        levels = fc.intervals[1].levels
        n_levels = length(levels)

        # Build quantile levels from confidence levels
        quantile_levels = Float64[]
        for level in levels
            # Convert confidence level to quantile levels
            # e.g., 0.95 confidence → 0.025 and 0.975 quantiles
            lower_q = (1.0 - level) / 2.0
            upper_q = 1.0 - lower_q
            push!(quantile_levels, lower_q)
            push!(quantile_levels, upper_q)
        end

        # Add median quantile if present
        has_median = fc.median !== nothing
        if has_median
            push!(quantile_levels, 0.5)
        end

        # Build quantile matrix: rows = horizons, cols = quantile levels
        n_quantiles = length(quantile_levels)
        quantiles_matrix = zeros(n_horizons, n_quantiles)

        for h in 1:n_horizons
            interval = fc.intervals[h]
            col_idx = 1
            # Fill lower and upper bounds for each level
            for i in 1:n_levels
                quantiles_matrix[h, col_idx] = interval.lower[i]
                quantiles_matrix[h, col_idx + 1] = interval.upper[i]
                col_idx += 2
            end
            # Fill median if present
            if has_median
                median_val = if isa(fc.median, AbstractArray)
                    fc.median[h]
                else
                    fc.median  # Scalar median
                end
                quantiles_matrix[h, col_idx] = Float64(median_val)
            end
        end

        # Sort columns by quantile level
        perm = sortperm(quantile_levels)
        quantiles = quantiles_matrix[:, perm]
        quantile_levels_sorted = quantile_levels[perm]
    end

    return (
        horizon = to_r_compatible(fc.horizon),
        mean = to_r_compatible(fc.mean),
        median = to_r_compatible(fc.median),
        truth = to_r_compatible(fc.truth),
        model_name = String(fc.model_name),
        quantiles = quantiles,
        quantile_levels = quantile_levels_sorted,
        trajectories = nothing,  # Skip trajectories for now
    )
end

"""
Convert interval forecast results to a NamedTuple for translation to R.
"""
function interval_result_to_r(res)
    # res is a tuple: (point, median, intervals, trajectories)
    return (
        point = to_r_compatible(res[1]),
        median = to_r_compatible(res[2]),
        intervals = res[3],
        trajectories = res[4],
    )
end

"""
Construct an `ARMAModel` with `k` seasonal harmonic waves and an optional
linear trend.

The keyword constructor in ForecastBaselines.jl hardcodes a single harmonic
pair, so this builds the mean function directly and passes it to the inner
constructor:

    μ(θ, t) = θ[1] + ∑_{j=1}^k [θ[2j] sin(2πjt/s) + θ[2j+1] cos(2πjt/s)]
              (+ θ[2k+2] t if trend)

With `k = 1` this reproduces the keyword constructor exactly.
"""
function arma_model(p::Int, q::Int, s::Int, k::Int, trend::Bool)
    if s <= 0 || k <= 0
        μ = trend ? (θ, t) -> θ[1] + θ[2] * t : (θ, t) -> θ[1]
        return ForecastBaselines.ARMAModel(p, q, μ, trend ? 2 : 1)
    end

    μDim = 1 + 2 * k + (trend ? 1 : 0)
    μ = function (θ, t)
        m = θ[1]
        for j in 1:k
            ω = 2 * π * j * t / s
            m += θ[2 * j] * sin(ω) + θ[2 * j + 1] * cos(ω)
        end
        if trend
            m += θ[2 * k + 2] * t
        end
        return m
    end
    return ForecastBaselines.ARMAModel(p, q, μ, μDim)
end
