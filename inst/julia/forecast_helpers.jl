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
Convert a ForecastBaselines.Forecast object to an R-compatible Dict.
"""
function forecast_to_r_dict(fc)
    result = Dict{String, Any}()

    # Convert each field, handling Nothing and 0-d arrays
    result["horizon"] = to_r_compatible(fc.horizon)
    result["mean"] = to_r_compatible(fc.mean)
    result["median"] = to_r_compatible(fc.median)
    result["truth"] = to_r_compatible(fc.truth)
    result["model_name"] = String(fc.model_name)

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
        result["quantiles"] = quantiles_matrix[:, perm]
        result["quantile_levels"] = quantile_levels[perm]
    else
        result["quantiles"] = nothing
        result["quantile_levels"] = nothing
    end

    # Skip trajectories for now
    result["trajectories"] = nothing

    return result
end

"""
Convert interval forecast results to R-compatible format.
"""
function interval_result_to_r_dict(res)
    result = Dict{String, Any}()

    # res is a tuple: (point, median, intervals, trajectories)
    result["point"] = to_r_compatible(res[1])
    result["median"] = to_r_compatible(res[2])
    result["intervals"] = res[3]  # Skip for now
    result["trajectories"] = res[4]  # Skip for now

    return result
end
