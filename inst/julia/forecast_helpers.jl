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
        # intervals is a Dict{Float64, ForecastInterval}
        # keys are levels (e.g., 0.95), values are ForecastInterval objects
        # ForecastInterval has .lower, .upper, .levels fields (all Vector{Float64})

        levels = sort(collect(keys(fc.intervals)))
        n_horizons = length(fc.horizon)

        # Build quantile levels and quantile matrix
        quantile_levels = Float64[]
        quantiles_list = Vector{Float64}[]

        for level in levels
            interval_obj = fc.intervals[level]
            # interval_obj is a ForecastInterval with .lower and .upper vectors
            # lower quantile = (1 - level) / 2
            # upper quantile = 1 - (1 - level) / 2
            lower_q = (1.0 - level) / 2.0
            upper_q = 1.0 - lower_q

            push!(quantile_levels, lower_q)
            push!(quantile_levels, upper_q)
            push!(quantiles_list, Vector{Float64}(interval_obj.lower))  # lower bounds
            push!(quantiles_list, Vector{Float64}(interval_obj.upper))  # upper bounds
        end

        # Add median if present
        if fc.median !== nothing
            push!(quantile_levels, 0.5)
            push!(quantiles_list, to_r_compatible(fc.median))
        end

        # Sort by quantile level and arrange as matrix
        perm = sortperm(quantile_levels)
        quantile_levels = quantile_levels[perm]
        quantiles_list = quantiles_list[perm]

        # Convert to matrix: rows = horizons, cols = quantile levels
        result["quantiles"] = hcat(quantiles_list...)  # Shape: (n_horizons, n_quantile_levels)
        result["quantile_levels"] = quantile_levels
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
