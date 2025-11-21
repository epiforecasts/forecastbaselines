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

    # Skip complex types (intervals, trajectories) for now
    # These would need recursive conversion if needed
    result["intervals"] = nothing
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
