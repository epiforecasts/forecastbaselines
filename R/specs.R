# Models, interval methods and transformations are held in R as descriptions
# of the Julia constructor call, built in Julia only when they are used.

# Internal: describe a ForecastBaselines.jl constructor call.
baseline_spec <- function(name, role, ...,
                          .fn = paste0("ForecastBaselines.", name)) {
  juliabridge::julia_spec(.fn, ..., .role = role)
}

#' Build the Julia object a model, interval method or transformation describes
#'
#' Models, interval methods and transformations are descriptions of a call to
#' a ForecastBaselines.jl constructor, built in Julia when they are used. This
#' builds one now, for inspecting the Julia object it stands for.
#'
#' @param x A model, interval method, transformation or [TemporalInfo()].
#'   Anything else, such as an object already in Julia, is returned as it is.
#'
#' @return A proxy for the Julia object, whose fields can be read with `$`.
#' @export
#'
#' @examples
#' \dontrun{
#' as_julia_object(LSDModel(s = 12))$w
#' }
as_julia_object <- function(x) {
  if (!inherits(x, "julia_spec")) {
    return(x)
  }
  check_setup()
  juliabridge::eval_julia(juliabridge::as_julia(x))
}
