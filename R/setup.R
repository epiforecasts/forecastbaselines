# Package environment to track Julia setup state
.fbr_env <- new.env(parent = emptyenv())

NULL

#' Setup Julia and load ForecastBaselines.jl
#'
#' Initialises Julia against the pinned Julia project the package ships in
#' `inst/julia`, and loads the bridge code. Idempotent: subsequent calls are
#' no-ops once setup has succeeded. Heavy lifting is delegated to
#' [juliabridge::julia_ready()] and [juliabridge::julia_load_bridge()].
#'
#' Note: this is *not* called from `.onLoad`. Eager Julia initialisation
#' on package load can crash R during attach when interacting with other
#' compiled backends (e.g. Stan). Call this function explicitly, or rely
#' on the lazy [check_setup()] guard inside the package's exported
#' functions.
#'
#' @param verbose Whether to print progress messages.
#'
#' @return Invisibly returns `TRUE` if setup succeeded.
#' @export
#'
#' @examples
#' \dontrun{
#' setup_ForecastBaselines()
#' }
setup_ForecastBaselines <- function(verbose = TRUE) {
  juliabridge::julia_ready(
    packages  = "ForecastBaselines",
    state_env = .fbr_env,
    project   = system.file("julia", package = "forecastbaselines"),
    install   = FALSE,
    verbose   = verbose
  )
  juliabridge::julia_load_bridge(
    package = "forecastbaselines",
    files   = "forecast_helpers.jl",
    verbose = verbose
  )
  if (verbose) message("forecastbaselines setup complete!")
  invisible(TRUE)
}

#' Check if Julia and ForecastBaselines.jl are set up
#'
#' Tests whether Julia is configured and ForecastBaselines.jl is accessible.
#'
#' @return TRUE if setup is complete, FALSE otherwise
#' @export
#'
#' @examples
#' \dontrun{
#' if (is_setup()) {
#'   # Use forecasting functions
#' } else {
#'   setup_ForecastBaselines()
#' }
#' }
is_setup <- function() {
  isTRUE(.fbr_env$ready) &&
    tryCatch(
      juliabridge::eval_julia("isdefined(Main, :ForecastBaselines)"),
      error = function(e) FALSE
    )
}

# Internal: lazy-init guard. Call from any function that needs Julia.
check_setup <- function() {
  juliabridge::ensure_julia(.fbr_env, setup_ForecastBaselines)
}

# Internal: JuliaConnectoR sends a length-one R vector as a Julia scalar. A
# one-dimensional array keeps it a Vector, which is what the Julia functions
# taking levels, horizons or values expect.
julia_vector <- function(x) {
  array(x, dim = length(x))
}
