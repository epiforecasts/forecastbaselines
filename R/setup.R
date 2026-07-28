# Package environment to track Julia setup state
.fbr_env <- new.env(parent = emptyenv())

NULL

#' Setup Julia and load ForecastBaselines.jl
#'
#' Initialises Julia, installs ForecastBaselines.jl if needed, and loads
#' the bridge code. Idempotent: subsequent calls are no-ops once setup
#' has succeeded. Heavy lifting is delegated to
#' [juliaready::julia_ready()] and [juliaready::julia_load_bridge()].
#'
#' Note: this is *not* called from `.onLoad`. Eager Julia initialisation
#' on package load can crash R during attach when interacting with other
#' compiled backends (e.g. Stan). Call this function explicitly, or rely
#' on the lazy [check_setup()] guard inside the package's exported
#' functions.
#'
#' @param install_package Whether to install ForecastBaselines.jl if not
#'   already installed.
#' @param verbose Whether to print progress messages.
#'
#' @return Invisibly returns `TRUE` if setup succeeded.
#' @export
#'
#' @examples
#' \dontrun{
#' setup_ForecastBaselines()
#' }
setup_ForecastBaselines <- function(install_package = TRUE,
                                    verbose = TRUE) {
  juliaready::julia_ready(
    packages  = "ForecastBaselines",
    state_env = .fbr_env,
    project   = system.file("julia", package = "forecastbaselines"),
    install   = FALSE,
    verbose   = verbose
  )
  juliaready::julia_load_bridge(
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
      juliaready::eval_julia("isdefined(Main, :ForecastBaselines)"),
      error = function(e) FALSE
    )
}

# Internal: lazy-init guard. Call from any function that needs Julia.
check_setup <- function() {
  juliaready::ensure_julia(.fbr_env, setup_ForecastBaselines)
}

# Internal: evaluate Julia code and fully translate the result into R.
# eval_julia() returns composite values (NamedTuples, Dicts, structs) as
# proxy environments; get_julia() converts them to plain R objects, so a
# NamedTuple arrives as a named list.
eval_julia_value <- function(code) {
  juliaready::get_julia(juliaready::eval_julia(code))
}
