# Package environment to store Julia connection
.fbr_env <- new.env(parent = emptyenv())

#' @importFrom JuliaCall julia_setup julia_eval julia_call julia_assign
.onLoad <- function(libname, pkgname) {
  # Try automatic setup
  tryCatch(
    {
      JuliaCall::julia_setup()
      JuliaCall::julia_eval("using ForecastBaselines")

      helper_file <- system.file(
        "julia", "forecast_helpers.jl",
        package = "forecastbaselines"
      )
      if (file.exists(helper_file)) {
        JuliaCall::julia_eval(
          sprintf('include("%s"); nothing', helper_file)
        )
      }

      .fbr_env$setup_ok <- TRUE
    },
    error = function(e) {
      .fbr_env$setup_ok <- FALSE
    }
  )
}

.onAttach <- function(libname, pkgname) {
  if (isTRUE(.fbr_env$setup_ok)) {
    packageStartupMessage(
      "forecastbaselines: Julia backend loaded successfully"
    )
  } else {
    packageStartupMessage(
      "forecastbaselines: R interface to ForecastBaselines.jl\n",
      "Julia setup incomplete. Run setup_ForecastBaselines() to configure."
    )
  }
}

#' Setup Julia and load ForecastBaselines.jl
#'
#' This function initializes Julia, installs ForecastBaselines.jl if needed,
#' and loads the package. Must be called before using any forecasting functions.
#'
#' @param JULIA_HOME Path to Julia installation (optional, will auto-detect
#'   if not provided)
#' @param install_package Whether to install ForecastBaselines.jl if not
#'   already installed
#' @param rebuild Whether to rebuild the Julia system image
#' @param verbose Whether to print verbose output during setup
#'
#' @return Invisibly returns TRUE if setup was successful
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic setup (auto-detect Julia)
#' setup_ForecastBaselines()
#'
#' # Specify Julia location
#' setup_ForecastBaselines(JULIA_HOME = "/usr/local/julia/bin")
#' }
setup_ForecastBaselines <- function(JULIA_HOME = NULL,
                                    install_package = TRUE,
                                    rebuild = FALSE,
                                    verbose = TRUE) {
  if (verbose) {
    message("Initializing Julia...")
  }

  # Setup Julia
  if (is.null(JULIA_HOME)) {
    .fbr_env$julia <- JuliaCall::julia_setup(rebuild = rebuild)
  } else {
    .fbr_env$julia <- JuliaCall::julia_setup(
      JULIA_HOME = JULIA_HOME,
      rebuild = rebuild
    )
  }

  if (verbose) {
    message("Julia initialized successfully")
  }

  # Install ForecastBaselines.jl if requested
  if (install_package) {
    if (verbose) {
      message("Checking ForecastBaselines.jl installation...")
    }

    # Try to load the package, install if it fails
    tryCatch(
      {
        JuliaCall::julia_eval("using ForecastBaselines")
        if (verbose) {
          message("ForecastBaselines.jl is already installed")
        }
      },
      error = function(e) {
        if (verbose) {
          message("Installing ForecastBaselines.jl...")
        }
        JuliaCall::julia_eval(
          paste0(
            'using Pkg; Pkg.add(url="',
            'https://github.com/ManuelStapper/ForecastBaselines.jl")'
          )
        )
        JuliaCall::julia_eval("using ForecastBaselines")
        if (verbose) {
          message("ForecastBaselines.jl installed successfully")
        }
      }
    )
  } else {
    # Just try to load
    JuliaCall::julia_eval("using ForecastBaselines")
  }

  # Load helper functions for R conversion
  if (verbose) {
    message("Loading R conversion helpers...")
  }
  helper_file <- system.file(
    "julia", "forecast_helpers.jl",
    package = "forecastbaselines"
  )
  if (file.exists(helper_file)) {
    JuliaCall::julia_eval(sprintf('include("%s"); nothing', helper_file))
  } else {
    warning(
      "Could not find forecast_helpers.jl - ",
      "some functions may not work correctly"
    )
  }

  if (verbose) {
    message("forecastbaselines setup complete!")
  }

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
  tryCatch(
    {
      # Check if Julia is available and ForecastBaselines is loaded
      JuliaCall::julia_eval("isdefined(Main, :ForecastBaselines)")
    },
    error = function(e) {
      # Julia not available at all
      FALSE
    }
  )
}

# Internal function to check setup and give helpful error message
check_setup <- function() {
  if (!is_setup()) {
    stop(
      "Julia and ForecastBaselines.jl are not set up. ",
      "Please run setup_ForecastBaselines() first.",
      call. = FALSE
    )
  }
}
