# Installation Guide for forecastbaselines

This guide will walk you through installing forecastbaselines and its dependencies.

## Prerequisites

### 1. Install Julia

forecastbaselines requires Julia (version 1.9 or higher) to be installed on your system.

#### On Linux/macOS:

Download and install Julia from [julialang.org](https://julialang.org/downloads/):

```bash
# Example for Linux (adjust version as needed)
wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.0-linux-x86_64.tar.gz
tar xzf julia-1.10.0-linux-x86_64.tar.gz
sudo mv julia-1.10.0 /opt/
sudo ln -s /opt/julia-1.10.0/bin/julia /usr/local/bin/julia
```

Verify installation:
```bash
julia --version
```

#### On Windows:

1. Download the Windows installer from [julialang.org](https://julialang.org/downloads/)
2. Run the installer
3. Add Julia to your PATH if not done automatically

Verify installation in Command Prompt:
```cmd
julia --version
```

### 2. Install R

forecastbaselines requires R version 3.5.0 or higher.

Download and install R from [CRAN](https://cran.r-project.org/).

### 3. Install RStudio (Optional but Recommended)

Download and install RStudio from [posit.co](https://posit.co/downloads/).

## Installing forecastbaselines

### Step 1: Install JuliaCall

Open R or RStudio and run:

```r
install.packages("JuliaCall")
```

### Step 2: Install forecastbaselines

#### Option A: Install from local directory

```r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install from local directory
devtools::install_local("/path/to/forecastbaselines")
```

Replace `/path/to/forecastbaselines` with the actual path to the package directory.

#### Option B: Install from GitHub (when available)

```r
devtools::install_github("ManuelStapper/ForecastBaselines.jl", subdir = "forecastbaselines")
```

### Step 3: Load and Setup

```r
library(forecastbaselines)

# Initialize Julia and install/load ForecastBaselines.jl
setup_ForecastBaselines()
```

This command will:
1. Detect and initialize Julia
2. Install ForecastBaselines.jl if not already installed
3. Load the ForecastBaselines.jl package
4. Make all forecasting functions available

**Note**: The first run may take a few minutes to install ForecastBaselines.jl and its dependencies.

## Verifying Installation

Run a simple test to verify everything is working:

```r
library(forecastbaselines)
setup_ForecastBaselines()

# Create simple data
data <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# Fit a simple model
model <- ConstantModel()
fitted <- fit_baseline(data, model)

# Generate forecast
fc <- point_forecast(fitted, horizon = 1:3)
print(fc)
```

If you see forecast values printed, installation was successful!

## Troubleshooting

### Problem: Julia not found

**Solution**: Specify Julia location explicitly:

```r
setup_ForecastBaselines(JULIA_HOME = "/path/to/julia/bin")
```

On Linux/macOS, find Julia location with:
```bash
which julia
```

On Windows, it's typically in:
```
C:/Users/YourUsername/AppData/Local/Programs/Julia-1.x.x/bin
```

### Problem: ForecastBaselines.jl installation fails

**Solution**: Install manually in Julia, then load in R:

1. Open Julia terminal
2. Run:
```julia
using Pkg
Pkg.add(url="https://github.com/ManuelStapper/ForecastBaselines.jl")
```

3. Then in R:
```r
setup_ForecastBaselines(install_package = FALSE)
```

### Problem: JuliaCall errors on load

**Solution**: Rebuild Julia system image:

```r
setup_ForecastBaselines(rebuild = TRUE)
```

### Problem: Permission errors on Linux/macOS

**Solution**: Install Julia packages to local directory:

In Julia, run:
```julia
ENV["JULIA_DEPOT_PATH"] = expanduser("~/.julia")
using Pkg
Pkg.add(url="https://github.com/ManuelStapper/ForecastBaselines.jl")
```

### Problem: Slow first run

This is normal! Julia uses Just-In-Time (JIT) compilation, so the first run of each function will be slower. Subsequent runs will be much faster.

### Problem: Package conflicts

If you have conflicts with other Julia packages, create a new Julia environment:

In Julia:
```julia
using Pkg
Pkg.activate("forecastenv")
Pkg.add(url="https://github.com/ManuelStapper/ForecastBaselines.jl")
```

Then in R:
```r
# Point to the custom environment
Sys.setenv(JULIA_PROJECT = "~/forecastenv")
setup_ForecastBaselines()
```

## Testing Installation

Run the example scripts to test functionality:

```r
# Test basic usage
source("examples/basic_usage.R")

# Test seasonal forecasting
source("examples/seasonal_forecasting.R")
```

## Getting Help

If you encounter issues not covered here:

1. Check the [README](README.md) for usage examples
2. Open an issue on [GitHub](https://github.com/ManuelStapper/ForecastBaselines.jl/issues)
3. Check the [JuliaCall documentation](https://non-contradiction.github.io/JuliaCall/)

## Next Steps

Once installation is complete:

1. Read the [README](README.md) for an overview of features
2. Run the examples in the `examples/` directory
3. Check the function documentation with `?function_name` in R
4. Start forecasting!

## System Requirements

- **R**: >= 3.5.0
- **Julia**: >= 1.9
- **Operating Systems**: Linux, macOS, Windows
- **RAM**: At least 2GB recommended (4GB+ for large datasets)
- **Disk Space**: ~500MB for Julia and packages

## Optional Dependencies

For enhanced functionality:

```r
# For better plotting
install.packages("ggplot2")

# For data manipulation
install.packages("dplyr")

# For working with dates
install.packages("lubridate")
```
