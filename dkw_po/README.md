# DKW Bounds for Treatment Effects

This folder contains organized code for computing DKW (Dvoretzky-Kiefer-Wolfowitz) bounds for treatment effect distributions.

## File Structure

- **generate.R**: Contains functions for generating data and computing initial DKW bounds on the ECDFs
  - `generate_true_sample()`: Simulates potential outcomes
  - `get_bounds()`: Computes initial DKW bounds for both treatment and control groups
  
- **bounds.R**: Contains functions for computing bounds on the treatment effect distribution
  - `ecdf_Y1_bound()`: Computes bounds for treatment group ECDF
  - `ecdf_Y0_bound()`: Computes bounds for control group ECDF
  - `marakov_bounds_t()`: Computes Makarov bounds for a specific quantile
  - `get_marakov_bounds_all_t()`: Computes bounds for all quantiles and plots results
  
- **main.R**: Main execution script that orchestrates the analysis

## Usage

To run the complete analysis:

```r
source("dkw_po/main.R")
```

Or run individual components:

```r
source("dkw_po/generate.R")
source("dkw_po/bounds.R")

# Generate data
result <- get_bounds()

# Compute bounds
bounds <- result$bounds
data <- result$data
obs_Y1 <- sort(result$obs_Y1)
obs_Y0 <- sort(result$obs_Y0)

diff_bounds <- get_marakov_bounds_all_t(bounds, obs_Y1, obs_Y0, data)
```

## Method

The code implements:
1. DKW bounds on empirical CDFs for treatment and control groups
2. Makarov bounds on the distribution of individual treatment effects
3. Comparison with the RIQITE Stephenson rank-based method



