# Statistical Modeling in R

This project implements statistical models using R, focusing on nonlinear models and state space models for time series data.

## Overview

The project analyzes the provided dataset (`WWdata.csv`) using various statistical modeling techniques:

1. **Nonlinear Models**:
   - Exponential model: Y ~ a * exp(b * X)
   - Power model: Y ~ a * X^b
   - Multi-predictor model: Y ~ a * X^b * exp(c * T) * W^d

2. **State Space Models**:
   - Basic local level model (DLM)
   - Structural time series model with regression components

## Requirements

The following R packages are required:
- tidyverse
- nlme
- nls2
- dlm
- forecast
- KFAS
- minpack.lm
- pacman (used to install other packages)

## Data

The dataset (`WWdata.csv`) contains time series data with the following variables:
- Year(t): Time period
- Y: Response variable
- X, T, W: Predictor variables

## Usage

1. Ensure R is installed on your system
2. Make sure the data file is in the `assets` directory
3. Run the script:

```bash
Rscript statistical_models.R
```

## Output

The script will:
1. Read and display summary statistics for the data
2. Create exploratory data visualizations in the `plots` directory
3. Fit nonlinear models and display results
4. Fit state space models and display results
5. Provide a framework for model comparison

## Files

- `statistical_models.R`: Main R script for statistical modeling
- `assets/WWdata.csv`: Data file
- `plots/`: Directory containing generated plots

## Model Details

### Nonlinear Models
The script fits several forms of nonlinear models to understand the relationship between the response variable (Y) and predictors (X, T, W). R-squared values are calculated to assess model fit.

### State Space Models
The state space models are implemented using both the dlm package (for basic local level models) and the KFAS package (for structural time series models with regression components). These models are useful for understanding the underlying structure of the time series and making forecasts.

## Model Comparison
Models are compared based on:
- Information criteria (AIC, BIC)
- Prediction accuracy (RMSE, MAE)
- Residual diagnostics
- Cross-validation results 