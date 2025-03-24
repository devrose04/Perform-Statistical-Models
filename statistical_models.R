#!/usr/bin/env Rscript

# Load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,    # Data manipulation and visualization
  nlme,         # Nonlinear mixed-effects models
  nls2,         # Nonlinear regression
  dlm,          # Dynamic linear models (state space)
  forecast,     # Time series forecasting
  KFAS,         # Kalman filter and smoother
  minpack.lm    # More robust nonlinear least squares
)

# Function to read and prepare data
read_data <- function() {
  # Read the CSV file
  data <- read.csv("assets/WWdata.csv")
  
  # Display structure of the data
  cat("Data structure:\n")
  print(str(data))
  
  # Display summary statistics
  cat("\nSummary statistics:\n")
  print(summary(data))
  
  # Return the cleaned data
  return(data)
}

# Function for exploratory data analysis
explore_data <- function(data) {
  # Create a directory for plots if it doesn't exist
  if (!dir.exists("plots")) dir.create("plots")
  
  # Time series plots
  png("plots/time_series.png", width = 800, height = 600)
  par(mfrow = c(2, 2))
  plot(data$`Year.t.`, data$Y, type = "l", col = "blue", 
       main = "Y over Time", xlab = "Year", ylab = "Y")
  plot(data$`Year.t.`, data$X, type = "l", col = "red", 
       main = "X over Time", xlab = "Year", ylab = "X")
  plot(data$`Year.t.`, data$T, type = "l", col = "green", 
       main = "T over Time", xlab = "Year", ylab = "T")
  plot(data$`Year.t.`, data$W, type = "l", col = "purple", 
       main = "W over Time", xlab = "Year", ylab = "W")
  dev.off()
  
  # Scatter plots to explore relationships
  png("plots/scatter_plots.png", width = 800, height = 600)
  par(mfrow = c(2, 2))
  plot(data$X, data$Y, main = "Y vs X", xlab = "X", ylab = "Y")
  plot(data$T, data$Y, main = "Y vs T", xlab = "T", ylab = "Y")
  plot(data$W, data$Y, main = "Y vs W", xlab = "W", ylab = "Y")
  plot(data$X, data$W, main = "W vs X", xlab = "X", ylab = "W")
  dev.off()
  
  cat("Exploratory plots saved to 'plots' directory.\n")
}

# Function to fit nonlinear models
fit_nonlinear_models <- function(data) {
  cat("\n--- Fitting Nonlinear Models ---\n")
  
  # 1. Exponential model: Y ~ a * exp(b * X)
  try({
    exp_model <- nls(Y ~ a * exp(b * X), 
                    data = data,
                    start = list(a = 1, b = 0.01),
                    control = nls.control(maxiter = 100),
                    algorithm = "port")
    
    cat("\nExponential Model (Y ~ a * exp(b * X)):\n")
    print(summary(exp_model))
    
    # Calculate R-squared
    y_pred <- predict(exp_model)
    sse <- sum((data$Y - y_pred)^2)
    sst <- sum((data$Y - mean(data$Y))^2)
    rsq <- 1 - sse/sst
    cat("\nR-squared:", rsq, "\n")
  }, silent = TRUE)
  
  # 2. Power model: Y ~ a * X^b
  try({
    power_model <- nls(Y ~ a * X^b, 
                      data = data,
                      start = list(a = 1000, b = -0.5),
                      control = nls.control(maxiter = 100),
                      algorithm = "port")
    
    cat("\nPower Model (Y ~ a * X^b):\n")
    print(summary(power_model))
    
    # Calculate R-squared
    y_pred <- predict(power_model)
    sse <- sum((data$Y - y_pred)^2)
    sst <- sum((data$Y - mean(data$Y))^2)
    rsq <- 1 - sse/sst
    cat("\nR-squared:", rsq, "\n")
  }, silent = TRUE)
  
  # 3. Multi-predictor model: Y ~ a * X^b * exp(c * T) * W^d
  try({
    multi_model <- nlsLM(Y ~ a * X^b * exp(c * T) * W^d, 
                       data = data,
                       start = list(a = 1000, b = -0.5, c = 0.1, d = 1),
                       control = nls.control(maxiter = 200))
    
    cat("\nMulti-predictor Model (Y ~ a * X^b * exp(c * T) * W^d):\n")
    print(summary(multi_model))
    
    # Calculate R-squared
    y_pred <- predict(multi_model)
    sse <- sum((data$Y - y_pred)^2)
    sst <- sum((data$Y - mean(data$Y))^2)
    rsq <- 1 - sse/sst
    cat("\nR-squared:", rsq, "\n")
    
    # Create diagnostic plots
    png("plots/nonlinear_model_diagnostics.png", width = 800, height = 600)
    par(mfrow = c(2, 2))
    plot(y_pred, data$Y, main = "Predicted vs Actual",
         xlab = "Predicted Y", ylab = "Actual Y")
    abline(0, 1, col = "red")
    
    residuals <- data$Y - y_pred
    plot(y_pred, residuals, main = "Residuals vs Predicted",
         xlab = "Predicted Y", ylab = "Residuals")
    abline(h = 0, col = "red")
    
    qqnorm(residuals)
    qqline(residuals, col = "red")
    
    hist(residuals, main = "Histogram of Residuals")
    dev.off()
    
    cat("Nonlinear model diagnostic plots saved to 'plots' directory.\n")
  }, silent = TRUE)
}

# Function to fit state space models
fit_state_space_models <- function(data) {
  cat("\n--- Fitting State Space Models ---\n")
  
  # Create a time series object
  y_ts <- ts(data$Y, start = min(data$`Year.t.`), end = max(data$`Year.t.`))
  
  # Basic local level model using dlm package
  try({
    mod_dlm <- dlmModPoly(order = 1, dV = 1, dW = 1)
    fit_dlm <- dlmMLE(y_ts, mod_dlm)
    mod_dlm_fit <- dlmModPoly(order = 1, dV = exp(fit_dlm$par[1]), 
                            dW = exp(fit_dlm$par[2]))
    
    # Kalman filter and smoother
    filter_dlm <- dlmFilter(y_ts, mod_dlm_fit)
    smooth_dlm <- dlmSmooth(filter_dlm)
    
    cat("\nBasic Local Level Model (DLM):\n")
    cat("Observation variance:", exp(fit_dlm$par[1]), "\n")
    cat("State variance:", exp(fit_dlm$par[2]), "\n")
    cat("Log-likelihood:", fit_dlm$value, "\n")
    
    # Create state space model plot
    png("plots/state_space_model.png", width = 800, height = 600)
    par(mfrow = c(1, 1))
    plot(y_ts, main = "Local Level State Space Model", 
         ylab = "Y", type = "o", pch = 20, col = "darkgray")
    lines(dropFirst(filter_dlm$m), col = "blue", lwd = 2)
    lines(dropFirst(smooth_dlm$s), col = "red", lwd = 2)
    legend("topright", legend = c("Observed", "Filtered", "Smoothed"),
           col = c("darkgray", "blue", "red"), lwd = c(1, 2, 2), pch = c(20, NA, NA))
    dev.off()
    
    cat("State space model plot saved to 'plots' directory.\n")
  }, silent = TRUE)
  
  # For time series with seasonality or trend components, use structural time series
  try({
    # Create a more complex model with KFAS package
    data$t <- seq_along(data$Y)
    X_reg <- cbind(data$X, data$T, data$W)
    colnames(X_reg) <- c("X", "T", "W")
    
    # Define the state space model with trend and regression components
    model_spec <- SSModel(Y ~ SSMtrend(1, Q = NA) + SSMregression(~X+T+W, 
                                                                data = as.data.frame(X_reg), Q = 0), 
                        data = data, H = NA)
    
    # Fit the model
    fit_kfas <- fitSSM(model_spec, inits = c(0, 0))
    model_fitted <- fit_kfas$model
    
    # Extract components
    comp <- KFS(model_fitted, smoothing = c("state", "mean"))
    
    cat("\nStructural Time Series Model (KFAS):\n")
    cat("Observation variance:", model_fitted$H, "\n")
    cat("State variance:", model_fitted$Q, "\n")
    cat("Log-likelihood:", logLik(model_fitted), "\n")
    
    # Create forecast from the state space model
    forecast_kfas <- predict(model_fitted, n.ahead = 5)
    
    # Plot the results
    png("plots/structural_time_series.png", width = 800, height = 600)
    plot(data$Y, type = "o", pch = 20, col = "darkgray",
         main = "Structural Time Series Model with Regression",
         ylab = "Y", xlab = "Time")
    lines(comp$alphahat[, 1], col = "blue", lwd = 2)
    
    # Add forecast if available
    if(!is.null(forecast_kfas)) {
      points(seq(length(data$Y) + 1, length(data$Y) + length(forecast_kfas$mean)), 
             forecast_kfas$mean, col = "red", pch = 20)
      lines(seq(length(data$Y) + 1, length(data$Y) + length(forecast_kfas$mean)), 
            forecast_kfas$mean, col = "red", lwd = 2)
    }
    
    legend("topright", legend = c("Observed", "Fitted", "Forecast"),
           col = c("darkgray", "blue", "red"), lwd = c(1, 2, 2), pch = c(20, NA, 20))
    dev.off()
    
    cat("Structural time series model plot saved to 'plots' directory.\n")
  }, silent = TRUE)
}

# Compare models
compare_models <- function() {
  cat("\n--- Model Comparison ---\n")
  cat("Model comparison could be based on:\n")
  cat("1. Information criteria (AIC, BIC)\n")
  cat("2. Prediction accuracy (RMSE, MAE)\n")
  cat("3. Residual diagnostics\n")
  cat("4. Cross-validation results\n\n")
  cat("To perform a formal comparison, select the best models from each category\n")
  cat("and evaluate them based on the criteria above.\n")
}

# Main function
main <- function() {
  cat("--- Statistical Models Analysis ---\n")
  
  # Read and prepare data
  data <- read_data()
  
  # Exploratory data analysis
  explore_data(data)
  
  # Fit nonlinear models
  fit_nonlinear_models(data)
  
  # Fit state space models
  fit_state_space_models(data)
  
  # Compare models
  compare_models()
  
  cat("\nAnalysis complete. Review the model outputs and plots for insights.\n")
}

# Run the analysis
main() 