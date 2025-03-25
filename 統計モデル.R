#!/usr/bin/env Rscript

# 必要なライブラリの読み込み
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,    # データ操作と可視化
  nlme,         # 非線形混合効果モデル
  nls2,         # 非線形回帰
  dlm,          # 動的線形モデル（状態空間モデル）
  forecast,     # 時系列予測
  KFAS,         # カルマンフィルターとスムーザー
  minpack.lm    # より堅牢な非線形最小二乗法
)

# データの読み込みと準備を行う関数
read_data <- function() {
  # CSVファイルの読み込み
  data <- read.csv("assets/WWdata.csv")
  
  # データ構造の表示
  cat("データ構造:\n")
  print(str(data))
  
  # 要約統計量の表示
  cat("\n要約統計量:\n")
  print(summary(data))
  
  # クリーニングされたデータを返す
  return(data)
}

# 探索的データ分析を行う関数
explore_data <- function(data) {
  # プロット用ディレクトリの作成
  if (!dir.exists("plots")) dir.create("plots")
  
  # 時系列プロット
  png("plots/time_series.png", width = 800, height = 600)
  par(mfrow = c(2, 2))
  plot(data$`Year.t.`, data$Y, type = "l", col = "blue", 
       main = "Yの時系列推移", xlab = "年", ylab = "Y")
  plot(data$`Year.t.`, data$X, type = "l", col = "red", 
       main = "Xの時系列推移", xlab = "年", ylab = "X")
  plot(data$`Year.t.`, data$T, type = "l", col = "green", 
       main = "Tの時系列推移", xlab = "年", ylab = "T")
  plot(data$`Year.t.`, data$W, type = "l", col = "purple", 
       main = "Wの時系列推移", xlab = "年", ylab = "W")
  dev.off()
  
  # 散布図による関係性の探索
  png("plots/scatter_plots.png", width = 800, height = 600)
  par(mfrow = c(2, 2))
  plot(data$X, data$Y, main = "YとXの散布図", xlab = "X", ylab = "Y")
  plot(data$T, data$Y, main = "YとTの散布図", xlab = "T", ylab = "Y")
  plot(data$W, data$Y, main = "YとWの散布図", xlab = "W", ylab = "Y")
  plot(data$X, data$W, main = "WとXの散布図", xlab = "X", ylab = "W")
  dev.off()
  
  cat("探索的プロットを'plots'ディレクトリに保存しました。\n")
}

# 非線形モデルのフィッティングを行う関数
fit_nonlinear_models <- function(data) {
  cat("\n--- 非線形モデルのフィッティング ---\n")
  
  # 1. 指数モデル: Y ~ a * exp(b * X)
  try({
    exp_model <- nls(Y ~ a * exp(b * X), 
                    data = data,
                    start = list(a = 1, b = 0.01),
                    control = nls.control(maxiter = 100),
                    algorithm = "port")
    
    cat("\n指数モデル (Y ~ a * exp(b * X)):\n")
    print(summary(exp_model))
    
    # R二乗の計算
    y_pred <- predict(exp_model)
    sse <- sum((data$Y - y_pred)^2)
    sst <- sum((data$Y - mean(data$Y))^2)
    rsq <- 1 - sse/sst
    cat("\nR二乗:", rsq, "\n")
  }, silent = TRUE)
  
  # 2. べき乗モデル: Y ~ a * X^b
  try({
    power_model <- nls(Y ~ a * X^b, 
                      data = data,
                      start = list(a = 1000, b = -0.5),
                      control = nls.control(maxiter = 100),
                      algorithm = "port")
    
    cat("\nべき乗モデル (Y ~ a * X^b):\n")
    print(summary(power_model))
    
    # R二乗の計算
    y_pred <- predict(power_model)
    sse <- sum((data$Y - y_pred)^2)
    sst <- sum((data$Y - mean(data$Y))^2)
    rsq <- 1 - sse/sst
    cat("\nR二乗:", rsq, "\n")
  }, silent = TRUE)
  
  # 3. 多変量モデル: Y ~ a * X^b * exp(c * T) * W^d
  try({
    multi_model <- nlsLM(Y ~ a * X^b * exp(c * T) * W^d, 
                       data = data,
                       start = list(a = 1000, b = -0.5, c = 0.1, d = 1),
                       control = nls.control(maxiter = 200))
    
    cat("\n多変量モデル (Y ~ a * X^b * exp(c * T) * W^d):\n")
    print(summary(multi_model))
    
    # R二乗の計算
    y_pred <- predict(multi_model)
    sse <- sum((data$Y - y_pred)^2)
    sst <- sum((data$Y - mean(data$Y))^2)
    rsq <- 1 - sse/sst
    cat("\nR二乗:", rsq, "\n")
    
    # 診断プロットの作成
    png("plots/nonlinear_model_diagnostics.png", width = 800, height = 600)
    par(mfrow = c(2, 2))
    plot(y_pred, data$Y, main = "予測値と実測値の比較",
         xlab = "予測値", ylab = "実測値")
    abline(0, 1, col = "red")
    
    residuals <- data$Y - y_pred
    plot(y_pred, residuals, main = "残差と予測値の散布図",
         xlab = "予測値", ylab = "残差")
    abline(h = 0, col = "red")
    
    qqnorm(residuals)
    qqline(residuals, col = "red")
    
    hist(residuals, main = "残差のヒストグラム")
    dev.off()
    
    cat("非線形モデルの診断プロットを'plots'ディレクトリに保存しました。\n")
  }, silent = TRUE)
}

# 状態空間モデルのフィッティングを行う関数
fit_state_space_models <- function(data) {
  cat("\n--- 状態空間モデルのフィッティング ---\n")
  
  # 時系列オブジェクトの作成
  y_ts <- ts(data$Y, start = min(data$`Year.t.`), end = max(data$`Year.t.`))
  
  # dlmパッケージを使用した基本的なローカルレベルモデル
  try({
    mod_dlm <- dlmModPoly(order = 1, dV = 1, dW = 1)
    fit_dlm <- dlmMLE(y_ts, mod_dlm)
    mod_dlm_fit <- dlmModPoly(order = 1, dV = exp(fit_dlm$par[1]), 
                            dW = exp(fit_dlm$par[2]))
    
    # カルマンフィルターとスムーザー
    filter_dlm <- dlmFilter(y_ts, mod_dlm_fit)
    smooth_dlm <- dlmSmooth(filter_dlm)
    
    cat("\n基本的なローカルレベルモデル (DLM):\n")
    cat("観測分散:", exp(fit_dlm$par[1]), "\n")
    cat("状態分散:", exp(fit_dlm$par[2]), "\n")
    cat("対数尤度:", fit_dlm$value, "\n")
    
    # 状態空間モデルのプロット作成
    png("plots/state_space_model.png", width = 800, height = 600)
    par(mfrow = c(1, 1))
    plot(y_ts, main = "ローカルレベル状態空間モデル", 
         ylab = "Y", type = "o", pch = 20, col = "darkgray")
    lines(dropFirst(filter_dlm$m), col = "blue", lwd = 2)
    lines(dropFirst(smooth_dlm$s), col = "red", lwd = 2)
    legend("topright", legend = c("観測値", "フィルター値", "スムーズ値"),
           col = c("darkgray", "blue", "red"), lwd = c(1, 2, 2), pch = c(20, NA, NA))
    dev.off()
    
    cat("状態空間モデルのプロットを'plots'ディレクトリに保存しました。\n")
  }, silent = TRUE)
  
  # 季節性やトレンド成分を含む時系列のための構造時系列モデル
  try({
    # KFASパッケージを使用したより複雑なモデル
    data$t <- seq_along(data$Y)
    X_reg <- cbind(data$X, data$T, data$W)
    colnames(X_reg) <- c("X", "T", "W")
    
    # トレンドと回帰成分を含む状態空間モデルの定義
    model_spec <- SSModel(Y ~ SSMtrend(1, Q = NA) + SSMregression(~X+T+W, 
                                                                data = as.data.frame(X_reg), Q = 0), 
                        data = data, H = NA)
    
    # モデルのフィッティング
    fit_kfas <- fitSSM(model_spec, inits = c(0, 0))
    model_fitted <- fit_kfas$model
    
    # 成分の抽出
    comp <- KFS(model_fitted, smoothing = c("state", "mean"))
    
    cat("\n構造時系列モデル (KFAS):\n")
    cat("観測分散:", model_fitted$H, "\n")
    cat("状態分散:", model_fitted$Q, "\n")
    cat("対数尤度:", logLik(model_fitted), "\n")
    
    # 状態空間モデルからの予測
    forecast_kfas <- predict(model_fitted, n.ahead = 5)
    
    # 結果のプロット
    png("plots/structural_time_series.png", width = 800, height = 600)
    plot(data$Y, type = "o", pch = 20, col = "darkgray",
         main = "回帰成分を含む構造時系列モデル",
         ylab = "Y", xlab = "時間")
    lines(comp$alphahat[, 1], col = "blue", lwd = 2)
    
    # 予測値の追加
    if(!is.null(forecast_kfas)) {
      points(seq(length(data$Y) + 1, length(data$Y) + length(forecast_kfas$mean)), 
             forecast_kfas$mean, col = "red", pch = 20)
      lines(seq(length(data$Y) + 1, length(data$Y) + length(forecast_kfas$mean)), 
            forecast_kfas$mean, col = "red", lwd = 2)
    }
    
    legend("topright", legend = c("観測値", "フィット値", "予測値"),
           col = c("darkgray", "blue", "red"), lwd = c(1, 2, 2), pch = c(20, NA, 20))
    dev.off()
    
    cat("構造時系列モデルのプロットを'plots'ディレクトリに保存しました。\n")
  }, silent = TRUE)
}

# モデルの比較を行う関数
compare_models <- function() {
  cat("\n--- モデルの比較 ---\n")
  cat("モデルの比較は以下の基準に基づいて行うことができます：\n")
  cat("1. 情報量基準（AIC、BIC）\n")
  cat("2. 予測精度（RMSE、MAE）\n")
  cat("3. 残差診断\n")
  cat("4. 交差検証結果\n\n")
  cat("正式な比較を行うには、各カテゴリーから最適なモデルを選択し、\n")
  cat("上記の基準に基づいて評価してください。\n")
}

# メイン関数
main <- function() {
  cat("--- 統計モデル分析 ---\n")
  
  # データの読み込みと準備
  data <- read_data()
  
  # 探索的データ分析
  explore_data(data)
  
  # 非線形モデルのフィッティング
  fit_nonlinear_models(data)
  
  # 状態空間モデルのフィッティング
  fit_state_space_models(data)
  
  # モデルの比較
  compare_models()
  
  cat("\n分析が完了しました。モデルの出力とプロットを確認して洞察を得てください。\n")
}

# 分析の実行
main() 