# -----------------------------------------------------------------
# 05_Machine_Learning
# XGBoost para predecir respuesta fenotípica bajo estrés hídrico
# -----------------------------------------------------------------

set.seed(911)

n <- 700
datos <- data.frame(
  evapotranspiracion = rnorm(n, 5.2, 1.1),
  humedad_suelo = runif(n, 8, 42),
  temperatura = rnorm(n, 30, 3.5),
  ndvi = runif(n, 0.2, 0.9),
  genotipo = factor(sample(paste0("G", 1:10), n, replace = TRUE))
)

datos$indice_resiliencia <- 60 +
  3.5 * datos$ndvi -
  1.8 * datos$evapotranspiracion +
  0.7 * datos$humedad_suelo -
  0.6 * datos$temperatura +
  ifelse(datos$genotipo %in% c("G3", "G7"), 4, -1) +
  rnorm(n, 0, 4)

# Matriz de diseño (one-hot encoding para factores)
X <- model.matrix(indice_resiliencia ~ . - 1, data = datos)
y <- datos$indice_resiliencia

id_train <- sample(seq_len(n), size = floor(0.8 * n))
X_train <- X[id_train, ]; y_train <- y[id_train]
X_test <- X[-id_train, ]; y_test <- y[-id_train]

if (requireNamespace("xgboost", quietly = TRUE)) {
  dtrain <- xgboost::xgb.DMatrix(data = X_train, label = y_train)
  dtest <- xgboost::xgb.DMatrix(data = X_test, label = y_test)

  params <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    eta = 0.08,
    max_depth = 5,
    subsample = 0.8,
    colsample_bytree = 0.8
  )

  xgb <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = 250,
    watchlist = list(train = dtrain, test = dtest),
    verbose = 0
  )

  pred <- predict(xgb, X_test)
  rmse <- sqrt(mean((pred - y_test)^2))

  cat("\n--- XGBoost ---\n")
  cat("RMSE (test):", round(rmse, 3), "\n")
} else {
  message("Paquete 'xgboost' no disponible. Instálalo con: install.packages('xgboost')")
}
