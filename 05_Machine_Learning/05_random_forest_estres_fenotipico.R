# --------------------------------------------------------------------
# 05_Machine_Learning
# Random Forest para predecir índice fenotípico bajo estrés térmico
# --------------------------------------------------------------------

set.seed(910)

n <- 500
datos <- data.frame(
  temperatura = rnorm(n, 32, 4),
  humedad_suelo = runif(n, 10, 45),
  clorofila = rnorm(n, 38, 7),
  conductancia = rnorm(n, 250, 55),
  genotipo = factor(sample(paste0("G", 1:8), n, replace = TRUE))
)

datos$respuesta_fenotipica <- 80 -
  1.1 * datos$temperatura +
  0.9 * datos$humedad_suelo +
  0.5 * datos$clorofila +
  ifelse(datos$genotipo %in% c("G2", "G5"), 6, 0) +
  rnorm(n, 0, 5)

id_train <- sample(seq_len(n), size = floor(0.8 * n))
train <- datos[id_train, ]
test <- datos[-id_train, ]

if (requireNamespace("randomForest", quietly = TRUE)) {
  rf <- randomForest::randomForest(
    respuesta_fenotipica ~ .,
    data = train,
    ntree = 400,
    mtry = 3,
    importance = TRUE
  )

  pred <- predict(rf, newdata = test)
  rmse <- sqrt(mean((pred - test$respuesta_fenotipica)^2))

  cat("\n--- Random Forest ---\n")
  cat("RMSE (test):", round(rmse, 3), "\n")
  print(importance(rf))
} else {
  message("Paquete 'randomForest' no disponible. Instálalo con: install.packages('randomForest')")
}
