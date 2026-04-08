# -----------------------------------------------------------------
# 05_Machine_Learning
# SVM para clasificar tolerancia fenotípica a estrés térmico
# -----------------------------------------------------------------

set.seed(912)

n <- 450
datos <- data.frame(
  temperatura_hoja = rnorm(n, 34, 2.8),
  eficiencia_fotosintetica = rnorm(n, 0.68, 0.08),
  contenido_relativo_agua = rnorm(n, 78, 9),
  prolina = rnorm(n, 22, 7)
)

score <- -35 +
  0.9 * datos$eficiencia_fotosintetica * 100 +
  0.15 * datos$contenido_relativo_agua -
  0.6 * datos$temperatura_hoja +
  0.1 * datos$prolina +
  rnorm(n, 0, 3)

datos$clase_tolerancia <- factor(ifelse(score > median(score), "Alta", "Baja"))

id_train <- sample(seq_len(n), size = floor(0.75 * n))
train <- datos[id_train, ]
test <- datos[-id_train, ]

if (requireNamespace("e1071", quietly = TRUE)) {
  svm_fit <- e1071::svm(
    clase_tolerancia ~ .,
    data = train,
    kernel = "radial",
    gamma = 0.08,
    cost = 4
  )

  pred <- predict(svm_fit, newdata = test)
  acc <- mean(pred == test$clase_tolerancia)

  cat("\n--- SVM (clasificación) ---\n")
  cat("Accuracy (test):", round(acc, 3), "\n")
  print(table(Predicho = pred, Real = test$clase_tolerancia))
} else {
  message("Paquete 'e1071' no disponible. Instálalo con: install.packages('e1071')")
}
