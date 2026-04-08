# --------------------------------------------------------------
# 04_Estadistica_Bayesiana
# ANOVA Bayesiano para rendimiento bajo estrés salino
# --------------------------------------------------------------

set.seed(777)

n <- 35
datos <- data.frame(
  tratamiento = factor(rep(c("Control", "Salinidad_media", "Salinidad_alta"), each = n)),
  rendimiento = c(
    rnorm(n, mean = 58, sd = 4),
    rnorm(n, mean = 49, sd = 5),
    rnorm(n, mean = 41, sd = 6)
  )
)

if (requireNamespace("BayesFactor", quietly = TRUE)) {
  cat("\n--- ANOVA Bayesiano (BayesFactor::anovaBF) ---\n")
  bf <- BayesFactor::anovaBF(rendimiento ~ tratamiento, data = datos)
  print(bf)

  # Comparación posterior por pares (opcional)
  cat("\n--- Comparaciones por pares bayesianas ---\n")
  print(BayesFactor::ttestBF(
    x = datos$rendimiento[datos$tratamiento == "Control"],
    y = datos$rendimiento[datos$tratamiento == "Salinidad_alta"]
  ))
} else {
  message("Paquete 'BayesFactor' no disponible. Instálalo con: install.packages('BayesFactor')")
}
