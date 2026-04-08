# 01_Exploracion_Supuestos/01_supuestos_modelo.R
# Ensayo simulado: estrés térmico en trigo

suppressPackageStartupMessages({
  library(car)
  library(performance)
  library(ggplot2)
})

set.seed(123)

datos <- data.frame(
  bloque = factor(rep(1:4, each = 15)),
  tratamiento = factor(rep(c("Control", "Moderado", "Severo"), times = 20)),
  rendimiento = c(
    rnorm(20, mean = 6.2, sd = 0.35),
    rnorm(20, mean = 5.5, sd = 0.45),
    rnorm(20, mean = 4.7, sd = 0.55)
  )
)

modelo <- lm(rendimiento ~ bloque + tratamiento, data = datos)

cat("=== Normalidad de residuos (Shapiro-Wilk) ===\n")
print(shapiro.test(residuals(modelo)))

cat("\n=== Homocedasticidad (Levene) ===\n")
print(car::leveneTest(rendimiento ~ tratamiento, data = datos))

cat("\n=== Diagnostico integral (performance) ===\n")
print(performance::check_model(modelo))

# Graficos base de diagnostico
par(mfrow = c(2, 2))
plot(modelo)
par(mfrow = c(1, 1))

# Grafico complementario
p <- ggplot(datos, aes(tratamiento, rendimiento, fill = tratamiento)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Rendimiento por nivel de estrés térmico",
    x = "Tratamiento",
    y = "Rendimiento (t/ha)"
  )
print(p)
