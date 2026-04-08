# 00_exploracion_y_supuestos.R
# Exploración inicial y chequeo de supuestos para modelos lineales.

# Paquetes necesarios
required_packages <- c("car", "lmtest", "performance", "ggplot2", "dplyr")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste(
      "Instala los siguientes paquetes antes de ejecutar este script:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(car)
library(lmtest)
library(performance)
library(ggplot2)
library(dplyr)

set.seed(123)

# Datos experimentales de ejemplo
n <- 90
datos <- data.frame(
  tratamiento = factor(rep(c("Control", "T1", "T2"), each = 30)),
  bloque = factor(rep(1:10, each = 9)),
  covariable = rnorm(n, mean = 50, sd = 10)
)

# Respuesta continua con efecto de tratamiento + covariable + ruido
efecto_tratamiento <- c(Control = 0, T1 = 2.5, T2 = 4)
datos$respuesta <- 20 +
  efecto_tratamiento[datos$tratamiento] +
  0.2 * datos$covariable +
  rnorm(n, sd = 2.5)

# Exploración rápida
str(datos)
summary(datos)

# Visualizaciones básicas
print(
  ggplot(datos, aes(x = tratamiento, y = respuesta, fill = tratamiento)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Respuesta por tratamiento")
)

print(
  ggplot(datos, aes(x = covariable, y = respuesta, color = tratamiento)) +
    geom_point(alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_minimal() +
    labs(title = "Relación respuesta-covariable")
)

# Modelo base para diagnósticos
modelo_lm <- lm(respuesta ~ tratamiento + covariable, data = datos)

# 1) Normalidad de residuos (Shapiro-Wilk)
shapiro_res <- shapiro.test(residuals(modelo_lm))
print(shapiro_res)

# Complemento visual de normalidad (recomendado para muestras moderadas/grandes)
invisible({
  qqnorm(residuals(modelo_lm))
  qqline(residuals(modelo_lm), col = "red", lwd = 2)
})

# 2) Homocedasticidad
# 2a. Levene (por factor de tratamiento)
levene_res <- leveneTest(respuesta ~ tratamiento, data = datos)
print(levene_res)

# 2b. Breusch-Pagan (sobre modelo lineal)
bp_res <- bptest(modelo_lm)
print(bp_res)

# 3) Detección de outliers
# 3a. Outliers por IQR en la respuesta
q1 <- quantile(datos$respuesta, 0.25)
q3 <- quantile(datos$respuesta, 0.75)
iqr <- IQR(datos$respuesta)
lim_inf <- q1 - 1.5 * iqr
lim_sup <- q3 + 1.5 * iqr
outliers_iqr <- datos %>%
  mutate(id = row_number()) %>%
  filter(respuesta < lim_inf | respuesta > lim_sup)

cat("\nOutliers por criterio IQR:\n")
print(outliers_iqr)

# 3b. Puntos influyentes por distancia de Cook
cooks <- cooks.distance(modelo_lm)
umbral_cook <- 4 / nrow(datos)
influyentes <- which(cooks > umbral_cook)
cat("\nObservaciones influyentes (Cook >", round(umbral_cook, 4), "):\n")
print(influyentes)

# 4) Diagnóstico integral con performance
# Incluye chequeo visual de residuos, linealidad y colinealidad (VIF)
print(check_model(modelo_lm))

# VIF explícito
print(check_collinearity(modelo_lm))
