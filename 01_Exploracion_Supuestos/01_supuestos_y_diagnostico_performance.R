# ---------------------------------------------------------------
# 01_Exploracion_Supuestos
# Pruebas de normalidad, homocedasticidad y diagnóstico de modelo
# Contexto: biomasa de plantas bajo estrés hídrico simulado
# ---------------------------------------------------------------

set.seed(123)

# Datos simulados: tres tratamientos de riego
n <- 40
datos <- data.frame(
  tratamiento = factor(rep(c("Control", "Deficit_moderado", "Deficit_severo"),
                           each = n)),
  biomasa = c(
    rnorm(n, mean = 42, sd = 4.5),
    rnorm(n, mean = 35, sd = 5.5),
    rnorm(n, mean = 29, sd = 6.0)
  )
)

# Modelo lineal simple para evaluación de supuestos
mod <- lm(biomasa ~ tratamiento, data = datos)

# 1) Normalidad de residuos
cat("\n--- Prueba de Shapiro-Wilk sobre residuos ---\n")
print(shapiro.test(residuals(mod)))

# 2) Homocedasticidad (Levene, si car está disponible)
if (requireNamespace("car", quietly = TRUE)) {
  cat("\n--- Prueba de Levene ---\n")
  print(car::leveneTest(biomasa ~ tratamiento, data = datos))
} else {
  message("Paquete 'car' no disponible. Instálalo para correr Levene: install.packages('car')")
}

# 3) Diagnóstico integral con performance
if (requireNamespace("performance", quietly = TRUE)) {
  cat("\n--- Diagnóstico global de supuestos (performance::check_model) ---\n")
  print(performance::check_model(mod))
} else {
  message("Paquete 'performance' no disponible. Instálalo con: install.packages('performance')")
}

# 4) Gráficos diagnósticos base de R
par(mfrow = c(2, 2))
plot(mod)
par(mfrow = c(1, 1))

# Este script se puede adaptar a otras variables (altura, área foliar, rendimiento).
