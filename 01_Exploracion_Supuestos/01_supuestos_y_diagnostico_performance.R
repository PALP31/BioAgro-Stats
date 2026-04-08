# ============================================================================
# 01_Exploracion_Supuestos
# Guía Exhaustiva de Evaluación de Supuestos Estadísticos
# Desde métodos clásicos hasta diagnóstico avanzado con easystats
# ============================================================================
# Autor: Claude Scholar
# Paquetes requeridos: performance, easystats, car, lmtest, DHARMa, ggplot2
# ============================================================================

library(tidyverse)
library(performance)
library(easystats)
library(car)
library(lmtest)
library(DHARMa)
library(patchwork)

set.seed(42)

# ============================================================================
# PARTE I: DATOS GAUSSIANOS - Rendimiento de Trigo Bajo Estrés Térmico
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE I: SUPUESTOS PARA MODELOS GAUSSIANOS\n")
cat("Contexto: Rendimiento de trigo bajo estrés térmico\n")
cat(rep("=", 80), "\n\n", sep = "")

# ----------------------------------------------------------------------------
# 1.1 SIMULACIÓN DE DATOS
# ----------------------------------------------------------------------------

n_por_grupo <- 50
grupos <- c("Control", "Estres_moderado", "Estres_severo")

datos_trigo <- tibble(
  grupo = factor(rep(grupos, each = n_por_grupo)),
  temperatura = c(
    rnorm(n_por_grupo, 28, 2),
    rnorm(n_por_grupo, 34, 2),
    rnorm(n_por_grupo, 40, 2)
  ),
  humedad = c(
    rnorm(n_por_grupo, 65, 5),
    rnorm(n_por_grupo, 55, 5),
    rnorm(n_por_grupo, 45, 5)
  ),
  clorofila = c(
    rnorm(n_por_grupo, 45, 6),
    rnorm(n_por_grupo, 38, 6),
    rnorm(n_por_grupo, 30, 6)
  ),
  biomasa = c(
    rnorm(n_por_grupo, 120, 12),
    rnorm(n_por_grupo, 95, 15),
    rnorm(n_por_grupo, 70, 18)
  ),
  rendimiento = c(
    rnorm(n_por_grupo, 85, 8),
    rnorm(n_por_grupo, 68, 10),
    rnorm(n_por_grupo, 50, 14)
  )
)

# Modelo lineal múltiple
modelo_gaussiano <- lm(rendimiento ~ grupo + temperatura + humedad + clorofila,
                       data = datos_trigo)

cat("\n--- Resumen del Modelo Gaussiano ---\n")
summary(modelo_gaussiano)

# ----------------------------------------------------------------------------
# 1.2 MÉTODOS CLÁSICOS - NIVEL BÁSICO
# ----------------------------------------------------------------------------

cat("\n", rep("-", 60), "\n", sep = "")
cat("1.2 MÉTODOS CLÁSICOS DE DIAGNÓSTICO\n")
cat(rep("-", 60), "\n\n", sep = "")

# ---- 1.2.1 Normalidad de Residuos (Shapiro-Wilk) ----
cat(">>> 1. Shapiro-Wilk: Normalidad de residuos\n")
residuos <- residuals(modelo_gaussiano)

shapiro_result <- shapiro.test(residuos)
print(shapiro_result)

cat("\nInterpretación: ")
if (shapiro_result$p.value > 0.05) {
  cat("No hay evidencia para rechazar la normalidad (p =",
      round(shapiro_result$p.value, 4), ")\n")
} else {
  cat("Los residuos se desvían significativamente de la normalidad (p <",
      round(shapiro_result$p.value, 4), ")\n")
}

# Gráfico Q-Q
qq_plot <- ggplot(data.frame(residuos = residuos), aes(sample = residuos)) +
  stat_qq(color = "#2E86AB", size = 3, alpha = 0.7) +
  stat_qq_line(color = "#E74C3C", linewidth = 1.2) +
  labs(title = "Q-Q Plot de Residuos",
       subtitle = "Evaluación visual de normalidad",
       x = "Cuantiles teóricos", y = "Cuantiles observados") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

print(qq_plot)
ggsave("01_Exploracion_Supuestos/01_qq_plot_gaussiano.png", qq_plot,
       width = 6, height = 5, dpi = 300)

# ---- 1.2.2 Homocedasticidad (Levene) ----
cat("\n>>> 2. Prueba de Levene: Homocedasticidad\n")
levene_result <- leveneTest(rendimiento ~ grupo, data = datos_trigo)
print(levene_result)

cat("\nInterpretación: ")
if (levene_result$`Pr(>F)`[1] > 0.05) {
  cat("Homocedasticidad confirmada (p =",
      round(levene_result$`Pr(>F)`[1], 4), ")\n")
} else {
  cat("Heterocedasticidad detectada (p <",
      round(levene_result$`Pr(>F)`[1], 4), ")\n")
}

# ---- 1.2.3 Homocedasticidad (Breusch-Pagan) ----
cat("\n>>> 3. Breusch-Pagan: Homocedasticidad (alternativa)\n")
bp_result <- bptest(modelo_gaussiano)
print(bp_result)

cat("\nInterpretación: ")
if (bp_result$p.value > 0.05) {
  cat("No hay evidencia de heterocedasticidad (p =",
      round(bp_result$p.value, 4), ")\n")
} else {
  cat("Heterocedasticidad detectada (p <",
      round(bp_result$p.value, 4), ")\n")
}

# ----------------------------------------------------------------------------
# 1.3 MÉTODOS AVANZADOS - ECOSISTEMA EASYSTATS
# ----------------------------------------------------------------------------

cat("\n", rep("-", 60), "\n", sep = "")
cat("1.3 DIAGNÓSTICO AVANZADO CON EASYSTATS/PERFORMANCE\n")
cat(rep("-", 60), "\n\n", sep = "")

# ---- 1.3.1 check_model: Diagnóstico Global ----
cat(">>> 4. check_model: Evaluación integral de supuestos\n")
cat("Generando 4 paneles: Normalidad, Homocedasticidad, Independencia, Influencia\n")

modelo_check <- check_model(modelo_gaussiano,
                            panel = TRUE,
                            theme = "bw",
                            colors = list(viridis = TRUE))

print(modelo_check)
ggsave("01_Exploracion_Supuestos/02_check_model_gaussiano.png",
       plot = modelo_check, width = 10, height = 8, dpi = 300)

# ---- 1.3.2 check_normality ----
cat("\n>>> 5. check_normality: Prueba formal de normalidad\n")
normality_check <- check_normality(modelo_gaussiano)
print(normality_check)

# ---- 1.3.3 check_homoscedasticity ----
cat("\n>>> 6. check_homoscedasticity: Homocedasticidad\n")
homo_check <- check_homoscedasticity(modelo_gaussiano)
print(homo_check)

# ---- 1.3.4 check_outliers: Detección de Valores Atípicos ----
cat("\n>>> 7. check_outliers: Detección de outliers\n")
outliers_check <- check_outliers(modelo_gaussiano,
                                 method = "all",
                                 threshold = NULL)
print(outliers_check)

# Gráfico de outliers
outliers_plot <- check_outliers(modelo_gaussiano, plot = TRUE)
print(outliers_plot)
ggsave("01_Exploracion_Supuestos/03_outliers_gaussiano.png",
       plot = outliers_plot, width = 7, height = 5, dpi = 300)

# ---- 1.3.5 check_collinearity: Multicolinealidad ----
cat("\n>>> 8. check_collinearity: VIF y tolerancia\n")
collin_check <- check_collinearity(modelo_gaussiano)
print(collin_check)

cat("\nInterpretación VIF:\n")
vif_vals <- as.data.frame(collin_check)$VIF
for (i in seq_along(vif_vals)) {
  vif <- vif_vals[i]
  predictor <- rownames(as.data.frame(collin_check))[i]
  if (vif < 5) {
    cat(sprintf("  %s: VIF = %.2f (✓ Sin multicolinealidad)\n", predictor, vif))
  } else if (vif < 10) {
    cat(sprintf("  %s: VIF = %.2f (⚠ Multicolinealidad moderada)\n", predictor, vif))
  } else {
    cat(sprintf("  %s: VIF = %.2f (✗ Multicolinealidad severa)\n", predictor, vif))
  }
}

# ---- 1.3.6 check_independence ----
cat("\n>>> 9. check_independence: Independencia de residuos (Durbin-Watson)\n")
indep_check <- check_independence(modelo_gaussiano)
print(indep_check)

# ---- 1.3.7 check_sphericity ----
cat("\n>>> 10. check_sphericity: Esfericidad\n")
sphericity_check <- check_sphericity(modelo_gaussiano)
print(sphericity_check)

# ============================================================================
# PARTE II: DATOS DE CONTEO - Edafofauna con Exceso de Ceros
# ============================================================================

cat("\n\n", rep("=", 80), "\n", sep = "")
cat("PARTE II: SUPUESTOS PARA MODELOS DE CONTEO (POISSON/BINOMIAL NEGATIVA)\n")
cat("Contexto: Conteo de edafofauna del suelo con exceso de ceros\n")
cat(rep("=", 80), "\n\n", sep = "")

# ----------------------------------------------------------------------------
# 2.1 SIMULACIÓN DE DATOS CON EXCESO DE CEROS
# ----------------------------------------------------------------------------

n_suelos <- 200

# Simular conteo con exceso de ceros (zero-inflated)
set.seed(123)
# Proceso de generación de ceros estructurales
prob_cero <- 0.4  # 40% de ceros estructurales
es_cero <- rbinom(n_suelos, 1, prob_cero)

# Proceso de conteo (Poisson con sobredispersión)
conteo_puro <- rnbinom(n_suelos, mu = 8, size = 2)

# Combinar ambos procesos
conteo_final <- ifelse(es_cero == 1, 0, conteo_puro)

datos_edafo <- tibble(
  sitio = paste0("S", 1:n_suelos),
  tipo_suelo = factor(sample(c("Arcilloso", "Franco", "Arenoso"),
                             n_suelos, replace = TRUE)),
  ph = rnorm(n_suelos, 6.5, 0.8),
  materia_organica = rnorm(n_suelos, 3.2, 1.1),
  humedad_suelo = rnorm(n_suelos, 25, 8),
  profundidad_cm = sample(c(5, 10, 15, 20), n_suelos, replace = TRUE),
  conteo_fauna = conteo_final
)

cat("\nDistribución del conteo de fauna:\n")
print(datos_edafo %>%
        count(conteo_fauna) %>%
        arrange(conteo_fauna) %>%
        head(15))

cat("\nPorcentaje de ceros:",
    round(mean(datos_edafo$conteo_fauna == 0) * 100, 1), "%\n")

# Modelo Poisson
modelo_poisson <- glm(conteo_fauna ~ tipo_suelo + ph + materia_organica +
                        humedad_suelo + profundidad_cm,
                      family = poisson(link = "log"),
                      data = datos_edafo)

# Modelo Binomial Negativa (si MASS está disponible)
if (requireNamespace("MASS", quietly = TRUE)) {
  modelo_nb <- MASS::glm.nb(conteo_fauna ~ tipo_suelo + ph + materia_organica +
                              humedad_suelo + profundidad_cm,
                            data = datos_edafo)
  cat("\nModelo Binomial Negativa ajustado exitosamente.\n")
} else {
  cat("\nAdvertencia: Paquete MASS no disponible. Solo se usará Poisson.\n")
  modelo_nb <- NULL
}

# ----------------------------------------------------------------------------
# 2.2 EVALUACIÓN DE SUPUESTOS ESPECÍFICOS PARA CONTEOS
# ----------------------------------------------------------------------------

cat("\n", rep("-", 60), "\n", sep = "")
cat("2.2 SUPUESTOS ESPECÍFICOS PARA MODELOS DE CONTEO\n")
cat(rep("-", 60), "\n\n", sep = "")

# ---- 2.2.1 check_overdispersion: Sobredispersión ----
cat(">>> 11. check_overdispersion: Evaluación de sobredispersión\n")
overdisp_check <- check_overdispersion(modelo_poisson)
print(overdisp_check)

cat("\nInterpretación: ")
if (overdisp_check$overdispersion) {
  cat("SOBERDISPERSIÓN DETECTADA (razón =",
      round(overdisp_check$ratio, 2),
      "). Se recomienda Binomial Negativa en lugar de Poisson.\n")
} else {
  cat("No hay sobredispersión significativa. Poisson es adecuado.\n")
}

# Comparar Poisson vs Binomial Negativa
if (!is.null(modelo_nb)) {
  cat("\n>>> Comparación de modelos:\n")
  cat("Poisson AIC:", AIC(modelo_poisson), "\n")
  cat("Binomial Negativa AIC:", AIC(modelo_nb), "\n")

  if (AIC(modelo_nb) < AIC(modelo_poisson)) {
    cat("→ Binomial Negativa es preferible (menor AIC)\n")
  } else {
    cat("→ Poisson es preferible (menor AIC)\n")
  }
}

# ---- 2.2.2 check_zeroinflation: Inflación de Ceros ----
cat("\n>>> 12. check_zeroinflation: Evaluación de inflación de ceros\n")
zeroinf_check <- check_zeroinflation(modelo_poisson)
print(zeroinf_check)

cat("\nInterpretación: ")
if (zeroinf_check$zeroinflation) {
  cat("INFLACIÓN DE CEROS DETECTADA (razón =",
      round(zeroinf_check$ratio, 2),
      "). Se recomienda modelo Zero-Inflated o Hurdle.\n")
} else {
  cat("No hay inflación de ceros significativa.\n")
}

# ---- 2.2.3 DHARMa: Diagnóstico de Residuos Simulados ----
cat("\n>>> 13. DHARMa: Residuos simulados para GLM\n")

# Crear residuos DHARMa
if (requireNamespace("DHARMa", quietly = TRUE)) {
  simulation_output <- simulateResiduals(fittedModel = modelo_poisson, n = 1000)

  # Prueba de uniformidad
  cat("\nPrueba de uniformidad de residuos DHARMa:\n")
  print(testUniformity(simulation_output))

  # Prueba de dispersión
  cat("\nPrueba de dispersión DHARMa:\n")
  print(testDispersion(simulation_output))

  # Prueba de inflación de ceros
  cat("\nPrueba de inflación de ceros DHARMa:\n")
  print(testZeroInflation(simulation_output))

  # Gráfico de residuos DHARMa
  dharma_plot <- plot(simulation_output,
                      main = "Residuos DHARMa - Modelo Poisson")
  ggsave("01_Exploracion_Supuestos/04_dharma_poisson.png",
         plot = last_plot(), width = 7, height = 5, dpi = 300)
}

# ----------------------------------------------------------------------------
# 2.3 DIAGNÓSTICO GLOBAL CON EASYSTATS
# ----------------------------------------------------------------------------

cat("\n", rep("-", 60), "\n", sep = "")
cat("2.4 DIAGNÓSTICO GLOBAL - MODELO DE CONTEO\n")
cat(rep("-", 60), "\n\n", sep = "")

# Usar modelo Binomial Negativa si está disponible
modelo_final <- if (!is.null(modelo_nb)) modelo_nb else modelo_poisson

cat(">>> 14. check_model para modelo de conteo\n")
modelo_conteo_check <- check_model(modelo_final,
                                   check = c("normality", "outliers",
                                            "homoscedasticity"))
print(modelo_conteo_check)
ggsave("01_Exploracion_Supuestos/05_check_model_conteo.png",
       plot = modelo_conteo_check, width = 10, height = 8, dpi = 300)

cat("\n>>> 15. check_outliers para modelo de conteo\n")
outliers_conteo <- check_outliers(modelo_final, method = "all")
print(outliers_conteo)

# ============================================================================
# PARTE III: RESUMEN Y RECOMENDACIONES
# ============================================================================

cat("\n\n", rep("=", 80), "\n", sep = "")
cat("RESUMEN DE SUPUESTOS EVALUADOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Tabla resumen
resumen <- tibble(
  Supuesto = c("Normalidad", "Homocedasticidad", "Independencia",
               "No-multicolinealidad", "No-outliers", "Sobredispersión",
               "Inflación de ceros"),
  Funcion = c("shapiro.test(), check_normality()",
              "leveneTest(), bptest(), check_homoscedasticity()",
              "check_independence(), Durbin-Watson",
              "check_collinearity(), VIF",
              "check_outliers()",
              "check_overdispersion()",
              "check_zeroinflation(), DHARMa testZeroInflation()"),
  Paquete = c("stats, performance",
              "car, lmtest, performance",
              "performance",
              "performance",
              "performance",
              "performance",
              "performance, DHARMa")
)

print(resumen)

cat("\n\n", rep("=", 80), "\n", sep = "")
cat("RECOMENDACIONES POR TIPO DE MODELO\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("
┌─────────────────────────────────────────────────────────────────────────┐
│ MODELOS GAUSSIANOS (LM, LMM)                                            │
├─────────────────────────────────────────────────────────────────────────┤
│ ✓ check_model() - Diagnóstico integral                                  │
│ ✓ check_normality() - Shapiro-Wilk + Q-Q plot                           │
│ ✓ check_homoscedasticity() - Breusch-Pagan / Levene                     │
│ ✓ check_collinearity() - VIF < 5 (ideal), < 10 (aceptable)              │
│ ✓ check_outliers() - Identificar observaciones influyentes              │
│ ✓ check_independence() - Durbin-Watson para autocorrelación             │
└─────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────┐
│ MODELOS DE CONTEO (POISSON, BINOMIAL NEGATIVA)                          │
├─────────────────────────────────────────────────────────────────────────┤
│ ✓ check_overdispersion() - Razón > 1 indica sobredispersión             │
│ ✓ check_zeroinflation() - Razón > 1 indica exceso de ceros              │
│ ✓ DHARMa::simulateResiduals() - Residuos simulados                      │
│ ✓ DHARMa::testDispersion() - Prueba formal de dispersión                │
│ ✓ DHARMa::testZeroInflation() - Prueba formal de ceros                  │
│                                                                         │
│ Si hay sobredispersión → Usar Binomial Negativa                         │
│ Si hay inflación de ceros → Usar Zero-Inflated o Hurdle                 │
└─────────────────────────────────────────────────────────────────────────┘
")

cat("\n", rep("=", 80), "\n", sep = "")
cat("FIN DEL DIAGNÓSTICO\n")
cat(rep("=", 80), "\n")
