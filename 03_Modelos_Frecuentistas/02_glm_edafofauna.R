# ============================================================================
# 02_glm_edafofauna.R
# GLM para datos de conteo: Poisson y Binomial Negativa
# Contexto: Abundancia de edafofauna en huertos frutales con bioestimulantes
# ============================================================================
# Diseño: Completamente al Azar
# Factores: Tipo de bioestimulante, Dosis, Tipo de suelo
# Variable respuesta: Abundancia de colémbolos (conteo)
# ============================================================================

library(tidyverse)
library(MASS)        # Para glm.nb (Binomial Negativa)
library(emmeans)
library(performance)
library(DHARMa)
library(car)

set.seed(123)

# ============================================================================
# 1. SIMULACIÓN DE DATOS
# ============================================================================

# Parámetros
n_rep <- 15  # Repeticiones por tratamiento
bioestimulantes <- c("Control", "Micorrizas", "Trichoderma", "Consorcio")
dosis <- c("Baja", "Alta")
tipos_suelo <- c("Franco", "Arcilloso", "Arenoso")

# Crear diseño factorial completo
datos_edafo <- expand.grid(
  Bioestimulante = factor(bioestimulantes),
  Dosis = factor(dosis),
  Tipo_suelo = factor(tipos_suelo),
  Repeticion = 1:n_rep
)

# Efectos de los tratamientos
efecto_bio <- c(
  "Control" = 0,
  "Micorrizas" = 0.4,
  "Trichoderma" = 0.3,
  "Consorcio" = 0.7
)

efecto_dosis <- c("Baja" = 0, "Alta" = 0.25)

efecto_suelo <- c(
  "Franco" = 0,
  "Arcilloso" = -0.2,
  "Arenoso" = -0.4
)

# Calcular lambda (media esperada) para cada observación
datos_edafo <- datos_edafo %>%
  mutate(
    efecto_b = efecto_bio[as.character(Bioestimulante)],
    efecto_d = efecto_dosis[as.character(Dosis)],
    efecto_s = efecto_suelo[as.character(Tipo_suelo)],
    log_lambda = 2.8 + efecto_b + efecto_d + efecto_s,
    lambda = exp(log_lambda)
  )

# Simular conteos con sobredispersión (usando Binomial Negativa)
# size = 2 indica sobredispersión moderada
datos_edafo <- datos_edafo %>%
  mutate(
    abundancia_colémbolos = rnbinom(n(), mu = lambda, size = 2)
  ) %>%
  select(-efecto_b, -efecto_d, -efecto_s, -log_lambda, -lambda)

cat("\n=== Estructura de los Datos ===\n")
str(datos_edafo)

cat("\nDistribución de la abundancia:\n")
datos_edafo %>%
  count(abundancia_colémbolos) %>%
  arrange(abundancia_colémbolos) %>%
  print(n = 20)

cat("\nResumen por tratamiento:\n")
datos_edafo %>%
  group_by(Bioestimulante, Dosis, Tipo_suelo) %>%
  summarise(
    n = n(),
    media = mean(abundancia_colémbolos),
    sd = sd(abundancia_colémbolos),
    var = var(abundancia_colémbolos),
    .groups = "drop"
  ) %>%
  print()

# ============================================================================
# 2. EVALUACIÓN PRELIMINAR: ¿POISSON O BINOMIAL NEGATIVA?
# ============================================================================

cat("\n=== Evaluación de Sobredispersión ===\n")

# Ajustar modelo Poisson primero
modelo_poisson <- glm(
  abundancia_colémbolos ~ Bioestimulante * Dosis * Tipo_suelo,
  family = poisson(link = "log"),
  data = datos_edafo
)

cat("\nModelo Poisson - Devianza:\n")
cat("Devianza:", deviance(modelo_poisson), "\n")
cat("Grados de libertad:", df.residual(modelo_poisson), "\n")
cat("Razón Devianza/GL:", round(deviance(modelo_poisson) / df.residual(modelo_poisson), 2), "\n")

# Test formal de sobredispersión
cat("\nTest de sobredispersión (performance):\n")
overdisp_check <- check_overdispersion(modelo_poisson)
print(overdisp_check)

# Si hay sobredispersión, usar Binomial Negativa
if (overdisp_check$overdispersion) {
  cat("\n→ SOBREDISPERSIÓN DETECTADA. Usando Binomial Negativa.\n")

  modelo_nb <- glm.nb(
    abundancia_colémbolos ~ Bioestimulante * Dosis * Tipo_suelo,
    data = datos_edafo
  )

  cat("\nParámetro theta (Binomial Negativa):", round(modelo_nb$theta, 2), "\n")
  cat("Theta pequeño = mayor sobredispersión\n")

  modelo_final <- modelo_nb
  familia_usada <- "Binomial Negativa"
} else {
  cat("\n→ No hay sobredispersión significativa. Usando Poisson.\n")
  modelo_final <- modelo_poisson
  familia_usada <- "Poisson"
}

cat("\nFamilia utilizada:", familia_usada, "\n")

# ============================================================================
# 3. RESUMEN DEL MODELO
# ============================================================================

cat("\n=== Resumen del Modelo GLM ===\n")
summary(modelo_final)

# Tabla ANOVA tipo III
cat("\nTabla ANOVA (Tipo III):\n")
Anova(modelo_final, type = "III", test.statistic = "Chisq")

# ============================================================================
# 4. EVALUACIÓN DE SUPUESTOS
# ============================================================================

cat("\n=== Evaluación de Supuestos ===\n")

# 4.1 Residuos de Devianza
cat("\n1. Residuos de Devianza:\n")
residuos_dev <- residuals(modelo_final, type = "deviance")

# Shapiro-Wilk para residuos
shapiro_test <- shapiro.test(residuos_dev)
cat("Shapiro-Wilk (normalidad de residuos): p =", round(shapiro_test$p.value, 4), "\n")

# 4.2 Gráficos de residuos
cat("\n2. Gráficos de residuos:\n")
par(mfrow = c(2, 2))
plot(modelo_final, main = "Residuos vs Ajustados")
plot(modelo_final, which = 2, main = "Q-Q Plot")
par(mfrow = c(1, 1))

# 4.3 DHARMa - Residuos simulados
cat("\n3. DHARMa - Residuos simulados:\n")
if (familia_usada == "Binomial Negativa") {
  # Para NB, usar simulación
  simulation_output <- simulateResiduals(fittedModel = modelo_final, n = 1000)

  cat("\nPrueba de uniformidad:\n")
  print(testUniformity(simulation_output))

  cat("\nPrueba de dispersión:\n")
  print(testDispersion(simulation_output))

  cat("\nPrueba de inflación de ceros:\n")
  print(testZeroInflation(simulation_output))

  # Guardar gráfico
  dharma_plot <- plot(simulation_output, main = "Residuos DHARMa - GLM NB")
  ggsave("03_Modelos_Frecuentistas/02_glm_dharma.png",
         plot = last_plot(), width = 7, height = 5, dpi = 300)
}

# 4.4 check_model de performance
cat("\n4. Diagnóstico con check_model:\n")
tryCatch({
  check_plot <- check_model(modelo_final, panel = TRUE)
  print(check_plot)
  ggsave("03_Modelos_Frecuentistas/02_glm_diagnostico.png",
         plot = check_plot, width = 10, height = 8, dpi = 300)
}, error = function(e) {
  cat("Nota: check_model puede tener limitaciones con GLM NB\n")
})

# ============================================================================
# 5. COMPARACIONES MÚLTIPLES CON EMMEANS
# ============================================================================

cat("\n=== Comparaciones Múltiples (emmeans) ===\n")

# 5.1 Efecto principal: Bioestimulante
cat("\n--- Efecto Principal: Bioestimulante ---\n")
emm_bio <- emmeans(modelo_final, ~ Bioestimulante, type = "response")

cat("\nMedias marginales (escala de respuesta):\n")
print(emm_bio)

cat("\nComparaciones pairwise (Tukey):\n")
pairs(emm_bio, adjust = "tukey", type = "response")

# 5.2 Efecto principal: Dosis
cat("\n--- Efecto Principal: Dosis ---\n")
emm_dosis <- emmeans(modelo_final, ~ Dosis, type = "response")
print(emm_dosis)
cat("\nComparación (Tukey):\n")
pairs(emm_dosis, adjust = "tukey", type = "response")

# 5.3 Efecto principal: Tipo de suelo
cat("\n--- Efecto Principal: Tipo de Suelo ---\n")
emm_suelo <- emmeans(modelo_final, ~ Tipo_suelo, type = "response")
print(emm_suelo)
cat("\nComparaciones pairwise (Tukey):\n")
pairs(emm_suelo, adjust = "tukey", type = "response")

# ============================================================================
# 6. INTERACCIONES
# ============================================================================

cat("\n=== Análisis de Interacciones ===\n")

# Verificar si las interacciones son significativas
cat("\n¿Interacciones significativas?\n")
Anova(modelo_final, type = "III", test.statistic = "Chisq")

# Si la interacción triple es significativa, analizarla
# Medias para Bioestimulante x Dosis
cat("\nMedias - Bioestimulante x Dosis:\n")
emm_bio_dosis <- emmeans(modelo_final, ~ Bioestimulante * Dosis, type = "response")
print(emm_bio_dosis)

# Comparaciones simples: Efecto de Bioestimulante DENTRO de cada Dosis
cat("\nEfecto de Bioestimulante DENTRO de cada Dosis:\n")
pairs(emm_bio_dosis, by = "Dosis", adjust = "tukey", type = "response")

# Efecto de Dosis DENTRO de cada Bioestimulante
cat("\nEfecto de Dosis DENTRO de cada Bioestimulante:\n")
pairs(emm_bio_dosis, by = "Bioestimulante", adjust = "tukey", type = "response")

# ============================================================================
# 7. VISUALIZACIÓN DE RESULTADOS
# ============================================================================

cat("\n=== Gráficos de Resultados ===\n")

# Gráfico de barras con letras de significancia
datos_summary <- datos_edafo %>%
  group_by(Bioestimulante, Dosis) %>%
  summarise(
    media = mean(abundancia_colémbolos),
    ee = sd(abundancia_colémbolos) / sqrt(n()),
    .groups = "drop"
  )

# Obtener letras de cld
cld_data <- as.data.frame(emm_bio_dosis)
cld_result <- cld(emm_bio_dosis, alpha = 0.05)
cld_data$letras <- cld_result$.groups

# Unir datos
plot_data <- left_join(datos_summary, cld_data %>% select(Bioestimulante, Dosis, letras),
                       by = c("Bioestimulante", "Dosis"))

plot_barras <- ggplot(plot_data, aes(x = Bioestimulante, y = media,
                                      fill = Bioestimulante)) +
  geom_col(alpha = 0.8, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = media - ee, ymax = media + ee),
                width = 0.3, position = position_dodge(0.7)) +
  geom_text(aes(label = letras, y = media + ee + 1),
            size = 5, fontface = "bold", position = position_dodge(0.7)) +
  facet_wrap(~ Dosis, labeller = labeller(Dosis = c("Baja" = "Dosis Baja", "Alta" = "Dosis Alta"))) +
  labs(title = "Abundancia de Colémbolos por Bioestimulante y Dosis",
       subtitle = "Medias ± EE. Letras diferentes indican diferencias significativas (p < 0.05)",
       x = "Bioestimulante",
       y = "Abundancia Media (colémbolos/m²)") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 10),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray80"),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  scale_fill_brewer(palette = "Greens")

print(plot_barras)
ggsave("03_Modelos_Frecuentistas/02_glm_bioestimulantes.png",
       plot = plot_barras, width = 9, height = 6, dpi = 300)

# Gráfico de interacción con puntos
plot_interaccion <- ggplot(datos_edafo, aes(x = Dosis, y = abundancia_colémbolos,
                                             color = Bioestimulante, group = Bioestimulante)) +
  stat_summary(fun = "mean", geom = "point", size = 4,
               position = position_dodge(0.3)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
               position = position_dodge(0.3)) +
  facet_wrap(~ Tipo_suelo, labeller = labeller(Tipo_suelo = c(
    "Franco" = "Suelo Franco",
    "Arcilloso" = "Suelo Arcilloso",
    "Arenoso" = "Suelo Arenoso"
  ))) +
  labs(title = "Efecto de Bioestimulante y Dosis por Tipo de Suelo",
       x = "Dosis Aplicada",
       y = "Abundancia de Colémbolos") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray80"),
    strip.text = element_text(face = "bold")
  ) +
  scale_color_brewer(palette = "Set1", name = "Bioestimulante")

print(plot_interaccion)
ggsave("03_Modelos_Frecuentistas/02_glm_interaccion.png",
       plot = plot_interaccion, width = 10, height = 6, dpi = 300)

# ============================================================================
# 8. RESUMEN FINAL
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("RESUMEN DEL ANÁLISIS GLM\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Diseño: Completamente al Azar factorial 4x2x3\n")
cat("Familia:", familia_usada, "\n")
cat("Función de enlace: log\n")
cat("Total observaciones:", nrow(datos_edafo), "\n\n")

cat("Efectos significativos (α = 0.05):\n")
anova_tab <- Anova(modelo_final, type = "III", test.statistic = "Chisq")
sig_effects <- rownames(anova_tab)[anova_tab$`Pr(>Chisq)` < 0.05]
if (length(sig_effects) > 0) {
  for (eff in sig_effects) {
    cat("  ✓", eff, "(p <", round(anova_tab[eff, "Pr(>Chisq)"], 4), ")\n")
  }
} else {
  cat("  Ningún efecto principal significativo\n")
}

cat("\nSupuestos:\n")
cat("- Sobredispersión:", ifelse(overdisp_check$overdispersion, "Detectada", "No detectada"), "\n")
cat("- Normalidad de residuos (Shapiro): p =", round(shapiro_test$p.value, 4), "\n")

cat("\n", rep("=", 70), "\n", sep = "")
cat("FIN DEL ANÁLISIS\n")
cat(rep("=", 70), "\n", sep = "")
