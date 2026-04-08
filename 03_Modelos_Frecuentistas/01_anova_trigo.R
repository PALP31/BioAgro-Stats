# ============================================================================
# 01_anova_trigo.R
# ANOVA de 2 vías con Bloques Completos al Azar (BCA)
# Contexto: Biomasa de trigo bajo estrés salino y térmico
# ============================================================================
# Diseño: Factorial 2x2 con bloques
# Factores: Estrés salino (NaCl) x Estrés térmico (Temperatura)
# Variable respuesta: Biomasa aérea (g)
# ============================================================================

library(tidyverse)
library(emmeans)
library(performance)
library(car)
library(lmtest)

set.seed(42)

# ============================================================================
# 1. SIMULACIÓN DE DATOS
# ============================================================================

# Parámetros del diseño
n_bloques <- 5
niveles_salinidad <- c("0_mM", "50_mM", "100_mM", "150_mM")
niveles_temperatura <- c("25C", "32C", "38C")

# Crear diseño factorial con bloques
datos_trigo <- expand.grid(
  Bloque = factor(paste0("B", 1:n_bloques)),
  Salinidad = factor(niveles_salinidad),
  Temperatura = factor(niveles_temperatura)
)

# Efectos de los tratamientos
efecto_salinidad <- c("0_mM" = 0, "50_mM" = -3, "100_mM" = -8, "150_mM" = -15)
efecto_temperatura <- c("25C" = 0, "32C" = -4, "38C" = -10)

# Interacción: el estrés combinado es más severo
interaccion <- expand.grid(Salinidad = niveles_salinidad, Temperatura = niveles_temperatura)
interaccion$efecto <- with(interaccion,
  ifelse(Salinidad == "150_mM" & Temperatura == "38C", -8,
  ifelse(Salinidad == "100_mM" & Temperatura == "38C", -5,
  ifelse(Salinidad == "150_mM" & Temperatura == "32C", -4, 0)))
)

# Efecto del bloque (variabilidad espacial en el invernadero)
efecto_bloque <- rnorm(n_bloques, 0, 1.5)
names(efecto_bloque) <- paste0("B", 1:n_bloques)

# Simular biomasa
datos_trigo <- datos_trigo %>%
  mutate(
    efecto_sal = efecto_salinidad[as.character(Salinidad)],
    efecto_temp = efecto_temperatura[as.character(Temperatura)],
    efecto_int = interaccion$efecto[match(paste(Salinidad, Temperatura),
                                          paste(interaccion$Salinidad, interaccion$Temperatura))],
    efecto_bloq = efecto_bloque[as.character(Bloque)],
    biomasa = 45 + efecto_sal + efecto_temp + efecto_int + efecto_bloq + rnorm(n(), 0, 2.5)
  ) %>%
  select(-efecto_sal, -efecto_temp, -efecto_int, -efecto_bloq)

cat("\n=== Estructura de los Datos ===\n")
str(datos_trigo)
cat("\nResumen de biomasa por tratamiento:\n")
datos_trigo %>%
  group_by(Salinidad, Temperatura) %>%
  summarise(
    n = n(),
    media = mean(biomasa),
    sd = sd(biomasa),
    .groups = "drop"
  ) %>%
  print()

# ============================================================================
# 2. ANÁLISIS DE VARIANZA (ANOVA)
# ============================================================================

cat("\n=== ANOVA de 2 vías con Bloques ===\n")

# Ajustar modelo lineal
modelo_anova <- lm(biomasa ~ Bloque + Salinidad * Temperatura, data = datos_trigo)

# Tabla ANOVA tipo III (más apropiada para diseños factoriales)
cat("\nTabla ANOVA (Tipo III):\n")
Anova(modelo_anova, type = "III")

# Resumen completo del modelo
cat("\nResumen del modelo:\n")
summary(modelo_anova)

# ============================================================================
# 3. EVALUACIÓN DE SUPUESTOS
# ============================================================================

cat("\n=== Evaluación de Supuestos ===\n")

# 3.1 Normalidad de residuos (Shapiro-Wilk)
cat("\n1. Shapiro-Wilk - Normalidad de residuos:\n")
shapiro_test <- shapiro.test(residuals(modelo_anova))
print(shapiro_test)

# 3.2 Homocedasticidad (Levene)
cat("\n2. Levene - Homocedasticidad:\n")
# Para interacción
datos_trigo$grupo <- interaction(datos_trigo$Salinidad, datos_trigo$Temperatura)
levene_test <- leveneTest(biomasa ~ grupo, data = datos_trigo)
print(levene_test)

# 3.3 Diagnóstico gráfico con performance
cat("\n3. Diagnóstico gráfico (check_model):\n")
check_plot <- check_model(modelo_anova, panel = TRUE)
print(check_plot)
ggsave("03_Modelos_Frecuentistas/01_anova_diagnostico.png",
       plot = check_plot, width = 10, height = 8, dpi = 300)

# ============================================================================
# 4. COMPARACIONES MÚLTIPLES CON EMMEANS
# ============================================================================

cat("\n=== Comparaciones Múltiples (emmeans) ===\n")

# Medias marginales estimadas para Salinidad
cat("\nMedias marginales - Factor Salinidad:\n")
emm_sal <- emmeans(modelo_anova, ~ Salinidad)
print(emm_sal)

cat("\nComparaciones pairwise - Salinidad (Tukey):\n")
pairs(emm_sal, adjust = "tukey")

# Medias marginales estimadas para Temperatura
cat("\nMedias marginales - Factor Temperatura:\n")
emm_temp <- emmeans(modelo_anova, ~ Temperatura)
print(emm_temp)

cat("\nComparaciones pairwise - Temperatura (Tukey):\n")
pairs(emm_temp, adjust = "tukey")

# ============================================================================
# 5. INTERACCIÓN SALINIDAD x TEMPERATURA
# ============================================================================

cat("\n=== Análisis de la Interacción ===\n")

# Medias para cada combinación
cat("\nMedias para cada combinación Salinidad x Temperatura:\n")
emm_interaccion <- emmeans(modelo_anova, ~ Salinidad * Temperatura)
print(emm_interaccion)

# Comparaciones dentro de cada nivel de Temperatura
cat("\nComparaciones de Salinidad DENTRO de cada Temperatura:\n")
pairs(emm_interaccion, by = "Temperatura", adjust = "tukey")

# Comparaciones de Temperatura DENTRO de cada nivel de Salinidad
cat("\nComparaciones de Temperatura DENTRO de cada Salinidad:\n")
pairs(emm_interaccion, by = "Salinidad", adjust = "tukey")

# Contrastar interacción específica: ¿El efecto de 150mM NaCl es diferente a 38C vs 25C?
cat("\nContraste de interacción específica:\n")
contrast(emm_interaccion,
         interaction = "pairwise",
         adjust = "tukey")

# ============================================================================
# 6. VISUALIZACIÓN DE RESULTADOS
# ============================================================================

cat("\n=== Gráficos de Resultados ===\n")

# Gráfico de interacción
plot_interaccion <- ggplot(datos_trigo, aes(x = Salinidad, y = biomasa,
                                             color = Temperatura, group = Temperatura)) +
  stat_summary(fun = "mean", geom = "point", size = 4, position = position_dodge(0.3)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
               position = position_dodge(0.3)) +
  labs(title = "Biomasa de Trigo bajo Estrés Salino y Térmico",
       subtitle = "Medias ± Error Estándar",
       x = "Nivel de Salinidad (NaCl)",
       y = "Biomasa Aérea (g)") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(plot_interaccion)
ggsave("03_Modelos_Frecuentistas/01_anova_interaccion.png",
       plot = plot_interaccion, width = 8, height = 6, dpi = 300)

# Gráfico de letras de significancia
letras_data <- as.data.frame(emm_sal)
letras_data$letras <- cld(emm_sal, alpha = 0.05)$`.groups`

plot_salinidad <- ggplot(datos_trigo, aes(x = Salinidad, y = biomasa, fill = Salinidad)) +
  stat_summary(fun = "mean", geom = "bar", alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_text(data = letras_data,
            aes(x = Salinidad, y = max(datos_trigo$biomasa) + 2, label = letras),
            size = 6, fontface = "bold") +
  labs(title = "Efecto del Estrés Salino",
       x = "Salinidad (mM NaCl)",
       y = "Biomasa Aérea (g)") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_fill_brewer(palette = "Greens")

print(plot_salinidad)
ggsave("03_Modelos_Frecuentistas/01_anova_salinidad.png",
       plot = plot_salinidad, width = 7, height = 5, dpi = 300)

# ============================================================================
# 7. RESUMEN FINAL
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("RESUMEN DEL ANÁLISIS\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Diseño: Bloques Completos al Azar (BCA) factorial 4x3\n")
cat("Bloques: ", n_bloques, "\n", sep = "")
cat("Tratamientos: ", length(niveles_salinidad), " niveles de salinidad x ",
    length(niveles_temperatura), " niveles de temperatura\n", sep = "")
cat("Total observaciones: ", nrow(datos_trigo), "\n\n", sep = "")

cat("Resultados principales:\n")
cat("- Efecto Salinidad: ",
    ifelse(Anova(modelo_anova, type = "III")$`Pr(>F)`[2] < 0.05, "SIGNIFICATIVO", "No significativo"), "\n", sep = "")
cat("- Efecto Temperatura: ",
    ifelse(Anova(modelo_anova, type = "III")$`Pr(>F)`[3] < 0.05, "SIGNIFICATIVO", "No significativo"), "\n", sep = "")
cat("- Interacción S x T: ",
    ifelse(Anova(modelo_anova, type = "III")$`Pr(>F)`[4] < 0.05, "SIGNIFICATIVA", "No significativa"), "\n\n", sep = "")

cat("Supuestos:\n")
cat("- Normalidad (Shapiro-Wilk): p =", round(shapiro_test$p.value, 4), "\n")
cat("- Homocedasticidad (Levene): p =", round(levene_test$`Pr(>F)`[1], 4), "\n")

cat("\n", rep("=", 70), "\n", sep = "")
cat("FIN DEL ANÁLISIS\n")
cat(rep("=", 70), "\n", sep = "")
