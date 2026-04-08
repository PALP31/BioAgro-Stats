# ============================================================================
# 03_lmm_trichoderma.R
# Modelos Mixtos (LMM y GLMM) para medidas repetidas en el tiempo
# Contexto: Ensayo de antagonismo in vitro de Trichoderma vs Botrytis
# ============================================================================
# Diseño: Medidas repetidas en placas Petri
# Factor fijo: Cepa de Trichoderma
# Factor temporal: Días de evaluación
# Efecto aleatorio: Placa Petri (medidas repetidas)
# Variable respuesta: Diámetro de colonia de Botrytis (mm)
# ============================================================================

library(tidyverse)
library(lme4)        # lmer, glmer
library(glmmTMB)     # GLMM con diversas familias
library(lmerTest)    # p-valores para lme4
library(emmeans)     # Comparaciones múltiples
library(performance) # Diagnóstico de modelos
library(DHARMa)      # Residuos simulados para GLMM
library(car)

set.seed(456)

# ============================================================================
# PARTE A: MODELO LINEAL MIXTO (LMM) - DATOS CONTINUOS
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE A: MODELO LINEAL MIXTO (LMM) - Diámetro de Botrytis\n")
cat(rep("=", 80), "\n\n", sep = "")

# ----------------------------------------------------------------------------
# 1. SIMULACIÓN DE DATOS LONGITUDINALES
# ----------------------------------------------------------------------------

# Parámetros del diseño
n_placas_por_cepa <- 8  # Repeticiones por cepa
cepas <- c("Control", "T_asperellum", "T_harzianum", "T_viride", "T_koningii")
dias_eval <- 1:7  # Evaluación durante 7 días

# Crear estructura de datos
datos_antagonismo <- expand.grid(
  Cepa = factor(cepas),
  Placa = factor(1:n_placas_por_cepa),
  Dia = dias_eval
)

# Asignar ID único a cada placa
datos_antagonismo <- datos_antagonismo %>%
  mutate(
    Placa_ID = factor(paste0("P", as.numeric(Cepa) - 1, "_", Placa))
  )

# Efectos fijos de las cepas (reducción del crecimiento de Botrytis)
efecto_cepa <- c(
  "Control" = 0,           # Sin antagonismo
  "T_asperellum" = -8,     # Fuerte antagonismo
  "T_harzianum" = -5,      # Moderado
  "T_viride" = -6,         # Moderado-fuerte
  "T_koningii" = -3        # Débil
)

# Efecto del día (crecimiento natural de Botrytis)
efecto_dia <- c(0, 5, 10, 14, 17, 19, 20)  # Crecimiento logístico simplificado

# Interacción Cepa x Día (las cepas efectivas retardan el crecimiento)
interaccion <- expand.grid(Cepa = cepas, Dia = dias_eval)
interaccion$efecto <- with(interaccion,
  ifelse(Cepa == "Control", 0,
  ifelse(Cepa == "T_asperellum" & Dia >= 4, -seq(0, 12, length.out = 7)[Dia],
  ifelse(Cepa == "T_harzianum" & Dia >= 5, -seq(0, 8, length.out = 7)[Dia],
  ifelse(Cepa == "T_viride" & Dia >= 4, -seq(0, 9, length.out = 7)[Dia],
  ifelse(Cepa == "T_koningii" & Dia >= 6, -seq(0, 4, length.out = 7)[Dia], 0)))))
)

# Efectos aleatorios por placa (variabilidad biológica)
n_total_placas <- length(cepas) * n_placas_por_cepa
efecto_placa <- rnorm(n_total_placas, 0, 2)
names(efecto_placa) <- levels(datos_antagonismo$Placa_ID)

# Simular diámetro de colonia de Botrytis
datos_antagonismo <- datos_antagonismo %>%
  mutate(
    efecto_c = efecto_cepa[as.character(Cepa)],
    efecto_d = efecto_dia[Dia],
    efecto_int = interaccion$efecto[match(paste(Cepa, Dia),
                                          paste(interaccion$Cepa, interaccion$Dia))],
    efecto_p = efecto_placa[as.character(Placa_ID)],
    diametro_botrytis = 10 + efecto_c + efecto_d + efecto_int + efecto_p + rnorm(n(), 0, 1.5)
  ) %>%
  mutate(
    diametro_botrytis = pmax(diametro_botrytis, 0)  # No negativo
  ) %>%
  select(-efecto_c, -efecto_d, -efecto_int, -efecto_p)

cat("\n=== Estructura de los Datos ===\n")
str(datos_antagonismo)

cat("\nResumen del diámetro de Botrytis por Cepa y Día:\n")
datos_antagonismo %>%
  group_by(Cepa, Dia) %>%
  summarise(
    n = n(),
    media = mean(diametro_botrytis),
    sd = sd(diametro_botrytis),
    .groups = "drop"
  ) %>%
  print(n = 35)

# ----------------------------------------------------------------------------
# 2. AJUSTE DEL MODELO LINEAL MIXTO
# ----------------------------------------------------------------------------

cat("\n=== Ajuste del Modelo Lineal Mixto (LMM) ===\n")

# Modelo con efectos fijos (Cepa, Día, interacción) y aleatorio (Placa_ID)
modelo_lmm <- lmer(
  diametro_botrytis ~ Cepa * Dia + (1 | Placa_ID),
  data = datos_antagonismo,
  REML = FALSE  # Usar ML para comparaciones de modelos
)

cat("\nResumen del modelo:\n")
summary(modelo_lmm)

cat("\nVarianza de efectos aleatorios:\n")
print(VarCorr(modelo_lmm))

# ANOVA tipo III con aproximación de Satterthwaite
cat("\nTabla ANOVA (Tipo III):\n")
anova_lmm <- Anova(modelo_lmm, type = "III", test.statistic = "F")
print(anova_lmm)

# ----------------------------------------------------------------------------
# 3. EVALUACIÓN DE SUPUESTOS DEL LMM
# ----------------------------------------------------------------------------

cat("\n=== Evaluación de Supuestos del LMM ===\n")

# 3.1 Normalidad de residuos
cat("\n1. Normalidad de residuos (Shapiro-Wilk):\n")
residuos_lmm <- residuals(modelo_lmm)
shapiro_lmm <- shapiro.test(residuos_lmm)
print(shapiro_lmm)

# 3.2 Homocedasticidad
cat("\n2. Homocedasticidad (gráfico):\n")
plot(fitted(modelo_lmm), residuos_lmm,
     xlab = "Valores ajustados", ylab = "Residuos",
     main = "Residuos vs Ajustados - LMM")
abline(h = 0, col = "red", lwd = 2)

# 3.3 check_model de performance
cat("\n3. Diagnóstico integral (check_model):\n")
check_lmm <- check_model(modelo_lmm, panel = TRUE)
print(check_lmm)
ggsave("03_Modelos_Frecuentistas/03_lmm_diagnostico.png",
       plot = check_lmm, width = 10, height = 8, dpi = 300)

# 3.4 check_collinearity
cat("\n4. Multicolinealidad:\n")
check_collinearity(modelo_lmm)

# 3.5 check_outliers
cat("\n5. Detección de outliers:\n")
check_outliers(modelo_lmm)

# ----------------------------------------------------------------------------
# 4. COMPARACIONES MÚLTIPLES CON EMMEANS
# ----------------------------------------------------------------------------

cat("\n=== Comparaciones Múltiples (emmeans) ===\n")

# 4.1 Efecto principal: Cepa
cat("\n--- Efecto Principal: Cepa ---\n")
emm_cepa <- emmeans(modelo_lmm, ~ Cepa)
print(emm_cepa)

cat("\nComparaciones pairwise (Tukey):\n")
pairs(emm_cepa, adjust = "tukey")

# 4.2 Efecto principal: Día
cat("\n--- Efecto Principal: Día ---\n")
emm_dia <- emmeans(modelo_lmm, ~ Dia)
print(emm_dia)

# 4.3 Interacción Cepa x Día
cat("\n--- Interacción Cepa x Día ---\n")
emm_interaccion <- emmeans(modelo_lmm, ~ Cepa | Dia)

# Comparaciones de cepas DENTRO de cada día
cat("\nComparaciones de Cepas DENTRO de cada Día:\n")
for (d in c(3, 5, 7)) {  # Días clave
  cat("\nDía", d, ":\n")
  print(pairs(emmeans(modelo_lmm, ~ Cepa | Dia, at = list(Dia = d)), adjust = "tukey"))
}

# Comparaciones de días DENTRO de cada cepa
cat("\nComparaciones de Días DENTRO de cada Cepa:\n")
for (c in c("Control", "T_asperellum", "T_harzianum")) {
  cat("\nCepa:", c, "\n")
  emm_dia_cepa <- emmeans(modelo_lmm, ~ Dia | Cepa, at = list(Cepa = c))
  # Contrastar día 7 vs día 1
  print(contrast(emm_dia_cepa, method = "trt.vs.ctrl", adjust = "tukey"))
}

# Letras de significancia para el día 7 (última evaluación)
cat("\nLetras de significancia - Día 7:\n")
cld_result <- cld(emmeans(modelo_lmm, ~ Cepa | Dia, at = list(Dia = 7)),
                  alpha = 0.05, Letters = letters)
print(cld_result)

# ----------------------------------------------------------------------------
# 5. VISUALIZACIÓN DE RESULTADOS - LMM
# ----------------------------------------------------------------------------

cat("\n=== Gráficos de Resultados - LMM ===\n")

# Gráfico de trayectorias de crecimiento
plot_trayectorias <- ggplot(datos_antagonismo, aes(x = Dia, y = diametro_botrytis,
                                                    color = Cepa, group = Placa_ID)) +
  geom_line(alpha = 0.3, linewidth = 0.8) +
  stat_summary(fun = "mean", geom = "point", size = 3,
               position = position_dodge(0.5)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3,
               position = position_dodge(0.5)) +
  labs(title = "Crecimiento de Botrytis bajo Antagonismo de Trichoderma",
       subtitle = "Medidas repetidas en placas Petri (n = 8 por cepa)",
       x = "Días de incubación",
       y = "Diámetro de colonia (mm)") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(palette = "Set1", name = "Cepa")

print(plot_trayectorias)
ggsave("03_Modelos_Frecuentistas/03_lmm_trayectorias.png",
       plot = plot_trayectorias, width = 9, height = 6, dpi = 300)

# Gráfico de barras con letras para el día 7
datos_dia7 <- datos_antagonismo %>% filter(Dia == 7)

letras_dia7 <- as.data.frame(cld_result)

plot_barras_dia7 <- ggplot(datos_dia7, aes(x = Cepa, y = diametro_botrytis, fill = Cepa)) +
  stat_summary(fun = "mean", geom = "bar", alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_text(data = letras_dia7,
            aes(x = Cepa, y = max(datos_dia7$diametro_botrytis) + 3, label = .groups),
            size = 6, fontface = "bold") +
  labs(title = "Diámetro de Botrytis - Día 7",
       subtitle = "Letras diferentes indican diferencias significativas (Tukey, p < 0.05)",
       x = "Cepa de Trichoderma",
       y = "Diámetro (mm)") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_fill_brewer(palette = "Set1")

print(plot_barras_dia7)
ggsave("03_Modelos_Frecuentistas/03_lmm_dia7.png",
       plot = plot_barras_dia7, width = 8, height = 6, dpi = 300)

# ============================================================================
# PARTE B: GLMM PARA DATOS DE CONTEO O PROPORCIÓN
# ============================================================================

cat("\n\n", rep("=", 80), "\n", sep = "")
cat("PARTE B: GLMM - Inhibición del Crecimiento (Datos de Proporción)\n")
cat(rep("=", 80), "\n\n", sep = "")

# ----------------------------------------------------------------------------
# 1. SIMULACIÓN DE DATOS DE PROPORCIÓN
# ----------------------------------------------------------------------------

# Simular datos de inhibición: % de inhibición del crecimiento
# Esto es una proporción (0-100%) que requiere GLMM binomial o beta

datos_inhibicion <- datos_antagonismo %>%
  mutate(
    # Calcular inhibición como proporción del control
    diametro_max_control = max(diametro_botrytis[Cepa == "Control"]),
    inhibicion = pmax(0, pmin(1, (diametro_max_control - diametro_botrytis) / diametro_max_control)),
    # Convertir a "éxitos" sobre total (para binomial)
    exitos = round(inhibicion * 100),
    total = 100
  )

cat("\nDistribución de la inhibición:\n")
datos_inhibicion %>%
  group_by(Cepa) %>%
  summarise(
    media_inhibicion = mean(inhibicion),
    sd = sd(inhibicion),
    .groups = "drop"
  ) %>%
  print()

# ----------------------------------------------------------------------------
# 2. GLMM BINOMIAL
# ----------------------------------------------------------------------------

cat("\n=== Ajuste del GLMM Binomial ===\n")

# Usar glmmTMB para GLMM binomial con efectos aleatorios
modelo_glmm <- glmmTMB(
  cbind(exitos, total - exitos) ~ Cepa * Dia + (1 | Placa_ID),
  family = binomial(link = "logit"),
  data = datos_inhibicion
)

cat("\nResumen del GLMM:\n")
summary(modelo_glmm)

# ANOVA tipo III
cat("\nTabla ANOVA (Tipo III):\n")
Anova(modelo_glmm, type = "III", test.statistic = "Chisq")

# ----------------------------------------------------------------------------
# 3. EVALUACIÓN DE SUPUESTOS - GLMM
# ----------------------------------------------------------------------------

cat("\n=== Evaluación de Supuestos - GLMM ===\n")

# DHARMa para GLMM
cat("\n1. DHARMa - Residuos simulados:\n")
simulation_glmm <- simulateResiduals(modelo_glmm, n = 1000)

cat("\nPrueba de uniformidad:\n")
print(testUniformity(simulation_glmm))

cat("\nPrueba de dispersión:\n")
print(testDispersion(simulation_glmm))

# Gráfico de residuos DHARMa
plot(simulation_glmm, main = "Residuos DHARMa - GLMM Binomial")
ggsave("03_Modelos_Frecuentistas/03_glmm_dharma.png",
       plot = last_plot(), width = 7, height = 5, dpi = 300)

# check_model
cat("\n2. Diagnóstico con check_model:\n")
check_glmm <- check_model(modelo_glmm, panel = TRUE)
print(check_glmm)
ggsave("03_Modelos_Frecuentistas/03_glmm_diagnostico.png",
       plot = check_glmm, width = 10, height = 8, dpi = 300)

# ----------------------------------------------------------------------------
# 4. COMPARACIONES MÚLTIPLES - GLMM
# ----------------------------------------------------------------------------

cat("\n=== Comparaciones Múltiples - GLMM ===\n")

# emmeans para GLMM (en escala de respuesta = probabilidad)
emm_cepa_glmm <- emmeans(modelo_glmm, ~ Cepa, type = "response")

cat("\nMedias marginales (escala de respuesta - probabilidad):\n")
print(emm_cepa_glmm)

cat("\nComparaciones pairwise (Tukey):\n")
pairs(emm_cepa_glmm, adjust = "tukey", type = "response")

# Comparaciones por día
cat("\nComparaciones por Día (Día 7):\n")
emm_dia7_glmm <- emmeans(modelo_glmm, ~ Cepa, at = list(Dia = 7), type = "response")
print(emm_dia7_glmm)
pairs(emm_dia7_glmm, adjust = "tukey", type = "response")

# ============================================================================
# RESUMEN FINAL
# ============================================================================

cat("\n\n", rep("=", 80), "\n", sep = "")
cat("RESUMEN DEL ANÁLISIS DE MODELOS MIXTOS\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("DISEÑO EXPERIMENTAL:\n")
cat("- Tipo: Medidas repetidas en placas Petri\n")
cat("- Factor fijo: Cepa de Trichoderma (5 niveles)\n")
cat("- Factor temporal: Días de evaluación (7 días)\n")
cat("- Efecto aleatorio: Placa Petri (medidas repetidas)\n")
cat("- Total observaciones:", nrow(datos_antagonismo), "\n")
cat("- Total placas:", length(unique(datos_antagonismo$Placa_ID)), "\n\n")

cat("RESULTADOS LMM:\n")
cat("- Efecto Cepa:",
    ifelse(anova_lmm$`Pr(>F)`[1] < 0.05, "SIGNIFICATIVO", "No significativo"),
    "(F =", round(anova_lmm$`F value`[1], 2), ", p <",
    ifelse(anova_lmm$`Pr(>F)`[1] < 0.001, "0.001", round(anova_lmm$`Pr(>F)`[1], 4)), ")\n", sep = "")
cat("- Efecto Día:",
    ifelse(anova_lmm$`Pr(>F)`[2] < 0.05, "SIGNIFICATIVO", "No significativo"), "\n", sep = "")
cat("- Interacción Cepa x Día:",
    ifelse(anova_lmm$`Pr(>F)`[3] < 0.05, "SIGNIFICATIVA", "No significativa"), "\n\n", sep = "")

cat("COMPONENTES DE VARIANZA:\n")
var_comp <- as.data.frame(VarCorr(modelo_lmm))
cat("- Varianza entre placas:", round(var_comp$vcov[1], 2), "\n")
cat("- Varianza residual:", round(var_comp$vcov[2], 2), "\n")
cat("- ICC (coeficiente de correlación intraclase):",
    round(var_comp$vcov[1] / (var_comp$vcov[1] + var_comp$vcov[2]), 3), "\n\n")

cat("SUPUESTOS:\n")
cat("- Normalidad de residuos (Shapiro): p =", round(shapiro_lmm$p.value, 4), "\n")
cat("- DHARMa GLMM: uniformidad p =",
    round(testUniformity(simulation_glmm)$p.value, 4), "\n")

cat("\n", rep("=", 80), "\n", sep = "")
cat("FIN DEL ANÁLISIS\n")
cat(rep("=", 80), "\n", sep = "")
