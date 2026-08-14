# ============================================================================
# 04_medidas_repetidas_tiempo.R — MEDIDAS REPETIDAS EN EL TIEMPO (LMM)
# Análisis Longitudinal de Respuestas Fisiológicas en Cultivos
# ============================================================================
# Problema Biológico: Medir la misma planta/maceta a través del tiempo
# genera autocorrelación temporal y viola el supuesto de independencia del ANOVA clásico.
# Solución: Modelos Lineales Mixtos con intercepto aleatorio por sujeto (1 | ID).
# ============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
if (requireNamespace("easyModels", quietly = TRUE)) library(easyModels)

set.seed(314)

# 1. SIMULACIÓN DE DATOS LONGITUDINALES (DINÁMICA DE CRECIMIENTO)
# Factores: Tratamiento (Control, Sequía), Tiempo (D0, D7, D14, D21, D28)
# Unidades: 20 plantas (10 por tratamiento) evaluadas en 5 fechas (100 observaciones)

n_plantas <- 20
tiempos <- c("D0", "D7", "D14", "D21", "D28")

datos_sujetos <- data.frame(
  ID = factor(paste0("P_", sprintf("%02d", 1:n_plantas))),
  Tratamiento = factor(rep(c("Control", "Sequia"), each = 10)),
  Vigor_Basal = rnorm(n_plantas, mean = 0, sd = 2.0)  # Variabilidad intrínseca de cada planta
)

datos_tiempo <- expand.grid(
  ID = datos_sujetos$ID,
  Tiempo = factor(tiempos, levels = tiempos)
) %>%
  left_join(datos_sujetos, by = "ID") %>%
  mutate(
    dias_num = as.numeric(gsub("D", "", as.character(Tiempo))),
    # Tasa de crecimiento diferencial según tratamiento
    tasa = ifelse(Tratamiento == "Control", 1.25, 0.65),
    altura_media = 15 + (tasa * dias_num) + Vigor_Basal,
    Altura = round(altura_media + rnorm(n(), mean = 0, sd = 1.2), 2)
  ) %>%
  arrange(ID, Tiempo)

print(head(datos_tiempo, 10))

# ============================================================================
# 2. EL ERROR CLÁSICO: ANOVA DE MEDIDAS INDEPENDIENTES (Pseudorreplicación)
# ============================================================================
cat("\n--- [1] ANOVA CLÁSICO (INCORRECTO: IGNORA EL EFECTO SUJETO) ---\n")
mod_naive <- aov(Altura ~ Tratamiento * Tiempo, data = datos_tiempo)
print(summary(mod_naive))
cat("-> Nota: Este modelo asume 100 plantas independientes cuando en realidad son 20 plantas medidas 5 veces.\n")

# ============================================================================
# 3. RUTA CORRECTA: MODELO LINEAL MIXTO (LMM con lme4)
# ============================================================================
cat("\n--- [2] MODELO LINEAL MIXTO CON INTERCEPTO ALEATORIO (1 | ID) ---\n")
mod_lmm <- lmer(Altura ~ Tratamiento * Tiempo + (1 | ID), data = datos_tiempo)
print(anova(mod_lmm))

# Comparación de tratamientos en cada punto temporal (Efectos Simples)
emm_tiempo <- emmeans(mod_lmm, ~ Tratamiento | Tiempo)
cat("\n--- [3] CONTRASTES DE TRATAMIENTO EN CADA DÍA (EMMEANS) ---\n")
print(pairs(emm_tiempo, adjust = "bonferroni"))

# ============================================================================
# 4. RUTA AUTOMATIZADA CON easyModels
# ============================================================================
if (requireNamespace("easyModels", quietly = TRUE)) {
  cat("\n--- [4] AJUSTE CON easyModels::analizar_medidas_repetidas ---\n")
  mod_easy_rep <- analizar_medidas_repetidas(
    datos = datos_tiempo,
    formula_fijos = Altura ~ Tratamiento * Tiempo,
    sujeto = "ID",
    diagnosticos = FALSE
  )
  
  # Gráfico de líneas longitudinales de publicación
  p_rep <- graficar_predichos(
    modelo = mod_easy_rep,
    predictor = "Tiempo",
    por = "Tratamiento",
    tipo_grafico = "lineas",
    paleta = "teal",
    titulo = "Dinámica de Crecimiento en Altura de Planta",
    eje_x = "Días de Evaluación",
    eje_y = "Altura Promedio (cm)"
  )
  print(p_rep)
}

# ============================================================================
# 5. GRÁFICO AVANZADO GGPLOT2: TRAYECTORIAS INDIVIDUALES + MEDIA
# ============================================================================
p_spaghetti <- ggplot(datos_tiempo, aes(x = dias_num, y = Altura, color = Tratamiento)) +
  # Líneas semitransparentes por cada planta individual (Spaghetti Plot)
  geom_line(aes(group = ID), alpha = 0.35, linewidth = 0.6) +
  # Tendencia promedio general
  stat_summary(fun = mean, geom = "line", linewidth = 1.6, aes(group = Tratamiento)) +
  stat_summary(fun = mean, geom = "point", size = 3.5, aes(group = Tratamiento)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 1.0, linewidth = 0.8) +
  scale_color_manual(values = c("Control" = "#00A88F", "Sequia" = "#E65100")) +
  scale_x_continuous(breaks = c(0, 7, 14, 21, 28), labels = tiempos) +
  labs(
    title = "Curvas Longitudinales de Altura (Spaghetti Plot + Media ± SE)",
    subtitle = "Las líneas delgadas representan plantas individuales seguidas en el tiempo",
    x = "Tiempo de Evaluación",
    y = "Altura de Planta (cm)",
    color = "Tratamiento"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_spaghetti)

cat("\n--- Script 04 Medidas Repetidas finalizado correctamente ---\n")
