# ============================================================================
# 03_lmm_trichoderma.R (LMM: ANTAGONISMO IN VITRO)
# Modelos Mixtos para Medidas Repetidas en el Tiempo
# ============================================================================
library(lme4)        # Motor de modelos mixtos
library(lmerTest)    # P-valores por Satterthwaite
library(tidyverse)
library(performance) # Diagnóstico de modelos mixtos
library(emmeans)     # Comparaciones temporales y por cepa

set.seed(321)

# 1. SIMULACIÓN DE DATOS: DIÁMETRO DE BOTRYTIS (mm)
# Factores: Trichoderma (Control, T_harzianum, T_viride), Tiempo (Día 1 al 7)
# Aleatorio: Placa_Petri (Cada placa tiene su propia variabilidad)
cepas <- c("Control", "T_harzianum", "T_viride")
placas <- paste0("P", 1:12)

datos_tricho <- expand.grid(
  day = factor(c(1, 3, 5, 7)),
  strain = factor(cepas),
  plate = factor(placas)
) %>%
  # Cada placa solo tiene una cepa (diseño anidado)
  mutate(strain = rep(factor(cepas), each = 4, length.out = n())) %>%
  mutate(
    # Crecimiento lineal base
    growth_base = as.numeric(as.character(day)) * case_when(
      strain == "Control" ~ 12,
      strain == "T_harzianum" ~ 4,  # Fuerte inhibición
      strain == "T_viride" ~ 7
    ),
    # Efecto aleatorio de placa (intercepto aleatorio)
    ef_plate = rnorm(length(placas), 0, 2)[as.numeric(plate)],
    diameter = growth_base + ef_plate + rnorm(n(), 0, 1.5)
  )

# 2. AJUSTE DEL MODELO MIXTO (LMM)
# Formula: Diámetro ~ Tiempo * Cepa + (1|Placa)
mod_lmm <- lmer(diameter ~ day * strain + (1|plate), data = datos_tricho)

cat("\n--- [1] TABLA ANOVA (TIPO III) DEL MODELO MIXTO ---\n")
print(anova(mod_lmm))

# 3. EVALUACIÓN DE SUPUESTOS
cat("\n--- [2] DIAGNÓSTICO DEL MODELO MIXTO ---\n")
print(check_model(mod_lmm))

# 4. COMPARACIONES MÚLTIPLES (POR DÍA)
cat("\n--- [3] EFECTO DE LAS CEPAS DENTRO DE CADA DÍA ---\n")
emm_lmm <- emmeans(mod_lmm, ~ strain | day)
print(pairs(emm_lmm, adjust = "tukey"))

# 5. GRÁFICO PROFESIONAL: CURVAS DE CRECIMIENTO
# Extraemos medias estimadas para graficar
plot_data <- as.data.frame(emmeans(mod_lmm, ~ strain * day))

ggplot(plot_data, aes(x = day, y = emmean, color = strain, group = strain)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2) +
  labs(title = "Crecimiento de Botrytis cinerea vs Trichoderma spp.",
       subtitle = "Curvas de crecimiento estimadas por LMM (Medias ± SE)",
       x = "Días después de inoculación", y = "Diámetro de colonia (mm)",
       color = "Cepa de Trichoderma") +
  theme_classic() +
  scale_color_brewer(palette = "Set1") +
  theme(plot.title = element_text(face = "bold"))

cat("\n--- Script 03 finalizado correctamente ---\n")
