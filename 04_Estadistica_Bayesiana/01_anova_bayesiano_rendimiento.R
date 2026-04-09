# ============================================================================
# 01_anova_bayesiano_rendimiento.R
# Análisis Bayesiano de Comparación de Grupos (ANOVA)
# ============================================================================
# En este script evaluamos el efecto del estrés biótico en el rendimiento.
# Introducimos métricas modernas: ROPE (Región de Equivalencia Práctica)
# y Probability of Direction (pd), superando el p-valor tradicional.
# ============================================================================

library(brms)        # Motor bayesiano
library(tidybayes)   # Manipulación de la posterior
library(bayestestR)  # Métricas ROPE y pd
library(tidyverse)
library(see)         # Visualización de bayestestR

set.seed(123)

# 1. SIMULACIÓN DE DATOS (Biomasa bajo estrés)
datos_biom <- data.frame(
  trata = factor(rep(c("Control", "Estres"), each = 20)),
  biomasa = c(rnorm(20, 45, 4), rnorm(20, 38, 5))
)

# 2. AJUSTE DEL MODELO BAYESIANO
# Un ANOVA es equivalente a un modelo lineal con predictores categóricos.
mod_rend <- brm(
  biomasa ~ trata, 
  data = datos_biom,
  family = gaussian(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 123
)

# 3. MÉTRICAS BAYESIANAS MODERNAS (bayestestR)
cat("\n--- MÉTRICAS DE INFERENCIA BAYESIANA ---\n")
post_summary <- describe_posterior(
  mod_rend, 
  rope_range = c(-1, 1), # Rango de nulidad práctica (±1g de biomasa)
  test = c("p_direction", "rope")
)
print(post_summary)

# EXPLICACIÓN:
# - Probability of Direction (pd): Probabilidad de que el efecto sea mayor/menor a 0.
#   Un pd > 97.5% suele considerarse evidencia fuerte (similar a p < 0.05).
# - ROPE: Si la distribución posterior cae mayoritariamente fuera de la ROPE,
#   el efecto es significativamente práctico, no solo estadístico.

# 4. VISUALIZACIÓN (tidybayes: stat_halfeye)
ggplot(datos_biom %>% add_epred_draws(mod_rend), aes(x = .epred, y = trata, fill = trata)) +
  stat_halfeye(alpha = 0.6) +
  geom_vline(xintercept = 45, linetype = "dashed", color = "red") +
  labs(title = "Distribución Posterior de Biomasa Estimada",
       subtitle = "Punto = Mediana, Línea = HDI 95%. Línea roja = Media teórica Control.",
       x = "Biomasa Estimada (g)", y = "Tratamiento") +
  theme_minimal() +
  scale_fill_flat()

cat("\n--- Script 01 finalizado correctamente ---\n")
