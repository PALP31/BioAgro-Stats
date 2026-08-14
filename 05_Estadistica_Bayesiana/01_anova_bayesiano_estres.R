# ============================================================================
# 01_anova_bayesiano_estres.R
# ANOVA Bayesiano de 2 Vías con brms y tidybayes
# ============================================================================
# Contexto: Biomasa de trigo bajo estrés salino (Control, 150mM)
# y aplicación de bioestimulantes (Control, Bio_A).
# ============================================================================

library(brms)        # Interfaz de alto nivel para Stan
library(tidybayes)   # Manipulación y visualización de la distribución posterior
library(tidyverse)   # Manipulación de datos y ggplot2

set.seed(123)

# 1. SIMULACIÓN DE DATOS (2 Factores: Salinidad x Bioestimulante)
n_rep <- 10
datos_sal <- expand.grid(
  sal = factor(c("Control", "Salino")),
  bio = factor(c("Control", "Bio_A")),
  rep = 1:n_rep
) %>%
  mutate(
    yield_base = case_when(
      sal == "Control" & bio == "Control" ~ 35,
      sal == "Salino" & bio == "Control" ~ 22,
      sal == "Control" & bio == "Bio_A" ~ 38,
      sal == "Salino" & bio == "Bio_A" ~ 30  # Mitigación del estrés
    ),
    biomasa = yield_base + rnorm(n(), 0, 2.5)
  )

# 2. AJUSTE DEL MODELO BAYESIANO (brms)
# Usamos priors predeterminados para este ejemplo
mod_bayes_sal <- brm(
  biomasa ~ sal * bio,
  data = datos_sal,
  family = gaussian(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 123
)

cat("\n--- Resumen del Modelo Bayesiano (brms) ---\n")
print(summary(mod_bayes_sal))

# 3. EXTRACCIÓN DE LA DISTRIBUCIÓN POSTERIOR (tidybayes)
cat("\n--- Extrayendo Muestras Posteriores y HDI 95% ---\n")

# Obtener medias estimadas por combinación (Epred)
pred_posteriores <- datos_sal %>%
  group_by(sal, bio) %>%
  add_epred_draws(mod_bayes_sal)

# Resumen con HDI (Highest Density Interval)
resumen_hdi <- pred_posteriores %>%
  group_by(sal, bio) %>%
  median_hdi(.epred, .width = 0.95)

print(resumen_hdi)

# 4. VISUALIZACIÓN PROFESIONAL (Half-Eye Plot)
# ¿Para qué sirve? Muestra la densidad (forma), la mediana (punto)
# y el intervalo de credibilidad (línea).

plot_halfeye <- pred_posteriores %>%
  ggplot(aes(x = .epred, y = sal, fill = bio)) +
  stat_halfeye(alpha = 0.6, position = position_dodge(width = 0.5)) +
  labs(title = "Distribución Posterior: Biomasa de Trigo",
       subtitle = "Diseño Bayesiano de 2 Vías (Salinidad x Bioestimulante). Punto = Mediana, Línea = HDI 95%",
       x = "Biomasa Estimada (g)", y = "Nivel de Salinidad",
       fill = "Bioestimulante") +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "#2E86AB", "Bio_A" = "#228B22"))

print(plot_halfeye)

cat("\n--- Interpretación Bayesiana ---\n")
cat("Si el intervalo HDI del 95% de la diferencia entre grupos NO contiene el 0,
podemos afirmar con un 95% de probabilidad que existe un efecto del bioestimulante.\n")
