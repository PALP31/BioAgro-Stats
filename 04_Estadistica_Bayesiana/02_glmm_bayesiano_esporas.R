# ============================================================================
# 02_glmm_bayesiano_esporas.R
# GLMM Bayesiano (Binomial Negativa) para Conteos en el Tiempo
# ============================================================================
# Contexto: Conteo de esporas (UFC/ml) de Trichoderma evaluadas en 4 tiempos
# (0, 7, 14, 21 días) usando la placa Petri como efecto aleatorio.
# ============================================================================

library(brms)        # Interfaz de alto nivel para Stan
library(tidyverse)
library(tidybayes)
library(performance) # Diagnóstico y pp_check visual

set.seed(321)

# 1. SIMULACIÓN DE DATOS (Conteos Sobredispersos y Medidas Repetidas)
n_placas <- 15
tiempos <- c(0, 7, 14, 21)

datos_esporas <- expand.grid(
  day = tiempos,
  plate = factor(paste0("P", 1:n_placas))
) %>%
  mutate(
    # Efecto lineal del tiempo (log-link)
    log_lambda = 2 + (day * 0.1) + rnorm(n_placas, 0, 0.4)[as.numeric(plate)],
    # Generar conteos con distribución Binomial Negativa (Maneja sobredispersión)
    spores = rnbinom(n(), mu = exp(log_lambda), size = 1.2)
  )

# 2. DEFINICIÓN DE PRIORS (Priors débilmente informativos)
# Especificamos priors normales para los coeficientes de regresión.
# Esto orienta al modelo pero permite que los datos hablen.

priors_esporas <- c(
  prior(normal(2, 5), class = "Intercept"),
  prior(normal(0, 2), class = "b"),        # Efecto del tiempo
  prior(exponential(1), class = "sd"),     # Varianza de efectos aleatorios
  prior(gamma(0.01, 0.01), class = "shape") # Parámetro de dispersión
)

# 3. AJUSTE DEL MODELO GLMM BAYESIANO
mod_bayes_esporas <- brm(
  spores ~ day + (1 | plate),
  data = datos_esporas,
  family = negbinomial(),  # Maneja mejor los ceros y la varianza que Poisson
  prior = priors_esporas,
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 321
)

# 4. VERIFICACIÓN DE CONVERGENCIA (MCMC Diagnostics)
cat("\n--- Verificando Convergencia (MCMC) ---\n")
# El Rhat debe ser < 1.05 en todos los parámetros.

cat("\nResumen del modelo:\n")
print(summary(mod_bayes_esporas))

# Visualización de trazas (Trace plots)
# plot(mod_bayes_esporas) # Abre gráficos de trazas en R

# 5. CHEQUEO PREDICTIVO POSTERIOR (pp_check)
# ¿Qué hace? Compara los datos observados con datos simulados por el modelo.
cat("\n--- Chequeo Predictivo Posterior (PP Check) ---\n")
plot_pp <- pp_check(mod_bayes_esporas, ndraws = 50) +
  labs(title = "PP Check: Modelo de Esporas (Trichoderma)",
       subtitle = "Densidad observada (línea gruesa) vs simulada (líneas finas)") +
  theme_minimal()

print(plot_pp)

# 6. VISUALIZACIÓN DE TENDENCIAS (Tidybayes)
plot_temporal <- datos_esporas %>%
  group_by(day) %>%
  add_epred_draws(mod_bayes_esporas, ndraws = 100) %>%
  ggplot(aes(x = day, y = spores)) +
  stat_lineribbon(aes(y = .epred), .width = c(.95, .80, .50), alpha = 0.5) +
  geom_point(alpha = 0.4) +
  scale_fill_brewer(palette = "Greens", name = "Credibilidad") +
  labs(title = "Proyección Bayesiana de Esporulación",
       subtitle = "Curva ajustada por GLMM NegBinomial con Incertidumbre del 95%",
       x = "Tiempo (Días)", y = "Conteo de Esporas (UFC/ml)") +
  theme_minimal()

print(plot_temporal)

cat("\n--- Script 02 finalizado correctamente ---\n")
