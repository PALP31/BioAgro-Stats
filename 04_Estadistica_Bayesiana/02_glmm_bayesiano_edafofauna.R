# ============================================================================
# 02_glmm_bayesiano_edafofauna.R
# Modelo Mixto Bayesiano (Binomial Negativa) para Edafofauna
# ============================================================================
# En este script abordamos el modelado de conteos con sobredispersión y 
# estructura jerárquica (Parcelas como efecto aleatorio).
# Incluye el uso de Priors informados y diagnóstico MCMC exhaustivo.
# ============================================================================

library(brms)        # Motor bayesiano
library(bayesplot)   # Gráficos diagnósticos (trace, density, etc.)
library(performance) # pp_check e indicadores de convergencia
library(tidyverse)

set.seed(321)

# 1. SIMULACIÓN DE DATOS (Conteos de insectos en parcelas)
n_parcelas <- 20
datos_fauna <- data.frame(
  parcela = factor(paste0("P", 1:n_parcelas)),
  tratamiento = factor(rep(c("Convencional", "Ecologico"), each = 10))
) %>%
  mutate(
    # Intercepto aleatorio por parcela (variabilidad del suelo)
    ef_parc = rnorm(n_parcelas, 0, 0.5)[as.numeric(parcela)],
    mu = exp(2.5 + (tratamiento == "Ecologico")*0.8 + ef_parc),
    abundancia = rnbinom(n(), mu = mu, size = 1.2)
  )

# 2. DECLARACIÓN DE PRIORS (Priors débilmente informativos)
# ¿Por qué usarlos? Para orientar la convergencia y evitar valores biológicamente 
# imposibles, manteniendo flexibilidad en el modelo.

priors_fauna <- c(
  prior(normal(2.5, 5), class = "Intercept"),
  prior(normal(0, 2), class = "b"),        # Efecto del tratamiento
  prior(cauchy(0, 1), class = "sd"),       # Varianza aleatoria (siempre positiva)
  prior(gamma(0.01, 0.01), class = "shape") # Parámetro de dispersión
)

# 3. AJUSTE DEL MODELO GLMM BAYESIANO
mod_fauna <- brm(
  abundancia ~ tratamiento + (1 | parcela),
  data = datos_fauna,
  family = negbinomial(),
  prior = priors_fauna,
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 123
)

# 4. VERIFICACIÓN DE CONVERGENCIA (MCMC Diagnostics)
cat("\n--- DIAGNÓSTICO DE CONVERGENCIA ---")
cat("\n[1] Indicador Rhat (debe ser < 1.05):\n")
print(summary(mod_fauna)$fixed)

cat("\n[2] Effective Sample Size (ESS):\n")
cat("Un ESS bajo sugiere problemas de autocorrelación en las cadenas.")

# Visualización de Trazas (Trace plots)
# ¿Cómo leerlo? Debe parecer 'césped' o 'una oruga difusa' bien mezclada.
mcmc_plot(mod_fauna, type = "trace", pars = "^b_")

# 5. POSTERIOR PREDICTIVE CHECK (pp_check)
# ¿Qué hace? Simula datos del modelo y los compara con los observados.
cat("\n--- POSTERIOR PREDICTIVE CHECK ---")
pp_check(mod_fauna, ndraws = 50) +
  labs(title = "PP-Check: Abundancia de Edafofauna",
       subtitle = "La densidad negra (observada) debe encajar con las líneas azules (simuladas).") +
  theme_minimal()

cat("\n--- Script 02 finalizado correctamente ---\n")
