# =============================================================================
# Script: 02_modelos_bayesianos.R
# Descripción: Plantilla para análisis bayesiano de datos experimentales
#              biológicos / agronómicos con brms (Stan backend).
#              Modelos: regresión lineal, ANOVA bayesiano, GLMM bayesiano.
# Autor: [nombre]
# Fecha: [fecha]
# =============================================================================

# ── 0. Paquetes ──────────────────────────────────────────────────────────────
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,   # manipulación y visualización
  brms,        # modelos bayesianos via Stan
  bayesplot,   # gráficos de diagnóstico MCMC
  tidybayes,   # extracción y visualización de distribuciones posteriores
  ggdist,      # distribuciones para ggplot2
  posterior,   # resumen de draws
  loo          # comparación de modelos (LOO-CV)
)

# ── 1. Datos de ejemplo ──────────────────────────────────────────────────────
set.seed(42)
datos <- data.frame(
  tratamiento = rep(c("Control", "Dosis_baja", "Dosis_alta"), each = 20),
  bloque      = factor(rep(1:4, times = 15)),
  rendimiento = c(
    rnorm(20, mean = 5.0, sd = 0.8),
    rnorm(20, mean = 6.2, sd = 0.9),
    rnorm(20, mean = 7.5, sd = 1.0)
  ),
  conteo      = rpois(60, lambda = rep(c(3, 6, 10), each = 20))
)
datos$tratamiento <- factor(datos$tratamiento,
                            levels = c("Control", "Dosis_baja", "Dosis_alta"))

# ── 2. Priors informativos ───────────────────────────────────────────────────
# Explorar priors por defecto para el modelo
get_prior(rendimiento ~ tratamiento + (1 | bloque), data = datos, family = gaussian())

# Definir priors propios (ajustar según conocimiento del dominio)
priors_lmm <- c(
  prior(normal(5, 3),  class = "Intercept"),
  prior(normal(0, 2),  class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "sigma")
)

# ── 3. Modelo lineal mixto bayesiano ─────────────────────────────────────────
cat("\n====== LMM Bayesiano ======\n")
modelo_bayes_lmm <- brm(
  formula  = rendimiento ~ tratamiento + (1 | bloque),
  data     = datos,
  family   = gaussian(),
  prior    = priors_lmm,
  chains   = 4,
  iter     = 4000,
  warmup   = 1000,
  cores    = parallel::detectCores() - 1,
  seed     = 42,
  file     = "resultados/modelo_bayes_lmm"  # caché del modelo compilado
)

# Resumen del modelo
summary(modelo_bayes_lmm)
print(modelo_bayes_lmm, digits = 3)

# ── 4. Diagnósticos MCMC ─────────────────────────────────────────────────────
# Trazas y distribuciones
plot(modelo_bayes_lmm)

# R-hat y tamaño efectivo de muestra
bayesplot::mcmc_rhat(brms::rhat(modelo_bayes_lmm))
bayesplot::mcmc_neff(brms::neff_ratio(modelo_bayes_lmm))

# Posterior predictive check
pp_check(modelo_bayes_lmm, ndraws = 100)
pp_check(modelo_bayes_lmm, type = "boxplot", ndraws = 50)

# ── 5. Inferencia bayesiana ───────────────────────────────────────────────────
# Distribuciones posteriores de los coeficientes
posterior_samples <- as_draws_df(modelo_bayes_lmm)

# Intervalos de credibilidad al 95% (HDI)
posterior::summarise_draws(
  posterior_samples,
  mean, median,
  ~posterior::quantile2(.x, probs = c(0.025, 0.975)),
  posterior::default_convergence_measures()
)

# Probabilidad de que el efecto de Dosis_alta > 0
mean(posterior_samples$b_tratamientoDosis_alta > 0)

# ── 6. Comparaciones entre tratamientos ──────────────────────────────────────
# Medias marginales posteriores
em_bayes <- emmeans::emmeans(modelo_bayes_lmm, ~ tratamiento)
summary(em_bayes, point.est = mean)

# Contrastes
emmeans::contrast(em_bayes, method = "pairwise")

# Con tidybayes
datos %>%
  tidybayes::add_epred_draws(modelo_bayes_lmm) %>%
  ggplot(aes(x = .epred, y = tratamiento, fill = tratamiento)) +
  ggdist::stat_halfeye(.width = c(0.66, 0.95)) +
  labs(
    title = "Distribución posterior de medias por tratamiento",
    x     = "Rendimiento esperado",
    y     = "Tratamiento"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

# ── 7. GLMM bayesiano (Poisson) ───────────────────────────────────────────────
cat("\n====== GLMM Bayesiano - Poisson ======\n")
priors_glmm <- c(
  prior(normal(0, 2.5), class = "Intercept"),
  prior(normal(0, 1),   class = "b"),
  prior(exponential(1), class = "sd")
)

modelo_bayes_glmm <- brm(
  formula  = conteo ~ tratamiento + (1 | bloque),
  data     = datos,
  family   = poisson(link = "log"),
  prior    = priors_glmm,
  chains   = 4,
  iter     = 4000,
  warmup   = 1000,
  cores    = parallel::detectCores() - 1,
  seed     = 42,
  file     = "resultados/modelo_bayes_glmm_poisson"
)

summary(modelo_bayes_glmm)
pp_check(modelo_bayes_glmm, ndraws = 100)

# ── 8. Comparación de modelos (LOO-CV) ────────────────────────────────────────
cat("\n====== Comparación de modelos ======\n")

# Modelo nulo (solo intercepto + bloque)
modelo_nulo <- brm(
  formula = rendimiento ~ 1 + (1 | bloque),
  data    = datos,
  family  = gaussian(),
  prior   = c(
    prior(normal(5, 3),   class = "Intercept"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, iter = 4000, warmup = 1000,
  cores  = parallel::detectCores() - 1,
  seed   = 42,
  file   = "resultados/modelo_bayes_nulo"
)

# LOO
loo_lmm  <- loo(modelo_bayes_lmm)
loo_nulo <- loo(modelo_nulo)
loo_compare(loo_lmm, loo_nulo)

# ── 9. Exportar resultados ────────────────────────────────────────────────────
# Descomentar para guardar
# saveRDS(modelo_bayes_lmm, "resultados/modelo_bayes_lmm.rds")
# write.csv(
#   as.data.frame(fixef(modelo_bayes_lmm)),
#   "resultados/efectos_fijos_bayes.csv"
# )
