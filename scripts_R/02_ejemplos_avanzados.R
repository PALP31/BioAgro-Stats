# =============================================================================
# 02_ejemplos_avanzados.R
# Ejemplos de modelos estadísticos avanzados:
#   - GAM  (Generalized Additive Model) usando mgcv
#   - Modelo Lineal Bayesiano simple usando brms
# =============================================================================

# --- Paquetes necesarios -----------------------------------------------------
library(mgcv)   # GAM
library(brms)   # Modelos Bayesianos (requiere RStan instalado)
# Alternativa Bayesiana: library(rstanarm)


# =============================================================================
# 1. GAM – MODELO ADITIVO GENERALIZADO (mgcv)
# =============================================================================

# --- Cargar / simular datos de ejemplo ---------------------------------------
set.seed(42)
n <- 200
datos_gam <- data.frame(
  x1 = seq(0, 4 * pi, length.out = n),
  x2 = rnorm(n, mean = 5, sd = 1.5),
  grupo = factor(sample(c("A", "B"), n, replace = TRUE))
)
datos_gam$y <- sin(datos_gam$x1) + 0.3 * datos_gam$x2 +
  ifelse(datos_gam$grupo == "B", 1.5, 0) + rnorm(n, sd = 0.5)

# --- Ajustar modelo ----------------------------------------------------------
# s()  = término de suavizado (spline cúbico por defecto)
# k    = número máximo de nodos (ajustar según datos)
modelo_gam <- mgcv::gam(
  y ~ s(x1, k = 10) + s(x2, k = 5) + grupo,
  family = gaussian(),
  method = "REML",     # Estimación por máxima verosimilitud restringida
  data   = datos_gam
)

# --- Revisar supuestos y resumen ---------------------------------------------
summary(modelo_gam)

# Verificar si k es suficiente (edf cercano a k-1 indica que k podría ser mayor)
mgcv::k.check(modelo_gam)

# Gráficos de los términos suavizados
par(mfrow = c(1, 2))
plot(modelo_gam, pages = 1, residuals = TRUE, pch = 20, cex = 0.4)
par(mfrow = c(1, 1))

# Gráficos de diagnóstico de residuos
mgcv::gam.check(modelo_gam)

# --- Ejemplo: GAM con familia Poisson (conteos) ------------------------------
datos_gam_pois <- datos_gam
datos_gam_pois$conteo <- rpois(n, lambda = exp(0.5 + sin(datos_gam$x1) / 2))

modelo_gam_pois <- mgcv::gam(
  conteo ~ s(x1, k = 10) + x2,
  family = poisson(link = "log"),
  method = "REML",
  data   = datos_gam_pois
)
summary(modelo_gam_pois)


# =============================================================================
# 2. MODELO LINEAL BAYESIANO SIMPLE (brms)
# =============================================================================
# NOTA: brms requiere RStan. La primera vez que se ejecuta compila el modelo
# en C++, lo que puede tardar varios minutos.
# Para instalación: install.packages("brms")

# --- Cargar / simular datos de ejemplo ---------------------------------------
set.seed(7)
n_bayes <- 100
datos_bayes <- data.frame(
  x = rnorm(n_bayes, mean = 0, sd = 1)
)
datos_bayes$y <- 2 + 1.5 * datos_bayes$x + rnorm(n_bayes, sd = 0.8)

# --- Especificar priors ------------------------------------------------------
priors_lm <- c(
  brms::prior(normal(0, 10),  class = "b"),          # Prior sobre coeficientes
  brms::prior(normal(0, 10),  class = "Intercept"),  # Prior sobre intercepto
  brms::prior(exponential(1), class = "sigma")       # Prior sobre desviación estándar
)

# --- Ajustar modelo ----------------------------------------------------------
modelo_bayes_lm <- brms::brm(
  formula  = y ~ x,
  data     = datos_bayes,
  family   = gaussian(),
  prior    = priors_lm,
  chains   = 4,         # Número de cadenas MCMC
  iter     = 2000,      # Iteraciones por cadena (incluyendo warmup)
  warmup   = 1000,      # Iteraciones de calentamiento
  cores    = 4,         # Núcleos en paralelo (ajustar según CPU disponible)
  seed     = 42,
  file     = "datos/modelo_bayes_lm"  # Guarda el modelo compilado
)

# --- Revisar supuestos y diagnósticos ----------------------------------------
print(modelo_bayes_lm)
summary(modelo_bayes_lm)

# Diagnóstico de convergencia (R-hat y tamaño efectivo de muestra)
# R-hat debe ser <= 1.01; n_eff debe ser > 400 (preferiblemente > 1000)
brms::rhat(modelo_bayes_lm)
brms::neff_ratio(modelo_bayes_lm)

# Gráficos de cadenas y distribuciones posteriores
plot(modelo_bayes_lm)

# Posterior predictive check (¿el modelo reproduce los datos observados?)
brms::pp_check(modelo_bayes_lm, ndraws = 100)

# Intervalos de credibilidad (95%)
brms::posterior_summary(modelo_bayes_lm)


# =============================================================================
# ALTERNATIVA: MODELO LINEAL BAYESIANO CON rstanarm
# =============================================================================
# Descomenta las líneas siguientes si prefieres usar rstanarm (más sencillo
# de instalar y sin necesidad de compilar modelos).

# library(rstanarm)
#
# modelo_rstanarm <- rstanarm::stan_glm(
#   formula = y ~ x,
#   data    = datos_bayes,
#   family  = gaussian(),
#   prior   = rstanarm::normal(0, 2.5, autoscale = TRUE),
#   prior_intercept = rstanarm::normal(0, 10, autoscale = TRUE),
#   chains  = 4,
#   iter    = 2000,
#   seed    = 42
# )
#
# summary(modelo_rstanarm)
# rstanarm::posterior_interval(modelo_rstanarm, prob = 0.95)
# rstanarm::pp_check(modelo_rstanarm)
