# =============================================================================
# 01_ejemplos_frecuentistas.R
# Ejemplos de modelos estadísticos frecuentistas:
#   - Modelo Lineal (LM)
#   - Modelo Lineal Generalizado (GLM) con distribuciones Poisson y Binomial
#   - Análisis de Varianza (ANOVA)
#   - Modelo Lineal Mixto Generalizado (GLMM)
# =============================================================================

# --- Paquetes necesarios -----------------------------------------------------
library(lme4)      # GLMM
library(car)       # Pruebas de diagnóstico y Anova tipo II/III

# =============================================================================
# 1. MODELO LINEAL (LM)
# =============================================================================

# --- Cargar / simular datos de ejemplo ---------------------------------------
set.seed(123)
n <- 100
datos_lm <- data.frame(
  x1 = rnorm(n, mean = 50, sd = 10),
  x2 = rnorm(n, mean = 5,  sd = 2),
  grupo = factor(sample(c("A", "B", "C"), n, replace = TRUE))
)
datos_lm$y <- 3 + 0.5 * datos_lm$x1 - 1.2 * datos_lm$x2 + rnorm(n, sd = 5)

# --- Ajustar modelo ----------------------------------------------------------
modelo_lm <- lm(y ~ x1 + x2 + grupo, data = datos_lm)

# --- Revisar supuestos -------------------------------------------------------
summary(modelo_lm)
par(mfrow = c(2, 2))
plot(modelo_lm)          # Gráficos de diagnóstico (residuos, Q-Q, etc.)
par(mfrow = c(1, 1))

# Prueba de normalidad de residuos
shapiro.test(residuals(modelo_lm))

# Prueba de homocedasticidad (Breusch-Pagan via car)
car::ncvTest(modelo_lm)


# =============================================================================
# 2. GLM – DISTRIBUCIÓN POISSON (datos de conteos)
# =============================================================================

# --- Cargar / simular datos de ejemplo ---------------------------------------
datos_pois <- data.frame(
  exposicion = runif(n, 0.5, 2),
  x1 = rnorm(n)
)
datos_pois$conteos <- rpois(n, lambda = exp(0.8 + 0.6 * datos_pois$x1))

# --- Ajustar modelo ----------------------------------------------------------
modelo_pois <- glm(conteos ~ x1 + offset(log(exposicion)),
                   family = poisson(link = "log"),
                   data   = datos_pois)

# --- Revisar supuestos -------------------------------------------------------
summary(modelo_pois)

# Sobredispersión: razón entre devianza residual y grados de libertad
cat("Razón devianza/gl:", modelo_pois$deviance / modelo_pois$df.residual, "\n")

# Gráficos de diagnóstico
par(mfrow = c(2, 2))
plot(modelo_pois)
par(mfrow = c(1, 1))


# =============================================================================
# 3. GLM – DISTRIBUCIÓN BINOMIAL (datos binarios / proporciones)
# =============================================================================

# --- Cargar / simular datos de ejemplo ---------------------------------------
datos_bin <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n)
)
prob_exito <- plogis(-0.5 + 1.2 * datos_bin$x1 - 0.8 * datos_bin$x2)
datos_bin$exito <- rbinom(n, size = 1, prob = prob_exito)

# --- Ajustar modelo ----------------------------------------------------------
modelo_bin <- glm(exito ~ x1 + x2,
                  family = binomial(link = "logit"),
                  data   = datos_bin)

# --- Revisar supuestos -------------------------------------------------------
summary(modelo_bin)

# Razón devianza/gl (bondad de ajuste)
cat("Razón devianza/gl:", modelo_bin$deviance / modelo_bin$df.residual, "\n")

# Gráficos de diagnóstico
par(mfrow = c(2, 2))
plot(modelo_bin)
par(mfrow = c(1, 1))


# =============================================================================
# 4. ANOVA (una y dos vías)
# =============================================================================

# --- Cargar / simular datos de ejemplo ---------------------------------------
datos_anova <- data.frame(
  tratamiento = factor(rep(c("Control", "T1", "T2"), each = 30)),
  bloque      = factor(rep(1:3, times = 30)),
  respuesta   = c(
    rnorm(30, mean = 10, sd = 2),
    rnorm(30, mean = 13, sd = 2),
    rnorm(30, mean = 11, sd = 2)
  )
)

# --- Ajustar modelo ----------------------------------------------------------
# ANOVA de una vía
modelo_aov1 <- aov(respuesta ~ tratamiento, data = datos_anova)

# ANOVA de dos vías
modelo_aov2 <- aov(respuesta ~ tratamiento + bloque, data = datos_anova)

# --- Revisar supuestos -------------------------------------------------------
summary(modelo_aov1)
summary(modelo_aov2)

# Normalidad de residuos
shapiro.test(residuals(modelo_aov1))

# Homogeneidad de varianzas (Levene)
car::leveneTest(respuesta ~ tratamiento, data = datos_anova)

# Comparaciones múltiples post-hoc (Tukey)
TukeyHSD(modelo_aov1)

# Gráficos de diagnóstico
par(mfrow = c(2, 2))
plot(modelo_aov1)
par(mfrow = c(1, 1))


# =============================================================================
# 5. GLMM – MODELO LINEAL MIXTO GENERALIZADO (lme4)
# =============================================================================

# --- Cargar / simular datos de ejemplo ---------------------------------------
n_suj   <- 20
n_obs   <- 10
datos_glmm <- data.frame(
  sujeto    = factor(rep(1:n_suj, each = n_obs)),
  tiempo    = rep(1:n_obs, times = n_suj),
  tratamiento = factor(rep(c("A", "B"), each = n_suj * n_obs / 2))
)
# Efectos aleatorios por sujeto
b0 <- rnorm(n_suj, mean = 0, sd = 1)[datos_glmm$sujeto]
datos_glmm$conteo <- rpois(
  nrow(datos_glmm),
  lambda = exp(1.5 + 0.2 * datos_glmm$tiempo + b0)
)

# --- Ajustar modelo ----------------------------------------------------------
# GLMM Poisson con intercepto aleatorio por sujeto
modelo_glmm <- lme4::glmer(
  conteo ~ tiempo + tratamiento + (1 | sujeto),
  family = poisson(link = "log"),
  data   = datos_glmm
)

# --- Revisar supuestos -------------------------------------------------------
summary(modelo_glmm)

# Verificar convergencia
if (!is.null(modelo_glmm@optinfo$conv$lme4)) {
  warning("Posibles problemas de convergencia.")
}

# Residuos de Pearson vs valores ajustados
plot(fitted(modelo_glmm), residuals(modelo_glmm, type = "pearson"),
     xlab = "Valores ajustados", ylab = "Residuos de Pearson",
     main = "GLMM – Residuos vs Ajustados")
abline(h = 0, lty = 2)
