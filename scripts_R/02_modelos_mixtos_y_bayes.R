# 02_modelos_mixtos_y_bayes.R
# Estructura base para LMM, GLMM y un modelo bayesiano simple con brms.

required_packages <- c("lme4", "lmerTest", "brms")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste(
      "Instala los siguientes paquetes antes de ejecutar este script:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(lme4)
library(lmerTest)
library(brms)

# Asume que `datos` incluye:
# respuesta (continua), tratamiento, covariable
# bloque (factor aleatorio), y opcionalmente respuesta_count / respuesta_bin

# ------------------------------------------------------------------
# 1) LMM (respuesta continua)
# ------------------------------------------------------------------
modelo_lmm <- lmer(
  respuesta ~ tratamiento + covariable + (1 | bloque),
  data = datos,
  REML = TRUE
)
summary(modelo_lmm)
anova(modelo_lmm)

# ------------------------------------------------------------------
# 2) GLMM Poisson (conteos)
# ------------------------------------------------------------------
if ("respuesta_count" %in% names(datos)) {
  modelo_glmm_pois <- glmer(
    respuesta_count ~ tratamiento + covariable + (1 | bloque),
    family = poisson(link = "log"),
    data = datos
  )
  print(summary(modelo_glmm_pois))
}

# ------------------------------------------------------------------
# 3) GLMM Binomial (0/1)
# ------------------------------------------------------------------
if ("respuesta_bin" %in% names(datos)) {
  modelo_glmm_bin <- glmer(
    respuesta_bin ~ tratamiento + covariable + (1 | bloque),
    family = binomial(link = "logit"),
    data = datos
  )
  print(summary(modelo_glmm_bin))
}

# ------------------------------------------------------------------
# 4) Modelo bayesiano simple con brms
# ------------------------------------------------------------------
# Nota: requiere backend de Stan configurado (cmdstanr o rstan).
priors_basicos <- c(
  set_prior("normal(0, 5)", class = "b"),
  set_prior("student_t(3, 0, 10)", class = "Intercept"),
  set_prior("exponential(1)", class = "sd")
)

modelo_bayes <- brm(
  formula = respuesta ~ tratamiento + covariable + (1 | bloque),
  data = datos,
  family = gaussian(),
  prior = priors_basicos,
  chains = 2,
  iter = 2000,
  warmup = 1000,
  cores = 2,
  seed = 123
)

print(summary(modelo_bayes))
plot(modelo_bayes)
pp_check(modelo_bayes)
