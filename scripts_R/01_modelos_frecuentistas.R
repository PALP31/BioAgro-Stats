# 01_modelos_frecuentistas.R
# Estructura base para modelos frecuentistas (LM, ANOVA, GLM).

# Asume que existe un data.frame `datos` con:
# - respuesta continua: respuesta
# - tratamiento (factor)
# - covariable (numérica)
# - opcional: respuesta_count (conteos), respuesta_bin (0/1)

required_packages <- c("emmeans")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    paste(
      "Instala los siguientes paquetes antes de ejecutar este script:",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(emmeans)

# ------------------------------------------------------------------
# 1) Modelo lineal (LM) - usar si residuos ~ normales y varianza homogénea
# ------------------------------------------------------------------
modelo_lm <- lm(respuesta ~ tratamiento + covariable, data = datos)
summary(modelo_lm)
anova(modelo_lm)

# Comparaciones post-hoc
emmeans_lm <- emmeans(modelo_lm, specs = ~ tratamiento)
pairs(emmeans_lm, adjust = "tukey")

# ------------------------------------------------------------------
# 2) ANOVA clásica - diseño simple por tratamiento
# ------------------------------------------------------------------
modelo_aov <- aov(respuesta ~ tratamiento, data = datos)
summary(modelo_aov)
TukeyHSD(modelo_aov)

# ------------------------------------------------------------------
# 3) GLM Poisson - para conteos
# ------------------------------------------------------------------
if ("respuesta_count" %in% names(datos)) {
  modelo_pois <- glm(
    respuesta_count ~ tratamiento + covariable,
    family = poisson(link = "log"),
    data = datos
  )
  print(summary(modelo_pois))

  # Chequeo simple de sobredispersión
  sobredisp <- sum(residuals(modelo_pois, type = "pearson")^2) / modelo_pois$df.residual
  cat("Sobredispersión (Poisson):", round(sobredisp, 3), "\n")
}

# ------------------------------------------------------------------
# 4) GLM Binomial - para respuesta 0/1
# ------------------------------------------------------------------
if ("respuesta_bin" %in% names(datos)) {
  modelo_bin <- glm(
    respuesta_bin ~ tratamiento + covariable,
    family = binomial(link = "logit"),
    data = datos
  )
  print(summary(modelo_bin))

  # Odds ratios
  or <- exp(cbind(OR = coef(modelo_bin), confint(modelo_bin)))
  print(or)
}
