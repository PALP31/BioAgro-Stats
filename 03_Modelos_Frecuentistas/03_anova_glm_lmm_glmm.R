# ------------------------------------------------------------------
# 03_Modelos_Frecuentistas
# ANOVA, GLM y LMM/GLMM
# Contexto: crecimiento vegetal y control biológico con Trichoderma
# ------------------------------------------------------------------

set.seed(456)

# =========================
# A) ANOVA (altura vegetal)
# =========================
n <- 30
d_anova <- data.frame(
  tratamiento = factor(rep(c("Control", "Bioestimulante_A", "Bioestimulante_B"), each = n)),
  altura_cm = c(
    rnorm(n, 35, 4),
    rnorm(n, 40, 4),
    rnorm(n, 43, 5)
  )
)

mod_anova <- aov(altura_cm ~ tratamiento, data = d_anova)
cat("\n--- ANOVA (altura vegetal) ---\n")
print(summary(mod_anova))

# =======================================
# B) GLM binomial (control biológico)
# =======================================
# Respuesta: éxito de antagonismo (1/0) de aislados microbianos
n_glm <- 180
d_glm <- data.frame(
  tratamiento = factor(sample(c("Control", "Trichoderma", "Bacillus"), n_glm, replace = TRUE)),
  temperatura = rnorm(n_glm, 26, 3)
)

linpred <- -2 +
  ifelse(d_glm$tratamiento == "Trichoderma", 1.3,
         ifelse(d_glm$tratamiento == "Bacillus", 0.7, 0)) +
  0.06 * d_glm$temperatura

prob <- 1 / (1 + exp(-linpred))
d_glm$antagonismo_exitoso <- rbinom(n_glm, 1, prob)

mod_glm <- glm(antagonismo_exitoso ~ tratamiento + temperatura,
               family = binomial(link = "logit"), data = d_glm)
cat("\n--- GLM binomial (antagonismo) ---\n")
print(summary(mod_glm))

# ========================================
# C) LMM / GLMM (datos jerárquicos)
# ========================================
# Parcelas repetidas dentro de bloques y sitios
if (requireNamespace("lme4", quietly = TRUE)) {
  n_b <- 6
  n_s <- 3
  n_rep <- 8

  d_mix <- expand.grid(
    bloque = factor(paste0("B", 1:n_b)),
    sitio = factor(paste0("Sitio", 1:n_s)),
    repeticion = 1:n_rep,
    tratamiento = factor(c("Control", "Trichoderma"))
  )

  eb <- rnorm(n_b, 0, 1.8); names(eb) <- levels(d_mix$bloque)
  es <- rnorm(n_s, 0, 2.2); names(es) <- levels(d_mix$sitio)

  d_mix$crecimiento <- 30 +
    ifelse(d_mix$tratamiento == "Trichoderma", 4.5, 0) +
    eb[as.character(d_mix$bloque)] +
    es[as.character(d_mix$sitio)] +
    rnorm(nrow(d_mix), 0, 3.0)

  mod_lmm <- lme4::lmer(crecimiento ~ tratamiento + (1 | bloque) + (1 | sitio),
                        data = d_mix)
  cat("\n--- LMM (crecimiento continuo) ---\n")
  print(summary(mod_lmm))

  # GLMM Poisson: recuento de colonias inhibidas
  eta <- 1.7 +
    ifelse(d_mix$tratamiento == "Trichoderma", 0.35, 0) +
    0.03 * d_mix$crecimiento +
    eb[as.character(d_mix$bloque)] / 8
  d_mix$colonias_inhibidas <- rpois(nrow(d_mix), lambda = exp(eta))

  mod_glmm <- lme4::glmer(colonias_inhibidas ~ tratamiento + crecimiento + (1 | bloque),
                          family = poisson(link = "log"), data = d_mix)
  cat("\n--- GLMM Poisson (colonias inhibidas) ---\n")
  print(summary(mod_glmm))
} else {
  message("Paquete 'lme4' no disponible. Instálalo con: install.packages('lme4')")
}
