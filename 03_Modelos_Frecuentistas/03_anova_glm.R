# 03_Modelos_Frecuentistas/03_anova_glm.R
# Ejemplos: ANOVA 1 y 2 vías, GLM Poisson y Binomial

set.seed(100)

# -------------------------------
# 1) ANOVA de una vía
# -------------------------------
d1 <- data.frame(
  tratamiento = factor(rep(c("Control", "Moderado", "Severo"), each = 25))
)
d1$rendimiento <- with(d1,
  c(rnorm(25, 6.1, 0.4), rnorm(25, 5.6, 0.45), rnorm(25, 4.9, 0.5))
)

anova1 <- aov(rendimiento ~ tratamiento, data = d1)
cat("=== ANOVA 1 via ===\n")
print(summary(anova1))

# -------------------------------
# 2) ANOVA de dos vías con interacción
# -------------------------------
d2 <- expand.grid(
  tratamiento = factor(c("Control", "Moderado", "Severo")),
  variedad = factor(c("V1", "V2")),
  rep = 1:15
)

base <- with(d2,
  6.5 -
    ifelse(tratamiento == "Moderado", 0.6, ifelse(tratamiento == "Severo", 1.3, 0)) +
    ifelse(variedad == "V2", 0.25, 0) +
    ifelse(tratamiento == "Severo" & variedad == "V2", -0.2, 0)
)

d2$rendimiento <- base + rnorm(nrow(d2), 0, 0.35)

anova2 <- aov(rendimiento ~ tratamiento * variedad, data = d2)
cat("\n=== ANOVA 2 vias ===\n")
print(summary(anova2))

# -------------------------------
# 3) GLM Poisson (conteos)
# -------------------------------
pois <- data.frame(
  tratamiento = factor(rep(c("Control", "Moderado", "Severo"), each = 30))
)
lambda <- c(Control = 2.5, Moderado = 4.0, Severo = 6.1)
pois$lesiones <- rpois(nrow(pois), lambda = lambda[pois$tratamiento])

glm_pois <- glm(lesiones ~ tratamiento, family = poisson(link = "log"), data = pois)
cat("\n=== GLM Poisson ===\n")
print(summary(glm_pois))

# -------------------------------
# 4) GLM Binomial (supervivencia)
# -------------------------------
bin <- data.frame(
  tratamiento = factor(rep(c("Control", "Moderado", "Severo"), each = 35))
)
p <- c(Control = 0.88, Moderado = 0.70, Severo = 0.52)
bin$supervive <- rbinom(nrow(bin), size = 1, prob = p[bin$tratamiento])

glm_bin <- glm(supervive ~ tratamiento, family = binomial(link = "logit"), data = bin)
cat("\n=== GLM Binomial ===\n")
print(summary(glm_bin))
