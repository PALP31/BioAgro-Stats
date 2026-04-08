# 04_Modelos_Avanzados/04_glmm_gam_bayes.R
# Estructuras base: GLMM, GAM y Bayesiano (brms)

suppressPackageStartupMessages({
  library(lme4)
  library(mgcv)
  library(brms)
})

set.seed(321)

# -------------------------------
# 1) GLMM: efecto fijo de tratamiento + aleatorios temporal/espacial
# -------------------------------
n_plot <- 24
tiempo <- 1:8

glmm_dat <- expand.grid(
  parcela = factor(1:n_plot),
  tiempo = tiempo
)
glmm_dat$tratamiento <- factor(ifelse(as.numeric(glmm_dat$parcela) <= n_plot / 2, "Control", "Severo"))
glmm_dat$bloque <- factor(rep(1:6, each = 4)[as.numeric(glmm_dat$parcela)])
glmm_dat$sitio <- factor(rep(1:3, length.out = n_plot)[as.numeric(glmm_dat$parcela)])

eta <- with(glmm_dat,
  1.6 +
    ifelse(tratamiento == "Severo", 0.35, 0) +
    0.08 * tiempo
)

# Conteo de eventos de dano
u_parcela <- rnorm(n_plot, 0, 0.25)
u_tiempo <- rnorm(length(tiempo), 0, 0.15)

glmm_dat$eventos <- rpois(
  nrow(glmm_dat),
  lambda = exp(eta + u_parcela[glmm_dat$parcela] + u_tiempo[glmm_dat$tiempo])
)

mod_glmm <- glmer(
  eventos ~ tratamiento + tiempo + (1 | bloque/parcela) + (1 | tiempo),
  family = poisson,
  data = glmm_dat
)

cat("=== GLMM ===\n")
print(summary(mod_glmm))

# -------------------------------
# 2) GAM: respuesta no lineal en el tiempo
# -------------------------------
gam_dat <- data.frame(
  tiempo = rep(seq(1, 60, by = 1), 2),
  tratamiento = factor(rep(c("Control", "Severo"), each = 60))
)

gam_dat$clorofila <- with(gam_dat,
  42 +
    ifelse(tratamiento == "Severo", -3.5, 0) +
    sin(tiempo / 8) * ifelse(tratamiento == "Severo", -4, -2) +
    rnorm(nrow(gam_dat), 0, 1.2)
)

mod_gam <- gam(clorofila ~ tratamiento + s(tiempo, by = tratamiento, k = 10), data = gam_dat)

cat("\n=== GAM ===\n")
print(summary(mod_gam))

# -------------------------------
# 3) Estructura Bayesiana basica (brms)
# -------------------------------
# Nota: este bloque muestra la estructura. Puede requerir tiempo de compilacion.
bayes_dat <- data.frame(
  bloque = factor(rep(1:5, each = 20)),
  tratamiento = factor(rep(c("Control", "Moderado", "Severo", "Severo"), length.out = 100))
)

mu <- with(bayes_dat,
  6.0 + ifelse(tratamiento == "Moderado", -0.5, ifelse(tratamiento == "Severo", -1.1, 0))
)
bayes_dat$rendimiento <- mu + rnorm(nrow(bayes_dat), 0, 0.5)

formula_brm <- bf(rendimiento ~ tratamiento + (1 | bloque))
priors_brm <- c(
  prior(normal(0, 2), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sigma")
)

# Ejecutar solo cuando se desee ajustar:
# mod_bayes <- brm(
#   formula = formula_brm,
#   data = bayes_dat,
#   family = gaussian(),
#   prior = priors_brm,
#   chains = 2, iter = 2000, cores = 2, seed = 321
# )

cat("\n=== Estructura Bayesiana (brms) preparada ===\n")
cat("Formula: rendimiento ~ tratamiento + (1 | bloque)\n")
