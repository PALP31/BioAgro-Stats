# ------------------------------------------------------------------------
# 04_Estadistica_Bayesiana
# Modelo mixto bayesiano avanzado con brms
# Contexto: biomasa en genotipos y sitios bajo estrés hídrico
# ------------------------------------------------------------------------

set.seed(778)

n_gen <- 20
n_sitios <- 4
n_rep <- 10

datos <- expand.grid(
  genotipo = factor(paste0("G", 1:n_gen)),
  sitio = factor(paste0("Sitio", 1:n_sitios)),
  rep = 1:n_rep,
  tratamiento = factor(c("Control", "Deficit_hidrico"))
)

ef_gen <- rnorm(n_gen, 0, 2.8); names(ef_gen) <- levels(datos$genotipo)
ef_sit <- rnorm(n_sitios, 0, 3.1); names(ef_sit) <- levels(datos$sitio)

# Simulación de biomasa continua con estructura jerárquica
datos$biomasa <- 45 +
  ifelse(datos$tratamiento == "Deficit_hidrico", -6.5, 0) +
  ef_gen[as.character(datos$genotipo)] +
  ef_sit[as.character(datos$sitio)] +
  rnorm(nrow(datos), 0, 3.5)

if (requireNamespace("brms", quietly = TRUE)) {
  # Priors débiles informativos, típicos en agronomía cuando hay poca información previa
  priors <- c(
    brms::prior(normal(0, 10), class = "b"),
    brms::prior(student_t(3, 0, 10), class = "Intercept"),
    brms::prior(exponential(1), class = "sd"),
    brms::prior(exponential(1), class = "sigma")
  )

  # Modelo mixto bayesiano
  fit <- brms::brm(
    biomasa ~ tratamiento + (1 | genotipo) + (1 | sitio),
    data = datos,
    family = gaussian(),
    prior = priors,
    chains = 2,
    iter = 2000,
    warmup = 1000,
    seed = 2026,
    cores = 2,
    refresh = 0
  )

  cat("\n--- Resumen del modelo mixto bayesiano (brms) ---\n")
  print(summary(fit))

  # Chequeo posterior predictivo
  print(brms::pp_check(fit))
} else {
  message("Paquete 'brms' no disponible. Instálalo con: install.packages('brms')")
  message("Nota: brms requiere CmdStan o rstan para compilar modelos bayesianos.")
}
