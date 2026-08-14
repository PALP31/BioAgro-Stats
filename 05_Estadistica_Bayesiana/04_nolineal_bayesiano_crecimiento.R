# ============================================================================
# 04_nolineal_bayesiano_crecimiento.R
# Modelo No Lineal Bayesiano (Curva de Crecimiento Logístico)
# ============================================================================
# En este script modelamos el diámetro de Botrytis cinerea a lo largo del 
# tiempo usando la fórmula de crecimiento logístico.
# Utilizamos la sintaxis bf() de brms para modelos no lineales (NL).
# ============================================================================

library(brms)        # Motor bayesiano para modelos NL
library(tidyverse)
library(performance) # Diagnóstico visual

set.seed(789)

# 1. SIMULACIÓN DE DATOS (Crecimiento de Botrytis en 10 días)
# Función logística: Diámetro = asym / (1 + exp((xmid - tiempo) / scal))
# asym = asíntota máxima; xmid = tiempo de crecimiento medio; scal = tasa.

n_rep <- 10
tiempos <- 1:10
datos_growth <- expand.grid(tiempo = tiempos, rep = 1:n_rep) %>%
  mutate(
    asym = 60, xmid = 5, scal = 1.2,
    diametro_base = asym / (1 + exp((xmid - tiempo) / scal)),
    diametro = diametro_base + rnorm(n(), 0, 2)
  )

# 2. ESPECIFICACIÓN DEL MODELO NO LINEAL (bf)
# Definimos los parámetros de la curva como funciones constantes (interceptos).

formula_nl <- bf(
  diametro ~ asym / (1 + exp((xmid - tiempo) / scal)),
  asym + xmid + scal ~ 1,
  nl = TRUE
)

# 3. PRIORS INFORMADOS PARA MODELOS NO LINEALES
# En modelos NL, los priors son cruciales para la convergencia.
# asym (asíntota) ~ Normal(60, 10); xmid ~ Normal(5, 2); scal ~ Normal(1, 0.5).

priors_nl <- c(
  prior(normal(60, 10), nlpar = "asym"),
  prior(normal(5, 2), nlpar = "xmid"),
  prior(normal(1, 1), nlpar = "scal")
)

# 4. AJUSTE DEL MODELO BAYESIANO NO LINEAL
mod_growth <- brm(
  formula_nl,
  data = datos_growth,
  family = gaussian(),
  prior = priors_nl,
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 123
)

# 5. VISUALIZACIÓN: CURVA DE PREDICCIÓN MEDIA (Credible Interval)
# conditional_effects() genera la curva ajustada con su banda de incertidumbre.

plot_growth <- conditional_effects(mod_growth)
plot(plot_growth, points = TRUE) # Graficamos con los puntos observados

cat("\n--- RESULTADOS DEL MODELO NO LINEAL ---")
print(summary(mod_growth))

# EXPLICACIÓN:
# La banda de incertidumbre sombreada representa el 95% de probabilidad de
# donde se encuentra la verdadera curva de crecimiento biológico. 
# Esto supera a las curvas polinomiales simples por su realismo biológico.

cat("\n--- Script 04 finalizado correctamente ---\n")
