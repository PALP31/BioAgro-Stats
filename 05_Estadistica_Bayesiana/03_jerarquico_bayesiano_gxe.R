# ============================================================================
# 03_jerarquico_bayesiano_gxe.R
# Modelo Multinivel Bayesiano para GxE y Efecto de "Shrinkage"
# ============================================================================
# En este script modelamos la interacción Genotipo x Ambiente en trigo.
# Explicamos el fenómeno de contracción parcial (shrinkage), que mejora
# drásticamente las predicciones genotípicas ante datos faltantes.
# ============================================================================

library(brms)        # Motor bayesiano
library(bayesplot)   # Caterpillar plot (mcmc_intervals)
library(tidybayes)
library(tidyverse)

set.seed(456)

# 1. SIMULACIÓN DE DATOS (Genotipos en ambientes óptimo y estresado)
n_gen <- 10
n_env <- 2
datos_gxe <- expand.grid(
  genotipo = factor(paste0("G", 1:n_gen)),
  ambiente = factor(c("Optimo", "Calor")),
  rep = 1:2
) %>%
  mutate(
    yield_base = case_when(
      ambiente == "Optimo" ~ 4500,
      ambiente == "Calor" ~ 3200
    ),
    ef_gen = rnorm(n_gen, 0, 400)[as.numeric(genotipo)],
    rendimiento = yield_base + ef_gen + rnorm(n(), 0, 150)
  )

# 2. AJUSTE DEL MODELO JERÁRQUICO
# Formula: Rendimiento ~ Ambiente + (1 | Genotipo)
# Modelamos los genotipos como efectos aleatorios bajo una distribución común.

mod_gxe <- brm(
  rendimiento ~ ambiente + (1 | genotipo),
  data = datos_gxe,
  family = gaussian(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 123
)

# 3. EXPLICACIÓN: SHRINKAGE (Contracción Parcial)
# En el enfoque bayesiano jerárquico, las estimaciones de genotipos 
# extremos (que rinden muy poco o demasiado debido al azar) se "contraen"
# hacia la media global. Esto reduce errores de selección y
# compensa genotipos que tienen pocas observaciones.

cat("\n--- ESTIMACIONES GENOTÍPICAS (EFECTOS ALEATORIOS) ---")
print(ranef(mod_gxe)$genotipo[,,1])

# 4. VISUALIZACIÓN: CATERPILLAR PLOT (mcmc_intervals)
# Muestra el ranking de genotipos con su incertidumbre bayesiana.
mcmc_plot(mod_gxe, type = "intervals", variable = "^r_genotipo", regex = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Ranking de Genotipos (Caterpillar Plot)",
       subtitle = "Efectos aleatorios de genotipos. Punto = Mediana posterior.",
       x = "Efecto sobre el Rendimiento (kg/ha)") +
  theme_minimal()

cat("\n--- Script 03 finalizado correctamente ---\n")
