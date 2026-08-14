# ============================================================================
# 01_rf_espectroscopia.R
# Predicción de Biomasa mediante Random Forest y Tidymodels
# ============================================================================
# Contexto: Uso de 20 índices espectrales (NDVI, GNDVI, etc.) para estimar
# la biomasa aérea en cultivos de trigo.
# ============================================================================

library(tidymodels)
library(vip)         # Importancia de variables
library(tidyverse)

set.seed(123)

# 1. SIMULACIÓN DE DATOS (20 Índices Espectrales + Biomasa)
n <- 200
indices <- matrix(runif(n * 20, 0.2, 0.9), ncol = 20)
colnames(indices) <- c("NDVI", "GNDVI", "EVI", "SAVI", "MSAVI", "NDRE", "CI_RedEdge", 
                       paste0("Idx_", 8:20))

datos_rf <- as_tibble(indices) %>%
  mutate(
    # Biomasa depende fuertemente de NDVI y GNDVI
    biomasa = 15 + (NDVI * 40) + (GNDVI * 25) + (NDRE * 15) + rnorm(n, 0, 5)
  )

# 2. PARTICIÓN DE DATOS (Tidymodels: initial_split)
set.seed(456)
split_datos <- initial_split(datos_rf, prop = 0.8, strata = biomasa)
entrenamiento <- training(split_datos)
testeo <- testing(split_datos)

# 3. RECETA DE PREPROCESAMIENTO (Recipes)
receta_rf <- recipe(biomasa ~ ., data = entrenamiento) %>%
  step_normalize(all_predictors())

# 4. ESPECIFICACIÓN DEL MODELO (Parsnip)
espec_rf <- rand_forest(trees = 500, mtry = 6, min_n = 5) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("regression")

# 5. WORKFLOW Y VALIDACIÓN CRUZADA
vfold <- vfold_cv(entrenamiento, v = 5)

workflow_rf <- workflow() %>%
  add_recipe(receta_rf) %>%
  add_model(espec_rf)

# Ajuste final
ajuste_rf <- workflow_rf %>%
  fit(data = entrenamiento)

# 6. EVALUACIÓN (Yardstick)
predicciones <- predict(ajuste_rf, testeo) %>%
  bind_cols(testeo)

metricas <- metrics(predicciones, truth = biomasa, estimate = .pred)
cat("\n--- Métricas de Desempeño (Test Set) ---\n")
print(metricas)

# 7. IMPORTANCIA DE VARIABLES (vip)
cat("\n--- Generando Gráfico de Importancia de Variables ---\n")
importancia_plot <- ajuste_rf %>%
  extract_fit_parsnip() %>%
  vip(num_features = 10, geom = "point") +
  theme_minimal() +
  labs(title = "Top 10 Índices Espectrales (Random Forest)")

print(importancia_plot)
