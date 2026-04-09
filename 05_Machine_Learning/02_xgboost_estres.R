# ============================================================================
# 02_xgboost_estres.R
# Clasificación de Estrés Hídrico mediante XGBoost y Tidymodels
# ============================================================================
# Contexto: Clasificación de 3 niveles de estrés (Bajo, Medio, Alto)
# basándose en variables climáticas y fisiológicas (VPD, Leaf_Temp, Stomatal_Cond).
# ============================================================================

library(tidymodels)
library(xgboost)
library(tidyverse)

set.seed(123)

# 1. SIMULACIÓN DE DATOS (Factores climáticos + Nivel de Estrés)
n <- 300
datos_estres <- tibble(
  vpd = runif(n, 1, 4),           # Déficit de Presión de Vapor
  temp_hoja = runif(n, 20, 38),   # Temperatura de la hoja
  gs = runif(n, 0.1, 0.5)         # Conductancia estomática
) %>%
  mutate(
    prob_bajo = 1 / (1 + exp(-(4 - vpd*1.2 + gs*10))),
    prob_alto = 1 / (1 + exp(-(vpd*1.8 - gs*8 - 3))),
    estres_cat = case_when(
      prob_bajo > 0.6 ~ "Bajo",
      prob_alto > 0.6 ~ "Alto",
      TRUE ~ "Medio"
    ) %>% factor(levels = c("Bajo", "Medio", "Alto"))
  ) %>%
  select(vpd, temp_hoja, gs, estres_cat)

# 2. PARTICIÓN DE DATOS (Tidymodels: initial_split)
set.seed(456)
split_estres <- initial_split(datos_estres, prop = 0.8, strata = estres_cat)
entrena <- training(split_estres)
testeo <- testing(split_estres)

# 3. RECETA Y ESPECIFICACIÓN (Parsnip: tune)
receta_xgb <- recipe(estres_cat ~ ., data = entrena) %>%
  step_normalize(all_predictors())

espec_xgb <- boost_tree(
  trees = tune(),
  mtry = tune(),
  learn_rate = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# 4. TUNEO DE HIPERPARÁMETROS (Tune Grid)
workflow_xgb <- workflow() %>%
  add_recipe(receta_xgb) %>%
  add_model(espec_xgb)

vfold_xgb <- vfold_cv(entrena, v = 5)

# Grilla simple para ejemplo
grilla_xgb <- grid_regular(
  trees(range = c(100, 500)),
  mtry(range = c(1, 3)),
  learn_rate(range = c(-2, -1)),
  levels = 3
)

set.seed(789)
tuneo_xgb <- tune_grid(
  workflow_xgb,
  resamples = vfold_xgb,
  grid = grilla_xgb
)

# 5. MEJOR MODELO Y AJUSTE FINAL
mejor_xgb <- select_best(tuneo_xgb, metric = "accuracy")

final_xgb <- finalize_workflow(workflow_xgb, mejor_xgb) %>%
  fit(data = entrena)

# 6. EVALUACIÓN Y MATRIZ DE CONFUSIÓN (Autoplot)
preds_xgb <- predict(final_xgb, testeo) %>%
  bind_cols(testeo)

cat("\n--- Reporte de Desempeño (XGBoost) ---\n")
print(accuracy(preds_xgb, truth = estres_cat, estimate = .pred))

# Matriz de Confusión Visual
conf_mat_plot <- conf_mat(preds_xgb, truth = estres_cat, estimate = .pred) %>%
  autoplot(type = "heatmap") +
  labs(title = "Matriz de Confusión (Niveles de Estrés Hídrico)") +
  theme_minimal()

print(conf_mat_plot)
