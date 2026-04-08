# --------------------------------------------------------------------
# 05_Machine_Learning
# Random Forest con tidymodels para predecir rendimiento bajo estrés térmico
# Integración de datos fenotípicos y respuestas multi-ómicas simuladas
# --------------------------------------------------------------------

library(tidymodels)
library(tidyverse)
library(themis)
library(vip)

set.seed(910)

# ============================================================================
# 1. SIMULACIÓN DE DATOS MULTI-ÓMICOS Y FENOTÍPICOS
# ============================================================================

n <- 500

# Datos fenotípicos ambientales
datos_feno <- tibble(
  temperatura = rnorm(n, 32, 4),
  humedad_suelo = runif(n, 10, 45),
  radiacion_solar = rnorm(n, 25, 6),
  velocidad_viento = rnorm(n, 12, 4)
)

# Datos transcriptómicos simulados (expresión génica)
datos_transcriptoma <- tibble(
  gen_HSP = rnorm(n, 50, 15),
  gen_DREB = rnorm(n, 40, 12),
  gen_NAC = rnorm(n, 35, 10),
  gen_bZIP = rnorm(n, 30, 8)
)

# Datos metabolómicos simulados
datos_metaboloma <- tibble(
  prolina = rnorm(n, 25, 8),
  glicina_betaina = rnorm(n, 18, 6),
  azucares_solubles = rnorm(n, 45, 12),
  poliaminas = rnorm(n, 12, 4)
)

# Genotipos
genotipos <- factor(sample(paste0("G", 1:10), n, replace = TRUE))

# Respuesta objetivo: rendimiento bajo estrés térmico
datos <- bind_cols(
  datos_feno,
  datos_transcriptoma,
  datos_metaboloma,
  tibble(genotipo = genotipos)
)

datos <- datos %>%
  mutate(
    rendimiento_estres = 100 -
      1.2 * temperatura +
      0.8 * humedad_suelo +
      0.4 * clorofila +
      0.3 * gen_HSP +
      0.25 * gen_DREB +
      0.2 * prolina +
      0.15 * glicina_betaina +
      ifelse(genotipo %in% c("G2", "G5", "G8"), 8, 0) +
      rnorm(n, 0, 6)
  ) %>%
  mutate(
    rendimiento_estres = pmax(rendimiento_estres, 0),
    genotipo = as.factor(genotipo)
  )

# ============================================================================
# 2. PREPROCESAMIENTO Y DIVISIÓN DE DATOS
# ============================================================================

set.seed(910)
split <- initial_split(datos, prop = 0.8, strata = genotipo)
train <- training(split)
test <- testing(split)

# ============================================================================
# 3. RECETA DE PREPROCESAMIENTO
# ============================================================================

rf_recipe <- recipe(rendimiento_estres ~ ., data = train) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_corr(all_numeric_predictors(), threshold = 0.9) %>%
  step_zv(all_predictors())

# ============================================================================
# 4. ESPECIFICACIÓN DEL MODELO Y AJUSTE DE HIPERPARÁMETROS
# ============================================================================

rf_spec <- rand_forest(
  trees = tune(),
  min_n = tune(),
  mtry = tune()
) %>%
  set_mode("regression") %>%
  set_engine("ranger", importance = "impurity")

# ============================================================================
# 5. VALIDACIÓN CRUZADA
# ============================================================================

set.seed(910)
folds <- vfold_cv(train, v = 5, strata = genotipo)

rf_grid <- grid_regular(
  trees = values(100, 500),
  min_n = values(5, 15),
  mtry = values(3, 10),
  levels = list(trees = 3, min_n = 2, mtry = 3)
)

rf_wf <- workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(rf_spec)

set.seed(910)
rf_tune <- tune_grid(
  rf_wf,
  resamples = folds,
  grid = rf_grid,
  metrics = metric_set(rmse, rsq, mae),
  control = control_grid(verbose = TRUE, save_pred = TRUE)
)

# ============================================================================
# 6. SELECCIÓN DE MEJORES HIPERPARÁMETROS
# ============================================================================

mejores_params <- select_best(rf_tune, "rmse")
cat("\n=== Mejores Hiperparámetros ===\n")
print(mejores_params)

# ============================================================================
# 7. ENTRENAMIENTO DEL MODELO FINAL
# ============================================================================

rf_final <- finalize_model(rf_spec, mejores_params) %>%
  fit(rendimiento_estres ~ ., data = train)

# ============================================================================
# 8. EVALUACIÓN EN TEST
# ============================================================================

pred_test <- predict(rf_final, test) %>%
  bind_cols(test)

metricas_test <- pred_test %>%
  metrics(truth = rendimiento_estres, estimate = .pred)

cat("\n=== Métricas en Test ===\n")
print(metricas_test)

# Gráfico de valores observados vs predichos
plot_obs_pred <- pred_test %>%
  ggplot(aes(x = rendimiento_estres, y = .pred)) +
  geom_point(alpha = 0.5, color = "#2E86AB", size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E74C3C") +
  labs(
    title = "Random Forest: Observado vs Predicho",
    subtitle = "Rendimiento bajo estrés térmico en trigo duro",
    x = "Rendimiento Observado",
    y = "Rendimiento Predicho"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("05_Machine_Learning/05_observed_vs_predicted.png", plot_obs_pred, width = 7, height = 5, dpi = 300)
print(plot_obs_pred)

# ============================================================================
# 9. IMPORTANCIA DE VARIABLES
# ============================================================================

importancia <- vip::vi(rf_final) %>%
  arrange(desc(Importance))

cat("\n=== Importancia de Variables (Top 10) ===\n")
print(head(importancia, 10))

plot_importancia <- importancia %>%
  head(15) %>%
  mutate(
    Variable = str_remove(Variable, "^gen_"),
    Variable = fct_reorder(Variable, Importance)
  ) %>%
  ggplot(aes(x = Importance, y = Variable, fill = Importance)) +
  geom_col(alpha = 0.8) +
  scale_fill_gradient(low = "#A8D5E2", high = "#2E86AB") +
  labs(
    title = "Importancia de Variables - Random Forest",
    subtitle = "Predictores del rendimiento bajo estrés térmico",
    x = "Importancia (Gini)",
    y = "Variable"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40"),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave("05_Machine_Learning/05_importancia_variables.png", plot_importancia, width = 8, height = 6, dpi = 300)
print(plot_importancia)

# ============================================================================
# 10. RESUMEN FINAL
# ============================================================================

cat("\n=== Resumen del Modelo ===\n")
cat("Modelo: Random Forest (tidymodels)\n")
cat("Datos: ", nrow(train), " entrenamiento, ", nrow(test), " test\n", sep = "")
cat("Hiperparámetros óptimos:\n")
cat("  - trees: ", mejores_params$trees, "\n", sep = "")
cat("  - min_n: ", mejores_params$min_n, "\n", sep = "")
cat("  - mtry: ", mejores_params$mtry, "\n", sep = "")
cat("RMSE en test: ", round(metricas_test$estimate[1], 3), "\n", sep = "")
cat("R² en test: ", round(metricas_test$estimate[2], 3), "\n", sep = "")
