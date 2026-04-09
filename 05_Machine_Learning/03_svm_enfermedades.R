# ============================================================================
# 03_svm_enfermedades.R
# Detección de Infecciones Fúngicas mediante SVM y Tidymodels
# ============================================================================
# Contexto: Detección de severidad de roya (Infectado vs Sano)
# usando fluorescencia de clorofila y metabolitos secundarios.
# ============================================================================

library(tidymodels)
library(tidyverse)

set.seed(123)

# 1. SIMULACIÓN DE DATOS (Variables Fisiológicas + Severidad)
n <- 200
datos_fongos <- tibble(
  fluo = runif(n, 0.4, 0.8),      # Fluorescencia Fv/Fm
  fenoles = runif(n, 10, 50),     # Fenoles totales
  carotenoides = runif(n, 2, 8)   # Contenido de carotenoides
) %>%
  mutate(
    z = 12 - fluo*25 + fenoles*0.2 - carotenoides*1.2,
    prob = 1 / (1 + exp(-z)),
    severidad = ifelse(prob > 0.5, "Infectado", "Sano") %>% 
      factor(levels = c("Sano", "Infectado"))
  )

# 2. PARTICIÓN DE DATOS (Tidymodels: initial_split)
set.seed(456)
split_fongos <- initial_split(datos_fongos, prop = 0.75, strata = severidad)
entrena <- training(split_fongos)
testeo <- testing(split_fongos)

# 3. RECETA Y ESPECIFICACIÓN DEL MODELO (Parsnip: svm_poly o svm_rbf)
receta_svm <- recipe(severidad ~ ., data = entrena) %>%
  step_normalize(all_predictors())

espec_svm <- svm_rbf(
  cost = 1,
  rbf_sigma = 0.1
) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# 4. WORKFLOW Y AJUSTE FINAL
workflow_svm <- workflow() %>%
  add_recipe(receta_svm) %>%
  add_model(espec_svm)

ajuste_svm <- workflow_svm %>%
  fit(data = entrena)

# 5. EVALUACIÓN: CURVA ROC Y AUC (Yardstick)
preds_svm <- predict(ajuste_svm, testeo, type = "prob") %>%
  bind_cols(testeo)

cat("\n--- Reporte de Desempeño (SVM Radial) ---\n")
print(roc_auc(preds_svm, truth = severidad, .pred_Sano))

# Visualización Curva ROC
roc_plot <- roc_curve(preds_svm, truth = severidad, .pred_Sano) %>%
  autoplot() +
  labs(title = "Curva ROC - Detección de Roya Fúngica (SVM Radial)") +
  theme_minimal()

print(roc_plot)
