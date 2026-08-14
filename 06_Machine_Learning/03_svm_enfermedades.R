# ============================================================================
# 03_svm_enfermedades.R (VERSIÓN EXPERTA CORREGIDA)
# Detección de Infecciones Fúngicas con SVM y Tidymodels
# ============================================================================

library(tidymodels)
library(tidyverse)
library(kernlab) # <- El motor matemático que faltaba para SVM

set.seed(123)

cat("\n", rep("=", 60), "\n")
cat("1. SIMULACIÓN DE DATOS LÍMPIOS\n")
cat(rep("=", 60), "\n")

n <- 200
datos_fongos <- tibble(
  fluo = runif(n, 0.4, 0.8),      # Fluorescencia Fv/Fm
  fenoles = runif(n, 10, 50),     # Fenoles totales
  carotenoides = runif(n, 2, 8)   # Contenido de carotenoides
) %>%
  mutate(
    z = 12 - fluo*25 + fenoles*0.2 - carotenoides*1.2,
    prob = 1 / (1 + exp(-z)),
    # Ponemos "Infectado" primero, porque es el "evento" que queremos detectar
    severidad = ifelse(prob > 0.5, "Infectado", "Sano") %>% 
      factor(levels = c("Infectado", "Sano"))
  ) %>%
  # ¡CORRECCIÓN CRÍTICA! Eliminamos z y prob para que el modelo no haga trampa
  dplyr::select(-z, -prob)

print(table(datos_fongos$severidad))

cat("\n", rep("=", 60), "\n")
cat("2. PARTICIÓN DE DATOS (Training / Testing)\n")
cat(rep("=", 60), "\n")

set.seed(456)
split_fongos <- initial_split(datos_fongos, prop = 0.75, strata = severidad)
entrena <- training(split_fongos)
testeo <- testing(split_fongos)

cat("\n", rep("=", 60), "\n")
cat("3. RECETA Y ESPECIFICACIÓN DEL MODELO (SVM Radial)\n")
cat(rep("=", 60), "\n")

receta_svm <- recipe(severidad ~ ., data = entrena) %>%
  step_normalize(all_predictors())

# Usamos SVM con un kernel radial (RBF)
espec_svm <- svm_rbf(cost = 1, rbf_sigma = 0.1) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

cat("Ajustando el modelo SVM...\n")
workflow_svm <- workflow() %>%
  add_recipe(receta_svm) %>%
  add_model(espec_svm)

ajuste_svm <- workflow_svm %>%
  fit(data = entrena)

cat("\n", rep("=", 60), "\n")
cat("4. EVALUACIÓN: CURVA ROC Y AUC (Yardstick)\n")
cat(rep("=", 60), "\n")

# Hacemos predicciones sobre los datos de prueba (testeo)
preds_svm <- predict(ajuste_svm, testeo, type = "prob") %>%
  bind_cols(testeo)

cat("\n--- Área Bajo la Curva (AUC) ---\n")
auc_val <- roc_auc(preds_svm, truth = severidad, .pred_Infectado)
print(auc_val)

# Visualización Curva ROC
roc_plot <- roc_curve(preds_svm, truth = severidad, .pred_Infectado) %>%
  autoplot() +
  labs(title = "Curva ROC - Detección de Roya Fúngica",
       subtitle = paste("Modelo SVM Radial | AUC =", round(auc_val$.estimate, 3)),
       x = "Tasa de Falsos Positivos (1 - Especificidad)",
       y = "Tasa de Verdaderos Positivos (Sensibilidad)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(roc_plot)
cat("\n--- Script SVM Finalizado. Gráfico generado en la pestaña 'Plots' ---\n")
