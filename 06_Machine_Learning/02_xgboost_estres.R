# ============================================================================
# 02_xgboost_estres.R (VERSIÓN EXPERTA CORREGIDA)
# Clasificación de Estrés Hídrico mediante XGBoost y Tidymodels
# ============================================================================
# Contexto: Clasificación de 3 niveles de estrés (Bajo, Medio, Alto)
# basándose en variables climáticas y fisiológicas (VPD, Leaf_Temp, gs).
# ============================================================================

library(tidymodels)
library(xgboost)
library(tidyverse)

set.seed(123)

cat("\n", rep("=", 60), "\n")
cat("\n", rep("=", 60), "\n")
cat("\n", rep("=", 60), "\n")
cat("1. SIMULACIÓN DE DATOS (Factores climáticos + Nivel de Estrés)\n")
cat(rep("=", 60), "\n")

n <- 300
datos_estres <- tibble(
  vpd = runif(n, 1, 4),           # Déficit de Presión de Vapor
  temp_hoja = runif(n, 20, 38),   # Temperatura de la hoja
  gs = runif(n, 0.1, 0.5)         # Conductancia estomática
) %>%
  mutate(
    # Creamos un índice continuo simulado con lógica biológica
    # (VPD alto + Temp alta - gs alta = Más estrés)
    indice_bruto = (vpd * 2) + (temp_hoja * 0.5) - (gs * 10),
    # Lo estandarizamos para que sea fácil cortarlo en 3 pedazos
    indice_estandar = scale(indice_bruto)[,1] + rnorm(n, 0, 0.5),
    
    # Forzamos una distribución equilibrada en 3 clases
    estres_cat = case_when(
      indice_estandar < -0.5 ~ "Bajo",
      indice_estandar > 0.5 ~ "Alto",
      TRUE ~ "Medio"
    ) %>% factor(levels = c("Bajo", "Medio", "Alto"))
  ) %>%
  dplyr::select(vpd, temp_hoja, gs, estres_cat) 

print(table(datos_estres$estres_cat))

cat("\n", rep("=", 60), "\n")
cat("2. PARTICIÓN DE DATOS (Training / Testing)\n")
cat(rep("=", 60), "\n")

set.seed(456)
# Dividimos 80% para entrenar, 20% para probar (estratificando por clase)
split_estres <- initial_split(datos_estres, prop = 0.8, strata = estres_cat)
entrena <- training(split_estres)
testeo <- testing(split_estres)

cat("Datos de entrenamiento:", nrow(entrena), "filas\n")
cat("Datos de prueba:", nrow(testeo), "filas\n")

cat("\n", rep("=", 60), "\n")
cat("3. RECETA Y ESPECIFICACIÓN DEL MODELO (XGBoost)\n")
cat(rep("=", 60), "\n")

# La 'receta' prepara los datos. Para XGBoost, centramos y escalamos los predictores.
receta_xgb <- recipe(estres_cat ~ ., data = entrena) %>%
  step_normalize(all_numeric_predictors())

# Especificamos el modelo: Dejamos profundidad del árbol y tasa de aprendizaje para "afinar" (tune)
modelo_xgb <- boost_tree(
  trees = 100,             # Número de árboles
  tree_depth = tune(),     # Lo afinaremos luego
  learn_rate = tune()      # Lo afinaremos luego
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# Unimos la receta y el modelo en un flujo de trabajo (workflow)
flujo_xgb <- workflow() %>%
  add_recipe(receta_xgb) %>%
  add_model(modelo_xgb)

cat("\n", rep("=", 60), "\n")
cat("4. TUNEO DE HIPERPARÁMETROS (Validación Cruzada)\n")
cat(rep("=", 60), "\n")

set.seed(789)
# Creamos 5 pliegues para la validación cruzada
pliegues <- vfold_cv(entrena, v = 5, strata = estres_cat)

# Creamos una grilla de posibles valores para afinar
grilla <- grid_regular(tree_depth(), learn_rate(), levels = 3)

cat("Buscando los mejores hiperparámetros (Esto puede tardar unos segundos)...\n")
tuneo_res <- tune_grid(
  flujo_xgb,
  resamples = pliegues,
  grid = grilla
)

cat("\n", rep("=", 60), "\n")
cat("5. EVALUACIÓN FINAL CON DATOS DE PRUEBA\n")
cat(rep("=", 60), "\n")

# Seleccionamos la mejor combinación basada en la exactitud (accuracy)
mejor_xgb <- select_best(tuneo_res, metric = "accuracy")
cat("Mejores parámetros encontrados:\n")
print(mejor_xgb)

# Finalizamos el flujo con esos mejores parámetros
flujo_final <- finalize_workflow(flujo_xgb, mejor_xgb)

# Entrenamos por última vez con todos los datos de entrenamiento
# y evaluamos inmediatamente con los datos de prueba
ajuste_final <- last_fit(flujo_final, split_estres)

# Extraemos las métricas finales (¡Estas son las que reportas en el paper!)
metricas_finales <- collect_metrics(ajuste_final)
print(metricas_finales)

cat("\n", rep("=", 60), "\n")
cat("6. VISUALIZACIÓN: MATRIZ DE CONFUSIÓN\n")
cat(rep("=", 60), "\n")

# Extraemos las predicciones del modelo sobre el conjunto de prueba
predicciones <- collect_predictions(ajuste_final)

# Creamos la Matriz de Confusión como un mapa de calor (Heatmap)
grafico_confusion <- conf_mat(predicciones, truth = estres_cat, estimate = .pred_class) %>%
  autoplot(type = "heatmap") +
  labs(title = "Matriz de Confusión: XGBoost",
       subtitle = "Rendimiento del modelo en el conjunto de prueba (Datos nunca vistos)",
       x = "Clase Predicha por la IA",
       y = "Clase Real (Verdad de Campo)") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(grafico_confusion)

cat("\n--- Script de XGBoost Finalizado. Gráfico generado en la pestaña 'Plots' ---\n")
