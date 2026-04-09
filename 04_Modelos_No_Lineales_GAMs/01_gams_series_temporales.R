# ============================================================================
# 01_gams_series_temporales.R (FOTOSÍNTESIS Y ESTRÉS TÉRMICO)
# Modelos Aditivos Generalizados (GAMs) con mgcv y gratia
# ============================================================================
# Los GAMs son la herramienta definitiva para curvas biológicas no lineales
# (como la respuesta térmica), ya que permiten que los datos definan la
# forma de la curva sin forzar una estructura lineal o polinómica.
# ============================================================================

library(mgcv)      # Motor principal para GAMs
library(ggplot2)   # Visualización de alta calidad
library(gratia)    # Diagnóstico y visualización moderna de GAMs (imprescindible)
library(dplyr)     # Manipulación de datos

set.seed(456)

# 1. SIMULACIÓN DE DATOS: CURVA TÉRMICA DE FOTOSÍNTESIS
# Escenario: 5 plantas medidas de 20°C a 40°C. 
# La fotosíntesis sube hasta un óptimo (~28°C) y cae drásticamente por estrés.

crear_curva <- function(temp, id) {
  # Curva biológica teórica (Campana asimétrica)
  base <- 30 * exp(-0.5 * ((temp - 28) / 5)^2) 
  # Añadir variabilidad individual por planta
  variabilidad <- rnorm(1, 0, 2) + (temp * rnorm(1, 0, 0.1))
  return(base + variabilidad + rnorm(length(temp), 0, 1))
}

temperaturas <- seq(20, 40, length.out = 50)
plantas <- paste0("Planta_", 1:5)

datos_gam <- expand.grid(temp = temperaturas, planta = factor(plantas)) %>%
  group_by(planta) %>%
  mutate(fotosintesis = crear_curva(temp, planta)) %>%
  ungroup()

# ============================================================================
# 2. GAM BÁSICO: CURVA TÉRMICA GENERAL
# ============================================================================
# Modelamos la tendencia promedio de todas las plantas.
# s(temp) indica un 'smooth' (suavizado) basado en splines.

mod_gam <- gam(fotosintesis ~ s(temp), data = datos_gam, method = "REML")

cat("\n--- [1] RESUMEN DEL GAM BÁSICO ---\n")
print(summary(mod_gam))

# ¿Cómo interpretar? 
# El 'edf' (grados de libertad efectivos) indica la complejidad: 
# edf=1 es lineal, edf>1 es no lineal. Aquí esperamos un edf alto (~4-6).

# Visualización con gratia
cat("\n--- Generando Visualización del Término Suavizado ---\n")
draw(mod_gam, residuals = TRUE) 

# ============================================================================
# 3. GAMM: CURVAS INDIVIDUALES (Factor-Smooth Interactions)
# ============================================================================
# Aquí usamos bs='fs' (Factor Smooth). 
# Esto permite que CADA PLANTA tenga su propia curva (intercepto + forma),
# tratándolas como efectos aleatorios. Es un "Mixed GAM".

mod_gamm <- gam(fotosintesis ~ s(temp, planta, bs = "fs"), 
                data = datos_gam, method = "REML")

cat("\n--- [2] RESUMEN DEL GAMM (CURVAS INDIVIDUALES) ---\n")
print(summary(mod_gamm))

# Visualización de la variabilidad entre plantas
draw(mod_gamm)

# ============================================================================
# 4. VISUALIZACIÓN COMPARATIVA (GGPLOT2)
# ============================================================================

# Predecir valores para graficar las curvas ajustadas
datos_pred <- datos_gam %>%
  mutate(pred = predict(mod_gamm))

ggplot(datos_gam, aes(x = temp, y = fotosintesis, color = planta)) +
  geom_point(alpha = 0.4) +
  geom_line(data = datos_pred, aes(y = pred), linewidth = 1) +
  labs(title = "Modelado GAM de la Tasa Fotosintética",
       subtitle = "Curvas individuales por planta bajo gradiente térmico (20-40 °C)",
       x = "Temperatura de Hoja (°C)",
       y = "Fotosíntesis Neta (μmol CO2 m-2 s-1)") +
  theme_minimal() +
  scale_color_viridis_d()

# ============================================================================
# 5. DIAGNÓSTICO DEL MODELO (APPRAISE)
# ============================================================================
# gratia::appraise() es el equivalente a performance::check_model() para GAMs.
# Verifica si los splines son suficientes (k-index) y la normalidad.

cat("\n--- [3] DIAGNÓSTICO TÉCNICO DEL GAM ---\n")
appraise(mod_gamm)

cat("\n--- Script de GAMs finalizado correctamente ---\n")
