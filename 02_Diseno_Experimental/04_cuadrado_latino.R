# ============================================================================
# 04_cuadrado_latino.R — DISEÑO EN CUADRADO LATINO (LSD)
# Control de Doble Gradiente Ambiental (Filas y Columnas) en Agronomía
# ============================================================================
# Aplicación: Ensayos de invernadero con gradiente de luz (filas) y ventilación (columnas),
# o experimentos de nutrición animal (periodo x animal).
# ============================================================================

library(tidyverse)
library(emmeans)
library(car)
library(lme4)
if (requireNamespace("easyModels", quietly = TRUE)) library(easyModels)

set.seed(42)

# 1. SIMULACIÓN DE DATOS (CUADRADO LATINO 4x4)
# Tratamientos: 4 dosis de biofertilizante (T0, T1, T2, T3)
# Bloqueo doble: 4 Filas (gradiente térmico) y 4 Columnas (gradiente lumínico)

filas <- factor(rep(1:4, each = 4))
columnas <- factor(rep(1:4, times = 4))

# Distribución balanceada de tratamientos en cuadrado latino
tratamientos <- factor(c(
  "T0", "T1", "T2", "T3",
  "T1", "T2", "T3", "T0",
  "T2", "T3", "T0", "T1",
  "T3", "T0", "T1", "T2"
))

# Simular rendimiento de grano (ton/ha)
rendimiento_base <- case_when(
  tratamientos == "T0" ~ 4.2,
  tratamientos == "T1" ~ 5.8,
  tratamientos == "T2" ~ 6.9,
  tratamientos == "T3" ~ 6.3
)

efecto_fila <- as.numeric(filas) * 0.3
efecto_col <- as.numeric(columnas) * 0.2
ruido <- rnorm(16, mean = 0, sd = 0.35)

datos_latino <- data.frame(
  Fila = filas,
  Columna = columnas,
  Tratamiento = tratamientos,
  Rendimiento = round(rendimiento_base + efecto_fila + efecto_col + ruido, 2)
)

print(head(datos_latino, 8))

# ============================================================================
# 2. RUTA PEDAGÓGICA (Clásica frecuentista con lm y aov)
# ============================================================================
cat("\n--- [1] ANOVA CLÁSICO DE CUADRADO LATINO (Efectos Fijos) ---\n")
mod_lm <- lm(Rendimiento ~ Fila + Columna + Tratamiento, data = datos_latino)
anova_latino <- anova(mod_lm)
print(anova_latino)

# Comparaciones de medias con Tukey
emm_lat <- emmeans(mod_lm, ~ Tratamiento)
pairs_lat <- pairs(emm_lat, adjust = "tukey")
cat("\n--- [2] COMPARACIONES MÚLTIPLES TUKEY ---\n")
print(pairs_lat)

# ============================================================================
# 3. RUTA MODERNA DE EFECTOS MIXTOS (LMM con easyModels)
# ============================================================================
cat("\n--- [3] AJUSTE RÁPIDO CON easyModels ---\n")
# easyModels ajusta filas y columnas como efectos aleatorios cruzados: (1|Fila) + (1|Columna)
if (requireNamespace("easyModels", quietly = TRUE)) {
  mod_easy <- analizar_latino(
    datos = datos_latino,
    formula_fijos = Rendimiento ~ Tratamiento,
    fila = "Fila",
    columna = "Columna",
    diagnosticos = FALSE
  )
  
  # Auditoría de supuestos
  verificar_supuestos(mod_easy)
  
  # Gráfico de barras de publicación con letras Tukey (CLD)
  p_lat <- graficar_predichos(
    modelo = mod_easy,
    predictor = "Tratamiento",
    tipo_grafico = "barras",
    mostrar_letras = TRUE,
    paleta = "teal",
    titulo = "Efecto de Biofertilizantes en Diseño Cuadrado Latino",
    eje_x = "Dosis de Biofertilizante",
    eje_y = "Rendimiento Estimado (ton/ha)"
  )
  print(p_lat)
}

# ============================================================================
# 4. GRÁFICO PERSONALIZADO GGPLOT2 (Mapa de calor del diseño)
# ============================================================================
p_mapa <- ggplot(datos_latino, aes(x = Columna, y = Fila, fill = Tratamiento)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = paste0(Tratamiento, "\n(", Rendimiento, ")")), 
            color = "white", fontface = "bold", size = 4) +
  scale_fill_manual(values = c("T0" = "#7F8C8D", "T1" = "#2980B9", "T2" = "#27AE60", "T3" = "#E67E22")) +
  labs(
    title = "Distribución Espacial en Cuadrado Latino (4x4)",
    subtitle = "Tratamiento y (Rendimiento en ton/ha) por celda",
    x = "Columna (Gradiente Lumínico)",
    y = "Fila (Gradiente Térmico)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

print(p_mapa)

cat("\n--- Script 04 Cuadrado Latino finalizado correctamente ---\n")
