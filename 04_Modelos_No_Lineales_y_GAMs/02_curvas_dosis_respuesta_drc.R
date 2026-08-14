# ============================================================================
# 02_curvas_dosis_respuesta_drc.R — MODELOS DOSIS-RESPUESTA NO LINEALES (drc)
# Estandarización de EC50 / ED50 en Fitopatología y Protección de Cultivos
# ============================================================================
# Aplicación: Ensayos de inhibición micelial in vitro, eficacia de herbicidas
# y bioestimulantes mediante modelos Log-Logísticos de 4 parámetros (LL.4).
# ============================================================================

library(tidyverse)
library(drc)

set.seed(999)

# 1. SIMULACIÓN DE BIOENSAYO DE INHIBICIÓN FÚNGICA IN VITRO
# Cepa A (Sensible) vs Cepa B (Tolerante)
# Dosis de fungicida (mg/L): 0, 0.05, 0.1, 0.5, 1.0, 5.0, 10, 50, 100
# 4 réplicas por dosis y cepa

dosis_vect <- c(0, 0.05, 0.1, 0.5, 1.0, 5.0, 10, 50, 100)
n_rep <- 4

simular_curva_ll4 <- function(cepa, ec50, pendiente, max_val = 100, min_val = 0) {
  df <- expand.grid(
    Dosis = dosis_vect,
    Replica = 1:n_rep
  ) %>%
    mutate(
      Cepa = cepa,
      # Ecuación Log-Logística de 4 Parámetros (LL.4)
      inhibicion_esperada = ifelse(
        Dosis == 0,
        min_val,
        min_val + (max_val - min_val) / (1 + exp(pendiente * (log(Dosis) - log(ec50))))
      ),
      Inhibicion = pmin(pmax(inhibicion_esperada + rnorm(n(), mean = 0, sd = 4.0), 0), 100)
    )
  return(df)
}

# Cepa Sensible (EC50 = 0.8 mg/L) vs Cepa Tolerante (EC50 = 12.5 mg/L)
datos_sensible <- simular_curva_ll4("Sensible", ec50 = 0.8, pendiente = 1.4)
datos_tolerante <- simular_curva_ll4("Tolerante", ec50 = 12.5, pendiente = 1.2)
datos_bioensayo <- bind_rows(datos_sensible, datos_tolerante)

print(head(datos_bioensayo, 10))

# ============================================================================
# 2. AJUSTE DE MODELO DOSIS-RESPUESTA CON drc (fct = LL.4())
# ============================================================================
cat("\n--- [1] AJUSTE DEL MODELO LOG-LOGÍSTICO DE 4 PARÁMETROS (LL.4) ---\n")
# LL.4: y = c + (d - c) / (1 + exp(b * (log(x) - log(e))))
# b: pendiente (slope), c: límite inferior, d: límite superior, e: EC50

mod_drc <- drm(
  Inhibicion ~ Dosis,
  curveid = Cepa,
  data = datos_bioensayo,
  fct = LL.4(names = c("Pendiente", "Inferior", "Superior", "EC50"))
)

print(summary(mod_drc))

# ============================================================================
# 3. CÁLCULO DE EC50 / ED50 Y PRUEBA DE FACTOR DE RESISTENCIA
# ============================================================================
cat("\n--- [2] ESTIMACIÓN DE CONCENTRACIÓN EFECTIVA AL 50% (EC50) ---\n")
ec50_estimado <- ED(mod_drc, c(50), interval = "delta")
print(ec50_estimado)

cat("\n--- [3] PRUEBA DE RATIO DE RESISTENCIA (EDcomp: Tolerante vs Sensible) ---\n")
ratio_resistencia <- EDcomp(mod_drc, c(50, 50), interval = "delta")
print(ratio_resistencia)

# ============================================================================
# 4. GRÁFICO CIENTÍFICO DE PUBLICACIÓN CON GGPLOT2
# ============================================================================
# Generar malla fina de predicción para trazar curvas sigmoidales suaves
dosis_grid <- exp(seq(log(0.01), log(120), length.out = 200))
pred_grid <- expand.grid(
  Dosis = dosis_grid,
  Cepa = unique(datos_bioensayo$Cepa)
)

pred_grid$Inhibicion_Pred <- predict(mod_drc, newdata = pred_grid)

# Extraer valores puntuales de EC50 para las líneas de referencia
ec50_vals <- data.frame(
  Cepa = c("Sensible", "Tolerante"),
  EC50 = c(ec50_estimado[1, 1], ec50_estimado[2, 1])
)

p_drc <- ggplot() +
  # Datos experimentales (puntos promedio ± SE)
  stat_summary(
    data = datos_bioensayo,
    aes(x = Dosis, y = Inhibicion, color = Cepa),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.1,
    linewidth = 0.7
  ) +
  stat_summary(
    data = datos_bioensayo,
    aes(x = Dosis, y = Inhibicion, color = Cepa),
    fun = mean,
    geom = "point",
    size = 3.5
  ) +
  # Curva ajustada no lineal
  geom_line(
    data = pred_grid,
    aes(x = Dosis, y = Inhibicion_Pred, color = Cepa),
    linewidth = 1.2
  ) +
  # Líneas de referencia para EC50 (50% de inhibición)
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_segment(
    data = ec50_vals,
    aes(x = EC50, xend = EC50, y = 0, yend = 50, color = Cepa),
    linetype = "dotted",
    linewidth = 0.8
  ) +
  scale_x_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")
  ) +
  scale_color_manual(values = c("Sensible" = "#00A88F", "Tolerante" = "#E65100")) +
  labs(
    title = "Curvas Dosis-Respuesta de Inhibición Fúngica",
    subtitle = paste0("Modelo Log-Logístico LL.4 | Factor de Resistencia = ", round(ratio_resistencia[1, 1], 1), "x"),
    x = "Concentración de Fungicida (mg/L, escala log)",
    y = "Inhibición del Crecimiento Micelial (%)",
    color = "Aislado Fúngico"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "grey92", linetype = "dashed")
  )

print(p_drc)

cat("\n--- Script 02 Dosis-Respuesta finalizado correctamente ---\n")
