# ============================================================================
# 03_pca_biplot_fisiologia.R — PCA BIPLOT MULTIVARIADO EN FISIOLOGÍA VEGETAL
# Reducción de Dimensionalidad y Correlaciones Fisiológicas en Cultivos
# ============================================================================
# Aplicación: Integración de variables fisiológicas (fotosíntesis, SPAD, MDA,
# prolina, estabilidad de membrana) bajo estrés térmico en trigo durum.
# ============================================================================

library(tidyverse)
library(patchwork)

set.seed(555)

# 1. SIMULACIÓN DE DATASET MULTIDIMENSIONAL (FISIOLOGÍA + RENDIMIENTO)
# Factores: 4 Genotipos (G1_Tolerante, G2_Sensible, G3_Intermedio, G4_Tolerante)
# Tratamientos: Control (22°C) vs Estrés Térmico (38°C)
# 6 réplicas biológicas por combinación (N = 48)

genotipos <- c("G1_Tolerante", "G2_Sensible", "G3_Intermedio", "G4_Tolerante")
tratamientos <- c("Control", "Calor_Z65")
n_rep <- 6

df_base <- expand.grid(
  Genotipo = factor(genotipos, levels = genotipos),
  Tratamiento = factor(tratamientos, levels = tratamientos),
  Replica = 1:n_rep
)

# Generar variables fisiológicas correlacionadas
datos_fisiologia <- df_base %>%
  mutate(
    # Fotosíntesis neta (An: umol CO2 / m2 s)
    Fotosintesis = case_when(
      Tratamiento == "Control" ~ rnorm(n(), mean = 24, sd = 1.5),
      Genotipo %in% c("G1_Tolerante", "G4_Tolerante") ~ rnorm(n(), mean = 19, sd = 1.2),
      TRUE ~ rnorm(n(), mean = 11, sd = 1.4)
    ),
    # Conductancia estomática (gs: mol H2O / m2 s)
    Conductancia = Fotosintesis * 0.015 + rnorm(n(), 0, 0.03),
    # Contenido de clorofila (SPAD)
    SPAD = Fotosintesis * 1.5 + rnorm(n(), 15, 2.0),
    # Daño oxidativo: Malondialdehído (MDA: nmol/g MF) - Negativamente correlacionado
    MDA = case_when(
      Tratamiento == "Control" ~ rnorm(n(), mean = 8, sd = 1.0),
      Genotipo %in% c("G1_Tolerante", "G4_Tolerante") ~ rnorm(n(), mean = 14, sd = 1.5),
      TRUE ~ rnorm(n(), mean = 26, sd = 2.2)
    ),
    # Osmoprotector: Prolina libre (umol/g MF)
    Prolina = case_when(
      Tratamiento == "Control" ~ rnorm(n(), mean = 3.5, sd = 0.5),
      Genotipo %in% c("G1_Tolerante", "G4_Tolerante") ~ rnorm(n(), mean = 16.5, sd = 1.8),
      TRUE ~ rnorm(n(), mean = 7.0, sd = 1.2)
    ),
    # Índice de estabilidad de membrana (MSI: %)
    MSI = pmin(pmax(85 - (MDA * 1.8) + rnorm(n(), 0, 2.5), 30), 95),
    # Rendimiento de grano (g/planta)
    Rendimiento = (Fotosintesis * 0.4) + (MSI * 0.15) - (MDA * 0.2) + rnorm(n(), 0, 0.8)
  )

vars_numericas <- c("Fotosintesis", "Conductancia", "SPAD", "MDA", "Prolina", "MSI", "Rendimiento")
print(head(datos_fisiologia[, c("Genotipo", "Tratamiento", vars_numericas)], 6))

# ============================================================================
# 2. CÁLCULO DE PCA (Estandarización Z-score con scale = TRUE)
# ============================================================================
cat("\n--- [1] EJECUCIÓN DEL ANÁLISIS DE COMPONENTES PRINCIPALES (PCA) ---\n")
matriz_datos <- datos_fisiologia[, vars_numericas]
pca_res <- prcomp(matriz_datos, scale. = TRUE, center = TRUE)

# Varianza explicada por componente
var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
pc1_var <- round(var_exp[1], 1)
pc2_var <- round(var_exp[2], 1)

cat("Varianza explicada por PC1:", pc1_var, "%\n")
cat("Varianza explicada por PC2:", pc2_var, "%\n")
cat("Varianza acumulada (PC1 + PC2):", round(pc1_var + pc2_var, 1), "%\n")

# ============================================================================
# 3. EXTRACCIÓN DE SCORES Y CARGAS (LOADINGS)
# ============================================================================
scores_df <- as.data.frame(pca_res$x) %>%
  bind_cols(datos_fisiologia[, c("Genotipo", "Tratamiento", "Replica")])

# Cargas factoriales (vectores)
loadings_df <- as.data.frame(pca_res$rotation) %>%
  mutate(Variable = rownames(pca_res$rotation))

# Factor de escala para que las flechas encajen armónicamente en el biplot
mult_escala <- 4.5
loadings_df <- loadings_df %>%
  mutate(
    PC1_scale = PC1 * mult_escala,
    PC2_scale = PC2 * mult_escala
  )

# ============================================================================
# 4. CONSTRUCCIÓN DEL PCA BIPLOT DE NIVEL DE PUBLICACIÓN
# ============================================================================
p_biplot <- ggplot() +
  # Elipses de concentración del 95% por Tratamiento
  stat_ellipse(
    data = scores_df,
    aes(x = PC1, y = PC2, fill = Tratamiento, color = Tratamiento),
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    linetype = "dashed"
  ) +
  # Puntos de muestras individuales
  geom_point(
    data = scores_df,
    aes(x = PC1, y = PC2, color = Tratamiento, shape = Genotipo),
    size = 3.8,
    alpha = 0.9
  ) +
  # Vectores de carga de variables
  geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0, xend = PC1_scale, yend = PC2_scale),
    arrow = grid::arrow(length = grid::unit(0.25, "cm")),
    color = "#2C3E50",
    linewidth = 0.9
  ) +
  # Etiquetas de variables con fondo blanco para legibilidad
  geom_label(
    data = loadings_df,
    aes(x = PC1_scale * 1.12, y = PC2_scale * 1.12, label = Variable),
    color = "#2C3E50",
    fontface = "bold",
    size = 3.6,
    fill = "white",
    linewidth = 0.2
  ) +
  # Ejes en cero
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  # Colores y formas
  scale_color_manual(values = c("Control" = "#00A88F", "Calor_Z65" = "#E65100")) +
  scale_fill_manual(values = c("Control" = "#00A88F", "Calor_Z65" = "#E65100")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(
    title = "PCA Biplot: Respuesta Fisiológica y Rendimiento bajo Estrés Térmico",
    subtitle = paste0("Trigo Candeal (Triticum durum) | Varianza Total Explicada: ", round(pc1_var + pc2_var, 1), "%"),
    x = paste0("PC1 (", pc1_var, "% de varianza)"),
    y = paste0("PC2 (", pc2_var, "% de varianza)"),
    color = "Condición Térmica",
    fill = "Condición Térmica",
    shape = "Genotipo"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "right",
    legend.box = "vertical"
  )

print(p_biplot)

cat("\n--- Script 03 PCA Biplot finalizado correctamente ---\n")
