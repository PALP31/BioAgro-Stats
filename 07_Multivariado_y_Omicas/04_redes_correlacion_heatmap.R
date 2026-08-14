# ============================================================================
# 04_redes_correlacion_heatmap.R — MAPAS DE CALOR Y MATRICES DE CORRELACIÓN
# Integración Multi-Atributo (Ionómica, Fisiología y Rendimiento en Cultivos)
# ============================================================================
# Aplicación: Visualización de relaciones inter-rasgo (K+, Na+, Prolina, SPAD,
# MDA, Rendimiento) mediante Clustering Jerárquico y Matrices de Correlación.
# ============================================================================

library(tidyverse)
library(corrplot)
library(pheatmap)

set.seed(777)

# 1. SIMULACIÓN DE DATOS (IONÓMICA + FISIOLOGÍA BAJO ESTRÉS SALINO/TÉRMICO)
# 40 muestras experimentales de trigo
n <- 40

# Variables ionómicas y fisiológicas correlacionadas
K_foliar <- rnorm(n, mean = 28, sd = 4.5)
Na_foliar <- pmax(50 - (K_foliar * 1.2) + rnorm(n, 0, 3.5), 2)
Ratio_K_Na <- round(K_foliar / Na_foliar, 2)
Fotosintesis <- (K_foliar * 0.6) - (Na_foliar * 0.25) + rnorm(n, 10, 1.5)
SPAD <- (Fotosintesis * 1.3) + rnorm(n, 12, 2.0)
MDA <- (Na_foliar * 0.7) - (K_foliar * 0.3) + rnorm(n, 15, 2.5)
Prolina <- (Na_foliar * 0.45) + rnorm(n, 5, 1.2)
Rendimiento <- (Fotosintesis * 0.5) + (Ratio_K_Na * 1.8) - (MDA * 0.15) + rnorm(n, 0, 0.8)

datos_multiomica <- data.frame(
  K_Foliar = round(K_foliar, 2),
  Na_Foliar = round(Na_foliar, 2),
  Ratio_K_Na = Ratio_K_Na,
  Fotosintesis = round(Fotosintesis, 2),
  SPAD = round(SPAD, 2),
  MDA = round(MDA, 2),
  Prolina = round(Prolina, 2),
  Rendimiento = round(Rendimiento, 2)
)

print(head(datos_multiomica, 6))

# ============================================================================
# 2. CÁLCULO DE MATRIZ DE CORRELACIÓN Y P-VALORES
# ============================================================================
cat("\n--- [1] CÁLCULO DE CORRELACIONES DE PEARSON Y SIGNIFICANCIA ---\n")
mat_cor <- cor(datos_multiomica, method = "pearson")

# Función interna para calcular matriz de p-valores
calcular_p_mat <- function(df) {
  cols <- ncol(df)
  nombres <- colnames(df)
  p_mat <- matrix(NA, cols, cols)
  colnames(p_mat) <- nombres
  rownames(p_mat) <- nombres
  for (i in 1:(cols - 1)) {
    for (j in (i + 1):cols) {
      test <- cor.test(df[[i]], df[[j]])
      p_mat[i, j] <- p_mat[j, i] <- test$p.value
    }
  }
  diag(p_mat) <- 0
  return(p_mat)
}

p_matrix <- calcular_p_mat(datos_multiomica)

# ============================================================================
# 3. GRÁFICO 1: CORRPLOT CON SIGNIFICANCIA ESTADÍSTICA (P < 0.05)
# ============================================================================
# Colores: Azul verdoso (positivo), Blanco (cero), Rojo coral (negativo)
col_palette <- colorRampPalette(c("#E74C3C", "#FFFFFF", "#00A88F"))(200)

corrplot(
  mat_cor,
  method = "circle",
  type = "upper",
  order = "hclust",
  p.mat = p_matrix,
  sig.level = 0.05,
  insig = "blank", # Deja en blanco las correlaciones no significativas
  col = col_palette,
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.9,
  title = "Matriz de Correlación: Ionómica vs Fisiología (p < 0.05)",
  mar = c(0, 0, 2, 0)
)

# ============================================================================
# 4. GRÁFICO 2: HEATMAP CON CLUSTERING JERÁRQUICO COMPLETO (pheatmap)
# ============================================================================
# Normalizar datos por columna (Z-score) para comparar rasgos en diferentes escalas
datos_escalados <- scale(datos_multiomica)

# Anotación de variables por categoría biológica
anotacion_vars <- data.frame(
  Categoria = c("Ionomica", "Ionomica", "Ionomica", "Fisiologia", "Fisiologia", "Bioquimica", "Bioquimica", "Agronomia"),
  row.names = colnames(datos_multiomica)
)

colores_anotacion <- list(
  Categoria = c(
    "Ionomica" = "#3498DB",
    "Fisiologia" = "#2ECC71",
    "Bioquimica" = "#E67E22",
    "Agronomia" = "#9B59B6"
  )
)

pheatmap(
  t(datos_escalados),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#2C3E50", "#ECF0F1", "#E74C3C"))(100),
  annotation_row = anotacion_vars,
  annotation_colors = colores_anotacion,
  main = "Mapa de Calor Jerárquico: Perfiles Multi-Atributo en Trigo",
  fontsize = 10,
  fontsize_row = 10,
  show_colnames = FALSE,
  border_color = NA
)

cat("\n--- Script 04 Redes y Heatmap finalizado correctamente ---\n")
