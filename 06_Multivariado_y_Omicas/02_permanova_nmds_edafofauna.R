# ============================================================================
# 02_permanova_nmds_edafofauna.R
# Análisis de Comunidades: PERMANOVA y NMDS
# Contexto: Comunidades de Edafofauna tratadas con Bioestimulantes
# ============================================================================
# Este script utiliza 'vegan' para evaluar diferencias en la composición de
# comunidades, validando la homogeneidad de dispersiones antes de PERMANOVA.
# ============================================================================

library(vegan)
library(tidyverse)
library(ggsci) # Paletas de colores profesionales
library(viridis)

set.seed(456)

# 1. SIMULACIÓN DE DATOS (20 Muestras x 30 Especies)
# 2 Tratamientos: Control vs Bioestimulante (n=10)
n_muestras <- 20
n_especies <- 30
tratamientos <- factor(rep(c("Control", "Bioestimulante"), each = 10))

# Simulación con Binomial Negativa (Sobredispersión natural en conteos)
# El bioestimulante altera la abundancia de ciertas especies (indices 1:10)
datos_comunidad <- matrix(nrow = n_muestras, ncol = n_especies)
colnames(datos_comunidad) <- paste0("Spp_", 1:n_especies)

# Control: Baja abundancia general
datos_comunidad[1:10, ] <- rnbinom(10 * n_especies, mu = 2, size = 1.2)
# Bioestimulante: Alta abundancia en especies clave 1:10
datos_comunidad[11:20, 1:10] <- rnbinom(10 * 10, mu = 15, size = 1.2)
# Especies comunes
datos_comunidad[11:20, 11:30] <- rnbinom(10 * 20, mu = 3, size = 1.2)

# 2. ANÁLISIS DE DISIMILITUD (Bray-Curtis)
# Matriz de distancia para datos de conteo
dist_matrix <- vegdist(datos_comunidad, method = "bray")

# 3. CRÍTICO: ANÁLISIS DE HOMOGENEIDAD DE DISPERSIONES (PERMDISP)
# Un PERMANOVA significativo puede ser causado por diferencias en dispersión,
# no solo en posición (composición). Debemos descartar esto primero.

cat("\n--- Análisis de Homogeneidad de Dispersiones (PERMDISP) ---\n")
disp_calc <- betadisper(dist_matrix, group = tratamientos)
perm_disp <- permutest(disp_calc, permutations = 999)

print(perm_disp)

if(perm_disp$tab$`Pr(>F)`[1] > 0.05) {
  cat("Homogeneidad de dispersiones CONFIRMADA (p > 0.05). \n")
} else {
  cat("⚠ ADVERTENCIA: Heterogeneidad de dispersiones detectada. Proceda con precaución.\n")
}

# 4. PERMANOVA (adonis2)
# ¿Difieren significativamente las comunidades entre tratamientos?
cat("\n--- Análisis de Varianza Multivariado (PERMANOVA) ---\n")
permanova_res <- adonis2(datos_comunidad ~ tratamientos, method = "bray", permutations = 999)
print(permanova_res)

# 5. NMDS (Non-metric Multidimensional Scaling)
# Reducción de dimensiones para visualización
cat("\n--- Calculando Escalamiento Multidimensional No Métrico (NMDS) ---\n")
nmds_calc <- metaMDS(datos_comunidad, distance = "bray", k = 2, trymax = 100)

# Extraer Scores para ggplot2
nmds_scores <- as.data.frame(scores(nmds_calc)$sites) %>%
  mutate(Tratamiento = tratamientos)

# 6. VISUALIZACIÓN NMDS DE ALTA CALIDAD (ggplot2)
cat("\n--- Generando NMDS Biplot ---\n")

plot_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Tratamiento)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(fill = Tratamiento), geom = "polygon", alpha = 0.2, linetype = "dashed") +
  annotate("text", x = max(nmds_scores$NMDS1)*0.7, y = min(nmds_scores$NMDS2)*0.9, 
           label = paste("Stress =", round(nmds_calc$stress, 3)), fontface = "italic") +
  scale_color_npg() +
  scale_fill_npg() +
  labs(title = "Composición de Edafofauna (NMDS Bray-Curtis)",
       subtitle = paste("PERMANOVA p-value:", permanova_res$`Pr(>F)`[1]),
       x = "NMDS 1", y = "NMDS 2") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

print(plot_nmds)

cat("\n--- Script Finalizado (PERMANOVA/NMDS) ---\n")
