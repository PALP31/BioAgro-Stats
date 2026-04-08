# ============================================================================
# 01_dbca_basico.R
# Diseño de Bloques Completos al Azar (DBCA/RCBD)
# Contexto: Dosis de bioestimulantes en huertos de manzanas
# ============================================================================
# Paquete principal: agricolae
# Diseño: DBCA clásico con 4 tratamientos y 5 bloques
# Variable respuesta: Rendimiento de fruta (kg/árbol)
# ============================================================================

# Cargar paquetes necesarios
library(agricolae)   # Paquete principal para diseño experimental
library(tidyverse)   # Manipulación de datos y gráficos
library(performance) # Diagnóstico de modelos
library(car)         # Pruebas de homocedasticidad

set.seed(789)

# ============================================================================
# 1. GENERACIÓN DEL DISEÑO EXPERIMENTAL CON AGRICOLAE
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 1: GENERACIÓN DEL DISEÑO DBCA\n")
cat(rep("=", 80), "\n\n", sep = "")

# Definir tratamientos: 4 dosis de bioestimulante
tratamientos <- c("Testigo", "Dosis_baja", "Dosis_media", "Dosis_alta")
n_bloques <- 5

# Usar agricolae para generar el diseño DBCA
# design.rcbd() crea la asignación aleatoria de tratamientos a bloques
design_dbca <- design.rcbd(
  trt = tratamientos,  # Tratamientos
  r = n_bloques,       # Número de repeticiones/bloques
  serie = 2,           # Serie para numeración de parcelas
  seed = 789          # Semilla para reproducibilidad
)

cat("\n=== Esquema del Diseño DBCA ===\n")
cat("Tratamientos:", length(tratamientos), "\n")
cat("Bloques:", n_bloques, "\n")
cat("Total de parcelas:", length(tratamientos) * n_bloques, "\n")

cat("\n=== Libro de campo (fieldbook) ===\n")
print(design_dbca$book)

cat("\n=== Matriz de diseño (tratamientos por bloque) ===\n")
# Crear matriz para visualizar el diseño
matriz_diseno <- matrix(design_dbca$book$plots,
                        nrow = n_bloques,
                        ncol = length(tratamientos),
                        byrow = TRUE)
colnames(matriz_diseno) <- tratamientos
rownames(matriz_diseno) <- paste0("Bloque_", 1:n_bloques)
print(matriz_diseno)

# ============================================================================
# 2. SIMULACIÓN DE DATOS DEL ENSAYO
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 2: SIMULACIÓN DE DATOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Extraer información del diseño para simular datos
datos <- design_dbca$book %>%
  mutate(
    bloque = factor(rep(paste0("Bloque_", 1:n_bloques), each = length(tratamientos))),
    tratamiento = factor(tratamiento, levels = 1:4, labels = tratamientos)
  )

# Efectos verdaderos de los tratamientos (desviación del testigo)
efectos_tratamiento <- c(
  "Testigo" = 0,        # Referencia
  "Dosis_baja" = 3.5,   # Mejora moderada
  "Dosis_media" = 6.2,  # Mejora óptima
  "Dosis_alta" = 4.8    # Mejora menor (ley de rendimientos decrecientes)
)

# Efectos de bloque (variabilidad del terreno, exposición solar, etc.)
efectos_bloque <- c(
  "Bloque_1" = -1.5,   # Zona baja, más húmeda
  "Bloque_2" = 0.8,    # Zona media
  "Bloque_3" = 1.2,    # Zona alta, mejor suelo
  "Bloque_4" = -0.5,   # Zona media-baja
  "Bloque_5" = 0.0     # Zona de referencia
)

# Simular rendimiento con efectos fijos + error aleatorio
datos <- datos %>%
  mutate(
    # Media base: 45 kg/árbol
    efecto_trat = efectos_tratamiento[as.character(tratamiento)],
    efecto_bloq = efectos_bloque[as.character(bloque)],
    error = rnorm(n(), 0, 2.0),  # Error experimental
    rendimiento = 45 + efecto_trat + efecto_bloq + error
  ) %>%
  # Limpiar columnas auxiliares
  select(-efecto_trat, -efecto_bloq, -error, -plots)

cat("\n=== Estructura de los datos ===\n")
str(datos)

cat("\n=== Resumen descriptivo por tratamiento ===\n")
datos %>%
  group_by(tratamiento) %>%
  summarise(
    n = n(),
    media = mean(rendimiento),
    sd = sd(rendimiento),
    min = min(rendimiento),
    max = max(rendimiento),
    .groups = "drop"
  ) %>%
  print()

cat("\n=== Resumen por bloque ===\n")
datos %>%
  group_by(bloque) %>%
  summarise(
    n = n(),
    media = mean(rendimiento),
    sd = sd(rendimiento),
    .groups = "drop"
  ) %>%
  print()

# ============================================================================
# 3. ANÁLISIS DE VARIANZA (ANOVA)
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 3: ANÁLISIS DE VARIANZA\n")
cat(rep("=", 80), "\n\n", sep = "")

# Ajustar modelo lineal para DBCA
# NOTA: En DBCA, bloque es efecto fijo para el análisis (aunque representa
# variabilidad aleatoria en la población de bloques)
modelo <- aov(rendimiento ~ tratamiento + bloque, data = datos)

cat("\n=== Tabla ANOVA ===\n")
# Usar summary() para obtener la tabla ANOVA completa
anova_table <- summary(modelo)
print(anova_table)

# Extraer información clave
cat("\n=== Interpretación de Resultados ===\n")
tratamiento_p <- anova_table[[1]]$`Pr(>F)`[1]
bloque_p <- anova_table[[1]]$`Pr(>F)`[2]

cat("Efecto del Tratamiento:\n")
cat("  F =", round(anova_table[[1]]$`F value`[1], 2),
    ", p =", round(tratamiento_p, 5), "\n")
if (tratamiento_p < 0.05) {
  cat("  → DIFERENCIAS SIGNIFICATIVAS entre tratamientos (α = 0.05)\n")
} else {
  cat("  → NO hay diferencias significativas entre tratamientos\n")
}

cat("\nEfecto del Bloque:\n")
cat("  F =", round(anova_table[[1]]$`F value`[2], 2),
    ", p =", round(bloque_p, 5), "\n")
if (bloque_p < 0.05) {
  cat("  → El bloqueo fue EFECTIVO para controlar variabilidad\n")
} else {
  cat("  → El bloqueo NO fue efectivo (quizás CRD hubiera sido suficiente)\n")
}

# ============================================================================
# 4. EVALUACIÓN DE SUPUESTOS DEL ANOVA
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 4: EVALUACIÓN DE SUPUESTOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Extraer residuos
residuos <- residuals(modelo)

# ---- 4.1 Normalidad de Residuos ----
cat("\n4.1 Normalidad de Residuos (Shapiro-Wilk)\n")
cat("-" , rep("-", 50), "\n", sep = "")

shapiro_test <- shapiro.test(residuos)
cat("Estadístico W =", round(shapiro_test$statistic, 4), "\n")
cat("p-value =", round(shapiro_test$p.value, 4), "\n")

if (shapiro_test$p.value > 0.05) {
  cat("→ Los residuos SIGUEN una distribución normal ✓\n")
} else {
  cat("→ Los residuos se DESVÍAN de la normalidad ⚠\n")
}

# Gráfico Q-Q
qq_plot <- ggplot(data.frame(residuos = residuos), aes(sample = residuos)) +
  stat_qq(color = "#2E86AB", size = 3, alpha = 0.7) +
  stat_qq_line(color = "#E74C3C", linewidth = 1.2) +
  labs(title = "Q-Q Plot de Residuos",
       subtitle = "Evaluación visual de normalidad",
       x = "Cuantiles teóricos", y = "Cuantiles observados") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

print(qq_plot)
ggsave("02_Diseno_Experimental/01_dbca_qq_plot.png", qq_plot,
       width = 6, height = 5, dpi = 300)

# ---- 4.2 Homocedasticidad (Levene) ----
cat("\n4.2 Homocedasticidad (Prueba de Levene)\n")
cat("-" , rep("-", 50), "\n", sep = "")

levene_test <- leveneTest(rendimiento ~ tratamiento, data = datos)
print(levene_test)

if (levene_test$`Pr(>F)`[1] > 0.05) {
  cat("→ Homocedasticidad CONFIRMADA ✓\n")
} else {
  cat("→ Heterocedasticidad DETECTADA ⚠\n")
}

# ---- 4.3 Independencia de Residuos ----
cat("\n4.3 Independencia de Residuos (Durbin-Watson)\n")
cat("-" , rep("-", 50), "\n", sep = "")

dw_test <- durbinWatsonTest(modelo)
cat("Estadístico Durbin-Watson =", round(dw_test$dw, 4), "\n")
cat("p-value =", round(dw_test$p, 4), "\n")

if (dw_test$p > 0.05) {
  cat("→ No hay autocorrelación significativa ✓\n")
} else {
  cat("→ Posible autocorrelación de residuos ⚠\n")
}

# ---- 4.4 Diagnóstico Gráfico ----
cat("\n4.4 Diagnóstico Gráfico (check_model)\n")
cat("-" , rep("-", 50), "\n", sep = "")

check_plot <- check_model(modelo, panel = TRUE)
print(check_plot)
ggsave("02_Diseno_Experimental/01_dbca_diagnostico.png",
       plot = check_plot, width = 10, height = 8, dpi = 300)

# ============================================================================
# 5. PRUEBAS DE COMPARACIONES MÚLTIPLES (POST-HOC)
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 5: COMPARACIONES MÚLTIPLES (POST-HOC)\n")
cat(rep("=", 80), "\n\n", sep = "")

# agricolae proporciona varias pruebas de comparaciones múltiples
# Extraer Media Cuadrática del Error y grados de libertad del ANOVA
ame <- anova_table[[1]]$`Mean Sq`[3]  # Media Cuadrática del Error
gl_error <- anova_table[[1]]$`Df`[3]  # Grados de libertad del error

cat("Media Cuadrática del Error (CME):", round(ame, 3), "\n")
cat("Grados de libertad del error:", gl_error, "\n\n")

# ---- 5.1 Prueba de Tukey (HSD) ----
cat("\n5.1 Prueba de Tukey (Honestly Significant Difference)\n")
cat("-" , rep("-", 50), "\n", sep = "")

tukey_result <- HSD.test(
  model = modelo,
  trt = "tratamiento",
  alpha = 0.05,
  group = TRUE,
  console = FALSE
)

cat("Diferencia Mínima Significativa (HSD):", round(tukey_result$HSD, 3), "\n\n")
cat("Medias por tratamiento con letras de agrupamiento:\n")
print(tukey_result$groups)

# Interpretación de letras:
cat("\nInterpretación: Tratamientos con la MISMA letra NO son\n")
cat("significativamente diferentes (p > 0.05)\n")

# ---- 5.2 Prueba de Fisher (LSD) ----
cat("\n5.2 Prueba de Fisher (Least Significant Difference)\n")
cat("-" , rep("-", 50), "\n", sep = "")

lsd_result <- LSD.test(
  model = modelo,
  trt = "tratamiento",
  alpha = 0.05,
  group = TRUE,
  console = FALSE
)

cat("Diferencia Mínima Significativa (LSD):", round(lsd_result$HSD, 3), "\n")
cat("(LSD es menos conservador que Tukey)\n\n")
print(lsd_result$groups)

# ---- 5.3 Prueba de Duncan ----
cat("\n5.3 Prueba de Duncan\n")
cat("-" , rep("-", 50), "\n", sep = "")

duncan_result <- duncan.test(
  model = modelo,
  trt = "tratamiento",
  alpha = 0.05,
  group = TRUE,
  console = FALSE
)

print(duncan_result$groups)

# ---- 5.4 Comparación con el Testigo (Dunnett) ----
cat("\n5.4 Prueba de Dunnett (comparación vs Testigo)\n")
cat("-" , rep("-", 50), "\n", sep = "")

dunnett_result <- dunnett.test(
  model = modelo,
  trt = "tratamiento",
  control = "Testigo",  # Especificar cuál es el tratamiento control
  alpha = 0.05
)

cat("Comparaciones contra el Testigo:\n")
print(dunnett_result$groups)

cat("\nInterpretación: asterisco (*) indica diferencia significativa\n")
cat("con respecto al Testigo (p < 0.05)\n")

# ============================================================================
# 6. ANÁLISIS DE CONTRASTES ORTOGONALES
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 6: CONTRASTES ORTOGONALES\n")
cat(rep("=", 80), "\n\n", sep = "")

# Los contrastes ortogonales permiten probar hipótesis específicas
# Definir matriz de contrastes:
# Contraste 1: Testigo vs Todos los bioestimulantes
# Contraste 2: Dosis baja vs Dosis media+alta
# Contraste 3: Dosis media vs Dosis alta

contrastes <- matrix(c(
  # Testigo  Baja  Media  Alta
     -3,       1,     1,    1,   # C1: Control vs tratados
      0,      -2,     1,    1,   # C2: Baja vs Media+Alta
      0,       0,    -1,    1    # C3: Media vs Alta
), nrow = 4, byrow = FALSE)

colnames(contrastes) <- c("Control_vs_Tratados", "Baja_vs_MediaAlta", "Media_vs_Alta")
rownames(contrastes) <- tratamientos

cat("Matriz de contrastes ortogonales:\n")
print(contrastes)

# Verificar ortogonalidad (suma de productos cruzados = 0)
cat("\nVerificación de ortogonalidad:\n")
cat("C1 · C2 =", sum(contrastes[,1] * contrastes[,2]), "\n")
cat("C1 · C3 =", sum(contrastes[,1] * contrastes[,3]), "\n")
cat("C2 · C3 =", sum(contrastes[,2] * contrastes[,3]), "\n")

# Aplicar contrastes
contrastes_result <- contras.model(
  model = modelo,
  factores = list(tratamiento = contrastes),
  alpha = 0.05
)

cat("\nResultados de contrastes:\n")
print(contrastes_result)

# ============================================================================
# 7. VISUALIZACIÓN DE RESULTADOS
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 7: VISUALIZACIÓN DE RESULTADOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Añadir letras de Tukey a los datos
letras <- tukey_result$groups %>%
  rownames_to_column("tratamiento") %>%
  as_tibble()

# Gráfico de barras con letras
plot_barras <- ggplot(datos, aes(x = tratamiento, y = rendimiento, fill = tratamiento)) +
  stat_summary(fun = "mean", geom = "bar", alpha = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3) +
  geom_text(data = letras,
            aes(x = tratamiento, y = max(datos$rendimiento) + 1.5, label = groups),
            size = 6, fontface = "bold") +
  labs(title = "Efecto de Bioestimulantes en Rendimiento de Manzanas",
       subtitle = "DBCA con 5 bloques. Letras diferentes = diferencia significativa (Tukey, p < 0.05)",
       x = "Tratamiento",
       y = "Rendimiento (kg/árbol)") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_fill_brewer(palette = "Greens", name = "Tratamiento")

print(plot_barras)
ggsave("02_Diseno_Experimental/01_dbca_barras.png",
       plot = plot_barras, width = 8, height = 6, dpi = 300)

# Gráfico de cajas por bloque
plot_cajas <- ggplot(datos, aes(x = bloque, y = rendimiento, fill = bloque)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Rendimiento por Bloque",
       subtitle = "Evaluación del efecto de bloqueo",
       x = "Bloque", y = "Rendimiento (kg/árbol)") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  scale_fill_brewer(palette = "Blues")

print(plot_cajas)
ggsave("02_Diseno_Experimental/01_dbca_bloques.png",
       plot = plot_cajas, width = 7, height = 5, dpi = 300)

# ============================================================================
# 8. RESUMEN FINAL Y RECOMENDACIONES
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("RESUMEN DEL ANÁLISIS DBCA\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("DISEÑO EXPERIMENTAL:\n")
cat("  Tipo: Bloques Completos al Azar (DBCA/RCBD)\n")
cat("  Tratamientos:", length(tratamientos), "\n")
cat("  Bloques:", n_bloques, "\n")
cat("  Total de observaciones:", nrow(datos), "\n\n")

cat("RESULTADOS DEL ANOVA:\n")
cat("  Efecto Tratamiento: F =", round(anova_table[[1]]$`F value`[1], 2),
    ", p =", format(tratamiento_p, scientific = FALSE), "\n")
cat("  Efecto Bloque: F =", round(anova_table[[1]]$`F value`[2], 2),
    ", p =", format(bloque_p, scientific = FALSE), "\n\n")

cat("SUPUESTOS:\n")
cat("  Normalidad (Shapiro-Wilk): p =", round(shapiro_test$p.value, 4),
    ifelse(shapiro_test$p.value > 0.05, "✓", "⚠"), "\n")
cat("  Homocedasticidad (Levene): p =", round(levene_test$`Pr(>F)`[1], 4),
    ifelse(levene_test$`Pr(>F)`[1] > 0.05, "✓", "⚠"), "\n")
cat("  Independencia (Durbin-Watson): p =", round(dw_test$p, 4),
    ifelse(dw_test$p > 0.05, "✓", "⚠"), "\n\n")

cat("MEJOR TRATAMIENTO SEGÚN TUKEY:\n")
mejor_trat <- letras %>%
  filter(groups == min(groups)) %>%
  pull(tratamiento) %>%
  first()
cat("  →", mejor_trat, "con mayor rendimiento\n\n")

cat("EFICIENCIA DEL BLOQUEO:\n")
# Calcular eficiencia relativa del bloqueo
vc_bloque <- anova_table[[1]]$`Mean Sq`[2]
vc_error <- ame
eficiencia <- (1 - vc_error/vc_bloque) * 100
cat("  Eficiencia relativa:", round(eficiencia, 1), "%\n")
if (eficiencia > 0) {
  cat("  → El bloqueo MEJORÓ la precisión del experimento\n")
} else {
  cat("  → El bloqueo NO mejoró la precisión (considerar DCA en futuros ensayos)\n")
}

# ============================================================================
# 9. EJERCICIO AVANZADO: ESTRÉS POR SALINIDAD Y PGPR (Heterocedasticidad)
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 9: EJERCICIO AVANZADO - ESTRÉS SALINO Y PGPR\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("ESCENARIO:\n")
cat("Se evalúan 5 cepas de PGPR (Rhizobium, Pseudomonas, Bacillus, Azospirillum,
Enterobacter) en un suelo con alta salinidad (CE > 4 dS/m). 
El objetivo es mitigar el estrés osmótico en plántulas de tomate.
DIFICULTAD: Los datos presentan varianza no constante (heterocedasticidad).\n\n")

# 1. Simulación de datos con Heterocedasticidad
set.seed(123)
cepas <- c("Control", "Rhizobium", "Pseudomonas", "Bacillus", "Azospirillum")
bloques <- paste0("B", 1:4)

datos_salino <- expand.grid(Cepa = cepas, Bloque = bloques) %>%
  mutate(
    # Efectos de cepa (mejora de biomasa bajo salinidad)
    efecto = case_when(
      Cepa == "Control" ~ 0,
      Cepa == "Rhizobium" ~ 2.5,
      Cepa == "Pseudomonas" ~ 5.8,
      Cepa == "Bacillus" ~ 4.2,
      Cepa == "Azospirillum" ~ 3.1
    ),
    # Heterocedasticidad: mayor biomasa -> mayor varianza
    sd_error = case_when(
      Cepa == "Control" ~ 0.5,
      Cepa == "Pseudomonas" ~ 2.5,  # Mucha más varianza
      TRUE ~ 1.2
    ),
    Biomasa = 15 + efecto + rnorm(n(), 0, sd_error)
  )

cat("=== Diagnóstico de Heterocedasticidad ===\n")
mod_basico <- aov(Biomasa ~ Cepa + Bloque, data = datos_salino)
levene_sal <- car::leveneTest(Biomasa ~ Cepa, data = datos_salino)
print(levene_sal)

if(levene_sal$`Pr(>F)`[1] < 0.05) {
  cat("\n⚠ ALERTA: Heterocedasticidad detectada. El ANOVA clásico no es óptimo.\n")
}

# 2. Solución Avanzada: Mínimos Cuadrados Generalizados (GLS)
library(nlme)
cat("\n=== Ajuste con GLS (Modelando la estructura de varianza) ===\n")
# Permitimos que cada cepa tenga su propia varianza residual
mod_gls <- gls(Biomasa ~ Cepa + Bloque, 
               data = datos_salino,
               weights = varIdent(form = ~ 1 | Cepa))

anova_gls <- anova(mod_gls)
print(anova_gls)

# 3. Comparaciones Múltiples Robustas
cat("\n=== Comparaciones de Medias con errores estándar ajustados ===\n")
emm_sal <- emmeans::emmeans(mod_gls, ~ Cepa)
pairs_sal <- pairs(emm_sal, adjust = "tukey")
print(pairs_sal)

cat("\nRECOMENDACIÓN FINAL:\n")
cat("Cuando el estrés abiótico induce variabilidad diferencial entre tratamientos,
utilice modelos GLS o transformaciones (log, raíz cuadrada) para estabilizar
la varianza y evitar falsos positivos o pérdida de potencia.\n")

cat("\n", rep("=", 80), "\n", sep = "")
cat("FIN DEL ANÁLISIS\n")
cat(rep("=", 80), "\n", sep = "")
