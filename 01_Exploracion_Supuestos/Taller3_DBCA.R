
# ==============================================================================
# TALLER 3: DISEÑO DE BLOQUES COMPLETAMENTE AL AZAR (DBCA)
# Curso: Aplicaciones Estadísticas en Agronomía
# ==============================================================================
#
# OBJETIVOS DE APRENDIZAJE:
#   1. Entender la lógica y ventajas del DBCA frente al DCA.
#   2. Conocer el modelo estadístico y los grados de libertad.
#   3. Verificar los supuestos del modelo (normalidad, homocedasticidad).
#   4. Interpretar la tabla ANOVA con bloques.
#   5. Aplicar comparaciones múltiples de Tukey.
#   6. Generar gráficos profesionales con ggplot2.
#
# ==============================================================================

# ==============================================================================
# SECCIÓN 0: INSTALACIÓN Y CARGA DE LIBRERÍAS
# ==============================================================================

# Instala los paquetes si no los tienes (ejecuta una sola vez)
paquetes <- c("ggplot2", "agricolae", "car", "dplyr", "tidyr", "ggpubr",
              "performance", "see", "multcompView", "emmeans", "scales")

for (p in paquetes) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

# Carga de librerías
library(ggplot2)       # Gráficos elegantes
library(agricolae)     # Pruebas post-hoc (Tukey, Duncan, LSD)
library(car)           # Prueba de Levene (homocedasticidad)
library(dplyr)         # Manipulación de datos
library(tidyr)         # Organización de datos
library(ggpubr)        # Paneles de gráficos
library(performance)   # Diagnóstico del modelo
library(see)           # Visualización de supuestos
library(multcompView)  # Letras de significancia
library(emmeans)       # Medias marginales estimadas
library(scales)        # Escalas para ggplot2

cat("✅ Librerías cargadas correctamente.\n")


# ==============================================================================
# SECCIÓN 1: FUNDAMENTOS TEÓRICOS DEL DBCA
# ==============================================================================
#
# ── ¿QUÉ ES EL DBCA? ──────────────────────────────────────────────────────────
#
# El Diseño de Bloques Completamente al Azar (DBCA) es una extensión del DCA
# que se utiliza cuando las unidades experimentales NO son homogéneas.
#
# PRINCIPIO BÁSICO: Si existe una fuente de variación CONOCIDA que puede afectar
# la variable de respuesta (ej. fertilidad del suelo, pendiente del campo,
# operador en laboratorio), la "controlamos" agrupando las unidades en BLOQUES.
# Dentro de cada bloque, los tratamientos se asignan de forma COMPLETAMENTE
# al azar.
#
# ── COMPARACIÓN DCA vs DBCA ───────────────────────────────────────────────────
#
# ┌──────────────────┬────────────────────────────┬────────────────────────────┐
# │ Criterio         │ DCA                        │ DBCA                       │
# ├──────────────────┼────────────────────────────┼────────────────────────────┤
# │ Unid. exp.       │ Homogéneas                 │ Heterogéneas               │
# │ Fuente variación │ Solo tratamientos          │ Tratamientos + Bloques     │
# │ Modelo           │ Y = µ + τ_i + ε_ij         │ Y = µ + τ_i + β_j + ε_ij  │
# │ Grados libertad  │ Error = N - k              │ Error = (k-1)(r-1)         │
# │ Sensibilidad     │ Menor (si hay heterog.)    │ Mayor (aísla efecto bloque)│
# │ Uso típico       │ Invernadero uniforme       │ Campo con gradiente        │
# └──────────────────┴────────────────────────────┴────────────────────────────┘
#
# ── MODELO ESTADÍSTICO ────────────────────────────────────────────────────────
#
#   Y_ij = µ + τ_i + β_j + ε_ij
#
#   Donde:
#   Y_ij  = observación del tratamiento i en el bloque j
#   µ     = media general (constante)
#   τ_i   = efecto del tratamiento i   (i = 1, 2, ..., k)
#   β_j   = efecto del bloque j        (j = 1, 2, ..., r)
#   ε_ij  = error aleatorio            ε ~ N(0, σ²)
#
# ── GRADOS DE LIBERTAD ────────────────────────────────────────────────────────
#
#   Con k = tratamientos y r = bloques (réplicas), N = k × r observaciones:
#
#   ┌──────────────────────┬──────────────────────────────────────────┐
#   │ Fuente de Variación  │ Grados de Libertad (GL)                  │
#   ├──────────────────────┼──────────────────────────────────────────┤
#   │ Tratamientos         │ k - 1                                    │
#   │ Bloques              │ r - 1                                    │
#   │ Error Experimental   │ (k-1)(r-1) = N - k - r + 1              │
#   │ Total                │ N - 1 = kr - 1                           │
#   └──────────────────────┴──────────────────────────────────────────┘
#
#   COMPARACIÓN DE GL DEL ERROR:
#   - En el DCA:  GL_Error = N - k  = kr - k = k(r-1)
#   - En el DBCA: GL_Error = (k-1)(r-1)
#
#   EJEMPLO PRÁCTICO (k=5 tratamientos, r=4 bloques, N=20):
#   - DCA:  GL_Error = 20 - 5 = 15
#   - DBCA: GL_Error = (5-1)(4-1) = 4 × 3 = 12
#
#   ⚠️ NOTA IMPORTANTE: En el DBCA el GL del Error es MENOR que en el DCA.
#   Esto puede reducir la potencia estadística. Sin embargo, si los bloques
#   explican variabilidad real, la ganancia en precisión (menor SC_Error)
#   compensa con creces la pérdida de GL. Si los bloques NO explican variación,
#   el DBCA es MENOS poderoso que el DCA.
#
# ── TABLA ANOVA DEL DBCA ─────────────────────────────────────────────────────
#
# ┌──────────────────┬──────────┬──────────┬──────────┬──────────┬───────────┐
# │ FV               │ GL       │ SC       │ CM       │ F calc   │ p-valor   │
# ├──────────────────┼──────────┼──────────┼──────────┼──────────┼───────────┤
# │ Tratamientos     │ k-1      │ SC_Trat  │ CM_Trat  │ CM_T/CME │ P(F>Fc)   │
# │ Bloques          │ r-1      │ SC_Blq   │ CM_Blq   │ CM_B/CME │ P(F>Fc)   │
# │ Error            │(k-1)(r-1)│ SC_Error │ CM_Error │ —        │ —         │
# │ Total            │ N-1      │ SC_Total │ —        │ —        │ —         │
# └──────────────────┴──────────┴──────────┴──────────┴──────────┴───────────┘
#
# ==============================================================================


# ==============================================================================
# SECCIÓN 2: CARGA Y EXPLORACIÓN DE DATOS
# ==============================================================================
#
# CONTEXTO DEL EXPERIMENTO:
# ─────────────────────────
# Un agrónomo desea evaluar el efecto de 5 tratamientos fungicidas (incluyendo
# un control sin aplicación) sobre el rendimiento de frijol (Phaseolus vulgaris)
# en kg/ha.
#
# PROBLEMA: El campo experimental tiene un GRADIENTE DE FERTILIDAD del suelo
# de Norte a Sur. Si se usa un DCA, este gradiente puede confundirse con el
# efecto de los fungicidas.
#
# SOLUCIÓN: Se divide el campo en 4 BLOQUES perpendiculares al gradiente.
# Dentro de cada bloque, los 5 fungicidas se asignan al azar.
#
# DISEÑO: k=5 tratamientos × r=4 bloques = N=20 unidades experimentales.

# ── Carga de datos ────────────────────────────────────────────────────────────

datos <- read.csv("datos_taller3_dbca.csv")

# Verificar estructura
cat("\n📊 ESTRUCTURA DEL DATASET:\n")
str(datos)

cat("\n📋 PRIMERAS FILAS:\n")
print(head(datos, 10))

# Convertir a factores (¡muy importante para el ANOVA!)
datos$Bloque    <- factor(datos$Bloque,
                          levels = c("Bloque_1_Norte", "Bloque_2_CentroN",
                                     "Bloque_3_CentroS", "Bloque_4_Sur"),
                          labels = c("B1: Norte", "B2: Centro-N",
                                     "B3: Centro-S", "B4: Sur"))
datos$Fungicida <- factor(datos$Fungicida,
                          levels = c("Control", "Fungicida_A", "Fungicida_B",
                                     "Fungicida_C", "Fungicida_D"))

cat("\n✅ Variables convertidas a factores.\n")
cat("\n📊 NÚMERO DE OBSERVACIONES POR TRATAMIENTO Y BLOQUE:\n")
print(table(datos$Fungicida, datos$Bloque))


# ── Estadísticas descriptivas ─────────────────────────────────────────────────

cat("\n📊 ESTADÍSTICAS POR TRATAMIENTO (Fungicida):\n")
resumen_trt <- datos %>%
  group_by(Fungicida) %>%
  summarise(
    n       = n(),
    Media   = round(mean(Rendimiento_kgha), 1),
    DE      = round(sd(Rendimiento_kgha), 1),
    Min     = round(min(Rendimiento_kgha), 1),
    Max     = round(max(Rendimiento_kgha), 1),
    CV_pct  = round(sd(Rendimiento_kgha) / mean(Rendimiento_kgha) * 100, 1),
    .groups = "drop"
  )
print(as.data.frame(resumen_trt))

cat("\n📊 ESTADÍSTICAS POR BLOQUE:\n")
resumen_blq <- datos %>%
  group_by(Bloque) %>%
  summarise(
    n     = n(),
    Media = round(mean(Rendimiento_kgha), 1),
    DE    = round(sd(Rendimiento_kgha), 1),
    .groups = "drop"
  )
print(as.data.frame(resumen_blq))


# ==============================================================================
# SECCIÓN 3: EXPLORACIÓN VISUAL (HISTOGRAMA, BOXPLOT, HEATMAP)
# ==============================================================================

# ── Gráfico 1: Histograma general de la distribución ─────────────────────────
g1_hist <- ggplot(datos, aes(x = Rendimiento_kgha)) +
  geom_histogram(aes(y = after_stat(density)), bins = 10,
                 fill = "#2E86AB", color = "white", alpha = 0.85) +
  geom_density(color = "#E84855", linewidth = 1.2, linetype = "dashed") +
  labs(
    title    = "Distribución del Rendimiento de Frijol",
    subtitle = "Todos los datos del experimento DBCA",
    x        = "Rendimiento (kg/ha)",
    y        = "Densidad"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  )
print(g1_hist)

# ── Gráfico 2: Boxplot por tratamiento con puntos individuales ────────────────
paleta_trt <- c("#6C757D", "#2E86AB", "#A8C686", "#F4A261", "#E84855")

g2_box <- ggplot(datos, aes(x = Fungicida, y = Rendimiento_kgha, fill = Fungicida)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.8, color = "black") +
  scale_fill_manual(values = paleta_trt) +
  labs(
    title    = "Rendimiento de Frijol por Fungicida",
    subtitle = "Diseño de Bloques Completamente al Azar (n = 4 bloques por tratamiento)",
    x        = "Tratamiento Fungicida",
    y        = "Rendimiento (kg/ha)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x   = element_text(angle = 15, hjust = 1)
  )
print(g2_box)

# ── Gráfico 3: Boxplot por bloque (para visualizar el efecto del bloque) ──────
paleta_blq <- c("#457B9D", "#1D3557", "#A8DADC", "#E63946")

g3_blq <- ggplot(datos, aes(x = Bloque, y = Rendimiento_kgha, fill = Bloque)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.8, color = "black") +
  scale_fill_manual(values = paleta_blq) +
  labs(
    title    = "Rendimiento de Frijol por Bloque",
    subtitle = "¿Existe gradiente de fertilidad entre bloques?",
    x        = "Bloque (Gradiente Norte-Sur)",
    y        = "Rendimiento (kg/ha)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none"
  )
print(g3_blq)

# ── Gráfico 4: Heatmap de datos (matriz Bloque × Fungicida) ──────────────────
# Útil para visualizar si hay alguna tendencia o celda inusual
g4_heat <- ggplot(datos, aes(x = Fungicida, y = Bloque, fill = Rendimiento_kgha)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Rendimiento_kgha), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(
    low      = "#457B9D",
    mid      = "#F1FAEE",
    high     = "#E63946",
    midpoint = mean(datos$Rendimiento_kgha),
    name     = "kg/ha"
  ) +
  labs(
    title    = "Mapa de Rendimiento: Fungicida × Bloque",
    subtitle = "Cada celda = 1 unidad experimental",
    x        = "Fungicida",
    y        = "Bloque"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.text.x   = element_text(angle = 15, hjust = 1)
  )
print(g4_heat)

# ── Gráfico 5: Perfil de medias (interacción visual Bloque × Fungicida) ───────
medias_celda <- datos %>%
  group_by(Bloque, Fungicida) %>%
  summarise(Media = mean(Rendimiento_kgha), .groups = "drop")

g5_perfil <- ggplot(medias_celda,
                    aes(x = Bloque, y = Media, group = Fungicida,
                        color = Fungicida)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(values = paleta_trt) +
  labs(
    title    = "Perfil de Medias: Rendimiento por Bloque y Fungicida",
    subtitle = "Líneas aproximadamente paralelas → Sin interacción (asunción del DBCA)",
    x        = "Bloque",
    y        = "Rendimiento Medio (kg/ha)",
    color    = "Fungicida"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  )
print(g5_perfil)

cat("\n💡 INTERPRETACIÓN GRÁFICO DE PERFILES:\n")
cat("   Si las líneas son aproximadamente PARALELAS: la asunción de no-interacción\n")
cat("   (aditividad) del DBCA se cumple. El efecto del fungicida es consistente\n")
cat("   en todos los bloques. Si se cruzan mucho, hay interacción y el DBCA\n")
cat("   básico no es adecuado (necesitaría modelos más complejos).\n")


# ==============================================================================
# SECCIÓN 4: AJUSTE DEL MODELO DBCA
# ==============================================================================

cat("\n\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("SECCIÓN 4: AJUSTE DEL MODELO DBCA\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

# FÓRMULA DEL MODELO: Rendimiento ~ Tratamiento + Bloque
# ¡El orden importa conceptualmente! Primero el tratamiento de interés,
# luego los factores de bloqueo.

modelo_dbca <- aov(Rendimiento_kgha ~ Fungicida + Bloque, data = datos)

cat("\n📐 GRADOS DE LIBERTAD DEL EXPERIMENTO:\n")
k <- nlevels(datos$Fungicida)  # Número de tratamientos
r <- nlevels(datos$Bloque)     # Número de bloques
N <- nrow(datos)               # Total de observaciones

cat(sprintf("   k (tratamientos) = %d\n", k))
cat(sprintf("   r (bloques)      = %d\n", r))
cat(sprintf("   N (total obs.)   = %d = k × r = %d × %d\n", N, k, r))
cat(sprintf("\n   GL Tratamientos = k - 1       = %d - 1 = %d\n", k, k-1))
cat(sprintf("   GL Bloques      = r - 1       = %d - 1 = %d\n", r, r-1))
cat(sprintf("   GL Error        = (k-1)(r-1)  = %d × %d = %d\n", k-1, r-1, (k-1)*(r-1)))
cat(sprintf("   GL Total        = N - 1       = %d - 1 = %d\n", N, N-1))

cat("\n   Verificación: GL_Trat + GL_Bloq + GL_Error =", (k-1), "+", (r-1), "+", (k-1)*(r-1),
    "=", (k-1)+(r-1)+(k-1)*(r-1), "=", N-1, "✅\n")

cat("\n   ⚖️  COMPARACIÓN DE GL_ERROR:\n")
cat(sprintf("   - En DCA:  GL_Error = N - k = %d - %d = %d\n", N, k, N-k))
cat(sprintf("   - En DBCA: GL_Error = (k-1)(r-1) = %d\n", (k-1)*(r-1)))
cat(sprintf("   → El DBCA tiene %d GL menos en el Error (transferidos a Bloques)\n",
            (N-k) - (k-1)*(r-1)))


# ==============================================================================
# SECCIÓN 5: VERIFICACIÓN DE SUPUESTOS
# ==============================================================================

cat("\n\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("SECCIÓN 5: VERIFICACIÓN DE SUPUESTOS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

cat("
El ANOVA asume que los residuos del modelo cumplen 3 supuestos:

  1. NORMALIDAD: Los residuos siguen una distribución Normal N(0, σ²)
     Prueba: Shapiro-Wilk  → H₀: Los residuos son normales
             Si p > 0.05: ✅ No se rechaza H₀ (normalidad ok)

  2. HOMOCEDASTICIDAD: Las varianzas son iguales en todos los tratamientos
     Prueba: Levene       → H₀: σ²₁ = σ²₂ = ... = σ²ₖ
             Si p > 0.05: ✅ No se rechaza H₀ (varianzas homogéneas)

  3. INDEPENDENCIA: Las observaciones son independientes
     Garantizada por el diseño (aleatorización dentro de bloques).

  ⚠️  IMPORTANTE: Los supuestos se evalúan sobre los RESIDUOS del modelo,
      NO sobre los datos crudos.
")

# Extraer residuos del modelo
residuos <- residuals(modelo_dbca)
ajustados <- fitted(modelo_dbca)

# ── Supuesto 1: Normalidad (Shapiro-Wilk) ─────────────────────────────────────
cat("─── PRUEBA 1: NORMALIDAD DE RESIDUOS (Shapiro-Wilk) ───────────────────\n")
sw_test <- shapiro.test(residuos)
print(sw_test)
cat(sprintf("\nEstadístico W = %.4f,  p-valor = %.4f\n", sw_test$statistic, sw_test$p.value))
if (sw_test$p.value > 0.05) {
  cat("✅ Conclusión: p > 0.05. NO se rechaza H₀. Los residuos son normales.\n")
  cat("   El supuesto de normalidad SE CUMPLE.\n")
} else {
  cat("⚠️  Conclusión: p ≤ 0.05. Se rechaza H₀. Los residuos NO son normales.\n")
  cat("   Considerar transformaciones (log, raíz cuadrada) o pruebas no paramétricas.\n")
}

# ── Supuesto 2: Homocedasticidad (Levene) ─────────────────────────────────────
cat("\n─── PRUEBA 2: HOMOCEDASTICIDAD (Levene) ───────────────────────────────\n")
lev_test <- leveneTest(Rendimiento_kgha ~ Fungicida, data = datos)
print(lev_test)
if (lev_test$`Pr(>F)`[1] > 0.05) {
  cat("✅ Conclusión: p > 0.05. NO se rechaza H₀. Las varianzas son homogéneas.\n")
  cat("   El supuesto de homocedasticidad SE CUMPLE.\n")
} else {
  cat("⚠️  Conclusión: p ≤ 0.05. Se rechaza H₀. Las varianzas NO son homogéneas.\n")
  cat("   Alternativas: Transformar datos, usar Welch ANOVA o prueba de Kruskal-Wallis.\n")
}

# ── Gráficos de diagnóstico de supuestos ─────────────────────────────────────
cat("\n📊 Generando gráficos de diagnóstico de residuos...\n")

# Q-Q plot de normalidad
gg_qq <- ggplot(data.frame(residuos = residuos), aes(sample = residuos)) +
  stat_qq(color = "#2E86AB", size = 2.5) +
  stat_qq_line(color = "#E84855", linewidth = 1.2) +
  labs(
    title    = "Q-Q Plot de Normalidad",
    subtitle = "Los puntos deben seguir la línea roja",
    x        = "Cuantiles Teóricos (Normal)",
    y        = "Cuantiles Observados (Residuos)"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"))

# Residuos vs ajustados (homocedasticidad)
gg_resid <- ggplot(data.frame(Ajustados = ajustados, Residuos = residuos),
                   aes(x = Ajustados, y = Residuos)) +
  geom_point(size = 3, color = "#2E86AB", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#E84855", linewidth = 1) +
  geom_smooth(method = "loess", se = TRUE, color = "#F4A261", linewidth = 1) +
  labs(
    title    = "Residuos vs Valores Ajustados",
    subtitle = "Patrón aleatorio → homocedasticidad cumplida",
    x        = "Valores Ajustados (Fitted)",
    y        = "Residuos"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"))

# Histograma de residuos
gg_hist_res <- ggplot(data.frame(residuos = residuos), aes(x = residuos)) +
  geom_histogram(aes(y = after_stat(density)), bins = 8,
                 fill = "#A8C686", color = "white", alpha = 0.85) +
  geom_density(color = "#E84855", linewidth = 1.2) +
  stat_function(fun = dnorm,
                args = list(mean = mean(residuos), sd = sd(residuos)),
                color = "#2E86AB", linewidth = 1.2, linetype = "dashed") +
  labs(
    title    = "Histograma de Residuos",
    subtitle = "Curva azul = distribución Normal esperada",
    x        = "Residuos",
    y        = "Densidad"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"))

# Panel combinado de supuestos
panel_supuestos <- ggarrange(gg_qq, gg_resid, gg_hist_res,
                             ncol = 3, nrow = 1,
                             labels = c("A", "B", "C"))
print(annotate_figure(panel_supuestos,
  top = text_grob("Panel de Diagnóstico de Supuestos - Modelo DBCA",
                  face = "bold", size = 14)))


# ==============================================================================
# SECCIÓN 6: ANÁLISIS DE VARIANZA (ANOVA)
# ==============================================================================

cat("\n\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("SECCIÓN 6: ANÁLISIS DE VARIANZA (ANOVA)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

cat("
HIPÓTESIS:

  Para TRATAMIENTOS (Fungicida):
    H₀: τ₁ = τ₂ = τ₃ = τ₄ = τ₅ = 0  (Ningún fungicida difiere del control)
    H₁: Al menos un τᵢ ≠ 0            (Al menos un fungicida tiene efecto distinto)

  Para BLOQUES:
    H₀: β₁ = β₂ = β₃ = β₄ = 0  (Los bloques no difieren entre sí)
    H₁: Al menos un βⱼ ≠ 0       (El gradiente de fertilidad existe)
    → Si los bloques son significativos: ¡el bloqueo fue ÚTIL y necesario!
")

tabla_anova <- summary(modelo_dbca)
cat("\n📊 TABLA ANOVA:\n")
print(tabla_anova)

# Extraer valores para interpretación
sc_trt   <- tabla_anova[[1]]$`Sum Sq`[1]
sc_blq   <- tabla_anova[[1]]$`Sum Sq`[2]
sc_err   <- tabla_anova[[1]]$`Sum Sq`[3]
sc_total <- sum(tabla_anova[[1]]$`Sum Sq`)
p_trt    <- tabla_anova[[1]]$`Pr(>F)`[1]
p_blq    <- tabla_anova[[1]]$`Pr(>F)`[2]

cat("\n📊 DESCOMPOSICIÓN DE LA VARIACIÓN TOTAL:\n")
cat(sprintf("   SC_Tratamientos = %8.1f  (%.1f%% de SC_Total)\n",
            sc_trt, sc_trt/sc_total*100))
cat(sprintf("   SC_Bloques      = %8.1f  (%.1f%% de SC_Total)\n",
            sc_blq, sc_blq/sc_total*100))
cat(sprintf("   SC_Error        = %8.1f  (%.1f%% de SC_Total)\n",
            sc_err, sc_err/sc_total*100))
cat(sprintf("   SC_Total        = %8.1f  (100%%)\n", sc_total))

cat("\n📊 CONCLUSIONES DEL ANOVA:\n")
if (p_trt < 0.05) {
  cat(sprintf("   ✅ TRATAMIENTOS (Fungicida): p = %.4f < 0.05\n", p_trt))
  cat("      → Se RECHAZA H₀. Existen diferencias significativas entre fungicidas.\n")
  cat("      → Se procede con comparaciones múltiples (Tukey).\n")
} else {
  cat(sprintf("   ⚪ TRATAMIENTOS (Fungicida): p = %.4f ≥ 0.05\n", p_trt))
  cat("      → NO se rechaza H₀. Los fungicidas no difieren significativamente.\n")
}

if (p_blq < 0.05) {
  cat(sprintf("\n   ✅ BLOQUES: p = %.4f < 0.05\n", p_blq))
  cat("      → El gradiente de fertilidad ES significativo.\n")
  cat("      → ¡El BLOQUEO FUE ÚTIL! Habría inflado el error en un DCA.\n")
} else {
  cat(sprintf("\n   ⚪ BLOQUES: p = %.4f ≥ 0.05\n", p_blq))
  cat("      → El efecto de bloques NO es significativo.\n")
  cat("      → El bloqueo no fue necesario en este caso.\n")
}

# Eficiencia relativa del DBCA vs DCA
cat("\n─── EFICIENCIA RELATIVA DEL BLOQUEO ───────────────────────────────────\n")
cm_err <- tabla_anova[[1]]$`Mean Sq`[3]
cm_blq <- tabla_anova[[1]]$`Mean Sq`[2]
gl_blq <- tabla_anova[[1]]$Df[2]
gl_err <- tabla_anova[[1]]$Df[3]

# Fórmula de eficiencia relativa (Cochran & Cox)
cm_err_dca <- (sc_blq + sc_err) / (gl_blq + gl_err)
ef_rel <- ((gl_blq + 1) * (gl_err + 3) * cm_err_dca) /
          ((gl_err + 1) * (gl_blq + 3) * cm_err) * 100

cat(sprintf("   CM_Error en DBCA  = %.2f\n", cm_err))
cat(sprintf("   CM_Error en DCA   = %.2f (estimado)\n", cm_err_dca))
cat(sprintf("   Eficiencia Relativa del DBCA = %.1f%%\n", ef_rel))
if (ef_rel > 100) {
  cat(sprintf("   → El DBCA es %.1f%% más eficiente que el DCA.\n", ef_rel - 100))
  cat("   → El bloqueo redujo el error experimental sustancialmente.\n")
} else {
  cat("   → El DCA habría sido más eficiente (bloques innecesarios).\n")
}


# ==============================================================================
# SECCIÓN 7: COMPARACIONES MÚLTIPLES - PRUEBA DE TUKEY
# ==============================================================================

cat("\n\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("SECCIÓN 7: COMPARACIONES MÚLTIPLES - PRUEBA DE TUKEY (HSD)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

cat("
FUNDAMENTO DE TUKEY:
   La Prueba de Tukey (HSD = Honest Significant Difference) controla la
   tasa de error tipo I de FAMILIA (α_familia = 0.05) al realizar TODAS
   las comparaciones por pares posibles de tratamientos.

   Calcula la Diferencia Mínima Significativa (HSD):
     HSD = q_α(k, GL_error) × √(CM_Error / r)

   Si |ȳᵢ - ȳⱼ| > HSD → La diferencia entre tratamientos i y j es
                           estadísticamente significativa.

   ASIGNACIÓN DE LETRAS: Tratamientos con la MISMA letra NO difieren
   significativamente. Tratamientos con letras DISTINTAS SÍ difieren.
")

# Prueba de Tukey
tukey_resultado <- HSD.test(modelo_dbca, "Fungicida", group = TRUE, console = FALSE)

cat("📊 GRUPOS DE TUKEY (letras de significancia):\n")
print(tukey_resultado$groups)

cat("\n📊 TABLA DE COMPARACIONES PAREADAS:\n")
tukey_pairs <- HSD.test(modelo_dbca, "Fungicida", group = FALSE, console = FALSE)
print(tukey_pairs$comparison)

cat("\n📊 ESTADÍSTICOS DE LA PRUEBA:\n")
cat(sprintf("   Valor crítico q (Tukey) = %.4f\n", tukey_resultado$statistics$Tukey))
cat(sprintf("   HSD (Diferencia Mínima) = %.2f kg/ha\n", tukey_resultado$statistics$MSD))
cat(sprintf("   CM_Error                = %.4f\n", tukey_resultado$statistics$MSerror))
cat(sprintf("   GL_Error                = %d\n", tukey_resultado$statistics$Df))


# ==============================================================================
# SECCIÓN 8: VISUALIZACIONES AVANZADAS CON GGPLOT2
# ==============================================================================

cat("\n\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("SECCIÓN 8: VISUALIZACIONES CON GGPLOT2\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

# Preparar datos de grupos para graficar
grupos_tukey <- tukey_resultado$groups
grupos_tukey$Fungicida <- rownames(grupos_tukey)
grupos_tukey <- grupos_tukey %>%
  rename(Media_Tukey = Rendimiento_kgha, Grupo_letra = groups) %>%
  arrange(desc(Media_Tukey))

# Agregar error estándar
se_datos <- datos %>%
  group_by(Fungicida) %>%
  summarise(
    Media = mean(Rendimiento_kgha),
    SE    = sd(Rendimiento_kgha) / sqrt(n()),
    .groups = "drop"
  )

grupos_plot <- left_join(grupos_tukey, se_datos, by = "Fungicida")
grupos_plot$Fungicida <- factor(grupos_plot$Fungicida,
                                levels = grupos_plot$Fungicida[order(grupos_plot$Media)])

# ── Gráfico 6: Barras con error estándar y letras Tukey ──────────────────────
g6_barras <- ggplot(grupos_plot, aes(x = Fungicida, y = Media, fill = Fungicida)) +
  geom_col(color = "black", width = 0.65, alpha = 0.88) +
  geom_errorbar(aes(ymin = Media - SE, ymax = Media + SE),
                width = 0.25, linewidth = 0.8) +
  geom_text(aes(label = Grupo_letra, y = Media + SE + 30),
            size = 5.5, fontface = "bold") +
  scale_fill_manual(values = paleta_trt) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = comma) +
  labs(
    title    = "Efecto de Fungicidas en el Rendimiento de Frijol",
    subtitle = "Barras de error: ±1 Error Estándar | Letras: Prueba de Tukey (α = 0.05)",
    x        = "Tratamiento Fungicida",
    y        = "Rendimiento Medio (kg/ha)",
    caption  = "Diseño de Bloques Completamente al Azar (DBCA), r = 4 bloques"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    plot.caption  = element_text(color = "gray50", hjust = 0),
    axis.text.x   = element_text(angle = 15, hjust = 1),
    legend.position = "none"
  )
print(g6_barras)

# ── Gráfico 7: Boxplot + Jitter + Letras Tukey (intermedio) ──────────────────
datos_plot <- left_join(datos, grupos_tukey[, c("Fungicida", "Grupo_letra")],
                        by = "Fungicida")
# Posición y de la letra (sobre el máximo de cada tratamiento)
max_por_trt <- datos %>%
  group_by(Fungicida) %>%
  summarise(y_letra = max(Rendimiento_kgha) + 40, .groups = "drop")

grupos_plot2 <- left_join(grupos_plot, max_por_trt, by = "Fungicida")

g7_box_tukey <- ggplot(datos_plot, aes(x = Fungicida, y = Rendimiento_kgha,
                                        fill = Fungicida)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.5,
               color = "gray30") +
  geom_jitter(aes(color = Bloque), width = 0.18, size = 3.5, alpha = 0.9) +
  geom_text(data = grupos_plot2,
            aes(x = Fungicida, y = y_letra, label = Grupo_letra),
            size = 6, fontface = "bold", color = "black") +
  scale_fill_manual(values = paleta_trt, name = "Fungicida") +
  scale_color_manual(values = paleta_blq, name = "Bloque") +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Rendimiento de Frijol por Fungicida",
    subtitle = "Puntos coloreados por bloque | Letras: Tukey HSD (α = 0.05)",
    x        = "Fungicida",
    y        = "Rendimiento (kg/ha)",
    caption  = "DBCA — 4 bloques × 5 tratamientos = 20 unidades experimentales"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    plot.caption  = element_text(color = "gray50"),
    axis.text.x   = element_text(angle = 15, hjust = 1),
    legend.position = "right"
  )
print(g7_box_tukey)

# ── Gráfico 8: Medias ± IC95% (interval plot profesional) ────────────────────
emm <- emmeans(modelo_dbca, ~ Fungicida)
emm_df <- as.data.frame(emm)
emm_df <- emm_df %>%
  left_join(grupos_tukey[, c("Fungicida", "Grupo_letra")], by = "Fungicida") %>%
  mutate(Fungicida = factor(Fungicida, levels = levels(datos$Fungicida)))

g8_ic <- ggplot(emm_df, aes(x = Fungicida, y = emmean, color = Fungicida)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.3, linewidth = 1.2) +
  geom_text(aes(label = Grupo_letra, y = upper.CL + 25),
            size = 5.5, fontface = "bold", color = "black") +
  geom_hline(yintercept = mean(datos$Rendimiento_kgha),
             linetype = "dashed", color = "gray50", linewidth = 0.8) +
  annotate("text", x = 0.6, y = mean(datos$Rendimiento_kgha) + 25,
           label = "Media\nglobal", size = 3.5, color = "gray50") +
  scale_color_manual(values = paleta_trt) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Medias Estimadas con IC 95% — Prueba de Tukey",
    subtitle = "Medias marginales estimadas (emmeans) — DBCA",
    x        = "Fungicida",
    y        = "Rendimiento Estimado (kg/ha)",
    caption  = "IC 95% calculado con el CM_Error del DBCA"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    plot.caption  = element_text(color = "gray50"),
    axis.text.x   = element_text(angle = 15, hjust = 1),
    legend.position = "none"
  )
print(g8_ic)

cat("\n✅ Todos los gráficos generados correctamente.\n")


# ==============================================================================
# SECCIÓN 9: SÍNTESIS Y RECOMENDACIONES AGRONÓMICAS
# ==============================================================================

cat("\n\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("SECCIÓN 9: SÍNTESIS Y RECOMENDACIONES\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

cat("
📋 RESUMEN DEL ANÁLISIS:

  1. DISEÑO: DBCA con k=5 tratamientos (fungicidas) y r=4 bloques
     (gradiente de fertilidad Norte-Sur).

  2. SUPUESTOS: Normalidad y homocedasticidad verificadas mediante
     pruebas de Shapiro-Wilk y Levene respectivamente.

  3. ANOVA: Se detectaron diferencias significativas entre fungicidas
     (p < 0.05) y el bloqueo fue estadísticamente significativo,
     confirmando que el gradiente de fertilidad existía.

  4. TUKEY: [Ver salida de grupos arriba]

  5. RECOMENDACIÓN AGRONÓMICA:
     Fungicida_D muestra el mayor rendimiento medio. Sin embargo, la
     decisión final debe considerar también:
     - Costo del fungicida por hectárea
     - Efecto residual en el suelo
     - Toxicidad para el ambiente y la salud
     - Disponibilidad en el mercado regional

  6. EFICIENCIA DEL BLOQUEO:
     El DBCA fue más eficiente que el DCA al aislar la variabilidad
     del gradiente de suelo, permitiendo una comparación más precisa
     entre fungicidas.
")

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")


# ==============================================================================
# SECCIÓN 10: 📝 TALLER PARA ESTUDIANTES
# ==============================================================================

cat("\n\n╔══════════════════════════════════════════════════════════╗\n")
cat("║          TALLER EVALUATIVO PARA ESTUDIANTES             ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

cat("
CONTEXTO:
─────────
Un investigador del Instituto Nacional de Innovación Agraria (INIA)
desea evaluar el efecto de 4 VARIEDADES DE MAÍZ (Zea mays) sobre el
rendimiento en grano (ton/ha).

PROBLEMA DE DISEÑO:
   Las fechas de siembra varían a lo largo del ciclo agrícola y las
   condiciones climáticas (humedad, temperatura) cambian semana a
   semana. Para controlar esta variabilidad, el investigador decide
   BLOQUEAR por semana de siembra.

DISEÑO:
   k = 4 variedades (tratamientos)
   r = 5 semanas de siembra (bloques)
   N = 4 × 5 = 20 unidades experimentales

BASE DE DATOS: tarea_taller3_dbca.csv
   Variables:
     - Bloque_Semana:  Semana de siembra (Bloque)
     - Variedad:       Variedad de maíz evaluada (Tratamiento)
     - Rendimiento_tha: Rendimiento en toneladas por hectárea

INSTRUCCIONES:
   Complete las 6 etapas del análisis estadístico:

   ETAPA 1: Carga y exploración inicial de datos
   ETAPA 2: Cálculo manual de los Grados de Libertad esperados
   ETAPA 3: Histograma y boxplots exploratorios
   ETAPA 4: Ajuste del modelo DBCA y verificación de supuestos
   ETAPA 5: ANOVA e interpretación de la tabla
   ETAPA 6: Prueba de Tukey y gráfico de resultados con letras

PREGUNTAS DE REFLEXIÓN:
   A. ¿Los bloques (semanas) resultaron significativos?
      ¿Fue útil el bloqueo por fecha de siembra?

   B. ¿Cuántos GL tiene el error en este diseño?
      Compare con el GL del error si se hubiera usado un DCA.

   C. ¿Qué variedad recomienda adoptar para la región?
      Argumente con el análisis estadístico y criterios agronómicos.

   D. ¿Se cumplen los supuestos del modelo? Explique las pruebas
      utilizadas y su interpretación.
")

# ─── ESPACIO PARA QUE LOS ESTUDIANTES ESCRIBAN SU CÓDIGO ───────────────────

cat("─── ETAPA 1: Carga y exploración de datos ──────────────────────────────\n")
# datos_tarea <- read.csv("tarea_taller3_dbca.csv")
# datos_tarea$Bloque_Semana <- factor(datos_tarea$Bloque_Semana)
# datos_tarea$Variedad      <- factor(datos_tarea$Variedad)
# str(datos_tarea)
# summary(datos_tarea)

cat("─── ETAPA 2: Grados de Libertad esperados ──────────────────────────────\n")
# k_tarea <- ___  # Número de variedades
# r_tarea <- ___  # Número de bloques (semanas)
# cat("GL Tratamientos =", k_tarea - 1)
# cat("GL Bloques      =", r_tarea - 1)
# cat("GL Error        =", (k_tarea - 1) * (r_tarea - 1))
# cat("GL Total        =", k_tarea * r_tarea - 1)

cat("─── ETAPA 3: Exploración gráfica ───────────────────────────────────────\n")
# ggplot(datos_tarea, aes(x = Rendimiento_tha)) +
#   geom_histogram(bins = 8, fill = "steelblue", color = "white") +
#   labs(title = "Distribución del Rendimiento de Maíz") +
#   theme_classic()
#
# ggplot(datos_tarea, aes(x = Variedad, y = Rendimiento_tha, fill = Variedad)) +
#   geom_boxplot(alpha = 0.7) +
#   geom_jitter(width = 0.2, size = 2.5) +
#   theme_classic()

cat("─── ETAPA 4: Ajuste del modelo y supuestos ─────────────────────────────\n")
# modelo_tarea <- aov(___ ~ ___ + ___, data = datos_tarea)
# residuos_tarea <- residuals(modelo_tarea)
# shapiro.test(residuos_tarea)
# leveneTest(Rendimiento_tha ~ Variedad, data = datos_tarea)

cat("─── ETAPA 5: Tabla ANOVA ───────────────────────────────────────────────\n")
# summary(modelo_tarea)
# # ¿Los bloques son significativos? ¿Los tratamientos son significativos?

cat("─── ETAPA 6: Tukey y gráfico final ─────────────────────────────────────\n")
# tukey_tarea <- HSD.test(modelo_tarea, "Variedad", group = TRUE)
# print(tukey_tarea$groups)
#
# # Gráfico con letras Tukey
# ggplot(...) + ...

cat("\n✏️  ¡Completa las etapas anteriores y responde las preguntas de reflexión!\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("FIN DEL TALLER 3 - DBCA\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")