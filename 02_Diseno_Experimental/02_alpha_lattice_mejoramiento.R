# ============================================================================
# 02_alpha_lattice_mejoramiento.R
# Diseño Alpha-Lattice para Ensayos de Mejoramiento Genético
# Contexto: Evaluación de 96 líneas de trigo bajo estrés hídrico
# ============================================================================
# Diseño: Alpha-Lattice (0,1) con bloques incompletos
# Paquete: lme4 para modelos mixtos
# Situación real: Ensayos con cientos de genotipos donde un DBCA sería
#                 ineficiente debido a la alta variabilidad del campo
# ============================================================================

library(tidyverse)
library(lme4)        # Modelos mixtos: lmer()
library(lmerTest)    # p-valores para modelos mixtos
library(emmeans)     # Comparaciones múltiples
library(performance) # Diagnóstico de modelos
library(agricolae)   # Generación del diseño alpha-lattice
library(car)         # Pruebas estadísticas

set.seed(456)

# ============================================================================
# 1. FUNDAMENTOS DEL DISEÑO ALPHA-LATTICE
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 1: FUNDAMENTOS DEL DISEÑO ALPHA-LATTICE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("
¿QUÉ ES UN DISEÑO ALPHA-LATTICE?
--------------------------------
El diseño Alpha-Lattice (Patterson & Williams, 1976) es un diseño de
bloques incompletos balanceados utilizado cuando:

1. Hay MUCHOS tratamientos (ej. 50-500 genotipos)
2. El tamaño de bloque completo sería demasiado grande
3. La variabilidad dentro del bloque sería alta

VENTAJAS vs DBCA:
- Mayor eficiencia en la estimación de efectos de tratamientos
- Mejor control de la variabilidad espacial
- Permite evaluar cientos de genotipos en un solo ensayo

ESTRUCTURA:
- Genotipos se agrupan en BLOQUES INCOMPLETOS
- Cada bloque incompleto contiene k genotipos (k < total)
- El diseño se replica r veces (repeticiones completas)
", sep = "")

# ============================================================================
# 2. GENERACIÓN DEL DISEÑO ALPHA-LATTICE CON AGRICOLAE
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 2: GENERACIÓN DEL DISEÑO ALPHA-LATTICE\n")
cat(rep("=", 80), "\n\n", sep = "")

# Parámetros del diseño
n_genotipos <- 96       # Número de líneas de trigo a evaluar
n_repeticiones <- 3     # Número de repeticiones (bloques completos)
k <- 12                 # Tamaño del bloque incompleto (genotipos por bloque)

# Verificar que el diseño es posible
if (n_genotipos %% k != 0) {
  warning("El número de genotipos debe ser múltiplo del tamaño de bloque")
}

n_bloques_incompletos <- n_genotipos / k  # 8 bloques incompletos por repetición

cat("PARÁMETROS DEL DISEÑO:\n")
cat("  Número de genotipos:", n_genotipos, "\n")
cat("  Repeticiones:", n_repeticiones, "\n")
cat("  Tamaño de bloque incompleto (k):", k, "\n")
cat("  Bloques incompletos por repetición:", n_bloques_incompletos, "\n")
cat("  Total de bloques incompletos:", n_bloques_incompletos * n_repeticiones, "\n")
cat("  Total de parcelas:", n_genotipos * n_repeticiones, "\n\n")

# Generar designios alpha-lattice con agricolae
# design.alpha() genera el diseño óptimo
design_alpha <- design.alpha(
  trt = paste0("L", sprintf("%03d", 1:n_genotipos)),  # Nombres de genotipos
  reps = n_repeticiones,
  k = k,
  seed = 456
)

cat("=== LIBRO DE CAMPO (PRIMERAS 20 FILAS) ===\n")
print(head(design_alpha$book, 20))

cat("\n=== ESTRUCTURA DE BLOQUES INCOMPLETOS ===\n")
cat("Total de observaciones:", nrow(design_alpha$book), "\n")
cat("Número de bloques incompletos únicos:",
    length(unique(design_alpha$book$block)), "\n")

# Verificar balance del diseño
cat("\n=== VERIFICACIÓN DEL BALANCE ===\n")
tabla_genotipos <- table(design_alpha$book$trt)
cat("Cada genotipo aparece", min(tabla_genotipos), "-", max(tabla_genotipos), "veces\n")
cat("¿Diseño balanceado?",
    ifelse(min(tabla_genotipos) == max(tabla_genotipos), "SÍ ✓", "NO ⚠"), "\n")

# ============================================================================
# 3. SIMULACIÓN DE DATOS DE ENSAYO DE CAMPO
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 3: SIMULACIÓN DE DATOS DE CAMPO\n")
cat(rep("=", 80), "\n\n", sep = "")

# Preparar datos para simulación
datos <- design_alpha$book %>%
  mutate(
    rep = factor(rep),
    block = factor(paste0(rep, "_", block)),  # ID único de bloque incompleto
    genotype = factor(trt)
  )

# Simular efectos genotípicos
# En mejoramiento, esperamos que algunos genotipos sean superiores
set.seed(456)
efectos_genotipo <- rnorm(n_genotipos, 0, 4)  # Varianza genotípica
names(efectos_genotipo) <- paste0("L", sprintf("%03d", 1:n_genotipos))

# Simular efectos de repetición (variabilidad entre replicaciones)
efectos_rep <- c("1" = 0, "2" = -1.2, "3" = 0.8)  # Diferencias entre reps

# Simular efectos de bloque incompleto (variabilidad dentro de cada rep)
n_bloques_total <- n_bloques_incompletos * n_repeticiones
efectos_bloque <- rnorm(n_bloques_total, 0, 2.5)
names(efectos_bloque) <- levels(datos$block)

# Simular rendimiento (kg/ha)
datos <- datos %>%
  mutate(
    efecto_gen = efectos_genotipo[as.character(genotype)],
    efecto_rep = efectos_rep[as.character(rep)],
    efecto_bloq = efectos_bloque[as.character(block)],
    error = rnorm(n(), 0, 3.0),
    # Rendimiento base: 4500 kg/ha + efectos
    rendimiento = 4500 + efecto_gen + efecto_rep + efecto_bloq + error
  ) %>%
  select(-efecto_gen, -efecto_rep, -efecto_bloq, -error)

cat("=== RESUMEN DESCRIPTIVO ===\n")
cat("Rendimiento medio:", round(mean(datos$rendimiento), 1), "kg/ha\n")
cat("Desviación estándar:", round(sd(datos$rendimiento), 1), "kg/ha\n")
cat("Rango:", round(min(datos$rendimiento), 1), "-",
    round(max(datos$rendimiento), 1), "kg/ha\n\n")

cat("=== DISTRIBUCIÓN POR REPETICIÓN ===\n")
datos %>%
  group_by(rep) %>%
  summarise(
    n = n(),
    media = mean(rendimiento),
    sd = sd(rendimiento),
    cv = sd / media * 100,
    .groups = "drop"
  ) %>%
  print()

# ============================================================================
# 4. ANÁLISIS CON MODELO MIXTO PARA ALPHA-LATTICE
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 4: ANÁLISIS CON MODELO MIXTO\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("MODELO MIXTO PARA ALPHA-LATTICE:\n")
cat("--------------------------------\n")
cat("y_ijk = μ + Rep_i + Gen_j + B(Rep)_ik + ε_ijk\n\n")
cat("Donde:\n")
cat("  Rep_i     = Efecto fijo de la repetición i\n")
cat("  Gen_j     = Efecto fijo/aleatorio del genotipo j\n")
cat("  B(Rep)_ik = Efecto ALEATORIO del bloque incompleto k dentro de rep i\n")
cat("  ε_ijk     = Error experimental\n\n")

cat("NOTA: Los bloques incompletos se modelan como ALEATORIOS porque:\n")
cat("1. Representan una muestra de la variabilidad espacial\n")
cat("2. Queremos hacer inferencia más allá de estos bloques específicos\n")
cat("3. Esto permite 'shrinkage' y mejores predicciones BLUP\n\n")

# ----------------------------------------------------------------------------
# 4.1 MODELO CON GENOTIPOS COMO EFECTO FIJO (para pruebas F)
# ----------------------------------------------------------------------------

cat("\n4.1 MODELO CON GENOTIPOS FIJOS (ANOVA tipo III)\n")
cat("-" , rep("-", 60), "\n", sep = "")

modelo_fijo <- lmer(
  rendimiento ~ rep + genotype + (1 | block),
  data = datos,
  REML = FALSE
)

cat("\nResumen del modelo:\n")
print(summary(modelo_fijo))

cat("\nTabla ANOVA (Tipo III con aproximación Satterthwaite):\n")
anova_fijo <- Anova(modelo_fijo, type = "III", test.statistic = "F")
print(anova_fijo)

# Interpretar resultados
cat("\n=== INTERPRETACIÓN ===\n")
genotipo_p <- anova_fijo$`Pr(>F)`[2]
cat("Efecto Genotipo: F =", round(anova_fijo$`F value`[2], 2),
    ", p <", ifelse(genotipo_p < 0.001, "0.001", round(genotipo_p, 4)), "\n")
if (genotipo_p < 0.05) {
  cat("→ DIFERENCIAS SIGNIFICATIVAS entre genotipos ✓\n")
  cat("→ Existe variabilidad genética aprovechable para selección\n")
} else {
  cat("→ NO hay diferencias significativas entre genotipos\n")
}

# ----------------------------------------------------------------------------
# 4.2 MODELO CON GENOTIPOS COMO EFECTO ALEATORIO (para BLUPs)
# ----------------------------------------------------------------------------

cat("\n4.2 MODELO CON GENOTIPOS ALEATORIOS (Predicción BLUP)\n")
cat("-" , rep("-", 60), "\n", sep = "")

modelo_aleatorio <- lmer(
  rendimiento ~ rep + (1 | genotype) + (1 | block),
  data = datos,
  REML = TRUE
)

cat("\nResumen del modelo:\n")
print(summary(modelo_aleatorio))

cat("\nComponentes de varianza:\n")
var_comp <- as.data.frame(VarCorr(modelo_aleatorio))
print(var_comp)

# Calcular parámetros genéticos importantes
var_gen <- var_comp$vcov[var_comp$grp == "genotype"]
var_bloque <- var_comp$vcov[var_comp$grp == "block"]
var_error <- var_comp$vcov[var_comp$grp == "Residual"]

cat("\n=== PARÁMETROS GENÉTICOS ===\n")

# Heredabilidad en sentido amplio (H²)
# H² = σ²g / (σ²g + σ²b/r + σ²e/rk)
r <- n_repeticiones
k <- n_bloques_incompletos
h2 <- var_gen / (var_gen + var_bloque/r + var_error/(r*k))
cat("Heredabilidad (H²):", round(h2, 3), "\n")
if (h2 > 0.5) {
  cat("  → Heredabilidad ALTA: buena respuesta a selección\n")
} else if (h2 > 0.3) {
  cat("  → Heredabilidad MODERADA: selección posible\n")
} else {
  cat("  → Heredabilidad BAJA: difícil mejorar por selección\n")
}

# Coeficiente de variación genotípica
cv_gen <- sqrt(var_gen) / mean(datos$rendimiento) * 100
cat("CV genotípico:", round(cv_gen, 2), "%\n")

# Eficiencia del diseño alpha-lattice
# Comparar varianza del bloque incompleto vs residual
eficiencia_bloqueo <- (1 - var_bloque / (var_bloque + var_error)) * 100
cat("Eficiencia de bloques incompletos:", round(eficiencia_bloqueo, 1), "%\n")

# ============================================================================
# 5. PREDICCIÓN DE GENOTIPOS SUPERIORES (BLUPs)
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 5: PREDICCIÓN DE GENOTIPOS CON BLUPs\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("¿QUÉ SON LOS BLUPs?\n")
cat("-------------------\n")
cat("BLUP = Best Linear Unbiased Predictor\n")
cat("- Son las mejores predicciones de efectos genotípicos\n")
cat("- Incorporan 'shrinkage': genotipos con menos datos se contraen hacia la media\n")
cat("- Ideales para rankear y seleccionar genotipos\n\n")

# Extraer BLUPs de genotipos
blups_genotipo <- ranef(modelo_aleatorio)$genotype %>%
  as_tibble(rownames = "genotype") %>%
  rename(blup = `(Intercept)`) %>%
  mutate(
    rendimiento_predicho = mean(datos$rendimiento) + blup,
    rank = rank(desc(rendimiento_predicho))
  ) %>%
  arrange(desc(rendimiento_predicho))

cat("=== TOP 15 GENOTIPOS SEGÚN BLUPs ===\n")
print(head(blups_genotipo, 15))

cat("\n=== PEORES 15 GENOTIPOS SEGÚN BLUPs ===\n")
print(tail(blups_genotipo, 15))

# Guardar rankings para selección
write_csv(blups_genotipo, "02_Diseno_Experimental/alpha_lattice_rankings.csv")
cat("\nRankings guardados en: 02_Diseno_Experimental/alpha_lattice_rankings.csv\n")

# ============================================================================
# 6. COMPARACIONES MÚLTIPLES CON EMMEANS
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 6: COMPARACIONES MÚLTIPLES (EMMEANS)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("NOTA: Con 96 genotipos, las comparaciones pairwise generan\n")
cat("4,560 comparaciones (96*95/2). Mostraremos solo el TOP 10.\n\n")

# Calcular emmeans para genotipos
emm_gen <- emmeans(modelo_fijo, ~ genotype)

# Obtener los 10 mejores genotipos según BLUPs
top_10 <- blups_genotipo %>% head(10) %>% pull(genotype)

cat("=== MEDIAS AJUSTADAS - TOP 10 GENOTIPOS ===\n")
emm_top10 <- emm_gen %>%
  as.data.frame() %>%
  filter(genotype %in% top_10) %>%
  arrange(desc(emmean))
print(emm_top10)

cat("\n=== COMPARACIONES PAIRWISE DENTRO DEL TOP 10 ===\n")
# Comparaciones solo entre los mejores genotipos
pairs(emm_gen, adjust = "tukey") %>%
  as.data.frame() %>%
  filter(str_detect(contrast, paste(top_10, collapse = "|"))) %>%
  head(20) %>%
  print()

# Letras de agrupamiento para el TOP 10
cat("\n=== LETRAS DE AGRUPAMIENTO (TOP 10) ===\n")
cld_top10 <- cld(emm_gen, alpha = 0.05, Letters = letters) %>%
  as.data.frame() %>%
  filter(genotype %in% top_10) %>%
  arrange(desc(emmean))
print(cld_top10)

# ============================================================================
# 7. EVALUACIÓN DE SUPUESTOS DEL MODELO MIXTO
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 7: EVALUACIÓN DE SUPUESTOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Extraer residuos
residuos <- residuals(modelo_aleatorio)
ajustados <- fitted(modelo_aleatorio)

# ---- 7.1 Normalidad de Residuos ----
cat("\n7.1 Normalidad de Residuos (Shapiro-Wilk)\n")
cat("-" , rep("-", 50), "\n", sep = "")

shapiro_test <- shapiro.test(residuos)
cat("Estadístico W =", round(shapiro_test$statistic, 4), "\n")
cat("p-value =", round(shapiro_test$p.value, 4), "\n")

if (shapiro_test$p.value > 0.05) {
  cat("→ Residuos siguen distribución normal ✓\n")
} else {
  cat("→ Desviación de normalidad detectada ⚠\n")
}

# ---- 7.2 Gráficos de Diagnóstico ----
cat("\n7.2 Diagnóstico Gráfico (check_model)\n")
cat("-" , rep("-", 50), "\n", sep = "")

check_mixto <- check_model(modelo_aleatorio, panel = TRUE)
print(check_mixto)
ggsave("02_Diseno_Experimental/02_alpha_lattice_diagnostico.png",
       plot = check_mixto, width = 10, height = 8, dpi = 300)

# ---- 7.3 Gráfico de Residuos vs Ajustados ----
residuo_plot <- ggplot(data.frame(ajustados = ajustados, residuos = residuos),
                       aes(x = ajustados, y = residuos)) +
  geom_point(alpha = 0.4, color = "#2E86AB", size = 2) +
  geom_hline(yintercept = 0, color = "#E74C3C", linetype = "dashed") +
  labs(title = "Residuos vs Valores Ajustados",
       subtitle = "Alpha-Lattice: 96 genotipos, 3 repeticiones",
       x = "Valores ajustados", y = "Residuos") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

print(residuo_plot)
ggsave("02_Diseno_Experimental/02_alpha_lattice_residuos.png",
       plot = residuo_plot, width = 7, height = 5, dpi = 300)

# ============================================================================
# 8. VISUALIZACIÓN DE RESULTADOS
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 8: VISUALIZACIÓN DE RESULTADOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Gráfico de los 20 mejores genotipos con barras de error
top_20 <- blups_genotipo %>% head(20)

# Calcular errores estándar de BLUPs
# SE(BLUP) ≈ sqrt(σ²e / r) para diseño balanceado
se_blup <- sqrt(var_error / n_repeticiones)

plot_top20 <- ggplot(top_20, aes(x = reorder(genotype, rendimiento_predicho),
                                  y = rendimiento_predicho)) +
  geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
  geom_errorbar(aes(ymin = rendimiento_predicho - se_blup,
                    ymax = rendimiento_predicho + se_blup),
                width = 0.3, color = "#E74C3C") +
  geom_hline(yintercept = mean(datos$rendimiento),
             linetype = "dashed", color = "gray40") +
  coord_flip() +
  labs(title = "Top 20 Genotipos de Trigo - Rendimiento BLUP",
       subtitle = "Línea punteada = media del ensayo (4500 kg/ha)",
       x = "Genotipo", y = "Rendimiento Predicho (kg/ha)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9),
        panel.grid.minor = element_blank())

print(plot_top20)
ggsave("02_Diseno_Experimental/02_alpha_lattice_top20.png",
       plot = plot_top20, width = 9, height = 7, dpi = 300)

# Histograma de rendimientos de todos los genotipos
hist_plot <- ggplot(blups_genotipo, aes(x = rendimiento_predicho)) +
  geom_histogram(binwidth = 50, fill = "#2E86AB", color = "white", alpha = 0.7) +
  geom_vline(xintercept = mean(datos$rendimiento),
             color = "#E74C3C", linewidth = 1.2, linetype = "dashed") +
  labs(title = "Distribución de Rendimientos de Genotipos",
       subtitle = paste("Media =", round(mean(datos$rendimiento), 0), "kg/ha"),
       x = "Rendimiento Predicho (kg/ha)", y = "Frecuencia") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

print(hist_plot)
ggsave("02_Diseno_Experimental/02_alpha_lattice_histograma.png",
       plot = hist_plot, width = 7, height = 5, dpi = 300)

# ============================================================================
# 9. RESUMEN FINAL Y RECOMENDACIONES DE SELECCIÓN
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("RESUMEN DEL ENSAYO ALPHA-LATTICE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("DISEÑO EXPERIMENTAL:\n")
cat("  Tipo: Alpha-Lattice (0,1)\n")
cat("  Genotipos:", n_genotipos, "\n")
cat("  Repeticiones:", n_repeticiones, "\n")
cat("  Tamaño de bloque incompleto:", k, "\n")
cat("  Total de parcelas:", nrow(datos), "\n\n")

cat("COMPONENTES DE VARIANZA:\n")
cat("  σ² Genotipo:", round(var_gen, 2), "\n")
cat("  σ² Bloque incompleto:", round(var_bloque, 2), "\n")
cat("  σ² Error:", round(var_error, 2), "\n\n")

cat("PARÁMETROS GENÉTICOS:\n")
cat("  Heredabilidad (H²):", round(h2, 3), "\n")
cat("  CV genotípico:", round(cv_gen, 2), "%\n")
cat("  Eficiencia de bloqueo:", round(eficiencia_bloqueo, 1), "%\n\n")

cat("GENOTIPOS SELECCIONADOS (TOP 5):\n")
for (i in 1:5) {
  cat("  ", i, ". ", blups_genotipo$genotype[i],
      " (", round(blups_genotipo$rendimiento_predicho[i], 0), " kg/ha)\n", sep = "")
}

cat("\nRECOMENDACIÓN:\n")
if (h2 > 0.5) {
  cat("  → Alta heredabilidad: seleccionar TOP 10-20% de genotipos\n")
  cat("  → Avanzar a ensayos multi-ambientales\n")
} else if (h2 > 0.3) {
  cat("  → Heredabilidad moderada: seleccionar TOP 5-10%\n")
  cat("  → Considerar selección indexada con otras variables\n")
} else {
  cat("  → Baja heredabilidad: aumentar repeticiones o usar más parcelas\n")
  cat("  → Considerar selección familiar en lugar de individual\n")
}

# ============================================================================
# 10. EJERCICIO AVANZADO: ENSAYO MULTI-AMBIENTAL (MET) Y ESTRÉS TÉRMICO
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 10: EJERCICIO AVANZADO - MET Y ESTRÉS TÉRMICO\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("ESCENARIO:\n")
cat("Evaluamos los mismos 96 genotipos en 2 localidades:
1. Localidad A (Ambiente Óptimo): 20-25°C
2. Localidad B (Ambiente de Estrés Térmico): 35-40°C durante floración.
El objetivo es identificar genotipos tolerantes al calor (GxE).\n\n")

# 1. Simulación de datos MET (Genotipo x Ambiente)
set.seed(321)
ambientes <- c("Optimo", "Estres_Termico")

datos_met <- expand.grid(
  genotype = levels(datos$genotype),
  env = ambientes,
  rep = factor(1:2) # 2 reps por sitio
) %>%
  mutate(
    # Efecto del ambiente
    efecto_env = ifelse(env == "Optimo", 0, -1200),
    # GxE: Algunos genotipos sufren más que otros
    tol_calor = rnorm(n_genotipos, 0, 300)[as.numeric(genotype)],
    efecto_gxe = ifelse(env == "Estres_Termico", tol_calor, 0),
    # Rendimiento
    rendimiento = 4500 + efecto_env + efecto_gxe + rnorm(n(), 0, 200)
  )

# 2. Análisis GxE con Modelo Mixto
cat("=== Análisis de Interacción Genotipo x Ambiente (GxE) ===\n")
mod_met <- lmer(rendimiento ~ env + (1|genotype) + (1|genotype:env), 
                data = datos_met)
print(summary(mod_met))

# 3. Estabilidad y Selección
cat("\n=== Componentes de Varianza GxE ===\n")
var_met <- as.data.frame(VarCorr(mod_met))
print(var_met)

cat("\nINTERPRETACIÓN GxE:\n")
cat("Si la varianza G:E es alta, la selección debe hacerse por ambiente.
Genotipos con BLUPs positivos en 'Estres_Termico' son candidatos para
programas de mejoramiento en zonas cálidas.\n")

# 4. Gráfico de Estabilidad (Reactividad al ambiente)
library(ggrepel)
df_estabilidad <- datos_met %>%
  group_by(genotype, env) %>%
  summarise(media = mean(rendimiento), .groups = "drop") %>%
  pivot_wider(names_from = env, values_from = media)

plot_met <- ggplot(df_estabilidad, aes(x = Optimo, y = Estres_Termico)) +
  geom_point(color = "#2E86AB", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_text_repel(data = head(df_estabilidad[order(df_estabilidad$Estres_Termico, decreasing = T),], 5),
                  aes(label = genotype), box.padding = 0.5) +
  labs(title = "Estabilidad de Genotipos: Óptimo vs Estrés Térmico",
       subtitle = "Genotipos sobre la línea roja rinden proporcionalmente mejor bajo estrés",
       x = "Rendimiento en Ambiente Óptimo",
       y = "Rendimiento en Estrés Térmico") +
  theme_minimal()

print(plot_met)
ggsave("02_Diseno_Experimental/02_alpha_lattice_met.png", plot_met, width = 8, height = 6)

cat("\n", rep("=", 80), "\n", sep = "")
cat("FIN DEL ANÁLISIS MET\n")
cat(rep("=", 80), "\n", sep = "")
