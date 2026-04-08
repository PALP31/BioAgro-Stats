# ============================================================================
# 03_parcelas_divididas_splitplot.R
# Diseño de Parcelas Divididas (Split-Plot)
# Contexto: Riego x Trichoderma en cultivo de tomate
# ============================================================================
# Diseño: Split-Plot en Bloques Completos al Azar
# Parcela Principal (Main Plot): Régimen de riego
# Subparcela (Subplot): Inoculación con Trichoderma
# Paquete: lme4 para modelos mixtos con estructura de error correcta
# ============================================================================

library(tidyverse)
library(lme4)        # Modelos mixtos: lmer(), glmer()
library(lmerTest)    # p-valores para modelos mixtos
library(emmeans)     # Comparaciones múltiples
library(performance) # Diagnóstico de modelos
library(car)         # ANOVA tipo III
library(agricolae)   # Diseño split-plot

set.seed(321)

# ============================================================================
# 1. FUNDAMENTOS DEL DISEÑO SPLIT-PLOT
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 1: FUNDAMENTOS DEL DISEÑO SPLIT-PLOT\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("
¿QUÉ ES UN DISEÑO DE PARCELAS DIVIDIDAS (SPLIT-PLOT)?
-----------------------------------------------------
El diseño split-plot se utiliza cuando:

1. Hay factores que requieren diferentes tamaños de unidad experimental
2. Un factor es difícil/costoso de aplicar (se aplica a parcelas grandes)
3. Otro factor es fácil de aplicar (se aplica a subparcelas pequeñas)
4. Se quiere estudiar la interacción entre ambos factores

ESTRUCTURA JERÁRQUICA:
----------------------
BLOQUES
  └── PARCELA PRINCIPAL (Main Plot)
        └── SUBPARCELA (Subplot)

TÉRMINOS DE ERROR:
------------------
En un diseño split-plot hay DOS errores experimentales:

1. Error A (Main Plot Error): Variabilidad entre parcelas principales
   - Se usa para probar el efecto del factor A (riego)
   - Mayor varianza porque las parcelas son más grandes

2. Error B (Subplot Error): Variabilidad entre subparcelas
   - Se usa para probar el factor B (Trichoderma) y la interacción A×B
   - Menor varianza, mayor precisión

MODELO ESTADÍSTICO:
-------------------
y_ijk = μ + Bloque_i + A_j + ErrorA_ij + B_k + (A×B)_jk + ErrorB_ijk

Donde:
  Bloque_i  = Efecto del bloque i (aleatorio)
  A_j       = Efecto del factor A (riego) - FIJO
  ErrorA_ij = Error de parcela principal (aleatorio)
  B_k       = Efecto del factor B (Trichoderma) - FIJO
  (A×B)_jk  = Interacción - FIJA
  ErrorB_ijk = Error de subparcela (aleatorio)
", sep = "")

# ============================================================================
# 2. GENERACIÓN DEL DISEÑO SPLIT-PLOT
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 2: GENERACIÓN DEL DISEÑO\n")
cat(rep("=", 80), "\n\n", sep = "")

# Parámetros del diseño
n_bloques <- 4
riego <- c("Control", "Deficit_moderado", "Deficit_severo")  # Factor A (Main Plot)
trichoderma <- c("Sin_inocular", "T_asperellum", "T_harzianum", "T_viride")  # Factor B (Subplot)

# Calcular total de unidades experimentales
n_parcelas_principales <- n_bloques * length(riego)
n_subparcelas <- n_parcelas_principales * length(trichoderma)

cat("PARÁMETROS DEL DISEÑO:\n")
cat("  Bloques:", n_bloques, "\n")
cat("  Factor A (Riego):", length(riego), "niveles\n")
cat("  Factor B (Trichoderma):", length(trichoderma), "niveles\n")
cat("  Parcelas principales:", n_parcelas_principales, "\n")
cat("  Subparcelas totales:", n_subparcelas, "\n\n")

# Generar diseño split-plot manualmente
# En campo, primero se asignan los riegos a parcelas grandes dentro de cada bloque
# Luego se asignan los Trichoderma a subparcelas dentro de cada parcela principal

datos <- expand.grid(
  Bloque = factor(paste0("B", 1:n_bloques)),
  Riego = factor(riego),
  Trichoderma = factor(trichoderma)
)

# Crear ID único para cada parcela principal (combinación Bloque x Riego)
datos <- datos %>%
  mutate(
    Parcela_ID = factor(paste0(Bloque, "_", Riego)),
    # Ordenar para simular la estructura jerárquica
    Row = as.numeric(Bloque) * 10 + as.numeric(Riego)
  ) %>%
  arrange(Row, Trichoderma) %>%
  select(-Row)

cat("=== ESTRUCTURA DEL DISEÑO (primeras 20 filas) ===\n")
print(head(datos, 20))

cat("\n=== VERIFICACIÓN DE LA ESTRUCTURA JERÁRQUICA ===\n")
cat("Número de parcelas principales únicas:", length(unique(datos$Parcela_ID)), "\n")
cat("Cada parcela principal contiene:",
    nrow(datos) / length(unique(datos$Parcela_ID)), "subparcelas\n")
cat("Total de observaciones:", nrow(datos), "\n")

# ============================================================================
# 3. SIMULACIÓN DE DATOS
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 3: SIMULACIÓN DE DATOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Efectos verdaderos de los tratamientos
efecto_riego <- c(
  "Control" = 0,
  "Deficit_moderado" = -8,    # Estrés hídrico moderado reduce rendimiento
  "Deficit_severo" = -18      # Estrés severo tiene mayor impacto
)

efecto_trichoderma <- c(
  "Sin_inocular" = 0,
  "T_asperellum" = 5.5,       # T. asperellum mejora rendimiento
  "T_harzianum" = 4.2,
  "T_viride" = 3.8
)

# Interacción: Trichoderma mitiga el efecto del estrés hídrico
interaccion <- expand.grid(Riego = riego, Trichoderma = trichoderma)
interaccion$efecto <- with(interaccion,
  ifelse(Riego == "Control", 0,
  ifelse(Trichoderma == "Sin_inocular", 0,
  # Trichoderma ayuda más bajo estrés
  ifelse(Riego == "Deficit_moderado", 2,
  ifelse(Riego == "Deficit_severo", 4, 0))))
)

# Efectos aleatorios
efecto_bloque <- rnorm(n_bloques, 0, 2.5)  # Variabilidad entre bloques
names(efecto_bloque) <- paste0("B", 1:n_bloques)

# Error de parcela principal (variabilidad entre parcelas grandes)
error_parcela <- rnorm(length(unique(datos$Parcela_ID)), 0, 3.0)
names(error_parcela) <- levels(datos$Parcela_ID)

# Datos simulados
datos <- datos %>%
  mutate(
    efecto_bloq = efecto_bloque[as.character(Bloque)],
    efecto_riego = efecto_riego[as.character(Riego)],
    efecto_trich = efecto_trichoderma[as.character(Trichoderma)],
    efecto_int = interaccion$efecto[match(paste(Riego, Trichoderma),
                                          paste(interaccion$Riego, interaccion$Trichoderma))],
    error_parcela = error_parcela[as.character(Parcela_ID)],
    error_subparcela = rnorm(n(), 0, 2.0),
    # Rendimiento base: 65 kg/parcela
    rendimiento = 65 + efecto_bloq + efecto_riego + error_parcela +
                  efecto_trich + efecto_int + error_subparcela
  ) %>%
  select(-efecto_bloq, -efecto_riego, -efecto_trich, -efecto_int,
         -error_parcela, -error_subparcela)

cat("=== RESUMEN DESCRIPTIVO ===\n")
cat("Rendimiento medio:", round(mean(datos$rendimiento), 2), "kg/parcela\n")
cat("Desviación estándar:", round(sd(datos$rendimiento), 2), "kg/parcela\n\n")

cat("=== MEDIAS POR FACTOR ===\n")
cat("\nPor Riego:\n")
datos %>%
  group_by(Riego) %>%
  summarise(n = n(), media = mean(rendimiento), sd = sd(rendimiento)) %>%
  print()

cat("\nPor Trichoderma:\n")
datos %>%
  group_by(Trichoderma) %>%
  summarise(n = n(), media = mean(rendimiento), sd = sd(rendimiento)) %>%
  print()

# ============================================================================
# 4. ANÁLISIS CON aov() - ENFOQUE CLÁSICO
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 4: ANÁLISIS CLÁSICO CON aov()\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("IMPORTANTE: Especificación correcta de términos de error\n")
cat("---------------------------------------------------------\n")
cat("En split-plot, debemos declarar DOS errores:\n\n")
cat("1. Error(Parcela_ID) para el factor A (Riego)\n")
cat("   - Parcela_ID = Bloque x Riego (unidad experimental del factor A)\n")
cat("   - Este error captura la variabilidad entre parcelas principales\n\n")
cat("2. Error residual para el factor B (Trichoderma) y la interacción\n")
cat("   - La subparcela es la unidad experimental del factor B\n")
cat("   - Mayor precisión (menor error) que el factor A\n\n")

# ----------------------------------------------------------------------------
# 4.1 MODELO aov() CON ERROR SPLIT-PLOT
# ----------------------------------------------------------------------------

cat("\n4.1 AJUSTE DEL MODELO aov() CON Error()\n")
cat("-" , rep("-", 60), "\n", sep = "")

# La sintaxis correcta en aov() para split-plot es:
# aov(y ~ A * B + Error(Bloque/A), data)
# Donde Bloque/A expande a: Bloque + Bloque:A (error de parcela principal)

modelo_aov <- aov(
  rendimiento ~ Riego * Trichoderma + Error(Parcela_ID),
  data = datos
)

cat("\n=== TABLA ANOVA SPLIT-PLOT ===\n")
print(summary(modelo_aov))

# Extraer las tablas ANOVA por estrato
cat("\n=== INTERPRETACIÓN POR ESTRATO ===\n")
cat("\n[STRATUM 1: Parcela_ID (Error de Main Plot)]\n")
cat("  Aquí se prueba el efecto del RIEGO\n")
cat("  El Error A es la variabilidad entre parcelas principales\n")

cat("\n[STRATUM 2: Within (Error de Subplot)]\n")
cat("  Aquí se prueba Trichoderma y la interacción Riego x Trichoderma\n")
cat("  El Error B es la variabilidad entre subparcelas\n")

# ----------------------------------------------------------------------------
# 4.2 EXTRAER COMPONENTES DE VARIANZA
# ----------------------------------------------------------------------------

cat("\n4.2 COMPONENTES DE VARIANZA POR ESTRATO\n")
cat("-" , rep("-", 60), "\n", sep = "")

# Extraer suma de cuadrados y grados de libertad
anova_riego <- summary(modelo_aov)[[1]]$statistics
anova_sub <- summary(modelo_aov)[[2]]$statistics

cat("\nANOVA del Factor RIEGO (Stratum Parcela_ID):\n")
print(anova_riego)

cat("\nANOVA de Trichoderma e Interacción (Stratum Within):\n")
print(anova_sub)

# Calcular cuadrados medios
cat("\n=== CUADRADOS MEDIOS ===\n")
cm_riego <- anova_riego[1, 2]  # MS Riego
cm_error_a <- anova_riego[2, 2]  # MS Error A
cat("CM Riego:", round(cm_riego, 2), "\n")
cat("CM Error A (Parcela):", round(cm_error_a, 2), "\n")
cat("F Riego =", round(cm_riego / cm_error_a, 2), "\n")

cm_trich <- anova_sub[1, 2]
cm_inter <- anova_sub[2, 2]
cm_error_b <- anova_sub[3, 2]
cat("\nCM Trichoderma:", round(cm_trich, 2), "\n")
cat("CM Interacción:", round(cm_inter, 2), "\n")
cat("CM Error B (Subparcela):", round(cm_error_b, 2), "\n")
cat("F Trichoderma =", round(cm_trich / cm_error_b, 2), "\n")
cat("F Interacción =", round(cm_inter / cm_error_b, 2), "\n")

# ============================================================================
# 5. ANÁLISIS CON MODELO MIXTO (lmer) - ENFOQUE MODERNO
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 5: ANÁLISIS CON MODELO MIXTO (lmer)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("VENTAJAS DEL ENFOQUE MIXTO:\n")
cat("---------------------------\n")
cat("1. Estimación más precisa de componentes de varianza\n")
cat("2. Manejo natural de datos desbalanceados\n")
cat("3. Predicciones BLUP para bloques y parcelas\n")
cat("4. Extensión directa a GLMM para datos no normales\n\n")

cat("ESPECIFICACIÓN DEL MODELO:\n")
cat("--------------------------\n")
cat("lmer(rendimiento ~ Riego * Trichoderma + (1|Bloque) + (1|Parcela_ID))\n\n")
cat("Donde:\n")
cat("  (1|Bloque)      = Efecto aleatorio del bloque\n")
cat("  (1|Parcela_ID)  = Efecto aleatorio de parcela principal (Error A)\n")
cat("  Residual        = Error de subparcela (Error B)\n\n")

# ----------------------------------------------------------------------------
# 5.1 AJUSTE DEL MODELO MIXTO
# ----------------------------------------------------------------------------

modelo_mixto <- lmer(
  rendimiento ~ Riego * Trichoderma + (1 | Bloque) + (1 | Parcela_ID),
  data = datos,
  REML = TRUE
)

cat("\n=== RESUMEN DEL MODELO MIXTO ===\n")
print(summary(modelo_mixto))

# ----------------------------------------------------------------------------
# 5.2 COMPONENTES DE VARIANZA
# ----------------------------------------------------------------------------

cat("\n5.2 COMPONENTES DE VARIANZA\n")
cat("-" , rep("-", 60), "\n", sep = "")

var_comp <- as.data.frame(VarCorr(modelo_mixto))
print(var_comp)

# Extraer varianzas
var_bloque <- var_comp$vcov[var_comp$grp == "Bloque"]
var_parcela <- var_comp$vcov[var_comp$grp == "Parcela_ID"]
var_residual <- var_comp$vcov[var_comp$grp == "Residual"]

cat("\nInterpretación:\n")
cat("  σ² Bloque:", round(var_bloque, 3), "(variabilidad entre bloques)\n")
cat("  σ² Parcela:", round(var_parcela, 3), "(Error A - main plot)\n")
cat("  σ² Residual:", round(var_residual, 3), "(Error B - subplot)\n")

# Calcular proporción de varianza
var_total <- var_bloque + var_parcela + var_residual
cat("\nProporción de varianza:\n")
cat("  Bloque:", round(var_bloque/var_total*100, 1), "%\n")
cat("  Parcela (Error A):", round(var_parcela/var_total*100, 1), "%\n")
cat("  Residual (Error B):", round(var_residual/var_total*100, 1), "%\n")

# ----------------------------------------------------------------------------
# 5.3 ANOVA TIPO III CON SATTERTHWAITE
# ----------------------------------------------------------------------------

cat("\n5.3 ANOVA TIPO III (lmerTest con Satterthwaite)\n")
cat("-" , rep("-", 60), "\n", sep = "")

anova_mixto <- Anova(modelo_mixto, type = "III", test.statistic = "F")
print(anova_mixto)

cat("\n=== COMPARACIÓN DE ENFOQUES ===\n")
cat("Los valores F y p del modelo mixto deberían ser similares\n")
cat("a los del aov() si el diseño es balanceado.\n")

# ============================================================================
# 6. COMPARACIONES MÚLTIPLES CON EMMEANS
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 6: COMPARACIONES MÚLTIPLES (EMMEANS)\n")
cat(rep("=", 80), "\n\n", sep = "")

# ----------------------------------------------------------------------------
# 6.1 EFECTO PRINCIPAL: RIEGO
# ----------------------------------------------------------------------------

cat("\n6.1 EFECTO DEL RIEGO (Factor A - Main Plot)\n")
cat("-" , rep("-", 60), "\n", sep = "")

emm_riego <- emmeans(modelo_mixto, ~ Riego)
cat("\nMedias marginales:\n")
print(emm_riego)

cat("\nComparaciones pairwise (Tukey):\n")
pairs(emm_riego, adjust = "tukey")

# ----------------------------------------------------------------------------
# 6.2 EFECTO PRINCIPAL: TRICHODERMA
# ----------------------------------------------------------------------------

cat("\n6.2 EFECTO DE TRICHODERMA (Factor B - Subplot)\n")
cat("-" , rep("-", 60), "\n", sep = "")

emm_trich <- emmeans(modelo_mixto, ~ Trichoderma)
cat("\nMedias marginales:\n")
print(emm_trich)

cat("\nComparaciones pairwise (Tukey):\n")
pairs(emm_trich, adjust = "tukey")

# ----------------------------------------------------------------------------
# 6.3 INTERACCIÓN RIEGO x TRICHODERMA
# ----------------------------------------------------------------------------

cat("\n6.3 INTERACCIÓN RIEGO x TRICHODERMA\n")
cat("-" , rep("-", 60), "\n", sep = "")

emm_inter <- emmeans(modelo_mixto, ~ Riego * Trichoderma)
cat("\nMedias de la interacción:\n")
print(emm_inter)

# Comparaciones simples: Efecto de Trichoderma DENTRO de cada nivel de Riego
cat("\nEfecto de Trichoderma DENTRO de cada Riego:\n")
for (r in riego) {
  cat("\n---", r, "---\n")
  emm_trich_riego <- emmeans(modelo_mixto, ~ Trichoderma | Riego, at = list(Riego = r))
  print(pairs(emm_trich_riego, adjust = "tukey"))
}

# Comparaciones simples: Efecto de Riego DENTRO de cada Trichoderma
cat("\nEfecto de Riego DENTRO de cada Trichoderma:\n")
for (t in trichoderma) {
  cat("\n---", t, "---\n")
  emm_riego_trich <- emmeans(modelo_mixto, ~ Riego | Trichoderma, at = list(Trichoderma = t))
  print(pairs(emm_riego_trich, adjust = "tukey"))
}

# Letras de agrupamiento
cat("\n=== LETRAS DE AGRUPAMIENTO (CLD) ===\n")
cld_result <- cld(emm_inter, alpha = 0.05, Letters = letters)
print(cld_result)

# ============================================================================
# 7. EVALUACIÓN DE SUPUESTOS
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 7: EVALUACIÓN DE SUPUESTOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Extraer residuos
residuos <- residuals(modelo_mixto)
ajustados <- fitted(modelo_mixto)

# ---- 7.1 Normalidad ----
cat("\n7.1 Normalidad de Residuos (Shapiro-Wilk)\n")
cat("-" , rep("-", 50), "\n", sep = "")

shapiro_test <- shapiro.test(residuos)
cat("W =", round(shapiro_test$statistic, 4), ", p =", round(shapiro_test$p.value, 4), "\n")
if (shapiro_test$p.value > 0.05) {
  cat("→ Normalidad CONFIRMADA ✓\n")
} else {
  cat("→ Desviación de normalidad ⚠\n")
}

# ---- 7.2 Homocedasticidad ----
cat("\n7.2 Homocedasticidad (gráfico)\n")
cat("-" , rep("-", 50), "\n", sep = "")

homog_plot <- ggplot(data.frame(ajustados = ajustados, residuos = residuos),
                     aes(x = ajustados, y = residuos)) +
  geom_point(alpha = 0.5, color = "#2E86AB", size = 3) +
  geom_hline(yintercept = 0, color = "#E74C3C", linetype = "dashed") +
  labs(title = "Residuos vs Ajustados",
       x = "Valores ajustados", y = "Residuos") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

print(homog_plot)
ggsave("02_Diseno_Experimental/03_splitplot_residuos.png",
       plot = homog_plot, width = 7, height = 5, dpi = 300)

# ---- 7.3 Diagnóstico integral ----
cat("\n7.3 Diagnóstico Integral (check_model)\n")
cat("-" , rep("-", 50), "\n", sep = "")

check_split <- check_model(modelo_mixto, panel = TRUE)
print(check_split)
ggsave("02_Diseno_Experimental/03_splitplot_diagnostico.png",
       plot = check_split, width = 10, height = 8, dpi = 300)

# ============================================================================
# 8. VISUALIZACIÓN DE RESULTADOS
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("PARTE 8: VISUALIZACIÓN DE RESULTADOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Datos para gráficos con medias y errores
datos_summary <- datos %>%
  group_by(Riego, Trichoderma) %>%
  summarise(
    media = mean(rendimiento),
    ee = sd(rendimiento) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Añadir letras de CLD
letras_df <- as.data.frame(cld_result) %>%
  select(Riego, Trichoderma, .groups)

datos_plot <- left_join(datos_summary, letras_df, by = c("Riego", "Trichoderma"))

# ----------------------------------------------------------------------------
# 8.1 GRÁFICO DE INTERACCIÓN
# ----------------------------------------------------------------------------

plot_interaccion <- ggplot(datos_plot, aes(x = Riego, y = media,
                                            color = Trichoderma, group = Trichoderma)) +
  geom_point(size = 4, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = media - ee, ymax = media + ee),
                width = 0.2, position = position_dodge(0.3)) +
  geom_line(position = position_dodge(0.3)) +
  labs(title = "Interacción Riego x Trichoderma",
       subtitle = "Rendimiento de tomate bajo estrés hídrico e inoculación",
       x = "Régimen de Riego",
       y = "Rendimiento (kg/parcela)") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(palette = "Set1", name = "Trichoderma")

print(plot_interaccion)
ggsave("02_Diseno_Experimental/03_splitplot_interaccion.png",
       plot = plot_interaccion, width = 9, height = 6, dpi = 300)

# ----------------------------------------------------------------------------
# 8.2 GRÁFICO DE BARRAS CON LETRAS
# ----------------------------------------------------------------------------

plot_barras <- ggplot(datos_plot, aes(x = Trichoderma, y = media, fill = Trichoderma)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(aes(ymin = media - ee, ymax = media + ee), width = 0.3) +
  geom_text(aes(label = .groups, y = media + ee + 0.5),
            size = 5, fontface = "bold") +
  facet_wrap(~ Riego, labeller = labeller(Riego = c(
    "Control" = "Riego Control",
    "Deficit_moderado" = "Déficit Moderado",
    "Deficit_severo" = "Déficit Severo"
  ))) +
  labs(title = "Efecto de Trichoderma por Régimen de Riego",
       subtitle = "Letras diferentes = diferencia significativa (Tukey, p < 0.05)",
       x = "Inoculación", y = "Rendimiento (kg/parcela)") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray80"),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  scale_fill_brewer(palette = "Greens")

print(plot_barras)
ggsave("02_Diseno_Experimental/03_splitplot_barras.png",
       plot = plot_barras, width = 10, height = 6, dpi = 300)

# ============================================================================
# 9. RESUMEN FINAL Y RECOMENDACIONES AGRONÓMICAS
# ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("RESUMEN DEL DISEÑO SPLIT-PLOT\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("DISEÑO EXPERIMENTAL:\n")
cat("  Tipo: Parcelas Divididas (Split-Plot) en DBCA\n")
cat("  Bloques:", n_bloques, "\n")
cat("  Factor A (Main Plot): Riego -", length(riego), "niveles\n")
cat("  Factor B (Subplot): Trichoderma -", length(trichoderma), "niveles\n")
cat("  Parcelas principales:", n_parcelas_principales, "\n")
cat("  Subparcelas totales:", n_subparcelas, "\n\n")

cat("COMPONENTES DE VARIANZA:\n")
cat("  σ² Bloque:", round(var_bloque, 3), "\n")
cat("  σ² Parcela (Error A):", round(var_parcela, 3), "(para probar Riego)\n")
cat("  σ² Residual (Error B):", round(var_residual, 3), "(para probar Trichoderma)\n\n")

cat("RESULTADOS DEL ANOVA:\n")
riego_p <- anova_mixto$`Pr(>F)`[1]
trich_p <- anova_mixto$`Pr(>F)`[2]
inter_p <- anova_mixto$`Pr(>F)`[3]

cat("  Riego: F =", round(anova_mixto$`F value`[1], 2),
    ", p <", ifelse(riego_p < 0.001, "0.001", round(riego_p, 4)), "\n")
cat("  Trichoderma: F =", round(anova_mixto$`F value`[2], 2),
    ", p <", ifelse(trich_p < 0.001, "0.001", round(trich_p, 4)), "\n")
cat("  Interacción: F =", round(anova_mixto$`F value`[3], 2),
    ", p <", ifelse(inter_p < 0.001, "0.001", round(inter_p, 4)), "\n\n")

cat("RECOMENDACIONES AGRONÓMICAS:\n")

# Encontrar el mejor tratamiento
mejor_trat <- datos_plot %>%
  arrange(desc(media)) %>%
  head(1)

cat("  Mejor combinación:", mejor_trat$Riego, "+", mejor_trat$Trichoderma, "\n")
cat("  Rendimiento esperado:", round(mejor_trat$media, 1), "kg/parcela\n\n")

if (inter_p < 0.05) {
  cat("  → Interacción SIGNIFICATIVA: el efecto de Trichoderma depende del riego\n")
  cat("  → Bajo déficit hídrico, la inoculación mitiga las pérdidas\n")
  cat("  → Recomendar T. asperellum para condiciones de estrés\n")
} else {
  cat("  → Interacción NO significativa: efectos son aditivos\n")
  cat("  → Se pueden recomendar tratamientos independientemente\n")
}

cat("\nPRECISIÓN EXPERIMENTAL:\n")
cv_error_b <- sqrt(var_residual) / mean(datos$rendimiento) * 100
cat("  CV (Error B):", round(cv_error_b, 2), "%\n")
if (cv_error_b < 10) {
  cat("  → Precisión EXCELENTE ✓\n")
} else if (cv_error_b < 20) {
  cat("  → Precisión ACEPTABLE ✓\n")
} else {
  cat("  → Precisión BAJA ⚠ considerar más repeticiones\n")
}

cat("\n", rep("=", 80), "\n", sep = "")
cat("FIN DEL ANÁLISIS\n")
cat(rep("=", 80), "\n", sep = "")
