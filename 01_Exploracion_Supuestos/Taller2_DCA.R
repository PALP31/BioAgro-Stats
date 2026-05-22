# ==============================================================================
# **Taller 2: Diseño Completamente al Azar (DCA) y Comparaciones Múltiples**
# Curso: Diseño Experimental y Aplicaciones Estadísticas
# Nivel: Universitario (Pregrado / Posgrado Inicial)
# ==============================================================================

# ==============================================================================
# 0. CONFIGURACIÓN DE ENTORNO Y CARGA SILENCIOSA DE PAQUETES
# ==============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))

preparar_entorno <- function(paquetes) {
  for (p in paquetes) {
    if (!require(p, character.only = TRUE)) {
      message(paste("Instalando paquete:", p, "... (esto puede tomar un momento)"))
      install.packages(p, dependencies = TRUE, quiet = TRUE)
      library(p, character.only = TRUE)
    }
  }
  message("¡Entorno R configurado con éxito!")
}

paquetes_requeridos <- c("ggplot2", "performance", "see", "agricolae", "DescTools", "car")
suppressMessages(preparar_entorno(paquetes_requeridos))


# ==============================================================================
# 1. CASO DE ESTUDIO PRÁCTICO: CONSERVACIÓN DE FRESAS
# ==============================================================================
# Contexto Agrícola/Alimentario:
# Evaluar la eficacia de 3 recubrimientos biodegradables sobre la vida de anaquel
# (días de conservación óptima) de fresas frente a un grupo Control.
# DCA balanceado con n=6 réplicas por tratamiento (N=24).

datos_fresas <- data.frame(
  Tratamiento = factor(rep(c("Control", "Almidon", "Gelatina", "Quitosano"), each = 6)),
  Dias = c(4.2, 5.1, 5.8, 4.5, 6.2, 4.2,   # Control (Media = 5.0)
           6.9, 7.6, 8.3, 7.0, 8.6, 7.2,   # Almidón (Media = 7.6)
           7.6, 8.9, 10.0, 8.2, 10.3, 9.0, # Gelatina (Media = 9.0)
           9.5, 11.2, 12.5, 10.1, 12.8, 11.1) # Quitosano (Media = 11.2)
)

# Mostrar resumen descriptivo
print("--- RESUMEN DE MEDIAS POR TRATAMIENTO ---")
print(aggregate(Dias ~ Tratamiento, data = datos_fresas, 
                function(x) c(Media = mean(x), DesvEst = sd(x))))


# ==============================================================================
# 2. VISUALIZACIÓN EXPLORATORIA DE TRES NIVELES
# ==============================================================================

# 2.1 Nivel Principiante: Boxplot nativo de R
boxplot(Dias ~ Tratamiento, data = datos_fresas, 
        main = 'Boxplot Principiante: Días por Tratamiento', 
        xlab = 'Tratamiento', ylab = 'Días de Conservación',
        col = 'lightblue', border = 'darkblue')

# 2.2 Nivel Intermedio: ggplot2 Boxplot Simple
ggplot(datos_fresas, aes(x = Tratamiento, y = Dias, fill = Tratamiento)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = 'Boxplot Intermedio: ggplot2 Simple', x = 'Tratamiento', y = 'Días') +
  theme(legend.position = "none")

# 2.3 Nivel Avanzado: ggplot2 Boxplot con Jitter Estilo Premium
ggplot(datos_fresas, aes(x = Tratamiento, y = Dias, fill = Tratamiento)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, color = "#2c3e50") +
  geom_jitter(width = 0.15, size = 3.5, aes(color = Tratamiento), alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Visualización Avanzada en DCA",
    subtitle = "Días de conservación de fresas según recubrimiento",
    x = "Película Biodegradable (Tratamiento)",
    y = "Días de Conservación Óptima"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(face = "italic", hjust = 0.5),
    legend.position = "none"
  )


# ==============================================================================
# 3. AJUSTE DE MODELO Y DIAGNÓSTICO DE SUPUESTOS
# ==============================================================================

# Ajustar modelo lineal general
modelo_fresas <- aov(Dias ~ Tratamiento, data = datos_fresas)

# Extraer residuos del modelo
residuos_fresas <- residuals(modelo_fresas)

# 3.1 Normalidad Visual (QQ-Plot)
qqnorm(residuos_fresas, main = "Gráfico QQ-Plot de los Residuos", col = "darkblue", pch = 19)
qqline(residuos_fresas, col = "red", lwd = 2)

# 3.2 Prueba de Normalidad de Shapiro-Wilk
# H0: Los residuos provienen de una distribución normal.
print("--- PRUEBA DE SHAPIRO-WILK ---")
print(shapiro.test(residuos_fresas))

# 3.3 Prueba de Homocedasticidad de Levene (Librería car)
# H0: Las varianzas son homogéneas entre tratamientos.
print("--- PRUEBA DE LEVENE ---")
print(leveneTest(Dias ~ Tratamiento, data = datos_fresas))

# 3.4 Diagnóstico Integral de Supuestos (Ecosistema performance)
check_model(modelo_fresas, check = c("normality", "homogeneity"))


# ==============================================================================
# 4. ANÁLISIS DE VARIANZA (ANOVA)
# ==============================================================================
print("--- TABLA ANOVA GLOBAL ---")
print(summary(modelo_fresas))


# ==============================================================================
# 5. COMPARACIONES MÚLTIPLES DE MEDIAS (POST-HOC)
# ==============================================================================

# 5.1 Prueba de Tukey HSD (Conservadora)
print("--- PRUEBA DE TUKEY (HSD) ---")
tukey_fresas <- HSD.test(modelo_fresas, "Tratamiento", group = TRUE)
print(tukey_fresas$groups)

# Graficar grupos con la visualización nativa de agricolae
plot(tukey_fresas, main = "Grupos de Comparación Rápida (agricolae)")

# 5.2 Gráfico Premium de Medias, EEM y Letras de Tukey
df_tukey <- data.frame(
  Tratamiento = rownames(tukey_fresas$means),
  Media = tukey_fresas$means$Dias,
  SD = tukey_fresas$means$std,
  Rep = tukey_fresas$means$r
)
df_tukey$EEM <- df_tukey$SD / sqrt(df_tukey$Rep)
df_tukey$Grupo <- tukey_fresas$groups[rownames(df_tukey), "groups"]

ggplot(df_tukey, aes(x = Tratamiento, y = Media, fill = Tratamiento)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7, width = 0.5) +
  geom_errorbar(aes(ymin = Media - EEM, ymax = Media + EEM), width = 0.15, size = 0.8) +
  geom_text(aes(label = Grupo, y = Media + EEM + 0.3), size = 6, fontface = "bold") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Comparación de Medias con Grupos de Tukey",
    subtitle = "Letras distintas indican diferencias significativas (p < 0.05). Barras representan el EEM.",
    x = "Tratamiento",
    y = "Media de Días de Conservación"
  ) +
  theme(legend.position = "none")

# 5.3 Métodos de comparación alternativos
print("--- PRUEBA DE DUNCAN ---")
print(duncan.test(modelo_fresas, "Tratamiento", group = TRUE)$groups)

print("--- PRUEBA LSD DE FISHER ---")
print(LSD.test(modelo_fresas, "Tratamiento", group = TRUE)$groups)

print("--- PRUEBA DE DUNNETT (VS. CONTROL) ---")
print(DunnettTest(Dias ~ Tratamiento, data = datos_fresas, control = "Control"))


# ==============================================================================
# ==============================================================================
# **TALLER EVALUATIVO (ESTUDIANTES)**
# ==============================================================================
# ==============================================================================

# ==============================================================================
# EJERCICIO 1: Inoculantes en Pino (Sintaxis Guiada)
# ==============================================================================
# Contexto Forestal:
# Evaluar el efecto de 3 inoculantes de micorrizas (M1, M2, M3) sobre la altura (cm)
# de plántulas de Pinus radiata frente a plántulas Control. (DCA, N=20).

datos_pino <- data.frame(
  Micorriza = factor(rep(c("Control", "M1", "M2", "M3"), each = 5)),
  Altura = c(12.5, 11.8, 13.1, 12.2, 11.9,  # Control
             15.2, 16.1, 14.8, 15.5, 15.9,  # M1
             18.5, 19.2, 17.8, 18.1, 18.9,  # M2
             14.1, 13.8, 14.5, 15.0, 14.2)  # M3
)

# [Guía de Completación]
# Completa los espacios vacíos reemplazando los ___ con los comandos adecuados:

# 1. Ajustar el modelo lineal en DCA
# modelo_pino <- aov(___ ~ ___, data = datos_pino)

# 2. Extraer residuos y probar Normalidad formalmente
# residuos_pino <- residuals(___)
# shapiro.test(___)

# 3. Probar Homocedasticidad formalmente con Levene
# leveneTest(___ ~ ___, data = datos_pino)

# 4. Desplegar e interpretar la tabla ANOVA
# summary(___)

# 5. Ejecutar la prueba de Tukey HSD
# tukey_pino <- HSD.test(___, "Micorriza", group = TRUE)
# print(tukey_pino$groups)


# ==============================================================================
# EJERCICIO 2: Toma de Decisiones en Raleo de Eucalipto
# ==============================================================================
# Contexto Forestal:
# Investigar el efecto de 3 intensidades de raleo (Control, Ligero, Fuerte) sobre
# el Diámetro a la Altura del Pecho (DAP en cm) de Eucalyptus globulus (DCA, N=18).
# Restricción comercial: Evitar a toda costa falsos positivos (Error Tipo I).

datos_eucalipto <- data.frame(
  Raleo = factor(rep(c("Control", "Ligero", "Fuerte"), each = 6)),
  DAP = c(15.2, 14.8, 15.6, 16.1, 14.9, 15.4,   # Control
          18.5, 19.1, 17.9, 18.7, 19.5, 18.2,   # Ligero
          22.1, 23.5, 21.8, 22.9, 24.0, 21.5)   # Fuerte
)

# [Escribe aquí tu análisis completo paso a paso siguiendo la Hoja de Ruta:]
# Paso 1: Visualización con ggplot2 (Boxplot + Jitter + Brewer Set2).
# Paso 2: Validación cuantitativa de supuestos (Shapiro y Levene).
# Paso 3: Análisis de Varianza (ANOVA).
# Paso 4: Selección de la prueba Post-Hoc adecuada (Justifica y ejecuta: LSD, Duncan o Tukey).
