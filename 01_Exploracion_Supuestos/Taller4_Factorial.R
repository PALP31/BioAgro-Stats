# ==============================================================================
# TALLER 4: ARREGLOS FACTORIALES EN DCA Y DBCA
# Curso: Aplicaciones Estadísticas en Agronomía
# ==============================================================================
#
# OBJETIVOS DE APRENDIZAJE:
#   1. Comprender qué es un Arreglo Factorial y cuándo aplicarlo.
#   2. Diferenciar entre un Arreglo Factorial en DCA y en DBCA.
#   3. Estudiar la pérdida de grados de libertad al bloquear (el "costo" estadístico).
#   4. Identificar la Unidad Experimental y trazar un croquis de campo.
#   5. Ajustar modelos factoriales con interacción en R y verificar sus supuestos.
#   6. Priorizar jerárquicamente la interpretación del término de la Interacción.
#   7. Realizar pruebas post-hoc de Tukey para la interacción.
#   8. Generar visualizaciones de calidad científica nativas y con ggplot2.
#   9. Desarrollar de manera autónoma un taller completo importando datos de Excel.
#
# ==============================================================================

# Evitar diálogos interactivos de selección de espejo CRAN
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Instalación y carga silenciosa de paquetes estadísticos
preparar_entorno <- function(paquetes) {
  for (p in paquetes) {
    if (!require(p, character.only = TRUE)) {
      message(paste("Instalando paquete:", p, "..."))
      install.packages(p, dependencies = TRUE, quiet = TRUE)
      library(p, character.only = TRUE)
    }
  }
  message("¡Entorno de R configurado con éxito!")
}

paquetes_requeridos <- c("ggplot2", "performance", "see", "agricolae", "DescTools", "readxl")
suppressMessages(preparar_entorno(paquetes_requeridos))


# ==============================================================================
# SECCIÓN 1: FUNDAMENTOS TEÓRICOS DE LOS DISEÑOS FACTORIALES
# ==============================================================================
#
# ── ¿QUÉ ES UN ARREGLO FACTORIAL? ─────────────────────────────────────────────
#
# Consiste en evaluar simultáneamente el efecto de dos o más factores independientes
# combinando todos sus niveles.
# Permite estudiar:
#   1. Efectos Principales: El impacto aislado de cada factor (Variedad, Riego).
#   2. Efecto de Interacción: Si el comportamiento de un factor cambia según el 
#      nivel del otro factor (Sinergismo o Antagonismo).
#
# ── ARREGLO FACTORIAL EN DCA vs. DBCA ─────────────────────────────────────────
#
# * Factorial en DCA:
#   Se usa si las unidades experimentales son homogéneas (ej. macetas uniformes).
#   Modelo: Y_ijk = µ + α_i + β_j + (αβ)_ij + ε_ijk
#
# * Factorial en DBCA:
#   Se usa si las unidades presentan variabilidad en una dirección (ej. pendiente, sombra).
#   Se crean r bloques y las combinaciones de factores se aleatorizan dentro de cada bloque.
#   Modelo: Y_ijk = µ + α_i + β_j + (αβ)_ij + γ_k + ε_ijk
#   Donde γ_k representa el efecto del Bloque k.
#
# ── LA PÉRDIDA DE GRADOS DE LIBERTAD: EL COSTO DEL BLOQUEO ────────────────────
#
# Al bloquear, "consumimos" grados de libertad del error experimental:
#
# +------------------------+-------------------+-------------------+
# | Fuente de Variación    | DCA (GL)          | DBCA (GL)         |
# +------------------------+-------------------+-------------------+
# | Factor A               | a - 1             | a - 1             |
# | Factor B               | b - 1             | b - 1             |
# | Interacción (A x B)    | (a - 1)(b - 1)    | (a - 1)(b - 1)    |
# | Bloques (r)            | —                 | r - 1             |
# | Error Residual         | ab(r - 1)         | (ab - 1)(r - 1)   |
# | Total                  | abr - 1           | abr - 1           |
# +------------------------+-------------------+-------------------+
#
# Si el terreno era en realidad homogéneo, bloquear innecesariamente reduce los GL
# del error, aumentando el Cuadrado Medio del Error y disminuyendo la potencia del
# experimento para detectar diferencias o interacciones reales.
#
# ── IDENTIFICACIÓN DE LA UNIDAD EXPERIMENTAL ──────────────────────────────────
#
# La Unidad Experimental es la menor porción de material a la que se le aplica una
# combinación de tratamientos de forma independiente y al azar.
# En agronomía, NO es la planta individual; es la parcela individual dentro de un
# invernadero/bloque que recibe la combinación de variedad de semilla y tipo de riego.
#
# ── CROQUIS DE DISTRIBUCIÓN EN CAMPO (FACTORIAL 2x2 EN DBCA) ──────────────────
#
# Evaluamos 2 Variedades (V1, V2) y 2 Riegos (R1: Goteo, R2: Aspersión) en 3 Bloques
# perpendiculares a la pendiente del terreno.
#
#                     [ DIRECCIÓN DE LA PENDIENTE: ZONA ALTA / SECA ]
# Norte  ========================================================================
#        BLOQUE 1 (Zona Alta)
#        +-------------------+-------------------+-------------------+-------------------+
#        |    V1R2 (T2)      |    V2R1 (T3)      |    V1R1 (T1)      |    V2R2 (T4)      |  <- Aleatorio
#        +-------------------+-------------------+-------------------+-------------------+
#        ========================================================================
#        BLOQUE 2 (Zona Media)
#        +-------------------+-------------------+-------------------+-------------------+
#        |    V2R1 (T3)      |    V1R1 (T1)      |    V2R2 (T4)      |    V1R2 (T2)      |  <- Aleatorio
#        +-------------------+-------------------+-------------------+-------------------+
#        ========================================================================
#        BLOQUE 3 (Zona Baja)
#        +-------------------+-------------------+-------------------+-------------------+
#        |    V1R1 (T1)      |    V2R2 (T4)      |    V1R2 (T2)      |    V2R1 (T3)      |  <- Aleatorio
#        +-------------------+-------------------+-------------------+-------------------+
# Sur    ========================================================================
#                     [ DIRECCIÓN DE LA PENDIENTE: ZONA BAJA / HÚMEDA ]
#
# ==============================================================================


# ==============================================================================
# SECCIÓN 2: IMPORTACIÓN DE DATOS DESDE EXCEL Y CONVERSIÓN DE FACTORES
# ==============================================================================

# Nombre del archivo de datos Excel
excel_file <- "datos_taller4_factorial.xlsx"

# Descargar desde GitHub si no existe (ej. en Google Colab o entorno sin archivos locales)
if (!file.exists(excel_file)) {
  message("Descargando datos_taller4_factorial.xlsx desde GitHub...")
  url_git <- "https://raw.githubusercontent.com/PALP31/BioAgro-Stats/main/01_Exploracion_Supuestos/datos_taller4_factorial.xlsx"
  tryCatch({
    download.file(url_git, destfile = excel_file, mode = "wb", quiet = TRUE)
    message("¡Archivo descargado exitosamente!")
  }, error = function(e) {
    message("⚠️ Error en descarga. Cargando datos de forma programática...")
  })
}

# Carga de datos
if (file.exists(excel_file)) {
  datos_fact <- read_excel(excel_file)
} else {
  # Fallback programático idéntico en caso de fallo de red
  datos_fact <- data.frame(
    Variedad = c("V1", "V1", "V2", "V2", "V1", "V1", "V2", "V2", "V1", "V1", "V2", "V2"),
    Riego = c("Goteo", "Aspersion", "Goteo", "Aspersion", "Goteo", "Aspersion", "Goteo", "Aspersion", "Goteo", "Aspersion", "Goteo", "Aspersion"),
    Bloque = c("Inv_1", "Inv_1", "Inv_1", "Inv_1", "Inv_2", "Inv_2", "Inv_2", "Inv_2", "Inv_3", "Inv_3", "Inv_3", "Inv_3"),
    Produccion = c(15, 10, 18, 12, 16, 11, 20, 13, 14, 9, 17, 11)
  )
  message("✅ Base de datos cargada vía fallback programático.")
}

# Conversión mandatoria a factores categóricos (vital para el ANOVA en R)
datos_fact$Variedad <- as.factor(datos_fact$Variedad)
datos_fact$Riego    <- as.factor(datos_fact$Riego)
datos_fact$Bloque   <- as.factor(datos_fact$Bloque)

cat("\n📊 ESTRUCTURA DE LA BASE DE DATOS IMPORTADA:\n")
str(datos_fact)

cat("\n📋 CONFIGURACIÓN DE LOS DATOS INICIALES:\n")
print(head(datos_fact, 12))


# ==============================================================================
# SECCIÓN 3: GRÁFICOS EXPLORATORIOS Y PERFILES DE INTERACCIÓN
# ==============================================================================

cat("\n📊 Generando gráfico exploratorio de perfiles de interacción...\n")

# Gráfico de perfil nativo en R
# Líneas no paralelas = indicio visual de Interacción significativa
interaction.plot(
  x.factor = datos_fact$Riego,
  trace.factor = datos_fact$Variedad,
  response = datos_fact$Produccion,
  type = "b",
  pch = c(19, 17),
  col = c("#2E86AB", "#E84855"),
  xlab = "Método de Riego",
  ylab = "Producción Promedio (t/ha)",
  legend = TRUE,
  main = "Gráfico de Perfil: Interacción Riego * Variedad"
)


# ==============================================================================
# SECCIÓN 4: AJUSTE DEL MODELO FACTORIAL DBCA Y TABLA ANOVA
# ==============================================================================

# Ajuste del modelo lineal con interacción (*) y bloques (+)
# En R, Variedad * Riego es equivalente a Variedad + Riego + Variedad:Riego
modelo_fact <- aov(Produccion ~ Variedad * Riego + Bloque, data = datos_fact)

cat("\n📊 ANÁLISIS DE VARIANZA (ANOVA) FACTORIAL EN DBCA:\n")
print(summary(modelo_fact))

# JERARQUÍA DE INTERPRETACIÓN:
# 1. Evaluar primero la Interacción Variedad:Riego.
# 2. Si p < 0.05 (Interacción significativa): El efecto del Riego depende de la Variedad.
#    NO debemos interpretar efectos principales de forma aislada.
# 3. Si p >= 0.05: Los factores actúan independientemente. Interpretar por separado.


# ==============================================================================
# SECCIÓN 5: VALIDACIÓN DE LOS SUPUESTOS DEL MODELO (RESIDUOS)
# ==============================================================================

cat("\n📊 Realizando diagnóstico gráfico de supuestos con performance...\n")
# 5.1 Diagnóstico gráfico de residuos automatizado
check_model(modelo_fact, check = c("normality", "homogeneity"))

# 5.2 Prueba Formal Cuantitativa de Normalidad: Shapiro-Wilk (H0: Residuos normales)
residuos_f <- residuals(modelo_fact)
cat("\n─── PRUEBA DE NORMALIDAD: SHAPIRO-WILK ─────────────────────────────────\n")
print(shapiro.test(residuos_f))

# 5.3 Prueba Formal de Homocedasticidad: Levene (H0: Varianzas iguales)
# En arreglos factoriales, se evalúa sobre los tratamientos combinados.
datos_fact$Tratamiento_Combinado <- interaction(datos_fact$Variedad, datos_fact$Riego)

cat("\n─── PRUEBA DE HOMOCEDASTICIDAD: LEVENE ─────────────────────────────────\n")
print(LeveneTest(Produccion ~ Tratamiento_Combinado, data = datos_fact))


# ==============================================================================
# SECCIÓN 6: COMPARACIONES MÚLTIPLES DE TUKEY POST-HOC
# ==============================================================================

cat("\n📊 Calculando comparaciones múltiples de Tukey HSD...\n")

# Tukey HSD para las combinaciones de factores utilizando agricolae
tukey_fact <- HSD.test(modelo_fact, c("Variedad", "Riego"), group = TRUE, console = TRUE)


# ==============================================================================
# SECCIÓN 7: VISUALIZACIÓN AVANZADA DE MEDIAS Y GRUPOS DE TUKEY
# ==============================================================================

cat("\n📊 Generando visualizaciones avanzadas de comparación de medias...\n")

# 7.1 Gráfica rápida nativa de la librería agricolae
plot(tukey_fact, variation = 'SE', col = 'skyblue', main = 'Tukey HSD: Combinación Variedad * Riego')

# 7.2 Visualización Premium con ggplot2, Medias, EEM y Letras de Tukey
# Construcción del dataframe para ggplot
df_tukey_f <- data.frame(
  Tratamiento = rownames(tukey_fact$means),
  Media = tukey_fact$means$Produccion,
  SD = tukey_fact$means$std,
  Rep = tukey_fact$means$r
)

# Cálculo de las barras de error (Error Estándar de la Media - EEM)
df_tukey_f$EEM <- df_tukey_f$SD / sqrt(df_tukey_f$Rep)
df_tukey_f$Grupo <- tukey_fact$groups[rownames(df_tukey_f), 'groups']

# Gráfico premium con ggplot2
g_premium <- ggplot(df_tukey_f, aes(x = Tratamiento, y = Media, fill = Tratamiento)) +
  geom_bar(stat = 'identity', color = 'black', alpha = 0.75, width = 0.5) +
  geom_errorbar(aes(ymin = Media - EEM, ymax = Media + EEM), width = 0.15, size = 0.8) +
  geom_text(aes(label = Grupo, y = Media + EEM + 0.4), size = 6, fontface = 'bold') +
  scale_fill_brewer(palette = 'Accent') +
  theme_minimal(base_size = 14) +
  labs(
    title = 'Rendimiento Promedio por Combinación (Variedad:Riego)',
    subtitle = 'Letras distintas indican diferencias significativas (Tukey HSD, p < 0.05). Barras representan EEM.',
    x = 'Tratamiento Combinado (Variedad : Riego)',
    y = 'Producción Promedio (t/ha)'
  ) +
  theme(legend.position = 'none')

print(g_premium)


# ==============================================================================
# SECCIÓN 8: CUESTIONARIO TEÓRICO DE ANÁLISIS
# ==============================================================================
#
# 1. ¿El término de interacción Variedad:Riego resultó estadísticamente significativo
#    en el ANOVA? ¿Qué implicaciones agronómicas y prácticas tiene para un agricultor?
#    *Respuesta:*
#
# 2. De acuerdo con la prueba de Tukey, ¿cuál es la mejor combinación de Variedad 
#    y Riego para maximizar la producción? ¿Hay diferencias estadísticas entre ellas?
#    *Respuesta:*
#
# 3. ¿La inclusión del factor 'Bloque' (Invernaderos) fue efectiva en este diseño?
#    Justifica tu respuesta observando el p-valor de Bloque en la tabla ANOVA.
#    *Respuesta:*
#
# 4. ¿Qué supuestos matemáticos de los residuos fueron validados cuantitativamente
#    y qué nos indican sus respectivos p-valores obtenidos en Shapiro-Wilk y Levene?
#    *Respuesta:*


# ==============================================================================
# SECCIÓN 9: 📝 ACTIVIDAD PRÁCTICA AUTÓNOMA (Taller de Aplicación)
# ==============================================================================
#
# CONTEXTO DE LA TAREA:
# Un agrónomo desea evaluar cómo interactúan 2 niveles del factor Semilla (S1, S2)
# y 3 niveles del factor Espaciamiento (10cm, 20cm, 30cm) sobre la variable de
# respuesta Altura de la planta en centímetros (Altura_cm).
#
# DISEÑO:
# El ensayo de campo tiene variabilidad de textura de suelo (gradiente sistemático).
# Se agrupan en 3 Bloques (Franco, Arcilloso, Arenoso). Este diseño es un:
# Arreglo Factorial 2x3 en Diseño de Bloques Completos al Azar (DBCA).
#
# DATOS: tarea_taller4_factorial.xlsx
#
# OBJETIVOS:
#   1. Importar el archivo usando read_excel().
#   2. Convertir variables categóricas a factores.
#   3. Construir un perfil de interacción preliminar.
#   4. Ajustar el ANOVA factorial con bloques e interpretar significancias.
#   5. Validar supuestos (diagnóstico gráfico, Shapiro-Wilk, Levene).
#   6. Ejecutar la prueba Tukey HSD de comparaciones múltiples de la interacción.
#   7. Dibujar gráficos Tukey (nativo y premium ggplot2).
#   8. Resolver el cuestionario agronómico final.

# ── Paso 1: Importación del Excel de la Tarea ──────────────────────────────────
excel_tarea <- "tarea_taller4_factorial.xlsx"

if (!file.exists(excel_tarea)) {
  message("Descargando tarea_taller4_factorial.xlsx desde GitHub...")
  url_tarea_git <- "https://raw.githubusercontent.com/PALP31/BioAgro-Stats/main/01_Exploracion_Supuestos/tarea_taller4_factorial.xlsx"
  tryCatch({
    download.file(url_tarea_git, destfile = excel_tarea, mode = "wb", quiet = TRUE)
    message("¡Archivo de la tarea descargado exitosamente!")
  }, error = function(e) {
    message("⚠️ Error en descarga. Generando archivo de tarea local alternativo (CSV) como contingencia...")
    df_contingencia <- data.frame(
      Semilla = c("S1", "S1", "S1", "S2", "S2", "S2", "S1", "S1", "S1", "S2", "S2", "S2", "S1", "S1", "S1", "S2", "S2", "S2"),
      Espaciamiento = c("10cm", "20cm", "30cm", "10cm", "20cm", "30cm", "10cm", "20cm", "30cm", "10cm", "20cm", "30cm", "10cm", "20cm", "30cm", "10cm", "20cm", "30cm"),
      Textura_Suelo = c("Franco", "Franco", "Franco", "Franco", "Franco", "Franco", "Arcilloso", "Arcilloso", "Arcilloso", "Arcilloso", "Arcilloso", "Arcilloso", "Arenoso", "Arenoso", "Arenoso", "Arenoso", "Arenoso", "Arenoso"),
      Altura_cm = c(10.5, 12.1, 14.5, 11.2, 13.0, 15.2, 9.8, 11.5, 13.8, 10.5, 12.3, 14.5, 8.5, 10.2, 12.5, 9.2, 11.0, 13.1)
    )
    write.csv(df_contingencia, "tarea_taller4_factorial.csv", row.names = FALSE)
  })
}

# [Escribe tu código para importar el archivo Excel de la tarea aquí]
# datos_tarea4 <- read_excel(excel_tarea)
# (Si usas el fallback offline: datos_tarea4 <- read.csv("tarea_taller4_factorial.csv"))
datos_tarea4 <- ___


# ── Paso 2: Conversión obligatoria a factores categóricos ──────────────────────
# [Escribe tu código para la conversión obligatoria a factores aquí]
datos_tarea4$Semilla       <- ___(datos_tarea4$Semilla)
datos_tarea4$Espaciamiento <- ___(datos_tarea4$Espaciamiento)
datos_tarea4$Textura_Suelo <- ___(datos_tarea4$Textura_Suelo)

# Comprobar estructura
str(datos_tarea4)


# ── Paso 3: Perfil de Interacción Visual ────────────────────────────────────────
# [Escribe tu código para construir el gráfico de perfiles de interacción aquí]
interaction.plot(
  x.factor = datos_tarea4$___,
  trace.factor = datos_tarea4$___,
  response = datos_tarea4$___,
  type = "b",
  pch = c(19, 17),
  col = c("blue", "red"),
  xlab = "Espaciamiento",
  ylab = "Altura Promedio (cm)",
  legend = TRUE,
  main = "Gráfico de Perfil: Interacción Semilla * Espaciamiento"
)


# ── Paso 4: Ajuste del Modelo Factorial en DBCA y ANOVA ───────────────────────
# [Escribe tu código para ajustar el modelo y desplegar la tabla ANOVA aquí]
modelo_tarea4 <- aov(Altura_cm ~ Semilla ___ Espaciamiento + Textura_Suelo, data = datos_tarea4)
summary(modelo_tarea4)


# ── Paso 5: Validación de los Supuestos del Modelo (Residuos) ─────────────────
# 5.1 Diagnóstico Gráfico Integral
# [Escribe tu código para el diagnóstico gráfico automatizado de residuos aquí]
check_model(modelo_tarea4, check = c("normality", "homogeneity"))

# 5.2 Prueba Formal de Normalidad (Shapiro-Wilk)
# [Escribe tu código para evaluar la normalidad mediante Shapiro-Wilk aquí]
residuos_tarea4 <- residuals(___)
shapiro.test(residuos_tarea4)

# 5.3 Prueba Formal de Homocedasticidad (Levene por interacción)
# [Escribe tu código para evaluar la homocedasticidad mediante LeveneTest aquí]
datos_tarea4$Tratamiento_Combinado <- interaction(datos_tarea4$Semilla, datos_tarea4$Espaciamiento)
LeveneTest(Altura_cm ~ Tratamiento_Combinado, data = ___)


# ── Paso 6: Comparación de Medias mediante Tukey HSD ──────────────────────────
# [Escribe tu código para ejecutar el test Tukey HSD aquí]
tukey_tarea4 <- HSD.test(modelo_tarea4, c("Semilla", "Espaciamiento"), group = TRUE, console = TRUE)


# ── Paso 7: Visualización Avanzada de Medias y Letras de Tukey ──────────────────
# 7.1 Gráfica Básica (Nativa)
# [Escribe tu código para generar la gráfica nativa de Tukey aquí]
plot(tukey_tarea4, variation = 'SE', col = 'lightgreen', main = 'Tukey HSD: Tarea Semilla * Espaciamiento')

# 7.2 Gráfica Premium con ggplot2
# [Escribe tu código para elaborar la visualización premium de ggplot2 aquí]
df_tukey_t4 <- data.frame(
  Tratamiento = rownames(tukey_tarea4$means),
  Media = tukey_tarea4$means$Altura_cm,
  SD = tukey_tarea4$means$std,
  Rep = tukey_tarea4$means$r
)

# Cálculo de las barras de error (Error Estándar de la Media)
df_tukey_t4$EEM <- df_tukey_t4$SD / sqrt(df_tukey_t4$Rep)
df_tukey_t4$Grupo <- tukey_tarea4$groups[rownames(df_tukey_t4), 'groups']

# Gráfico premium con ggplot2
ggplot(df_tukey_t4, aes(x = Tratamiento, y = Media, fill = Tratamiento)) +
  geom_bar(stat = 'identity', color = 'black', alpha = 0.75, width = 0.5) +
  geom_errorbar(aes(ymin = Media - EEM, ymax = Media + EEM), width = 0.15, size = 0.8) +
  geom_text(aes(label = Grupo, y = Media + EEM + 0.3), size = 6, fontface = 'bold') +
  scale_fill_brewer(palette = 'Spectral') +
  theme_minimal(base_size = 14) +
  labs(
    title = 'Altura Promedio por Combinación de Tratamientos (Semilla:Espaciamiento)',
    subtitle = 'Letras distintas indican diferencias significativas (Tukey HSD, p < 0.05). Barras representan EEM.',
    x = 'Tratamiento Combinado (Semilla : Espaciamiento)',
    y = 'Altura Promedio (cm)'
  ) +
  theme(legend.position = 'none')


# ── Paso 8: Cuestionario Final de Conclusiones Técnicas ────────────────────────
#
# 1. ¿Existe una interacción estadísticamente significativa entre la Semilla y el
#    Espaciamiento en la altura de la planta (p < 0.05)? ¿Qué interpretación le das?
#    *Respuesta:*
#
# 2. De acuerdo al agrupamiento del test de Tukey, ¿cuál o cuáles combinaciones
#    de Semilla y Espaciamiento producen plantas significativamente más altas?
#    *Respuesta:*
#
# 3. ¿La textura del suelo (Bloques) resultó significativa en el ANOVA? ¿Fue adecuado
#    bloquear o se habría obtenido más potencia con un DCA simple?
#    *Respuesta:*
#
# 4. ¿Los supuestos de normalidad y homocedasticidad se cumplieron de manera
#    satisfactoria? Argumenta citando los p-valores obtenidos.
#    *Respuesta:*