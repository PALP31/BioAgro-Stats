import json

# Define the cells of the updated Jupyter Notebook
cells = [
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "# **Taller 2: Diseño Completamente al Azar (DCA) y Comparaciones Múltiples**\n",
            "\n",
            "**Curso:** Diseño Experimental y Aplicaciones Estadísticas  \n",
            "**Nivel:** Universitario (Pregrado / Posgrado Inicial)  \n",
            "**Instructor:** Experto en Diseño Experimental y Ciencia de Datos  \n",
            "\n",
            "---"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 🚀 **Instrucciones para Google Colab (R Kernel)**\n",
            "\n",
            "Si estás ejecutando este taller desde **Google Colab**, ten en cuenta lo siguiente:\n",
            "\n",
            "1. **Entorno en la Nube:** Colab ejecuta el cuaderno en un servidor temporal. Cada vez que abras el cuaderno, necesitarás instalar los paquetes de R requeridos.\n",
            "2. **Configuración Automática:** La celda de abajo instalará de forma silenciosa y automática todas las librerías necesarias (`ggplot2`, `performance`, `see`, `agricolae`, `DescTools`, `car`) sin inundar la pantalla con registros de descarga.\n",
            "3. **Sin Carga de Archivos Externos:** Todos los datos de este taller se generan directamente en el código R. No necesitas subir archivos `.csv` o `.xlsx` adicionales. ¡Todo está listo para ejecutarse!\n",
            "\n",
            "*(Para ejecutar una celda, haz clic en el botón de reproducción ▶️ a la izquierda de la celda o presiona `Shift + Enter`)*"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# 0. CONFIGURACIÓN DE ENTORNO EN GOOGLE COLAB / LOCAL (LIBRERÍAS)\n",
            "# ==============================================================================\n",
            "\n",
            "# Evitar diálogos interactivos de selección de espejo CRAN\n",
            "options(repos = c(CRAN = \"https://cloud.r-project.org\"))\n",
            "\n",
            "# Instalación y carga silenciosa de paquetes estadísticos\n",
            "preparar_entorno <- function(paquetes) {\n",
            "  for (p in paquetes) {\n",
            "    if (!require(p, character.only = TRUE)) {\n",
            "      message(paste(\"Instalando paquete:\", p, \"...\"))\n",
            "      install.packages(p, dependencies = TRUE, quiet = TRUE)\n",
            "      library(p, character.only = TRUE)\n",
            "    }\n",
            "  }\n",
            "  message(\"¡Entorno de R configurado con éxito!\")\n",
            "}\n",
            "\n",
            "paquetes_requeridos <- c(\"ggplot2\", \"performance\", \"see\", \"agricolae\", \"DescTools\", \"car\")\n",
            "suppressMessages(preparar_entorno(paquetes_requeridos))\n"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 1. Fundamentos Teóricos del DCA\n",
            "\n",
            "### ¿Por qué y cuándo utilizar un DCA?\n",
            "El **Diseño Completamente al Azar (DCA)** es el diseño experimental más simple, flexible y de más amplia aplicación cuando las condiciones externas son ideales. Se fundamenta en la asignación aleatoria de los tratamientos a las unidades experimentales, sin ningún tipo de restricción de bloqueo.\n",
            "\n",
            "**¿Cuándo se utiliza?**\n",
            "* Cuando las **unidades experimentales son altamente homogéneas** (ej. macetas de la misma variedad en un invernadero climatizado, placas de Petri en una incubadora, animales de la misma camada, edad y peso, o lotes uniformes de suelo en laboratorio).\n",
            "* Cuando el número de tratamientos o réplicas es reducido, ya que maximiza los grados de libertad del error.\n",
            "\n",
            "### Comparación conceptual: DCA vs. DBCA (Bloques)\n",
            "* **DCA (Sin restricciones):** Asume que la única fuente sistemática de variación son los **Tratamientos**. Cualquier otra variabilidad se agrupa bajo el término del error aleatorio.  \n",
            "  $$\\text{Modelo matemático del DCA:} \\quad Y_{ij} = \\mu + \\tau_i + \\varepsilon_{ij}$$\n",
            "  Donde $Y_{ij}$ es la observación en la unidad experimental $j$ del tratamiento $i$, $\\mu$ es la media general, $\\tau_i$ es el efecto del tratamiento $i$, y $\\varepsilon_{ij}$ es el error experimental aleatorio, que debe cumplir con $\\varepsilon_{ij} \\sim N(0, \\sigma^2)$ independientes.\n",
            "\n",
            "* **DBCA (Con restricción de aleatorización):** Si se detecta o sospecha la existencia de un gradiente de variabilidad no deseado pero predecible (por ejemplo, diferencias de luz, humedad o fertilidad en un terreno), la homogeneidad se rompe. Para controlarlo, agrupamos las unidades en **Bloques** homogéneos. La aleatorización se realiza *dentro* de cada bloque.  \n",
            "  $$\\text{Modelo matemático del DBCA:} \\quad Y_{ij} = \\mu + \\tau_i + \\beta_j + \\varepsilon_{ij}$$\n",
            "  Donde $\\beta_j$ es el efecto del bloque $j$. Al aislar la variación de los bloques, el error residual $\\varepsilon_{ij}$ disminuye sustancialmente, aumentando la precisión y potencia de la prueba ANOVA."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 2. Anatomía de la Tabla ANOVA para DCA\n",
            "\n",
            "El Análisis de Varianza (ANOVA) descompone la variabilidad total de los datos en dos fuentes principales: la variabilidad **Entre Tratamientos** (el efecto del factor estudiado) y la variabilidad **Dentro de Tratamientos** (el Error Residual o ruido aleatorio).\n",
            "\n",
            "| Fuente de Variación (SV) | Grados de Libertad (GL) | Suma de Cuadrados (SC) | Cuadrado Medio (CM) | $F$ Calculado ($F_{\\text{calc}}$) | Valor $p$ ($p\\text{-value}$) |\n",
            "| :--- | :--- | :--- | :--- | :--- | :--- |\n",
            "| **Tratamientos (Entre)** | $k - 1$ | $SC_{\\text{Trat}}$ | $CM_{\\text{Trat}} = \\frac{SC_{\\text{Trat}}}{k-1}$ | $F = \\frac{CM_{\\text{Trat}}}{CM_{\\text{Error}}}$ | $P(F_{(k-1, N-k)} \\ge F)$ |\n",
            "| **Error (Dentro)** | $N - k$ | $SC_{\\text{Error}}$ | $CM_{\\text{Error}} = \\frac{SC_{\\text{Error}}}{N-k}$ | | |\n",
            "| **Total** | $N - 1$ | $SC_{\\text{Total}}$ | | | |\n",
            "\n",
            "### ¿Cómo interpretar cada columna? (Conceptos Clave):\n",
            "* **Grados de Libertad (GL):** Es el número de piezas de información independiente. Para tratamientos es $k-1$ (donde $k$ es el número de grupos). Para el error es $N-k$ (donde $N$ es el número total de unidades).\n",
            "* **Suma de Cuadrados (SC):** Representa la dispersión o variabilidad. La $SC_{\\text{Trat}}$ mide qué tan lejos están las medias de los tratamientos respecto a la media general. La $SC_{\\text{Error}}$ mide la desviación de cada punto individual respecto a la media de su propio grupo.\n",
            "* **Cuadrados Medios (CM):** Son varianzas estimadas. Se obtienen dividiendo las sumas de cuadrados entre sus respectivos grados de libertad.\n",
            "* **$F$ Calculado ($F_{\\text{calc}}$):** Es el cociente entre el $CM_{\\text{Trat}}$ y el $CM_{\\text{Error}}$. Si es cercano a 1, la variabilidad entre tratamientos es similar al ruido de fondo. Si es significativamente mayor que 1, indica que los tratamientos causan cambios reales.\n",
            "* **Valor $p$:** La probabilidad de que las diferencias observadas se deban puramente al azar. Si $p < 0.05$, rechazamos la hipótesis nula ($H_0$: $\\mu_1 = \\mu_2 = ... = \\mu_k$)."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 3. Teoría de Comparaciones Múltiples (Post-Hoc)\n",
            "\n",
            "El ANOVA es una prueba *ómnibus*: nos dice si hay alguna diferencia significativa en el modelo global, pero no especifica qué grupos específicos difieren entre sí. Si rechazamos $H_0$, debemos realizar pruebas post-hoc:\n",
            "\n",
            "1. **Tukey (HSD - Honest Significant Difference):**\n",
            "   * **Enfoque:** Conservador.\n",
            "   * **Rigor:** Controla estrictamente la tasa de error por familia (Family-wise Error Rate, FWER), asegurando que la probabilidad de cometer al menos un Error Tipo I (falso positivo) en todo el conjunto de comparaciones simultáneas sea exactamente del $5\\%$. Excelente para comparaciones del tipo \"todos contra todos\".\n",
            "\n",
            "2. **LSD de Fisher (Least Significant Difference):**\n",
            "   * **Enfoque:** Liberal.\n",
            "   * **Rigor:** No controla el error acumulado por comparaciones múltiples. La probabilidad de cometer falsos positivos crece exponencialmente al aumentar el número de tratamientos. Solo debe emplearse si el ANOVA global resulta altamente significativo y el número de tratamientos es muy bajo (3 o 4).\n",
            "\n",
            "3. **Dunnett:**\n",
            "   * **Enfoque:** Altamente específico.\n",
            "   * **Rigor:** Compara todos los tratamientos activos **exclusivamente contra un tratamiento control o testigo**. Al reducir significativamente el número de comparaciones de $k(k-1)/2$ a solo $k-1$, ofrece un poder estadístico excelente para detectar diferencias frente a un testigo.\n",
            "\n",
            "4. **Duncan:**\n",
            "   * **Enfoque:** Rango Múltiple (Intermedio/Liberal).\n",
            "   * **Rigor:** Ajusta el nivel crítico según la distancia de rango entre las medias comparadas. Declara diferencias con mayor facilidad que Tukey, pero tiene un riesgo de falso positivo moderadamente alto."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 4. Caso de Estudio Práctico: Conservación de Fresas\n",
            "\n",
            "**Contexto Agrícola/Alimentario:**\n",
            "Un grupo de investigadores agrícolas desea probar la eficacia de 3 recubrimientos biodegradables sobre la vida de anaquel (en días de conservación óptima) de fresas (*Fragaria ananassa*) cosechadas en condiciones uniformes.\n",
            "\n",
            "* **Tratamientos (Factor):** Tipo de Recubrimiento (Control, Almidón, Gelatina, Quitosano).\n",
            "* **Variable de Respuesta:** Días de conservación sin daño fúngico o ablandamiento térmico.\n",
            "* **Diseño:** DCA balanceado con $n=6$ réplicas por tratamiento ($N=24$ unidades experimentales en condiciones controladas de laboratorio)."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# 2. CARGA Y ANÁLISIS DESCRIPTIVO DE DATOS\n",
            "# ==============================================================================\n",
            "\n",
            "# Creación del data.frame con valores decimales de alta sensibilidad\n",
            "datos_fresas <- data.frame(\n",
            "  Tratamiento = factor(rep(c(\"Control\", \"Almidon\", \"Gelatina\", \"Quitosano\"), each = 6)),\n",
            "  Dias = c(4.2, 5.1, 5.8, 4.5, 6.2, 4.2,   # Control (Media = 5.0)\n",
            "           6.9, 7.6, 8.3, 7.0, 8.6, 7.2,   # Almidón (Media = 7.6)\n",
            "           7.6, 8.9, 10.0, 8.2, 10.3, 9.0, # Gelatina (Media = 9.0)\n",
            "           9.5, 11.2, 12.5, 10.1, 12.8, 11.1) # Quitosano (Media = 11.2)\n",
            ")\n",
            "\n",
            "# Inspección preliminar de datos\n",
            "print(datos_fresas)\n",
            "\n",
            "# Estadísticas descriptivas de los días de conservación por tratamiento\n",
            "resumen_fresas <- aggregate(Dias ~ Tratamiento, data = datos_fresas, \n",
            "                            function(x) c(Media = mean(x), DesvEst = sd(x)))\n",
            "print(resumen_fresas)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### ⚠️ **La Trampa del Factor (The Factor Trap)**\n",
            "\n",
            "En R, cuando analizamos variables categóricas (como nuestros tratamientos `Control`, `Almidón`, etc.), siempre debemos declararlas explícitamente como **factores** usando la función `factor()`.  \n",
            "\n",
            "Si no lo haces (por ejemplo, si tus tratamientos se llaman `1`, `2`, `3` y los dejas cargados como números simples), R no entenderá que son categorías cualitativas independientes. En su lugar, el comando `aov()` los tratará como una variable continua e intentará hacer una regresión lineal (calculando un único grado de libertad), lo que producirá un **ANOVA completamente erróneo**. ¡Siempre verifica con `is.factor()`!"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### Visualización Exploratoria de Datos (Tres Niveles)\n",
            "\n",
            "Para adaptarnos a estudiantes de todos los niveles, implementamos tres formas progresivas de graficar los datos:\n",
            "1. **Nivel Principiante (Gráfico Base de R):** La forma más directa y simple sin necesidad de instalar librerías adicionales.\n",
            "2. **Nivel Intermedio (ggplot2 Boxplot Simple):** Una introducción a la potencia de ggplot2 y el sistema de capas.\n",
            "3. **Nivel Avanzado (Boxplot + Jitter + Estética Premium):** Para crear reportes científicos profesionales, superponiendo los datos individuales reales."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# 3. VISUALIZACIÓN EXPLORATORIA DE TRES NIVELES\n",
            "# ==============================================================================\n",
            "\n",
            "# 3.1 Nivel Principiante: Boxplot nativo de R\n",
            "boxplot(Dias ~ Tratamiento, data = datos_fresas, \n",
            "        main = 'Boxplot Principiante: Días por Tratamiento', \n",
            "        xlab = 'Tratamiento', ylab = 'Días de Conservación',\n",
            "        col = 'lightblue', border = 'darkblue')\n",
            "\n",
            "# 3.2 Nivel Intermedio: ggplot2 Boxplot Simple\n",
            "ggplot(datos_fresas, aes(x = Tratamiento, y = Dias, fill = Tratamiento)) +\n",
            "  geom_boxplot() +\n",
            "  theme_minimal() +\n",
            "  labs(title = 'Boxplot Intermedio: ggplot2 Simple', x = 'Tratamiento', y = 'Días') +\n",
            "  theme(legend.position = \"none\")\n",
            "\n",
            "# 3.3 Nivel Avanzado: ggplot2 Boxplot con Jitter Estilo Premium\n",
            "ggplot(datos_fresas, aes(x = Tratamiento, y = Dias, fill = Tratamiento)) +\n",
            "  geom_boxplot(alpha = 0.4, outlier.shape = NA, color = \"#2c3e50\") +\n",
            "  geom_jitter(width = 0.15, size = 3.5, aes(color = Tratamiento), alpha = 0.8) +\n",
            "  scale_fill_brewer(palette = \"Set2\") +\n",
            "  scale_color_brewer(palette = \"Set2\") +\n",
            "  theme_minimal(base_size = 14) +\n",
            "  labs(\n",
            "    title = \"Visualización Avanzada en DCA\",\n",
            "    subtitle = \"Días de conservación de fresas según recubrimiento\",\n",
            "    x = \"Película Biodegradable (Tratamiento)\",\n",
            "    y = \"Días de Conservación Óptima\"\n",
            "  ) +\n",
            "  theme(\n",
            "    plot.title = element_text(face = \"bold\", hjust = 0.5),\n",
            "    plot.subtitle = element_text(face = \"italic\", hjust = 0.5),\n",
            "    legend.position = \"none\"\n",
            "  )"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### 🔍 **¿Qué son los Residuos? (Concepto Clave)**\n",
            "\n",
            "Cuando ajustamos un modelo estadístico, cada observación real ($Y_{ij}$) se compone de una parte explicada por el tratamiento y una parte que no podemos explicar (el ruido aleatorio o error experimental). A este ruido lo llamamos **Residuo** ($\\varepsilon_{ij}$):\n",
            "\n",
            "$$\\text{Residuo} = \\text{Valor Observado} - \\text{Valor Predicho (Media de su Grupo)}$$\n",
            "\n",
            "Para que las conclusiones de nuestro ANOVA sean matemáticamente válidas y confiables, este ruido experimental debe comportarse de forma \"civilizada\":\n",
            "* Debe distribuirse de forma simétrica (como una campana de Gauss) -> **Supuesto de Normalidad**.\n",
            "* Su variabilidad (el ancho de su dispersión) debe ser constante en todos los tratamientos -> **Supuesto de Homocedasticidad**.\n",
            "* Las mediciones deben ser independientes unas de otras -> **Supuesto de Independencia**."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### Validación de Supuestos del Modelo\n",
            "\n",
            "Ajustamos el modelo lineal con la función `aov()` de R. Antes de interpretar el ANOVA, debemos diagnosticar si los residuos cumplen estrictamente con los supuestos:\n",
            "\n",
            "1. **Normalidad (QQ-plot y Shapiro-Wilk):** Evaluamos visualmente con el gráfico Cuantil-Cuantil (QQ-Plot), donde los puntos deben ajustarse a la línea diagonal de referencia. Validamos formalmente mediante la prueba cuantitativa de *Shapiro-Wilk* (buscamos $p > 0.05$).\n",
            "2. **Homocedasticidad (Prueba de Levene / Bartlett):** Evaluamos la igualdad de varianzas usando la robusta prueba de *Levene* (de la librería `car`) o *Bartlett*. Si $p > 0.05$, confirmamos que no hay diferencias significativas de varianza entre grupos.\n",
            "3. **Paquete `performance`:** Utilizaremos la función `check_model` para obtener un diagnóstico automatizado completo."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# 4. AJUSTE DE MODELO Y DIAGNÓSTICO DE SUPUESTOS COMPLETO\n",
            "# ==============================================================================\n",
            "\n",
            "# Ajuste del modelo\n",
            "modelo_fresas <- aov(Dias ~ Tratamiento, data = datos_fresas)\n",
            "\n",
            "# Extracción formal de residuos\n",
            "residuos_fresas <- residuals(modelo_fresas)\n",
            "\n",
            "# 4.1 QQ-Plot de los Residuos (Normalidad Visual)\n",
            "qqnorm(residuos_fresas, main = \"Gráfico QQ-Plot de los Residuos\", col = \"darkblue\", pch = 19)\n",
            "qqline(residuos_fresas, col = \"red\", lwd = 2)\n",
            "\n",
            "# 4.2 Prueba de Normalidad de Shapiro-Wilk\n",
            "# H0: Los residuos provienen de una distribución normal.\n",
            "prueba_shapiro <- shapiro.test(residuos_fresas)\n",
            "print(prueba_shapiro)\n",
            "\n",
            "# 4.3 Prueba de Homocedasticidad de Levene (Librería car)\n",
            "# H0: Las varianzas son homogéneas entre los tratamientos.\n",
            "prueba_levene <- leveneTest(Dias ~ Tratamiento, data = datos_fresas)\n",
            "print(prueba_levene)\n",
            "\n",
            "# 4.4 Diagnóstico automatizado integral con 'performance'\n",
            "check_model(modelo_fresas, check = c(\"normality\", \"homogeneity\"))"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### Interpretación del ANOVA\n",
            "Una vez confirmados todos los supuestos (QQ-Plot alineado, Levene y Shapiro con $p > 0.05$), desplegamos la tabla ANOVA con la función `summary()` para verificar si el efecto global es estadísticamente significativo."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# 5. ANÁLISIS DE VARIANZA (ANOVA)\n",
            "# ==============================================================================\n",
            "\n",
            "summary(modelo_fresas)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 5. Comparaciones Múltiples de Medias (Post-Hoc) y Graficación de Letras\n",
            "\n",
            "Cuando el ANOVA indica diferencias significativas, procedemos a realizar las comparaciones múltiples de medias. \n",
            "\n",
            "### Visualizar las Letras de Grupos de Comparación\n",
            "Para que los estudiantes entiendan visualmente qué significan las letras, la mejor práctica en ciencia de datos es graficar la media de cada grupo con sus respectivas **barras de error estándar** e imprimir la letra de los grupos del test directamente sobre las barras. Esto permite observar el solapamiento de los rangos estadísticos de un solo vistazo."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# 6. APLICACIÓN DE PRUEBAS POST-HOC Y VISUALIZACIÓN DE COMPARACIONES\n",
            "# ==============================================================================\n",
            "\n",
            "# 6.1 Prueba de Tukey HSD (Conservadora, estándar de oro)\n",
            "print(\"--- PRUEBA DE TUKEY (HSD) ---\")\n",
            "tukey_fresas <- HSD.test(modelo_fresas, \"Tratamiento\", group = TRUE)\n",
            "print(tukey_fresas$groups)\n",
            "\n",
            "# 6.2 Forma nativa rápida de graficar en agricolae\n",
            "plot(tukey_fresas, main = \"Grupos de Comparación Rápida (agricolae)\")\n",
            "\n",
            "# 6.3 Gráfico Científico Personalizado en ggplot2 con Medias, Barras de Error y Letras de Tukey\n",
            "# Extraemos las medias, desviaciones estándar y los grupos asignados por Tukey\n",
            "df_tukey <- data.frame(\n",
            "  Tratamiento = rownames(tukey_fresas$means),\n",
            "  Media = tukey_fresas$means$Dias,\n",
            "  SD = tukey_fresas$means$std,\n",
            "  Rep = tukey_fresas$means$r\n",
            ")\n",
            "# Calculamos el Error Estándar de la Media (EEM = SD / sqrt(n))\n",
            "df_tukey$EEM <- df_tukey$SD / sqrt(df_tukey$Rep)\n",
            "\n",
            "# Agregamos los grupos estadísticos correspondientes\n",
            "df_tukey$Grupo <- tukey_fresas$groups[rownames(df_tukey), \"groups\"]\n",
            "\n",
            "# Graficamos las Medias con Barras de Error y Letras\n",
            "ggplot(df_tukey, aes(x = Tratamiento, y = Media, fill = Tratamiento)) +\n",
            "  geom_bar(stat = \"identity\", color = \"black\", alpha = 0.7, width = 0.5) +\n",
            "  geom_errorbar(aes(ymin = Media - EEM, ymax = Media + EEM), width = 0.15, size = 0.8) +\n",
            "  geom_text(aes(label = Grupo, y = Media + EEM + 0.3), size = 6, fontface = \"bold\") +\n",
            "  scale_fill_brewer(palette = \"Set1\") +\n",
            "  theme_minimal(base_size = 14) +\n",
            "  labs(\n",
            "    title = \"Comparación de Medias con Grupos de Tukey\",\n",
            "    subtitle = \"Letras distintas indican diferencias significativas (p < 0.05). Barras representan el EEM.\",\n",
            "    x = \"Tratamiento\",\n",
            "    y = \"Media de Días de Conservación\"\n",
            "  ) +\n",
            "  theme(legend.position = \"none\")\n",
            "\n",
            "# 6.4 Otras pruebas para fines didácticos\n",
            "print(\"--- PRUEBA DE DUNCAN ---\")\n",
            "duncan_fresas <- duncan.test(modelo_fresas, \"Tratamiento\", group = TRUE)\n",
            "print(duncan_fresas$groups)\n",
            "\n",
            "print(\"--- PRUEBA LSD DE FISHER ---\")\n",
            "lsd_fresas <- LSD.test(modelo_fresas, \"Tratamiento\", group = TRUE)\n",
            "print(lsd_fresas$groups)\n",
            "\n",
            "print(\"--- PRUEBA DE DUNNETT (VS. CONTROL) ---\")\n",
            "dunnett_fresas <- DunnettTest(Dias ~ Tratamiento, data = datos_fresas, control = \"Control\")\n",
            "print(dunnett_fresas)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### 💡 **Guía Pedagógica: ¿Cómo interpretar las letras de comparación?**\n",
            "\n",
            "Las letras mostradas al lado de cada tratamiento en los resultados anteriores (agrupaciones de Tukey/Duncan/LSD) son una forma compacta y visual de resumir miles de pruebas de hipótesis simultáneas:\n",
            "\n",
            "* **Regla de oro:** Si dos tratamientos **comparten al menos una letra**, la diferencia entre sus medias **NO** es estadísticamente significativa ($p \\ge 0.05$).\n",
            "* **Si tienen letras totalmente diferentes:** La diferencia entre sus medias **SÍ** es estadísticamente significativa ($p < 0.05$).\n",
            "\n",
            "**Ejemplo práctico:**\n",
            "* Si el Quitosano tiene la letra `a` y el Almidón la letra `b`, significa que el Quitosano supera al Almidón de manera estadísticamente significativa.\n",
            "* Si la Gelatina y el Control comparten la letra `c`, se asume que ambos se comportan de forma estadísticamente idéntica.\n",
            "* Si un tratamiento hipotético tuviera la etiqueta `ab`, significa que no difiere estadísticamente de los que tienen `a`, ni tampoco de los que tienen `b`."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "# **TALLER EVALUATIVO**\n",
            "\n",
            "**Nombre del Estudiante:** ___________________________  \n",
            "**Fecha de Entrega:** ___________________________  \n",
            "\n",
            "---"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## EJERCICIO 1: Inoculantes en Pino (Sintaxis Guiada)\n",
            "\n",
            "**Contexto Forestal:**  \n",
            "Un vivero comercial evalúa el efecto de 3 inoculantes de micorrizas (`M1`, `M2`, `M3`) sobre la altura final (cm) de plántulas de *Pinus radiata* tras 6 meses en condiciones homogéneas de vivero, frente a plántulas Control (sin inoculación).\n",
            "\n",
            "**Instrucciones:** Completa los espacios vacíos `___` en el código a continuación para ajustar el modelo, validar supuestos de forma cuantitativa, realizar el ANOVA y ejecutar Tukey."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# EJERCICIO 1: CREACIÓN DE LOS DATOS DE PINO\n",
            "# ==============================================================================\n",
            "\n",
            "datos_pino <- data.frame(\n",
            "  Micorriza = factor(rep(c(\"Control\", \"M1\", \"M2\", \"M3\"), each = 5)),\n",
            "  Altura = c(12.5, 11.8, 13.1, 12.2, 11.9,  # Control\n",
            "             15.2, 16.1, 14.8, 15.5, 15.9,  # M1\n",
            "             18.5, 19.2, 17.8, 18.1, 18.9,  # M2\n",
            "             14.1, 13.8, 14.5, 15.0, 14.2)  # M3\n",
            ")\n",
            "\n",
            "head(datos_pino)"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# 1. Ajustar el modelo lineal en DCA\n",
            "# Completa con la variable de respuesta y el factor\n",
            "modelo_pino <- aov(___ ~ ___, data = datos_pino)\n",
            "\n",
            "# 2. Extraer residuos y probar Normalidad formalmente\n",
            "residuos_pino <- residuals(___)\n",
            "shapiro.test(___)\n",
            "\n",
            "# 3. Probar Homocedasticidad formalmente con la prueba de Levene (Librería car)\n",
            "leveneTest(___ ~ ___, data = datos_pino)\n",
            "\n",
            "# 4. Desplegar e interpretar la tabla ANOVA\n",
            "summary(___)\n",
            "\n",
            "# 5. Ejecutar la prueba de Tukey HSD para el modelo de pino\n",
            "tukey_pino <- HSD.test(___, \"Micorriza\", group = TRUE)\n",
            "print(tukey_pino$groups)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## EJERCICIO 2: Toma de Decisiones en Raleo de Eucalipto\n",
            "\n",
            "**Contexto Forestal:**  \n",
            "Un silvicultor industrial está investigando el efecto de 3 intensidades de raleo (`Control` sin raleo, `Ligero`, `Fuerte`) sobre el Diámetro a la Altura del Pecho (DAP en cm) de *Eucalyptus globulus* después de 5 años.\n",
            "\n",
            "**Restricción comercial crítica:**  \n",
            "La empresa realizará una fuerte inversión basada en estos resultados. Si recomiendan raleos innecesarios o agresivos debido a un \"falso positivo\" (Error Tipo I), las pérdidas operativas serán millonarias. Por lo tanto, el departamento técnico exige el análisis estadístico más conservador y riguroso posible.\n",
            "\n",
            "**Hoja de ruta obligatoria para tu análisis (4 Pasos):**\n",
            "1. **Visualización:** Genera un gráfico adecuado (diagrama de cajas y puntos jitter con `ggplot2`) para explorar la distribución de los tratamientos.\n",
            "2. **Supuestos:** Valida formalmente la normalidad (Shapiro-Wilk) y homocedasticidad (Levene) de los residuos de forma cuantitativa.\n",
            "3. **ANOVA:** Ajusta y evalúa si existen diferencias significativas globales en el DAP.\n",
            "4. **Post-Hoc Adecuado:** Elige y justifica de manera conceptual la prueba de comparaciones múltiples más idónea considerando la restricción comercial de evitar falsos positivos (Error Tipo I)."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# EJERCICIO 2: CARGA DE LOS DATOS DE EUCALIPTO\n",
            "# ==============================================================================\n",
            "\n",
            "datos_eucalipto <- data.frame(\n",
            "  Raleo = factor(rep(c(\"Control\", \"Ligero\", \"Fuerte\"), each = 6)),\n",
            "  DAP = c(15.2, 14.8, 15.6, 16.1, 14.9, 15.4,   # Control\n",
            "          18.5, 19.1, 17.9, 18.7, 19.5, 18.2,   # Ligero\n",
            "          22.1, 23.5, 21.8, 22.9, 24.0, 21.5)   # Fuerte\n",
            ")\n",
            "\n",
            "head(datos_eucalipto)"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe aquí tu análisis completo paso a paso (Visualización, Supuestos, ANOVA y Post-Hoc adecuado)]\n",
            "\n"
        ]
    }
]

# Configuración del archivo de notebook .ipynb con kernel de R
notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {
            "display_name": "R",
            "language": "R",
            "name": "ir"
        },
        "language_info": {
            "name": "R"
        }
    },
    "nbformat": 4,
    "nbformat_minor": 4
}

# Escribir el archivo final estructurado
for name in ["Taller_02_DCA_Comparaciones.ipynb", "Taller2_DCA.ipynb"]:
    with open(name, "w", encoding="utf-8") as f:
        json.dump(notebook, f, indent=2, ensure_ascii=False)

print("¡Notebooks 'Taller_02_DCA_Comparaciones.ipynb' y 'Taller2_DCA.ipynb' regenerados con éxito!")
