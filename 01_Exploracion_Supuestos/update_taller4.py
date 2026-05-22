import json

# Define the cells of the updated Jupyter Notebook for Taller 4
cells = [
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "# **Taller 4: Diseños Factoriales en DCA y DBCA (Arreglos Factoriales)**\n",
            "\n",
            "**Curso:** Diseño Experimental y Aplicaciones Estadísticas  \n",
            "**Nivel:** Universitario (Pregrado)\n",
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
            "2. **Configuración Automática:** La celda de abajo instalará de forma silenciosa y automática todas las librerías necesarias (`ggplot2`, `performance`, `see`, `agricolae`, `DescTools`, `readxl`) sin inundar la pantalla con registros de descarga.\n",
            "3. **Fácil Acceso a Datos:** El primer ejercicio descarga de forma remota el conjunto de datos de Excel desde GitHub (o usa un fallback local programático). ¡Todo es 100% interactivo y no requiere subir archivos manualmente!\n",
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
            "paquetes_requeridos <- c(\"ggplot2\", \"performance\", \"see\", \"agricolae\", \"DescTools\", \"readxl\")\n",
            "suppressMessages(preparar_entorno(paquetes_requeridos))\n"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 1. Fundamentos Teóricos de los Diseños Factoriales\n",
            "\n",
            "### ¿Qué es un Arreglo Factorial y cuándo se debe usar?\n",
            "En el diseño experimental clásico, evaluar un solo factor a la vez (como estudiar el riego por un lado y la variedad de semilla por otro) suele ser ineficiente e incompleto. Un **Arreglo Factorial** no es un diseño de aleatorización en sí mismo, sino una forma de estructurar los **tratamientos**. Consiste en combinar todos los niveles de dos o más factores independientes en un solo experimento, permitiendo evaluarlos simultáneamente.\n",
            "\n",
            "Se debe utilizar cuando queremos responder dos tipos de interrogantes de forma simultánea:\n",
            "1. **Efectos Principales:** ¿Cuál es el impacto aislado del Factor A (por ejemplo, Variedad) y el Factor B (por ejemplo, Riego) sobre la variable de respuesta?\n",
            "2. **Efecto de Interacción (A $\\times$ B):** ¿El comportamiento de los niveles del Factor A varía de forma significativa según los niveles del Factor B? (Sinergismo o Antagonismo). Si el efecto de un factor cambia según el nivel del otro, se dice que existe **Interacción**.\n",
            "\n",
            "### Diferencia Crítica: Arreglo Factorial en DCA vs. DBCA\n",
            "\n",
            "*   **Factorial en DCA (Diseño Completamente al Azar):**\n",
            "    Se utiliza cuando las unidades experimentales (ej. macetas en invernadero, bandejas homogéneas) son uniformes. Las combinaciones de los tratamientos ($a \\times b$) se asignan enteramente al azar a las unidades.\n",
            "    *   **Modelo Aditivo Lineal (DCA Factorial):**\n",
            "        $$Y_{ijk} = \\mu + \\alpha_i + \\beta_j + (\\alpha\\beta)_{ij} + \\varepsilon_{ijk}$$\n",
            "        Donde:\n",
            "        * $\\mu$: Media general del experimento.\n",
            "        * $\\alpha_i$: Efecto principal del Factor A (nivel $i$).\n",
            "        * $\\beta_j$: Efecto principal del Factor B (nivel $j$).\n",
            "        * $(\\alpha\\beta)_{ij}$: Efecto de la Interacción entre el nivel $i$ de A y el nivel $j$ de B.\n",
            "        * $\\varepsilon_{ijk}$: Error experimental residual $\\varepsilon_{ijk} \\sim N(0, \\sigma^2)$ independientes.\n",
            "\n",
            "*   **Factorial en DBCA (Diseño de Bloques Completos al Azar):**\n",
            "    Se utiliza cuando las unidades de campo presentan variabilidad en una dirección espacial (ej. ladera con pendiente, gradiente de sombra). Dividimos el campo en $r$ bloques homogéneos. Cada bloque debe ser lo suficientemente grande para albergar **todas las combinaciones de tratamientos** (completitud). La aleatorización de las combinaciones se realiza de forma independiente *dentro* de cada bloque.\n",
            "    *   **Modelo Aditivo Lineal (DBCA Factorial):**\n",
            "        $$Y_{ijk} = \\mu + \\alpha_i + \\beta_j + (\\alpha\\beta)_{ij} + \\gamma_k + \\varepsilon_{ijk}$$\n",
            "        Donde se añade:\n",
            "        * $\\gamma_k$: Efecto del Bloque $k$. Se asume que no existe interacción entre los bloques y los tratamientos combinados.\n",
            "\n",
            "---\n",
            "\n",
            "### La Pérdida de Grados de Libertad: El \"Costo\" del Bloqueo en Diseños Factoriales\n",
            "Al igual que en el diseño de bloques simple, el bloqueo en arreglos factoriales \"consume\" Grados de Libertad (GL) de la Suma de Cuadrados del Error. Veamos la distribución para un Factorial con dos factores (A y B) con $a$ y $b$ niveles y $r$ repeticiones/bloques:\n",
            "\n",
            "| Fuente de Variación | Factorial en DCA (GL) | Factorial en DBCA (GL) | Ejemplo Numérico ($a=2, b=2, r=3$) |\n",
            "| :--- | :--- | :--- | :--- |\n",
            "| **Factor A (a)** | $a - 1$ | $a - 1$ | $1$ |\n",
            "| **Factor B (b)** | $b - 1$ | $b - 1$ | $1$ |\n",
            "| **Interacción (A $\\times$ B)** | $(a - 1)(b - 1)$ | $(a - 1)(b - 1)$ | $1$ |\n",
            "| **Bloques (r)** | — | $r - 1$ | $2$ (Solo en DBCA) |\n",
            "| **Error Residual** | $ab(r - 1)$ | $(ab - 1)(r - 1)$ | DCA: $8$ GL vs. DBCA: $6$ GL |\n",
            "| **Total (N - 1)** | $abr - 1$ | $abr - 1$ | $11$ |\n",
            "\n",
            "**¿Por qué es esto crucial?**\n",
            "Si bloqueamos de forma innecesaria (sin un gradiente real de variación en el campo), habremos restado $r - 1$ grados de libertad al error. Esto disminuye la precisión de la varianza residual, infla el valor crítico de la distribución $F$ y **reduce la potencia del experimento** para detectar diferencias reales entre tratamientos o para detectar una interacción significativa.\n",
            "\n",
            "---\n",
            "\n",
            "### Identificación de la Unidad Experimental\n",
            "En diseño experimental, la **Unidad Experimental** es la fracción mínima de material experimental a la que se le aplica de forma independiente y aleatoria un tratamiento (en este caso, una combinación de factores).\n",
            "*   **En nuestro ejercicio principal (Variedad y Riego):** La unidad experimental **NO** es la planta individual de maíz. La unidad experimental es **la parcela individual de campo (sub-parcela)** en la que se implementa una combinación específica de Variedad y Riego dentro de un bloque determinado. La producción de la parcela es cosechada, pesada y reportada."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### Esquema de un Arreglo Factorial $2 \\times 2$ en DBCA (Croquis de Campo)\n",
            "Consideremos un ensayo agrícola con **2 Variedades de Maíz** ($V_1$, $V_2$) y **2 Frecuencias de Riego** ($R_1$: Goteo, $R_2$: Aspersión) sobre un terreno con una pendiente (gradiente de fertilidad de arriba a abajo).\n",
            "\n",
            "Tenemos $2 \\times 2 = 4$ Tratamientos combinados:\n",
            "*   **T1:** $V_1R_1$ (Variedad 1 con Goteo)\n",
            "*   **T2:** $V_1R_2$ (Variedad 1 con Aspersión)\n",
            "*   **T3:** $V_2R_1$ (Variedad 2 con Goteo)\n",
            "*   **T4:** $V_2R_2$ (Variedad 2 con Aspersión)\n",
            "\n",
            "Establecemos **3 Bloques** perpendiculares al gradiente. Cada bloque contiene las 4 combinaciones distribuidas enteramente al azar:\n",
            "\n",
            "```\n",
            "                    [ DIRECCIÓN DE LA PENDIENTE: ZONA ALTA / SECA ]\n",
            "Norte  ========================================================================\n",
            "       BLOQUE 1 (Zona Alta - Suelo menos profundo)\n",
            "       +-------------------+-------------------+-------------------+-------------------+\n",
            "       |    V1R2 (T2)      |    V2R1 (T3)      |    V1R1 (T1)      |    V2R2 (T4)      |  <- Aleatorio\n",
            "       +-------------------+-------------------+-------------------+-------------------+\n",
            "       ========================================================================\n",
            "       BLOQUE 2 (Zona Media)\n",
            "       +-------------------+-------------------+-------------------+-------------------+\n",
            "       |    V2R1 (T3)      |    V1R1 (T1)      |    V2R2 (T4)      |    V1R2 (T2)      |  <- Aleatorio\n",
            "       +-------------------+-------------------+-------------------+-------------------+\n",
            "       ========================================================================\n",
            "       BLOQUE 3 (Zona Baja - Suelo más profundo y húmedo)\n",
            "       +-------------------+-------------------+-------------------+-------------------+\n",
            "       |    V1R1 (T1)      |    V2R2 (T4)      |    V1R2 (T2)      |    V2R1 (T3)      |  <- Aleatorio\n",
            "       +-------------------+-------------------+-------------------+-------------------+\n",
            "Sur    ========================================================================\n",
            "                    [ DIRECCIÓN DE LA PENDIENTE: ZONA BAJA / HÚMEDA ]\n",
            "```\n",
            "\n",
            "*   **Nota de Aleatorización:** Las 4 combinaciones aparecen exactamente una vez dentro de cada bloque, pero su distribución espacial es aleatoria y única por bloque, lo cual controla el gradiente sistemático."
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 2. Configuración del Entorno de R\n",
            "\n",
            "Instalaremos y cargaremos la suite de paquetes indispensables para realizar análisis factoriales, importación de Excel, diagnóstico de supuestos y visualización científica.\n",
            "\n",
            "*(Nota: Si estás usando **Google Colab**, la celda autoloader del inicio ya instaló y cargó todo de forma silenciosa. Si trabajas de forma **local**, puedes ejecutar la celda de abajo para configurar tu entorno sin interrupciones interactivas)*"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Evitar diálogos interactivos de CRAN\n",
            "options(repos = c(CRAN = \"https://cloud.r-project.org\"))\n",
            "\n",
            "# Carga/Instalación no interactiva de librerías\n",
            "paquetes <- c(\"ggplot2\", \"performance\", \"see\", \"DescTools\", \"agricolae\", \"readxl\")\n",
            "for (p in paquetes) {\n",
            "  if (!require(p, character.only = TRUE)) {\n",
            "    install.packages(p, dependencies = TRUE, quiet = TRUE)\n",
            "    library(p, character.only = TRUE)\n",
            "  }\n",
            "}\n",
            "message(\"¡Entorno de R listo!\")"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 3. Importación de Datos de Excel y Estructuración en R\n",
            "\n",
            "### Carga de Datos desde Hojas de Cálculo (.xlsx)\n",
            "En el ámbito profesional y académico, los datos crudos se registran usualmente en Microsoft Excel. En R, empleamos el paquete `readxl` y su función `read_excel()` para importar las hojas de datos sin necesidad de convertirlas previamente a CSV.\n",
            "\n",
            "### La Conversión Mandatoria a Factores en Diseños Factoriales\n",
            "En un análisis factorial, es imperativo que las variables de clasificación (`Variedad`, `Riego`, `Bloque`) se conviertan de forma obligatoria a factores utilizando `as.factor()`. Si dejamos los factores como caracteres o números continuos, R interpretará incorrectamente el modelo como una regresión continua y no computará la descomposición de varianzas del ANOVA por tratamiento y bloque."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# 2. IMPORTACIÓN Y ESTRUCTURACIÓN DE DATOS DESDE EXCEL / GITHUB\n",
            "# ==============================================================================\n",
            "\n",
            "excel_file <- \"datos_taller4_factorial.xlsx\"\n",
            "\n",
            "# Descargar desde GitHub si no existe (ej. en Google Colab)\n",
            "if (!file.exists(excel_file)) {\n",
            "  message(\"Descargando datos_taller4_factorial.xlsx desde GitHub...\")\n",
            "  url_git <- \"https://raw.githubusercontent.com/PALP31/BioAgro-Stats/main/01_Exploracion_Supuestos/datos_taller4_factorial.xlsx\"\n",
            "  tryCatch({\n",
            "    download.file(url_git, destfile = excel_file, mode = \"wb\", quiet = TRUE)\n",
            "    message(\"¡Archivo descargado exitosamente!\")\n",
            "  }, error = function(e) {\n",
            "    message(\"⚠️ Error en descarga. Cargando datos de forma programática...\")\n",
            "  })\n",
            "}\n",
            "\n",
            "# Carga de datos\n",
            "if (file.exists(excel_file)) {\n",
            "  datos_fact <- read_excel(excel_file)\n",
            "} else {\n",
            "  # Fallback programático idéntico en caso de fallo de red\n",
            "  datos_fact <- data.frame(\n",
            "    Variedad = c(\"V1\", \"V1\", \"V2\", \"V2\", \"V1\", \"V1\", \"V2\", \"V2\", \"V1\", \"V1\", \"V2\", \"V2\"),\n",
            "    Riego = c(\"Goteo\", \"Aspersion\", \"Goteo\", \"Aspersion\", \"Goteo\", \"Aspersion\", \"Goteo\", \"Aspersion\", \"Goteo\", \"Aspersion\", \"Goteo\", \"Aspersion\"),\n",
            "    Bloque = c(\"Inv_1\", \"Inv_1\", \"Inv_1\", \"Inv_1\", \"Inv_2\", \"Inv_2\", \"Inv_2\", \"Inv_2\", \"Inv_3\", \"Inv_3\", \"Inv_3\", \"Inv_3\"),\n",
            "    Produccion = c(15, 10, 18, 12, 16, 11, 20, 13, 14, 9, 17, 11)\n",
            "  )\n",
            "  message(\"✅ Base de datos cargada vía fallback programático.\")\n",
            "}\n",
            "\n",
            "# Conversión mandatoria a factores categóricos\n",
            "datos_fact$Variedad <- as.factor(datos_fact$Variedad)\n",
            "datos_fact$Riego    <- as.factor(datos_fact$Riego)\n",
            "datos_fact$Bloque   <- as.factor(datos_fact$Bloque)\n",
            "\n",
            "# Visualización de la estructura y primeros datos\n",
            "head(datos_fact)\n",
            "str(datos_fact)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 4. Gráficos Exploratorios y Perfiles de Interacción\n",
            "\n",
            "Antes de ajustar el ANOVA, realizamos una exploración visual para comprender la naturaleza de los datos. La herramienta clave en los diseños factoriales es el **Gráfico de Interacción (Perfiles de Medias)**:\n",
            "*   Si las líneas que conectan las respuestas promedio de los niveles de un factor a través de los niveles del otro son **paralelas**, indica que los factores actúan de forma **aditiva independiente** (no hay interacción significativa).\n",
            "*   Si las líneas **se cruzan o muestran pendientes opuestas**, es un claro indicio visual de que existe una **Interacción significativa** (el efecto del riego depende de la variedad utilizada)."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# 3. GRÁFICOS EXPLORATORIOS DE INTERACCIÓN\n",
            "# ==============================================================================\n",
            "\n",
            "# Gráfico de perfil nativo en R\n",
            "interaction.plot(\n",
            "  x.factor = datos_fact$Riego,\n",
            "  trace.factor = datos_fact$Variedad,\n",
            "  response = datos_fact$Produccion,\n",
            "  type = \"b\",\n",
            "  pch = c(19, 17),\n",
            "  col = c(\"blue\", \"red\"),\n",
            "  xlab = \"Método de Riego\",\n",
            "  ylab = \"Producción Promedio (t/ha)\",\n",
            "  legend = TRUE,\n",
            "  main = \"Gráfico de Perfil: Interacción Riego * Variedad\"\n",
            ")"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 5. Ajuste del Modelo e Interpretación del ANOVA\n",
            "\n",
            "### La Especificación de la Interacción en R\n",
            "En R, el operador **`*`** representa la inclusión de los efectos principales y su interacción:\n",
            "`Variedad * Riego` equivale a escribir: `Variedad + Riego + Variedad:Riego` (donde `:` representa la interacción pura).\n",
            "Como nuestro diseño se aleatorizó en bloques completos al azar, agregamos el factor de control ambiental mediante `+ Bloque`. El modelo queda:\n",
            "`Produccion ~ Variedad * Riego + Bloque`\n",
            "\n",
            "### Jerarquía de la Interpretación en Diseños Factoriales\n",
            "**Regla de Oro:** Siempre evaluamos primero la significancia del término de la **Interacción (Variedad:Riego)**:\n",
            "1.  **Si la Interacción es Significativa ($p < 0.05$):** Concluimos que el efecto del Riego cambia según la Variedad. **No** debemos interpretar los efectos principales de forma aislada, ya que estaríamos dando conclusiones incompletas o erróneas. El análisis post-hoc de Tukey debe realizarse sobre las combinaciones de tratamientos.\n",
            "2.  **Si la Interacción NO es Significativa ($p > 0.05$):** Los factores son independientes. Procedemos a interpretar los efectos principales (Variedad y Riego por separado) como si fueran experimentos independientes."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# 4. AJUSTE DEL MODELO LINEAL FACTORIAL DBCA Y ANOVA\n",
            "# ==============================================================================\n",
            "\n",
            "# Ajuste del modelo lineal con interacción y bloques\n",
            "modelo_fact <- aov(Produccion ~ Variedad * Riego + Bloque, data = datos_fact)\n",
            "\n",
            "# Visualización de la tabla del Análisis de Varianza (ANOVA)\n",
            "summary(modelo_fact)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 6. Validación de los Supuestos del Modelo (Residuos)\n",
            "\n",
            "Antes de validar diferencias de medias, debemos comprobar que los residuos cumplan con los supuestos matemáticos de la inferencia lineal:\n",
            "\n",
            "1.  **Normalidad (Shapiro-Wilk):** Evaluamos si los residuos provienen de una distribución normal ($H_0$: Residuos normales).\n",
            "2.  **Homocedasticidad (Prueba de Levene):** En arreglos factoriales, el supuesto de homogeneidad de varianzas se verifica sobre los **tratamientos combinados** (la interacción de los factores), no sobre cada factor por separado. Usamos la función `LeveneTest` de `DescTools`."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# 6.1 Diagnóstico gráfico de residuos automatizado con performance\n",
            "check_model(modelo_fact, check = c(\"normality\", \"homogeneity\"))"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# 6.2 Prueba Formal Cuantitativa de Normalidad: Shapiro-Wilk (H0: Residuos normales)\n",
            "residuos_f <- residuals(modelo_fact)\n",
            "shapiro.test(residuos_f)"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# 6.3 Prueba Formal Cuantitativa de Homocedasticidad: Levene (H0: Varianzas iguales)\n",
            "# Creamos la interacción de tratamientos para Levene\n",
            "datos_fact$Tratamiento_Combinado <- interaction(datos_fact$Variedad, datos_fact$Riego)\n",
            "\n",
            "# Prueba formal de Levene\n",
            "LeveneTest(Produccion ~ Tratamiento_Combinado, data = datos_fact)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 7. Pruebas de Comparaciones Múltiples de Medias (Post-Hoc)\n",
            "\n",
            "Si la interacción `Variedad:Riego` es significativa, realizamos la comparación múltiple para evaluar la combinación de niveles mediante la prueba de **Tukey HSD**.\n",
            "\n",
            "Utilizando el paquete `agricolae`, especificamos un vector con los dos factores en el argumento correspondiente: `c(\"Variedad\", \"Riego\")`. Esto generará un análisis de agrupamiento para las combinaciones de tratamientos evaluando cuál combinación maximiza la respuesta."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Comparación múltiple de medias con Tukey HSD para la Interacción\n",
            "tukey_fact <- HSD.test(modelo_fact, c(\"Variedad\", \"Riego\"), group = TRUE, console = TRUE)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 8. Visualización Avanzada de Comparaciones de Medias y Letras\n",
            "\n",
            "La comunicación visual científica del arreglo factorial exige mostrar las medias y sus barras de error estándar para cada una de las 4 combinaciones, acompañadas por las letras de Tukey correspondientes.\n",
            "\n",
            "Dividimos esta visualización en dos celdas:\n",
            "1.  **Gráfica Básica (Nativa):** Utilizando el método `plot()` provisto por la librería `agricolae`.\n",
            "2.  **Gráfica Premium (ggplot2):** Creando un gráfico de barras de alta calidad editorial con barras de error estándar de la media (EEM) y letras de Tukey."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# 8.1 Gráfica rápida nativa de agricolae\n",
            "plot(tukey_fact, variation = 'SE', col = 'skyblue', main = 'Tukey HSD: Combinación Variedad * Riego')"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# 8.2 Visualización Premium con ggplot2, Medias, EEM y Letras de Tukey\n",
            "# Construcción del dataframe para ggplot\n",
            "df_tukey_f <- data.frame(\n",
            "  Tratamiento = rownames(tukey_fact$means),\n",
            "  Media = tukey_fact$means$Produccion,\n",
            "  SD = tukey_fact$means$std,\n",
            "  Rep = tukey_fact$means$r\n",
            ")\n",
            "\n",
            "# Cálculo de las barras de error (Error Estándar de la Media - EEM)\n",
            "df_tukey_f$EEM <- df_tukey_f$SD / sqrt(df_tukey_f$Rep)\n",
            "df_tukey_f$Grupo <- tukey_fact$groups[rownames(df_tukey_f), 'groups']\n",
            "\n",
            "# Gráfico premium con ggplot2\n",
            "ggplot(df_tukey_f, aes(x = Tratamiento, y = Media, fill = Tratamiento)) +\n",
            "  geom_bar(stat = 'identity', color = 'black', alpha = 0.75, width = 0.5) +\n",
            "  geom_errorbar(aes(ymin = Media - EEM, ymax = Media + EEM), width = 0.15, size = 0.8) +\n",
            "  geom_text(aes(label = Grupo, y = Media + EEM + 0.4), size = 6, fontface = 'bold') +\n",
            "  scale_fill_brewer(palette = 'Accent') +\n",
            "  theme_minimal(base_size = 14) +\n",
            "  labs(\n",
            "    title = 'Rendimiento Promedio por Combinación de Tratamientos (Variedad:Riego)',\n",
            "    subtitle = 'Letras distintas indican diferencias significativas (Tukey HSD, p < 0.05). Barras representan EEM.',\n",
            "    x = 'Tratamiento Combinado (Variedad : Riego)',\n",
            "    y = 'Producción Promedio (t/ha)'\n",
            "  ) +\n",
            "  theme(legend.position = 'none')"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 9. Cuestionario de Análisis del Experimento\n",
            "\n",
            "Responde de manera detallada y analítica las siguientes preguntas basadas en tu estudio factorial:\n",
            "\n",
            "1. **¿El término de interacción `Variedad:Riego` resultó estadísticamente significativo en el ANOVA? ¿Qué implicaciones agronómicas y prácticas tiene este resultado para un agricultor?**\n",
            "   \n",
            "   *Escribe tu respuesta aquí:*\n",
            "   \n",
            "2. **De acuerdo con la prueba de Tukey, ¿cuál es la mejor combinación de Variedad y Riego para maximizar la producción? ¿Hay diferencias estadísticas entre la mejor combinación y la segunda mejor?**\n",
            "   \n",
            "   *Escribe tu respuesta aquí:*\n",
            "   \n",
            "3. **¿La inclusión del factor 'Bloque' (Invernaderos) fue efectiva en este diseño experimental? Justifica tu respuesta observando el valor de p del factor Bloque en el ANOVA.**\n",
            "   \n",
            "   *Escribe tu respuesta aquí:*\n",
            "   \n",
            "4. **¿Qué supuestos matemáticos de los residuos fueron validados y qué nos indican los p-valores obtenidos en Shapiro-Wilk y Levene?**\n",
            "   \n",
            "   *Escribe tu respuesta aquí:*"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "# **ACTIVIDAD PRÁCTICA AUTÓNOMA (Taller de Aplicación)**\n",
            "\n",
            "**Contexto Agrícola de la Tarea:**  \n",
            "Un agrónomo desea evaluar cómo interactúan **2 niveles de factor Semilla** (`S1`, `S2`) y **3 niveles de factor Espaciamiento** entre plantas (`10cm`, `20cm`, `30cm`) sobre la variable de respuesta **Altura de la planta** en centímetros (`Altura_cm`). \n",
            "\n",
            "Debido a que el ensayo de campo se estableció en un terreno con variabilidad sistemática de la textura del suelo, el investigador agrupó las parcelas experimentales en **3 Bloques** basados en la textura predominante del suelo (`Franco`, `Arcilloso`, `Arenoso`). Este ensayo corresponde a un **Arreglo Factorial $2 \\times 3$ conducido bajo un Diseño de Bloques Completos al Azar (DBCA)**.\n",
            "\n",
            "**Conceptos de Preparación:**\n",
            "*   **Unidad Experimental:** Es la micro-parcela o surco de cultivo a la que se le aplica de forma independiente y aleatoria una de las 6 combinaciones de tratamiento (ej: Semilla S1 con Espaciamiento 10cm) dentro de un bloque de suelo determinado.\n",
            "*   **Distribución del Gradiente en el Campo:** Los 3 Bloques se organizan horizontalmente sobre cada franja de textura homogénea. Dentro de cada bloque se siembran las 6 combinaciones aleatorizadas.\n",
            "\n",
            "### Esquema del Croquis de Campo ($2 \\times 3$ Factorial en DBCA):\n",
            "```\n",
            "                 [ DIRECCIÓN DE LA TEXTURA DEL SUELO ]\n",
            "============================================================================\n",
            "BLOQUE 1 (Textura Suelo: Franco)\n",
            "+-----------+-----------+-----------+-----------+-----------+-----------+\n",
            "|  S1-10cm  |  S2-30cm  |  S1-20cm  |  S2-10cm  |  S1-30cm  |  S2-20cm  |  <- Aleatorio\n",
            "+-----------+-----------+-----------+-----------+-----------+-----------+\n",
            "============================================================================\n",
            "BLOQUE 2 (Textura Suelo: Arcilloso)\n",
            "+-----------+-----------+-----------+-----------+-----------+-----------+\n",
            "|  S2-20cm  |  S1-10cm  |  S2-10cm  |  S1-30cm  |  S2-30cm  |  S1-20cm  |  <- Aleatorio\n",
            "+-----------+-----------+-----------+-----------+-----------+-----------+\n",
            "============================================================================\n",
            "BLOQUE 3 (Textura Suelo: Arenoso)\n",
            "+-----------+-----------+-----------+-----------+-----------+-----------+\n",
            "|  S1-30cm  |  S1-20cm  |  S2-20cm  |  S2-30cm  |  S1-10cm  |  S2-10cm  |  <- Aleatorio\n",
            "+-----------+-----------+-----------+-----------+-----------+-----------+\n",
            "============================================================================\n",
            "```\n",
            "\n",
            "### Objetivos del Estudiante:\n",
            "1.  **Importación desde Excel:** Cargar el archivo de datos real `tarea_taller4_factorial.xlsx` usando `readxl`.\n",
            "2.  **Preparación de Datos:** Convertir las columnas categóricas a tipo factor.\n",
            "3.  **Visualización Exploratoria:** Generar un boxplot o gráfico exploratorio inicial.\n",
            "4.  **Ajuste del ANOVA Factorial:** Ajustar el modelo incluyendo la interacción `Semilla * Espaciamiento` más los bloques. Desplegar la tabla ANOVA.\n",
            "5.  **Validación de los Supuestos:**\n",
            "    *   Diagnóstico gráfico con `check_model()`.\n",
            "    *   Prueba cuantitativa de Shapiro-Wilk para los residuos.\n",
            "    *   Prueba cuantitativa de Levene para los tratamientos combinados.\n",
            "6.  **Prueba de Tukey Post-Hoc:** Ejecutar la prueba de Tukey HSD sobre las combinaciones de factores para agrupar las medias.\n",
            "7.  **Visualización de Comparación de Medias:**\n",
            "    *   Gráfica rápida nativa con `plot()`.\n",
            "    *   Gráfica de barras premium con `ggplot2` mostrando EEM y letras de Tukey.\n",
            "8.  **Cuestionario de Interpretación:** Responder a las interrogantes agronómicas de los agricultores."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# ==============================================================================\n",
            "# DOCENTE: DESCARGA AUTOMÁTICA DEL ARCHIVO DE LA TAREA PARA GOOGLE COLAB\n",
            "# ==============================================================================\n",
            "\n",
            "excel_tarea <- \"tarea_taller4_factorial.xlsx\"\n",
            "\n",
            "# Si el archivo no existe localmente (ej. en Google Colab), lo descargamos de GitHub\n",
            "if (!file.exists(excel_tarea)) {\n",
            "  message(\"Descargando tarea_taller4_factorial.xlsx desde GitHub...\")\n",
            "  url_tarea_git <- \"https://raw.githubusercontent.com/PALP31/BioAgro-Stats/main/01_Exploracion_Supuestos/tarea_taller4_factorial.xlsx\"\n",
            "  tryCatch({\n",
            "    download.file(url_tarea_git, destfile = excel_tarea, mode = \"wb\", quiet = TRUE)\n",
            "    message(\"¡Archivo de la tarea descargado exitosamente!\")\n",
            "  }, error = function(e) {\n",
            "    message(\"⚠️ Error en descarga. Generando archivo de tarea local alternativo (CSV) como contingencia...\")\n",
            "    df_contingencia <- data.frame(\n",
            "      Semilla = c(\"S1\", \"S1\", \"S1\", \"S2\", \"S2\", \"S2\", \"S1\", \"S1\", \"S1\", \"S2\", \"S2\", \"S2\", \"S1\", \"S1\", \"S1\", \"S2\", \"S2\", \"S2\"),\n",
            "      Espaciamiento = c(\"10cm\", \"20cm\", \"30cm\", \"10cm\", \"20cm\", \"30cm\", \"10cm\", \"20cm\", \"30cm\", \"10cm\", \"20cm\", \"30cm\", \"10cm\", \"20cm\", \"30cm\", \"10cm\", \"20cm\", \"30cm\"),\n",
            "      Textura_Suelo = c(\"Franco\", \"Franco\", \"Franco\", \"Franco\", \"Franco\", \"Franco\", \"Arcilloso\", \"Arcilloso\", \"Arcilloso\", \"Arcilloso\", \"Arcilloso\", \"Arcilloso\", \"Arenoso\", \"Arenoso\", \"Arenoso\", \"Arenoso\", \"Arenoso\", \"Arenoso\"),\n",
            "      Altura_cm = c(10.5, 12.1, 14.5, 11.2, 13.0, 15.2, 9.8, 11.5, 13.8, 10.5, 12.3, 14.5, 8.5, 10.2, 12.5, 9.2, 11.0, 13.1)\n",
            "    )\n",
            "    write.csv(df_contingencia, \"tarea_taller4_factorial.csv\", row.names = FALSE)\n",
            "  })\n",
            "}"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### **Paso 1: Importación del Excel de la Tarea**\n",
            "\n",
            "**Instrucciones:** Escribe la sentencia de R adecuada para importar el archivo `tarea_taller4_factorial.xlsx` usando la función `read_excel()`. Asigna los datos a la variable `datos_tarea4` y visualiza la estructura con `str()` y las primeras líneas con `head()`.\n",
            "\n",
            "*(Tip de contingencia: Si la descarga de red falló, puedes leer 'tarea_taller4_factorial.csv' usando la función `read.csv('tarea_taller4_factorial.csv')`)*"
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe tu código para importar el archivo Excel de la tarea aquí]\n",
            "datos_tarea4 <- ___(excel_tarea)\n",
            "\n",
            "# Visualizar la estructura y las primeras líneas de los datos\n",
            "str(datos_tarea4)\n",
            "head(datos_tarea4)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### **Paso 2: Conversión a Factores**\n",
            "\n",
            "**Instrucciones:** Asegura que las columnas `Semilla` (Factor A), `Espaciamiento` (Factor B) y `Textura_Suelo` (Bloques) sean convertidas a tipo factor de R de manera mandatoria utilizando `as.factor()`. Comprueba los cambios con `str()`."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe tu código para la conversión obligatoria a factores aquí]\n",
            "datos_tarea4$Semilla       <- ___(datos_tarea4$Semilla)\n",
            "datos_tarea4$Espaciamiento <- ___(datos_tarea4$Espaciamiento)\n",
            "datos_tarea4$Textura_Suelo <- ___(datos_tarea4$Textura_Suelo)\n",
            "\n",
            "# Comprobar estructura\n",
            "str(datos_tarea4)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### **Paso 3: Perfil de Interacción Visual**\n",
            "\n",
            "**Instrucciones:** Construye el gráfico de interacción (`interaction.plot`) para evaluar de forma preliminar si las pendientes de la altura de la planta a través de los diferentes espaciamientos difieren según el tipo de semilla. Dibuja las líneas usando `type = 'b'`."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe tu código para construir el gráfico de perfiles de interacción aquí]\n",
            "interaction.plot(\n",
            "  x.factor = datos_tarea4$___,\n",
            "  trace.factor = datos_tarea4$___,\n",
            "  response = datos_tarea4$___,\n",
            "  type = \"b\",\n",
            "  pch = c(19, 17),\n",
            "  col = c(\"blue\", \"red\"),\n",
            "  xlab = \"Espaciamiento\",\n",
            "  ylab = \"Altura Promedio (cm)\",\n",
            "  legend = TRUE,\n",
            "  main = \"Gráfico de Perfil: Interacción Semilla * Espaciamiento\"\n",
            ")"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### **Paso 4: Ajuste del Modelo Factorial en DBCA y ANOVA**\n",
            "\n",
            "**Instrucciones:** Ajusta el modelo factorial lineal aditivo en DBCA incluyendo la interacción de `Semilla` y `Espaciamiento`, y el factor de bloqueo `Textura_Suelo`. Asigna el resultado a `modelo_tarea4` y despliega el ANOVA con `summary()`."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe tu código para ajustar el modelo y desplegar la tabla ANOVA aquí]\n",
            "modelo_tarea4 <- aov(Altura_cm ~ Semilla ___ Espaciamiento + Textura_Suelo, data = datos_tarea4)\n",
            "summary(modelo_tarea4)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### **Paso 5: Validación de los Supuestos del Modelo (Residuos)**\n",
            "\n",
            "Verifica por separado los supuestos teóricos del modelo sobre los residuos de `modelo_tarea4`:\n",
            "\n",
            "#### **5.1 Diagnóstico Gráfico Integral**\n",
            "**Instrucciones:** Utiliza `check_model()` de la librería `performance` para inspeccionar la normalidad y la homogeneidad de varianza gráficamente."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe tu código para el diagnóstico gráfico automatizado de residuos aquí]\n",
            "check_model(modelo_tarea4, check = c(\"normality\", \"homogeneity\"))"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "#### **5.2 Prueba Formal de Normalidad**\n",
            "**Instrucciones:** Extrae los residuos de `modelo_tarea4` y evalúa cuantitativamente la normalidad mediante la prueba de `shapiro.test()`."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe tu código para evaluar la normalidad mediante Shapiro-Wilk aquí]\n",
            "residuos_tarea4 <- residuals(___)\n",
            "shapiro.test(residuos_tarea4)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "#### **5.3 Prueba Formal de Homocedasticidad**\n",
            "**Instrucciones:** Crea la interacción categórica de tratamientos y evalúa si existe homogeneidad de varianzas mediante la prueba de `LeveneTest()` (de la librería `DescTools`)."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe tu código para evaluar la homocedasticidad mediante LeveneTest aquí]\n",
            "datos_tarea4$Tratamiento_Combinado <- interaction(datos_tarea4$Semilla, datos_tarea4$Espaciamiento)\n",
            "LeveneTest(Altura_cm ~ Tratamiento_Combinado, data = ___)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### **Paso 6: Comparación de Medias mediante Tukey HSD**\n",
            "\n",
            "**Instrucciones:** Si la interacción `Semilla:Espaciamiento` resultó significativa en la tabla ANOVA, ejecuta la prueba post-hoc de Tukey HSD para las combinaciones de factores utilizando la función `HSD.test()` de `agricolae`. Guarda el resultado en `tukey_tarea4` e imprime el agrupamiento en consola."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe tu código para ejecutar el test Tukey HSD aquí]\n",
            "tukey_tarea4 <- HSD.test(modelo_tarea4, c(\"Semilla\", \"Espaciamiento\"), group = TRUE, console = TRUE)"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### **Paso 7: Visualización Avanzada de Medias y Letras de Tukey**\n",
            "\n",
            "#### **7.1 Gráfica Básica (Nativa)**\n",
            "**Instrucciones:** Dibuja la gráfica rápida nativa aplicando la función `plot()` sobre el objeto `tukey_tarea4`."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe tu código para generar la gráfica nativa de Tukey aquí]\n",
            "plot(tukey_tarea4, variation = 'SE', col = 'lightgreen', main = 'Tukey HSD: Tarea Semilla * Espaciamiento')"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "#### **7.2 Gráfica Premium con ggplot2**\n",
            "**Instrucciones:** Extrae los promedios y errores estándar de `tukey_tarea4` y construye un gráfico de barras premium con ggplot2. Añade las barras de error estándar de la media (EEM) y dibuja las letras de Tukey sobre cada barra."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# [Escribe tu código para elaborar la visualización premium de ggplot2 aquí]\n",
            "df_tukey_t4 <- data.frame(\n",
            "  Tratamiento = rownames(tukey_tarea4$means),\n",
            "  Media = tukey_tarea4$means$Altura_cm,\n",
            "  SD = tukey_tarea4$means$std,\n",
            "  Rep = tukey_tarea4$means$r\n",
            ")\n",
            "\n",
            "# Cálculo de las barras de error (Error Estándar de la Media)\n",
            "df_tukey_t4$EEM <- df_tukey_t4$SD / sqrt(df_tukey_t4$Rep)\n",
            "df_tukey_t4$Grupo <- tukey_tarea4$groups[rownames(df_tukey_t4), 'groups']\n",
            "\n",
            "# Gráfico premium con ggplot2\n",
            "ggplot(df_tukey_t4, aes(x = Tratamiento, y = Media, fill = Tratamiento)) +\n",
            "  geom_bar(stat = 'identity', color = 'black', alpha = 0.75, width = 0.5) +\n",
            "  geom_errorbar(aes(ymin = Media - EEM, ymax = Media + EEM), width = 0.15, size = 0.8) +\n",
            "  geom_text(aes(label = Grupo, y = Media + EEM + 0.3), size = 6, fontface = 'bold') +\n",
            "  scale_fill_brewer(palette = 'Spectral') +\n",
            "  theme_minimal(base_size = 14) +\n",
            "  labs(\n",
            "    title = 'Altura Promedio por Combinación de Tratamientos (Semilla:Espaciamiento)',\n",
            "    subtitle = 'Letras distintas indican diferencias significativas (Tukey HSD, p < 0.05). Barras representan EEM.',\n",
            "    x = 'Tratamiento Combinado (Semilla : Espaciamiento)',\n",
            "    y = 'Altura Promedio (cm)'\n",
            "  ) +\n",
            "  theme(legend.position = 'none')"
        ]
    },
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "### **Paso 8: Cuestionario Final de Conclusiones Técnicas**\n",
            "\n",
            "Responde a las siguientes preguntas analizando críticamente los resultados del ensayo de espaciamientos y semillas:\n",
            "\n",
            "1. **¿Existe una interacción estadísticamente significativa entre la Semilla y el Espaciamiento en la altura de la planta ($p < 0.05$)? ¿Qué interpretación física o agronómica le das a esta interacción?**\n",
            "   \n",
            "   *Escribe tu respuesta aquí:*\n",
            "   \n",
            "2. **De acuerdo al agrupamiento del test de Tukey, ¿cuál o cuáles combinaciones de Semilla y Espaciamiento producen plantas significativamente más altas? Justifica basándote en las letras de rango.**\n",
            "   \n",
            "   *Escribe tu respuesta aquí:*\n",
            "   \n",
            "3. **¿La textura del suelo (Bloques) resultó significativa en el ANOVA? ¿Fue adecuado bloquear el terreno experimental por textura o se habría obtenido más potencia estadística usando un DCA simple?**\n",
            "   \n",
            "   *Escribe tu respuesta aquí:*\n",
            "   \n",
            "4. **¿Los supuestos de normalidad y homocedasticidad se cumplieron de manera satisfactoria para los residuos de este modelo? Argumenta tu respuesta citando los valores de p obtenidos en Shapiro-Wilk y Levene.**\n",
            "   \n",
            "   *Escribe tu respuesta aquí:*"
        ]
    }
]

# Configure the JSON notebook structure with R kernel metadata
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

# Write the final .ipynb files
for name in ["Taller_04_Factorial.ipynb", "Taller4_Factorial.ipynb"]:
    with open(name, "w", encoding="utf-8") as f:
        json.dump(notebook, f, indent=2, ensure_ascii=False)

print("¡Notebooks 'Taller_04_Factorial.ipynb' y 'Taller4_Factorial.ipynb' creados con éxito!")
