# 🧠 Memoria de Progreso: Talleres de Aplicaciones Estadísticas

Esta memoria sirve como punto de control (checkpoint) para que el asistente de IA retome exactamente el trabajo donde se quedó, evitando repetir explicaciones, autoloader de paquetes o la depuración de sintaxis de R en Google Colab.

---

## 📌 Estado del Proyecto
*   **Repositorio GitHub:** [PALP31/BioAgro-Stats](https://github.com/PALP31/BioAgro-Stats)
*   **Rama Activa:** `main`
*   **Directorio de Trabajo local:** `/Users/paullopez/Documents/Aplicaciones_estadistica_talleres`
*   **Carpeta de Destino:** `01_Exploracion_Supuestos/` (aloja todos los entregables de los talleres para ser consumidos directamente en Google Colab).

---

## 🛠️ Talleres Completados (1 a 4)

Hemos finalizado y sincronizado con éxito la suite de **4 Talleres prácticos** para los estudiantes. Todos cuentan con:
1.  **Compatibilidad 100% con Google Colab (Kernel de R):** Se incluye el badge oficial e instrucciones para configurar R en la nube.
2.  **Autocargador Silencioso de CRAN (Silent Autoloader):** Instala de forma transparente y silenciosa los paquetes requeridos (`performance`, `DHARMa`, `agricolae`, `emmeans`, `readxl`, `ggplot2`, etc.) sin interrumpir al estudiante ni solicitar selección de espejos (mirrors).
3.  **Descarga Automática de Datos (GitHub Fallback):** Si se ejecuta en la nube (Colab), descarga los datos directamente de GitHub. Si se ejecuta en local, usa la ruta relativa del archivo local.
4.  **Ejercicios Rellenables (`___`):** Bloques de código con guiones bajos para fomentar el aprendizaje activo.

### Detalle de los Talleres:
*   **Taller 1: Introducción a R y Quarto** (`Taller_01_Introduccion_R.qmd`, `Taller_01_Introduccion_R.ipynb`, `Taller_01_Introduccion_R.R`)
    *   *Conceptos clave:* Sintaxis básica, vectores, dataframes, filtros y visualización exploratoria.
*   **Taller 2: Diseño Completamente al Azar (DCA)** (`Taller2_DCA.qmd`, `Taller2_DCA.ipynb`, `Taller2_DCA.R`, `Taller_02_DCA_Comparaciones.ipynb`, `update_taller2.py`)
    *   *Conceptos clave:* ANOVA de una vía, la trampa del factor (factor labels), análisis visual de residuos, supuestos de normalidad/homocedasticidad y letras de comparación múltiple (Tukey).
    *   *Dataset unificado:* Altura de plantas micorrizadas (4 niveles: Control, M1, M2, M3).
*   **Taller 3: Diseño en Bloques Completos al Azar (DBCA)** (`Taller3_DBCA.qmd`, `Taller3_DBCA.ipynb`, `Taller3_DBCA.R`, `Taller_03_DBCA_Importacion.ipynb`, `update_taller3.py`, datos, tareas)
    *   *Conceptos clave:* Control local (bloques), importación avanzada Excel/CSV, ANOVA bidireccional y pruebas post-hoc con letras compactas.
*   **Taller 4: Diseño Factorial** (`Taller4_Factorial.qmd`, `Taller4_Factorial.ipynb`, `Taller4_Factorial.R`, `Taller_04_Factorial.ipynb`, `update_taller4.py`, datos, tareas)
    *   *Conceptos clave:* Factores cruzados (A x B), efectos principales, término de interacción, perfiles de interacción (interaction plots) y Tukey multifactorial en R.

---

## 📈 Estándares de Diseño Implementados

*   **Identificadores Únicos y Metadatos:** Todos los cuadernos Jupyter cuentan con metadatos explícitos del kernel de R (`display_name: R`, `language: R`, `name: ir`) para que Colab inicie en modo R de inmediato.
*   **Explicaciones Visuales y Pedagógicas:**
    *   Explicación interactiva de cómo leer las letras de Tukey.
    *   Comparativas claras entre tratamientos para el entendimiento del estudiante.
*   **Scripts Python de Generación:** Los scripts `update_taller2.py`, `update_taller3.py`, `update_taller4.py` permiten reconstruir o modificar los `.ipynb` de forma programática.

---

## 🔮 Siguiente Paso (Taller 5 y posteriores)
Cuando solicites hacer un nuevo taller (por ejemplo, **Taller 5: Modelos Lineales Generalizados (GLMs)** o **Modelos Mixtos (LMMs)**):
1.  **No repitas** el código del autoloader ni la descarga de datos desde GitHub; copia directamente la estructura de la función `preparar_entorno` y el bloque `GitHub Fallback` utilizado en los talleres anteriores.
2.  Mantén el formato bilingüe/español con tonos académicos motivadores y pedagógicos.
3.  Utiliza los scripts generadores de Jupyter en Python (`update_tallerX.py`) para mantener sincronizados el archivo `.ipynb` y `.qmd`.

*Memoria guardada el 22 de mayo de 2026.*
