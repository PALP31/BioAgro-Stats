# 🧬 BioAgro-Stats: Repositorio de Análisis Estadístico Avanzado para Biotecnología y Agronomía

Este repositorio centraliza una colección robusta de scripts en **R** diseñados para el análisis riguroso de datos experimentales en ciencias biológicas, agrónomicas y biotecnológicas. Desde la validación exhaustiva de supuestos hasta el modelado bayesiano y machine learning, cada módulo está optimizado para garantizar la integridad científica y la reproducibilidad de los resultados.

---

## 📂 Contenido del Repositorio

| Módulo | Descripción | Herramientas Clave |
| :--- | :--- | :--- |
| **01. Exploración y Supuestos** | Diagnóstico experto de modelos. Validación de normalidad, homocedasticidad e independencia mediante residuos simulados y análisis visual de alto nivel. | `performance`, `DHARMa`, `easystats` |
| **02. Diseño Experimental** | Implementación de diseños clásicos y complejos: DBCA, Alpha-Lattice para mejoramiento genético y Parcelas Divididas (Split-plot). | `agricolae`, `emmeans`, `multcomp` |
| **03. Modelos Frecuentistas** | Modelado clásico mediante ANOVA, Modelos Lineales Generalizados (GLM) para conteos y Modelos Mixtos (LMM) para estructuras jerárquicas. | `lme4`, `lmerTest`, `car` |
| **04. Estadística Bayesiana** | Transición al paradigma bayesiano para modelos complejos: ANOVA bayesiano, GLMMs, interacción GxE y modelos de crecimiento no lineales. | `brms`, `tidybayes`, `bayesplot` |
| **05. Machine Learning** | Modelos predictivos aplicados: Random Forest para espectroscopía, XGBoost para estrés abiótico y SVM para clasificación de enfermedades. | `tidymodels`, `xgboost`, `vip` |

---

## 🚀 Otros Módulos Incluidos

- **04. Modelos No Lineales (GAMs):** Análisis de series temporales y curvas de respuesta dinámica mediante `mgcv` y `gratia`.
- **06. Multivariado y Ómicas:** Análisis de comunidades (PERMANOVA/NMDS) y algoritmos para datos de alta dimensión (sPLS-DA) con `vegan` y `mixOmics`.
- **Reportes Quarto:** Plantillas de alta calidad para la generación automática de informes técnicos y manuscritos científicos.

---

## 🛠️ Requisitos y Dependencias

Para ejecutar estos scripts, se recomienda tener instalada la última versión de **R** y **RStudio**. Los paquetes principales se agrupan por funcionalidad:

### 📦 Manipulación y Visualización
- `tidyverse` (dplyr, ggplot2, tidyr, purrr)
- `patchwork` (Combinación de gráficos)
- `ggsci` & `viridis` (Paletas científicas)

### 📊 Modelado Estadístico
- **Base:** `lme4`, `lmerTest`, `mgcv`, `nlme`, `MASS`
- **Post-hoc:** `emmeans`, `multcomp`, `multcompView`
- **Diseño:** `agricolae`

### 🕯️ Ecosistema Bayesiano
- `brms` (Interfaz para Stan)
- `tidybayes`, `bayesplot`, `bayestestR`

### 🤖 Machine Learning
- `tidymodels`, `xgboost`, `kernlab`, `vip`

### 🧪 Diagnóstico (Imprescindible)
- `performance` (Parte de `easystats`)
- `DHARMa` (Residuos simulados)

---

## 📖 Cómo usar este repositorio
| Paso | Instrucción |
| :---: | :--- |
| 1️⃣ | Clona el repositorio: `git clone https://github.com/PALP31/BioAgro-Stats.git` |
| 2️⃣ | Abre en **RStudio** el archivo `.Rproj` del proyecto para conservar rutas relativas y una sesión de trabajo reproducible. |
| 3️⃣ | Navega a la carpeta del módulo que desees explorar y ejecuta los scripts `.R` según tu flujo de análisis. |
| 4️⃣ | Usa los datos de ejemplo incluidos (o generados en los scripts) para practicar y adaptar los análisis a tus propios experimentos. |

---
**Desarrollado para el análisis de datos de alto impacto en Ciencias Agrarias.** 🌾🧪
