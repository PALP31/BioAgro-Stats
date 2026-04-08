# Scripts-estadistica- 📊

Repositorio de análisis estadístico para datos experimentales biológicos y agronómicos. Incluye scripts de ejemplo en **R** y **Python**, plantillas de reportes en **Quarto** y modelos estadísticos frecuentistas y bayesianos listos para adaptar a nuevos experimentos.

---

## 📁 Estructura del repositorio

```
Scripts-estadistica-/
├── R/
│   ├── 01_modelos_frecuentistas.R   # ANOVA, LM, LMM, GLM, GAM, GLMM
│   └── 02_modelos_bayesianos.R      # Análisis bayesiano con brms + Stan
├── Python/
│   ├── 01_modelos_estadisticos.py   # ANOVA, LM, LMM, GLM, GAM (statsmodels + pygam)
│   └── 02_modelos_bayesianos.py     # Análisis bayesiano con PyMC + ArviZ
├── quarto_templates/
│   ├── reporte_analisis.qmd         # Plantilla de reporte frecuentista
│   └── reporte_bayesiano.qmd        # Plantilla de reporte bayesiano
├── .gitignore                       # R, Python y Quarto
└── README.md
```

---

## 🔢 Modelos incluidos

| Modelo | R | Python |
|---|---|---|
| ANOVA (una vía / factorial) | `aov()` | `statsmodels` OLS + `anova_lm` |
| Modelo Lineal (LM) | `lm()` | `smf.ols()` |
| Modelo Lineal Mixto (LMM) | `lme4::lmer()` | `smf.mixedlm()` |
| GLM – Poisson | `glm(..., family = poisson)` | `smf.glm(..., family=Poisson)` |
| GLM – Gamma | `glm(..., family = Gamma)` | `smf.glm(..., family=Gamma)` |
| GAM | `mgcv::gam()` | `pygam.LinearGAM` |
| GLMM | `lme4::glmer()` | `pymer4` / R recomendado |
| Bayesiana (LM / LMM / GLMM) | `brms::brm()` | `pymc.Model` |

---

## ⚙️ Requisitos

### R

Instalar `pacman` (gestor de paquetes) y los paquetes se instalan automáticamente al ejecutar los scripts:

```r
install.packages("pacman")
```

Paquetes principales: `tidyverse`, `lme4`, `lmerTest`, `emmeans`, `mgcv`, `car`, `performance`, `brms`, `bayesplot`, `tidybayes`, `ggdist`, `posterior`, `loo`.

> **Stan / brms:** Sigue las instrucciones de instalación en <https://paul-buerkner.github.io/brms/>.

### Python

Crear un entorno virtual e instalar dependencias:

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install numpy pandas matplotlib seaborn statsmodels pygam pymc arviz
```

### Quarto

Instalar Quarto CLI desde <https://quarto.org/docs/get-started/>.  
Renderizar una plantilla:

```bash
quarto render quarto_templates/reporte_analisis.qmd
```

---

## 🚀 Uso rápido

### Scripts de R

```r
# 1. Abre R o RStudio
# 2. Ajusta la sección "Datos" para cargar tu propio archivo
# 3. Ejecuta el script completo o sección por sección
source("R/01_modelos_frecuentistas.R")
```

### Scripts de Python

```bash
python Python/01_modelos_estadisticos.py
```

### Reportes Quarto

1. Copia la plantilla que necesites desde `quarto_templates/`.
2. Ajusta los datos, título, autor y parámetros del modelo.
3. Renderiza:

```bash
quarto render quarto_templates/reporte_analisis.qmd --to html
quarto render quarto_templates/reporte_bayesiano.qmd --to pdf
```

---

## 📋 Flujo de trabajo sugerido

```
1. Preparar datos (CSV / Excel)
   ↓
2. Exploración y estadística descriptiva
   ↓
3. Seleccionar modelo apropiado (ver tabla)
   ↓
4. Verificar supuestos del modelo
   ↓
5. Ajustar e interpretar resultados
   ↓
6. Generar reporte con plantilla Quarto
```

---

## 🤝 Contribuciones

¡Las contribuciones son bienvenidas! Por favor:

1. Haz un *fork* del repositorio.
2. Crea una rama descriptiva: `git checkout -b feat/nuevo-modelo`.
3. Realiza tus cambios y abre un *pull request*.

---

## 📄 Licencia

Este repositorio se distribuye bajo la licencia **MIT**. Consulta el archivo `LICENSE` para más detalles.
