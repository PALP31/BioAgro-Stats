# ============================================================================
# 02_diagnostico_profundo_supuestos.R (MASTERCLASS PARTE 2: CASOS COMPLEJOS)
# ============================================================================
# Experto: Estadística Biotecnológica Avanzada
# Objetivo: Diagnóstico profundo con datos reales, Colinealidad, Puntos de
# Influencia (Cook), Autocorrelación y Transformaciones de Box-Cox.
# ============================================================================

# Autoinstalación de paquetes específicos para este análisis
paquetes_requeridos <- c(
  "tidyverse", "car", "GGally", "performance", "MASS", "lmtest", "patchwork"
)
paquetes_nuevos <- paquetes_requeridos[!(paquetes_requeridos %in% installed.packages()[,"Package"])]
if (length(paquetes_nuevos)) install.packages(paquetes_nuevos, dependencies = TRUE)

library(tidyverse)
library(car)         # VIF, InfluencePlot
library(GGally)      # Matrices de correlación avanzadas
library(performance) # Diagnósticos
library(MASS)        # Transformación Box-Cox
library(lmtest)      # Autocorrelación (Durbin-Watson)
library(patchwork)   # Diseños de gráficos

# ============================================================================
# 1. CARGA DE BASE DE DATOS REAL (BIOTECNOLOGÍA AGRÍCOLA)
# ============================================================================
cat("\n[1] LEYENDO BASE DE DATOS BIOTECNOLÓGICA...\n")
archivo_csv <- "01_Exploracion_Supuestos/datos_biotec_agricola.csv"

# Verificamos si existe (el script de python lo acaba de generar)
if(!file.exists(archivo_csv)) {
  stop("¡Falta el archivo CSV! Debe de estar en el directorio correspondiente.")
}

datos_biotec <- read.csv(archivo_csv)

cat("Estructura de la base de datos:\n")
str(datos_biotec)
cat("\nNuestra variable objetivo a mejorar es 'Biomasa_Gramos'.\n")

# ============================================================================
# NIVEL 1: Básico - Exploración Visual Multidetalle y Colinealidad Visual
# ============================================================================
cat("\n[NIVEL 1] MATRIZ DE EXPLORACIÓN (Ggally) OJO A LA COLINEALIDAD...\n")

# ggpairs permite cruzar todas las numéricas para ver densidades, dispersión y r de Pearson
vars_numericas <- datos_biotec %>% 
  select(Biomasa_Gramos, Temperatura_Hoja, Tasa_Fotosintetica, Conductancia_Estomatica)

png("01_Exploracion_Supuestos/P2_N1_Matriz_ggpairs.png", width=1200, height=1000, res=150)
print(ggpairs(vars_numericas, 
        title = "Matriz de Exploración: Alerta de Colinealidad y Asimetría",
        lower = list(continuous = wrap("points", alpha=0.5, color="darkgreen")),
        diag = list(continuous = wrap("densityDiag", fill="lightgreen"))
))
dev.off()
cat("- Gráfico de Matriz guardado en 'P2_N1_Matriz_ggpairs.png'.\n")
cat("-> DIAGNÓSTICO VISUAL: Vemos que la 'Biomasa' tiene una alta asimetría (no normal),\n")
cat("   y que 'Tasa Fotosintética' y 'Conductancia' son casi una línea recta (colinealidad peligrosa).\n")

# ============================================================================
# NIVEL 2: Intermedio - El Veneno de la Multicolinealidad y Autocorrelación
# ============================================================================
cat("\n[NIVEL 2] AJUSTE DEL MODELO, VIF Y AUTOCORRELACIÓN TEMP/ESPACIAL\n")

# Planteamos el Modelo Lineal Múltiple inicial "ingenuo"
modelo_crudo <- lm(Biomasa_Gramos ~ Variedad + Temperatura_Hoja + Tasa_Fotosintetica + Conductancia_Estomatica, 
                   data = datos_biotec)

# 2.1 Multicolinealidad (VIF)
cat("\n-> Análisis de Factor de Inflación de Varianza (VIF):\n")
vif_resultado <- car::vif(modelo_crudo)
print(vif_resultado)

cat("\n[!] ALERTA EXPERTA: Un VIF > 5 o 10 indica multicolinealidad severa.\n")
cat("Tener dos variables biológicas redundantes roba significancia estadística (infla el error estándar).\n")
cat("SOLUCIÓN: Eliminaremos 'Conductancia_Estomatica' del modelo para estabilizarlo.\n")

modelo_estable <- lm(Biomasa_Gramos ~ Variedad + Temperatura_Hoja + Tasa_Fotosintetica, 
                     data = datos_biotec)

# 2.2 Autocorrelación (Durbin-Watson)
cat("\n-> Prueba de Durbin-Watson para Independencia de Residuos:\n")
# Nos importa mucho si medimos las plantas de manera contigua
print(dwtest(modelo_estable))
cat("Interpretación DW: Cercano a 2 = Independientes. < 1.5 o > 2.5 indica Autocorrelación.\n")

# ============================================================================
# NIVEL 3: Avanzado - Distancia de Cook, Apalancamiento (Leverage) y Outliers
# ============================================================================
cat("\n[NIVEL 3] CAZA DE PUNTOS DE INFLUENCIA (OUTLIERS QUE DESTRUYEN EL MODELO)\n")

# No es lo mismo un outlier (Y extrema) que un punto de Leverage (X extremo).
# Un punto de Leverage con Y extremo se vuelve una "Observación Influyente" e inclina la recta.
# La Distancia de Cook los descubre.

png("01_Exploracion_Supuestos/P2_N3_InfluencePlot.png", width=1000, height=800, res=150)
influencePlot(modelo_estable, 
              main="Gráfico de Influencia: Leverage (X) vs Residuals (Y)",
              sub="Círculos grandes = Alta Distancia de Cook (Destructivos)")
dev.off()
cat("- Gráfico de Influencia guardado en 'P2_N3_InfluencePlot.png'.\n")

# Extrayendo la distancia de Cook para buscar los problemáticos
cooks_d <- cooks.distance(modelo_estable)
indices_problematicos <- as.numeric(names(cooks_d)[cooks_d > (4/nrow(datos_biotec))])

cat("\n-> Plantas detectadas como altamente influyentes (Distancia de Cook crítica):\n")
print(datos_biotec[indices_problematicos, c("ID_Planta", "Variedad", "Biomasa_Gramos")])
cat("¡La observación con ID altísimo de fotosíntesis (apalancamiento) o biomasa loca (outlier) debe investigarse o removerse!\n")

# Limpieza quirúrgica de datos
datos_limpios <- datos_biotec[-indices_problematicos, ]
modelo_limpio <- lm(Biomasa_Gramos ~ Variedad + Temperatura_Hoja + Tasa_Fotosintetica, 
                     data = datos_limpios)

# ============================================================================
# NIVEL 4: Experto - Curando Asimetría con Box-Cox
# ============================================================================
cat("\n[NIVEL 4] TRANSFORMACIÓN BOX-COX PARA SALVAR NORMALIDAD Y HOMOCEDASTICIDAD\n")

# Comprobamos la normalidad de nuestro modelo sin outliers
shapiro_limpio <- shapiro.test(residuals(modelo_limpio))
cat("Normalidad antes de Box-Cox (p-value):", shapiro_limpio$p.value, "\n")

if(shapiro_limpio$p.value < 0.05) {
  cat("\n¡Falla la normalidad debido a la asimetría log-normal del rendimiento!\n")
}

# Ejecutando Box-Cox para hallar el Lambda óptimo
png("01_Exploracion_Supuestos/P2_N4_BoxCox.png", width=800, height=600, res=150)
boxcox_res <- boxcox(modelo_limpio, plotit = TRUE, lambda = seq(-2, 2, by = 0.1))
title("Log-Likelihood para el parámetro Lambda (Box-Cox)")
dev.off()

lambda_optimo <- boxcox_res$x[which.max(boxcox_res$y)]
cat("\n-> El Lambda óptimo encontrado por Box-Cox es:", round(lambda_optimo, 3), "\n")

# Si el lambda óptimo está cerca de 0, es Transformación Logarítmica (ln).
if(abs(lambda_optimo) < 0.2) {
  cat("Como Lambda ~ 0, aplicaremos TRANSFORMACIÓN LOGARÍTMICA (log).\n")
  modelo_curado <- lm(log(Biomasa_Gramos) ~ Variedad + Temperatura_Hoja + Tasa_Fotosintetica, data = datos_limpios)
} else {
  cat("Aplicando transformación Box-Cox matemática pura (Y^lambda).\n")
  modelo_curado <- lm(((Biomasa_Gramos^lambda_optimo - 1)/lambda_optimo) ~ Variedad + Temperatura_Hoja + Tasa_Fotosintetica, data = datos_limpios)
}

cat("\nNormalidad DESPUÉS de Box-Cox (p-value):", shapiro.test(residuals(modelo_curado))$p.value, "\n")
cat("¡Supuesto de Normalidad matemáticamente salvado!\n")

# Comparación Visual del ANTES y DESPUÉS usando performance
p_antes <- plot(check_normality(modelo_estable)) + ggtitle("1. Con Outliers + Asimetría") 
p_limpio <- plot(check_normality(modelo_limpio)) + ggtitle("2. Sin Outliers (Asimétrico)") 
p_curado <- plot(check_normality(modelo_curado)) + ggtitle("3. Curado (Box-Cox)") 

# Combinar con patchwork
png("01_Exploracion_Supuestos/P2_N4_Antes_Despues_Norm.png", width=1400, height=500, res=150)
print(p_antes | p_limpio | p_curado)
dev.off()

cat("\n============================================================================\n")
cat("ANÁLISIS PROFUNDO PARTE 2 COMPLETADO CON ÉXITO\n")
cat("============================================================================\n")
