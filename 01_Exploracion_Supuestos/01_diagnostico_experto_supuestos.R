# ============================================================================
# 01_diagnostico_experto_supuestos.R (MASTERCLASS DE SUPUESTOS)
# ============================================================================
# Experto: Estadística y Diseño Experimental
# Objetivo: Guía progresiva (Sencillo -> Difícil) para la verificación de
# supuestos estadísticos en investigación agrícola, biotecnológica y de campo.
# ============================================================================

# Autoinstalación de paquetes faltantes
paquetes_requeridos <- c(
  "tidyverse", "performance", "see", "car", "lmtest", 
  "patchwork", "DHARMa", "nortest", "corrplot", "fitdistrplus", "lme4"
)
paquetes_nuevos <- paquetes_requeridos[!(paquetes_requeridos %in% installed.packages()[,"Package"])]
if (length(paquetes_nuevos)) install.packages(paquetes_nuevos, dependencies = TRUE)

# Cargar librerías
library(tidyverse)
library(performance) # Diagnósticos modernos y visuales
library(see)         # Visualización de ecosistema easystats
library(car)         # Pruebas clásicas de supuestos
library(lmtest)      # Pruebas heterocedasticidad y autocorrelación
library(patchwork)   # Combinación de gráficos complejos
library(DHARMa)      # Diagnóstico avanzado de residuos simulados
library(nortest)     # Evaluación exhaustiva de normalidad
library(corrplot)    # Correlogramas
library(fitdistrplus)# Ajuste experto de distribuciones estadísticas
library(lme4)        # Modelamiento mixto

set.seed(42) # Reproducibilidad científica

# Creación de directorio si no existe (para guardar los gráficos)
if(!dir.exists("01_Exploracion_Supuestos")) {
  dir.create("01_Exploracion_Supuestos")
}

# ============================================================================
# 1. SIMULACIÓN PARAMETRIZADA (Enfoque Agrícola/Biotecnológico)
# ============================================================================
cat("\n[1] GENERANDO MATRICES DE DATOS SIMULADOS...\n")

# Caso A: Laboratorio (Enzimología y Multi-rasgo Continuo)
# Se incluye correlación y outliers
datos_lab <- data.frame(
  genotipo = factor(rep(c("Silvestre", "Crispr"), each = 30)),
  clorofila = rnorm(60, mean = 45, sd = 4),
  nitrogeno_foliar = rnorm(60, mean = 2.5, sd = 0.3)
) %>%
  mutate(
    # Modelo generador con ruido normal
    act_enzimatica = 15 + 0.6 * clorofila + 3 * nitrogeno_foliar + 
                     ifelse(genotipo == "Crispr", 8, 0) + rnorm(60, 0, 3)
  )
# Inducir dos valores atípicos (outliers de pipeteo)
datos_lab$act_enzimatica[c(10, 50)] <- c(80, 5) 

# Caso B: Invernadero (Biomasa con Heterocedasticidad)
# Varianza se desestabiliza a mayor estrés
datos_inv <- data.frame(
  tratamiento = factor(rep(c("T0_Control", "T1_Sequia", "T2_Salinidad"), each = 20)),
  biomasa = c(
    rnorm(20, mean = 60, sd = 3),   # Homogéneo
    rnorm(20, mean = 35, sd = 9),   # Varianza Alta
    rnorm(20, mean = 25, sd = 15)   # Varianza Extrema
  )
)

# Caso C: Campo (Incidencia de Fitopatógenos - Conteo)
# Exceso de ceros y sobredispersión inherente a la naturaleza epifitiológica
datos_campo <- data.frame(
  manejo = factor(rep(c("Testigo", "Fungicida_Sintetico", "Trichoderma"), each = 40)),
  lesiones_hoja = c(rpois(40, lambda = 12), rpois(40, lambda = 1), rpois(40, lambda = 3))
)
# Agregamos exceso de ceros artificial (ej: plantas que escaparon del inóculo)
datos_campo$lesiones_hoja[sample(1:120, 25)] <- 0 

# ============================================================================
# NIVEL 1: Sencillo - Exploración Visual, Correlaciones y Outliers
# ============================================================================
cat("\n[NIVEL 1] INICIO DE ANÁLISIS SENCILLO (EXPLORATORIO)\n")

# 1.1 Matriz de Correlación
vars_num_lab <- datos_lab %>% select(clorofila, nitrogeno_foliar, act_enzimatica)
matriz_corr <- cor(vars_num_lab)

cat("-> Matriz de Correlación (Laboratorio):\n")
print(matriz_corr)

png("01_Exploracion_Supuestos/N1_correlaciones_lab.png", width=800, height=800, res=150)
corrplot(matriz_corr, method = "color", type = "upper", addCoef.col = "black",
         tl.col = "darkblue", tl.srt = 45, diag = FALSE,
         title = "Correlación Traits Biotecnológicos", mar=c(0,0,2,0))
dev.off()

# 1.2 Histogramas y Boxplots para Outliers
p_hist <- ggplot(datos_lab, aes(x = act_enzimatica, fill = genotipo)) +
  geom_histogram(bins = 20, alpha = 0.6, color = "black", position = "identity") +
  theme_minimal() + labs(title = "Distribución Enzimática", y = "Frecuencia")

p_box <- ggplot(datos_lab, aes(x = genotipo, y = act_enzimatica, fill = genotipo)) +
  geom_boxplot(outlier.colour = "red", outlier.size = 3) +
  theme_minimal() + labs(title = "Detección de Outliers (Puntos Rojos)")

ggsave("01_Exploracion_Supuestos/N1_exploracion_visual.png", 
       plot = (p_hist / p_box), width = 8, height = 7, dpi = 300)

cat("- Gráficos Nivel 1 guardados exitosamente.\n")

# ============================================================================
# NIVEL 2: Intermedio - Testeo Formal Clásico de Modelos Lineales (LM)
# ============================================================================
cat("\n[NIVEL 2] SUPUESTOS CLÁSICOS (HOMOCEDASTICIDAD Y NORMALIDAD)\n")

# Ajuste Modelo Lineal Invernadero
modelo_lm_inv <- lm(biomasa ~ tratamiento, data = datos_inv)
residuos_inv <- residuals(modelo_lm_inv)

# 2.1 Pruebas de Normalidad
cat("\n-> Tests de Normalidad Clásicos:\n")
cat("Shapiro-Wilk:", round(shapiro.test(residuos_inv)$p.value, 4), "\n")
cat("Anderson-Darling:", round(ad.test(residuos_inv)$p.value, 4), "\n")
cat("Lilliefors (K-S):", round(lillie.test(residuos_inv)$p.value, 4), "\n")

# 2.2 Pruebas de Homogeneidad de Varianzas
cat("\n-> Pruebas de Homocedasticidad:\n")
cat("Bartlett (Sensible a normalidad): p =", 
    round(bartlett.test(biomasa ~ tratamiento, data = datos_inv)$p.value, 4), "\n")
cat("Levene (Robusto): p =", 
    round(leveneTest(biomasa ~ tratamiento, data = datos_inv)$"Pr(>F)"[1], 4), "\n")

# 2.3 Diagnóstico visual clásico de R base
png("01_Exploracion_Supuestos/N2_diagnostico_lm_clasico.png", width=800, height=800, res=150)
par(mfrow = c(2,2))
plot(modelo_lm_inv, main = "Diagnóstico Base LM")
par(mfrow = c(1,1))
dev.off()
cat("- Gráficos Nivel 2 guardados.\n")

# ============================================================================
# NIVEL 3: Avanzado - Verificación Easystats y Naturaleza de Distribución
# ============================================================================
cat("\n[NIVEL 3] ANÁLISIS DE PERFORMANCE Y FITDISTRPLUS\n")

# 3.1 Verificación empírica de distribuciones (Cullen & Frey)
cat("-> Evaluando Cullen & Frey para la variable discreta 'Lesiones en Hoja'\n")
png("01_Exploracion_Supuestos/N3_CullenFrey_Conteo.png", width=800, height=600, res=150)
fitdistrplus::descdist(datos_campo$lesiones_hoja, discrete = TRUE, boot = 500)
dev.off()

# 3.2 Visualización Compleja mediante Easystats (Performance)
# Esto evalúa normalidad, colinealidad, homogeneidad y outliers de forma integral
reporte_modelo <- performance::check_model(modelo_lm_inv, panel = FALSE)

# Organizamos el panel personalizado
panel_easystats <- plot(reporte_modelo[[1]]) + plot(reporte_modelo[[2]]) + 
                   plot(reporte_modelo[[3]]) + plot(reporte_modelo[[4]]) +
                   plot_layout(ncol = 2)

ggsave("01_Exploracion_Supuestos/N3_Performance_Easystats.png", panel_easystats, 
       width = 12, height = 10, dpi = 300)
cat("- Panel avanzado Easystats guardado exitosamente.\n")

# ============================================================================
# NIVEL 4: Difícil - GLMs de Conteo, Exceso de Ceros y validación DHARMa
# ============================================================================
cat("\n[NIVEL 4] DATOS DE CAMPO, ZERO-INFLATION Y DHARMa\n")
cat("Contexto: Modelemos el conteo de lesiones asumiendo familia Poisson.\n")

modelo_glm_poisson <- glm(lesiones_hoja ~ manejo, family = poisson, data = datos_campo)

# 4.1 Chequeo rápido Easystats
cat("\n-> Evaluando Sobredispersión (Easystats):\n")
print(check_overdispersion(modelo_glm_poisson))
cat("\n-> Evaluando Inflación de Ceros (Easystats):\n")
print(check_zeroinflation(modelo_glm_poisson))

# 4.2 El estándar de Oro para GLM: DHARMa
# DHARMa simula respuestas del modelo ajustado y compara la distribución empírica
# de los residuos simulados mediante Q-Q plots estandarizados uniformes.
residuos_dharma <- simulateResiduals(fittedModel = modelo_glm_poisson, plot = FALSE)

png("01_Exploracion_Supuestos/N4_DHARMa_Validacion.png", width=1000, height=500, res=150)
plot(residuos_dharma)
dev.off()

# 4.3 Pruebas formales DHARMa en Consola
cat("\n-> Ejecutando pruebas rigurosas DHARMa en consola:\n")
cat("- Test de Dispersión Ajustado DHARMa:\n")
testDispersion(residuos_dharma)

cat("\n- Test de Inflación de Ceros DHARMa:\n")
testZeroInflation(residuos_dharma)

cat("\n============================================================================\n")
cat("EXPLICACIÓN EXPERTA FINAL:\n")
cat("Si DHARMa revela sobredispersión o desviación en el Q-Q plot de uniformidad,\n")
cat("el modelo Poisson no es adecuado. Dependiendo del contexto biotecnológico,\n")
cat("deberás transitar hacia distribuciones Binomial Negativa (MASS::glm.nb),\n")
cat("o Modelos Inflados con Ceros (Zero-Inflated / Hurdle usando glmmTMB o pscl).\n")
cat("La validación de supuestos es el corazón de la estadística aplicada al agro.\n")
cat("============================================================================\n")

# ============================================================================
# NIVEL 5: SUPER AVANZADO - DISEÑO EXPERIMENTAL EN INVERNADERO (LMM)
# ============================================================================
cat("\n[NIVEL 5] MÓDULO DE INVERNADERO (LME4 Y AGRICOLAE)...\n")

# Para este análisis debes instalar lmerTest, agricolae, emmeans
if (!require(lmerTest)) install.packages("lmerTest")
if (!require(agricolae)) install.packages("agricolae")
if (!require(emmeans)) install.packages("emmeans")

library(lme4)
library(lmerTest)
library(agricolae)
library(emmeans)

# Simulación de ensayo de Invernadero (Medidas Repetidas o Bloques)
# 7 Genotipos, 2 Tratamientos (Control vs Calor), 3 Bloques (Mesones del invernadero)
datos_invernadero <- expand.grid(
  genotipo = paste0("G", 1:7),
  tratamiento = c("Control", "Calor"),
  bloque = as.factor(1:3) # Mesón 1, 2 y 3
) %>%
  mutate(
    # Efecto aleatorio del bloque (ej. el mesón 3 es más caluroso)
    efecto_bloque = case_when(bloque == "1" ~ -2, bloque == "2" ~ 0, bloque == "3" ~ 3),
    # Rendimiento simulado
    rendimiento = 50 - (ifelse(tratamiento == "Calor", 15, 0)) + efecto_bloque + rnorm(n(), 0, 2)
  )

# ----------------------------------------------------------------------------
# Enfoque A: DBCA Clásico (Bloques Fijos)
# ----------------------------------------------------------------------------
cat("\n--- Enfoque A: DBCA (ANOVA Clásico) ---\n")
mod_dbca <- lm(rendimiento ~ genotipo * tratamiento + bloque, data = datos_invernadero)
print(anova(mod_dbca))

# Post-Hoc con Agricolae (Test de Tukey clásico en agronomía)
tukey_dbca <- HSD.test(mod_dbca, trt = c("genotipo", "tratamiento"), console = FALSE)
cat("\nGrupos de significancia (Agricolae):\n")
print(head(tukey_dbca$groups))

# ----------------------------------------------------------------------------
# Enfoque B: Modelo Lineal Mixto (LMM con lme4)
# ----------------------------------------------------------------------------
cat("\n--- Enfoque B: Modelo Mixto (Efectos Aleatorios) ---\n")
# '(1 | bloque)' le dice al modelo: "permite que cada bloque tenga su propio intercepto"
mod_mixto <- lmer(rendimiento ~ genotipo * tratamiento + (1 | bloque), data = datos_invernadero)

print(summary(mod_mixto))

# Medias Marginales Estimadas (EMMeans) - El estándar moderno para LMM
medias_mixto <- emmeans(mod_mixto, ~ tratamiento | genotipo)
cat("\nComparación de tratamientos dentro de cada genotipo (EMMeans):\n")
print(pairs(medias_mixto))

cat("\n--- Módulo de Invernadero Finalizado con Éxito ---\n")
