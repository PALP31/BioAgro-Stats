# ============================================================================
# 01_diagnostico_experto_supuestos.R (MASTERCLASS DE SUPUESTOS)
# ============================================================================
# Experto: Estadística y Diseño Experimental
# Objetivo: Guía progresiva (Sencillo -> Difícil) para la verificación de
# supuestos estadísticos en investigación agrícola, biotecnológica y de campo.
# ============================================================================

# Verificación de paquetes requeridos
paquetes_requeridos <- c(
  "tidyverse", "performance", "see", "car",
  "patchwork", "DHARMa", "nortest", "corrplot", "fitdistrplus", "lme4",
  "lmerTest", "agricolae", "emmeans"
)
paquetes_faltantes <- paquetes_requeridos[!(paquetes_requeridos %in% installed.packages()[, "Package"])]
if (length(paquetes_faltantes)) {
  stop(
    paste0(
      "Faltan paquetes requeridos para ejecutar este script: ",
      paste(paquetes_faltantes, collapse = ", "),
      ".\nInstálelos manualmente antes de continuar, por ejemplo con:\n",
      "install.packages(c(",
      paste(sprintf('"%s"', paquetes_faltantes), collapse = ", "),
      "), dependencies = TRUE)"
    ),
    call. = FALSE
  )
}

# Directorio único de salida para figuras y archivos generados
output_dir <- "01_Exploracion_Supuestos"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cargar librerías
library(tidyverse)
library(performance) # Diagnósticos modernos y visuales
library(see)         # Visualización de ecosistema easystats
library(car)         # Pruebas clásicas de supuestos
library(patchwork)   # Combinación de gráficos complejos
library(DHARMa)      # Diagnóstico avanzado de residuos simulados
library(nortest)     # Evaluación exhaustiva de normalidad
library(corrplot)    # Correlogramas
library(fitdistrplus)# Ajuste experto de distribuciones estadísticas
library(lme4)        # Modelamiento mixto
library(lmerTest)    # P-valores para modelos mixtos
library(agricolae)   # Pruebas post-hoc en agronomía (HSD.test)
library(emmeans)     # Medias marginales estimadas

set.seed(42) # Reproducibilidad científica

# Creación de directorio de salida si no existe (para guardar los gráficos)
# Evita duplicar "01_Exploracion_Supuestos" cuando el script se ejecuta
# con el working directory ya situado dentro de esa carpeta.
current_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
output_dir <- if (basename(current_dir) == "01_Exploracion_Supuestos") {
  current_dir
} else {
  file.path(current_dir, "01_Exploracion_Supuestos")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# 1. SIMULACIÓN PARAMETRIZADA (Enfoque Agrícola/Biotecnológico)
# ============================================================================
cat("\n[1] GENERANDO MATRICES DE DATOS SIMULADOS...\n")

# Caso A: Laboratorio (Enzimología y Multi-rasgo Continuo)
datos_lab <- data.frame(
  genotipo = factor(rep(c("Silvestre", "Crispr"), each = 30)),
  clorofila = rnorm(60, mean = 45, sd = 4),
  nitrogeno_foliar = rnorm(60, mean = 2.5, sd = 0.3)
) %>%
  mutate(
    act_enzimatica = 15 + 0.6 * clorofila + 3 * nitrogeno_foliar + 
                     ifelse(genotipo == "Crispr", 8, 0) + rnorm(60, 0, 3)
  )
datos_lab$act_enzimatica[c(10, 50)] <- c(80, 5) # Outliers severos

# Caso B: Invernadero (Biomasa con Heterocedasticidad)
datos_inv <- data.frame(
  tratamiento = factor(rep(c("T0_Control", "T1_Sequia", "T2_Salinidad"), each = 20)),
  biomasa = c(
    rnorm(20, mean = 60, sd = 3),   
    rnorm(20, mean = 35, sd = 9),   
    rnorm(20, mean = 25, sd = 15)   
  )
)

# Caso C: Campo (Incidencia de Fitopatógenos - Conteo)
datos_campo <- data.frame(
  manejo = factor(rep(c("Testigo", "Fungicida_Sintetico", "Trichoderma"), each = 40)),
  lesiones_hoja = c(rpois(40, lambda = 12), rpois(40, lambda = 1), rpois(40, lambda = 3))
)
datos_campo$lesiones_hoja[sample(1:120, 25)] <- 0 # Exceso de ceros

# ============================================================================
# NIVEL 1: Sencillo - Exploración Visual, Correlaciones y Outliers
# ============================================================================
cat("\n[NIVEL 1] INICIO DE ANÁLISIS SENCILLO (EXPLORATORIO)\n")

vars_num_lab <- datos_lab %>% select(clorofila, nitrogeno_foliar, act_enzimatica)
matriz_corr <- cor(vars_num_lab)
cat("-> Matriz de Correlación (Laboratorio):\n")
print(matriz_corr)

png(file.path(output_dir, "N1_correlaciones_lab.png"), width=800, height=800, res=150)
corrplot(matriz_corr, method = "color", type = "upper", addCoef.col = "black",
         tl.col = "darkblue", tl.srt = 45, diag = FALSE,
         title = "Correlación Traits Biotecnológicos", mar=c(0,0,2,0))
dev.off()

p_hist <- ggplot(datos_lab, aes(x = act_enzimatica, fill = genotipo)) +
  geom_histogram(bins = 20, alpha = 0.6, color = "black", position = "identity") +
  theme_minimal() + labs(title = "Distribución Enzimática", y = "Frecuencia")

p_box <- ggplot(datos_lab, aes(x = genotipo, y = act_enzimatica, fill = genotipo)) +
  geom_boxplot(outlier.colour = "red", outlier.size = 3) +
  theme_minimal() + labs(title = "Detección de Outliers (Puntos Rojos)")

ggsave(file.path(output_dir, "N1_exploracion_visual.png"), plot = (p_hist / p_box), width = 8, height = 7, dpi = 300)
cat("- Gráficos Nivel 1 guardados exitosamente.\n")

# ============================================================================
# NIVEL 2: Intermedio - Testeo Formal Clásico de Modelos Lineales (LM)
# ============================================================================
cat("\n[NIVEL 2] SUPUESTOS CLÁSICOS (HOMOCEDASTICIDAD Y NORMALIDAD)\n")

modelo_lm_inv <- lm(biomasa ~ tratamiento, data = datos_inv)
residuos_inv <- residuals(modelo_lm_inv)

cat("\n-> Tests de Normalidad Clásicos:\n")
cat("Shapiro-Wilk (p):", round(shapiro.test(residuos_inv)$p.value, 4), "\n")
cat("Anderson-Darling (p):", round(ad.test(residuos_inv)$p.value, 4), "\n")
cat("Lilliefors K-S (p):", round(lillie.test(residuos_inv)$p.value, 4), "\n")

cat("\n-> Pruebas de Homocedasticidad:\n")
cat("Bartlett (Sensible a normalidad): p =", round(bartlett.test(biomasa ~ tratamiento, data = datos_inv)$p.value, 4), "\n")
cat("Levene (Robusto): p =", round(leveneTest(biomasa ~ tratamiento, data = datos_inv)$"Pr(>F)"[1], 4), "\n")

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

cat("-> Evaluando Cullen & Frey para 'Lesiones en Hoja'\n")
png("01_Exploracion_Supuestos/N3_CullenFrey_Conteo.png", width=800, height=600, res=150)
fitdistrplus::descdist(datos_campo$lesiones_hoja, discrete = TRUE, boot = 500)
dev.off()

reporte_modelo <- performance::check_model(modelo_lm_inv, panel = FALSE)
panel_easystats <- reporte_modelo[[1]] + reporte_modelo[[2]] + 
                   reporte_modelo[[3]] + reporte_modelo[[4]] +
                   plot_layout(ncol = 2)

ggsave("01_Exploracion_Supuestos/N3_Performance_Easystats.png", panel_easystats, width = 12, height = 10, dpi = 300)
cat("- Panel avanzado Easystats guardado exitosamente.\n")

# ============================================================================
# NIVEL 4: Difícil - GLMs de Conteo, Exceso de Ceros y validación DHARMa
# ============================================================================
cat("\n[NIVEL 4] DATOS DE CAMPO, ZERO-INFLATION Y DHARMa\n")
cat("Contexto: Modelemos el conteo de lesiones asumiendo familia Poisson.\n")

modelo_glm_poisson <- glm(lesiones_hoja ~ manejo, family = poisson, data = datos_campo)

cat("\n-> Evaluando Sobredispersión (Easystats):\n")
print(check_overdispersion(modelo_glm_poisson))
cat("\n-> Evaluando Inflación de Ceros (Easystats):\n")
print(check_zeroinflation(modelo_glm_poisson))

residuos_dharma <- simulateResiduals(fittedModel = modelo_glm_poisson, plot = FALSE)
png("01_Exploracion_Supuestos/N4_DHARMa_Validacion.png", width=1000, height=500, res=150)
plot(residuos_dharma)
dev.off()

cat("\n-> Pruebas formales DHARMa:\n")
testDispersion(residuos_dharma)
testZeroInflation(residuos_dharma)
cat("============================================================================\n")

# ============================================================================
# NIVEL 5: SUPER AVANZADO - DISEÑO EXPERIMENTAL EN INVERNADERO (LMM)
# ============================================================================
cat("\n[NIVEL 5] MÓDULO DE INVERNADERO (LME4 Y AGRICOLAE)...\n")

# Simulación de ensayo de Invernadero (Bloques / Medidas Repetidas)
datos_invernadero <- expand.grid(
  genotipo = paste0("G", 1:7),
  tratamiento = c("Control", "Calor"),
  bloque = as.factor(1:3) 
) %>%
  mutate(
    efecto_bloque = case_when(bloque == "1" ~ -2, bloque == "2" ~ 0, bloque == "3" ~ 3),
    rendimiento = 50 - (ifelse(tratamiento == "Calor", 15, 0)) + efecto_bloque + rnorm(n(), 0, 2)
  )

# --- DBCA Clásico (Bloques Fijos) ---
cat("\n--- Enfoque A: DBCA (ANOVA Clásico) ---\n")
mod_dbca <- lm(rendimiento ~ genotipo * tratamiento + bloque, data = datos_invernadero)
print(anova(mod_dbca))

tukey_dbca <- HSD.test(mod_dbca, trt = c("genotipo", "tratamiento"), console = FALSE)
cat("\nGrupos de significancia (Agricolae):\n")
print(head(tukey_dbca$groups))

# --- Modelo Lineal Mixto (Efectos Aleatorios) ---
cat("\n--- Enfoque B: Modelo Mixto con lme4 ---\n")
mod_mixto <- lmer(rendimiento ~ genotipo * tratamiento + (1 | bloque), data = datos_invernadero)
print(summary(mod_mixto))

cat("\nComparación post-hoc (EMMeans):\n")
medias_mixto <- emmeans(mod_mixto, ~ tratamiento | genotipo)
print(pairs(medias_mixto))

cat("\n--- MASTERCLASS Y MÓDULO DE INVERNADERO FINALIZADOS CON ÉXITO ---\n")
