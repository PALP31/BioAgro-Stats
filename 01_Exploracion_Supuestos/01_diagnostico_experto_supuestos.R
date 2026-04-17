# ============================================================================
# 01_diagnostico_experto_supuestos.R (GUÍA DEFINITIVA - ACTUALIZADA)
# Evaluación de Supuestos: Estadística Clásica y Moderna (Easystats)
# ============================================================================
# Este script es una herramienta integral para validar modelos antes de su
# publicación. Cubre datos gaussianos (Biomasa) y de conteo (Edafofauna).
# ============================================================================

library(tidyverse)
library(performance) # El corazón de easystats para diagnósticos
library(see)         # Visualización de easystats
library(car)         # Pruebas clásicas (Levene)
library(lmtest)      # Breusch-Pagan, Durbin-Watson
library(patchwork)   # Combinación de gráficos para papers
library(DHARMa)      # Residuos simulados para modelos de conteo
library(nortest)     # Anderson-Darling para normalidad

set.seed(123)

# ============================================================================
# 1. SIMULACIÓN DE DATOS MEJORADA
# ============================================================================

# Caso A: Biomasa de Trigo (Datos Continuos - Gaussianos)
# Agregamos 'temperatura' como covariable continua
datos_trigo <- data.frame(
  trata = factor(rep(c("Control", "Estres"), each = 15)),
  temperatura = runif(30, 15, 35),
  biom = NA
)

# Simulación con efecto de temperatura y trata + outlier
datos_trigo <- datos_trigo %>%
  mutate(
    yield_base = ifelse(trata == "Control", 35, 22),
    biom = yield_base + (temperatura * 0.2) + rnorm(n(), 0, 3)
  )
datos_trigo$biom[5] <- 55 # Outlier artificial

# Caso B: Conteo de Insectos (Poisson con Exceso de Ceros y Sobredispersión)
datos_insectos <- data.frame(
  trata = factor(rep(c("Bosque", "Cultivo"), each = 20)),
  conteo = c(rpois(20, 1.5), rpois(20, 15))
)
# Inducir exceso de ceros y sobredispersión
datos_insectos$conteo[1:10] <- 0
datos_insectos$conteo[30:40] <- datos_insectos$conteo[30:40] * 3

# ============================================================================
# 2. PARTE 1: MODELO GAUSSAINO (BIOMASA)
# ============================================================================

cat("\n", rep("=", 60), "\n")
cat("CASO 1: MODELO LINEAL MULTIPLE (BIOMASA ~ TRATA + TEMP)\n")
cat(rep("=", 60), "\n")

mod_gauss <- lm(biom ~ trata + temperatura, data = datos_trigo)
residuos <- residuals(mod_gauss)

# ============================================================================
# PARTE 2: DIAGNÓSTICO VISUAL Y CLÁSICO (Paso a Paso)
# ============================================================================

cat("\n--- [A] DIAGNÓSTICO VISUAL BÁSICO ---\n")

# a) Histograma de Residuos
hist(residuos, main = "Histograma de Residuos", xlab = "Residuos", col = "skyblue")

# b) Q-Q Plot Clásico
qqnorm(residuos, main = "Q-Q Plot de Residuos")
qqline(residuos, col = "red", lwd = 2)

# c) Dispersión: Residuos vs Temperatura (Variable Continua)
plot(datos_trigo$temperatura, residuos, 
     main = "Residuos vs Temperatura", 
     xlab = "Temperatura (°C)", ylab = "Residuos", pch = 19)
abline(h = 0, lty = 2, col = "blue")

cat("\n--- [B] PRUEBAS FORMALES CLÁSICAS ---\n")

# Shapiro-Wilk: Evalúa la Normalidad de los residuos
# H0: Los residuos siguen una distribución normal
shapiro_test <- shapiro.test(residuos)
cat("- Shapiro-Wilk (Normalidad): p =", round(shapiro_test$p.value, 5), "\n")

# LeveneTest: Evalúa Homocedasticidad para el factor 'trata'
# H0: Las varianzas entre grupos son iguales
levene_test <- leveneTest(biom ~ trata, data = datos_trigo)
cat("- Levene Test (Varianza entre grupos): p =", round(levene_test$`Pr(>F)`[1], 5), "\n")

# Breusch-Pagan: Evalúa Heterocedasticidad contra la variable continua 'temperatura'
# H0: La varianza de los residuos es constante (Homocedasticidad)
bp_test <- bptest(mod_gauss)
cat("- Breusch-Pagan (Heterocedasticidad vs Temp): p =", round(bp_test$p.value, 5), "\n")

# ============================================================================
# PARTE 3: DIAGNÓSTICO MODERNO (performance / easystats)
# ============================================================================

cat("\n", rep("=", 60), "\n")
cat("DIAGNÓSTICO MODERNO (EASYSTATS)\n")
cat(rep("=", 60), "\n")

# Detección de outliers
outliers <- check_outliers(mod_gauss)
print(outliers)

# Chequeo integral (Genera múltiples gráficos diagnósticos)
diag_gauss <- check_model(mod_gauss, panel = FALSE)

# ============================================================================
# PARTE 4: MODELOS DE CONTEO (INSECTOS)
# ============================================================================

cat("\n", rep("=", 60), "\n")
cat("CASO 2: MODELO DE CONTEO (INSECTOS - POISSON)\n")
cat(rep("=", 60), "\n")

mod_count <- glm(conteo ~ trata, family = poisson, data = datos_insectos)

# Sobredispersión y Exceso de Ceros (easystats)
cat("\n[3] DIAGNÓSTICO DE CONTEO:\n")
overdisp <- check_overdispersion(mod_count)
print(overdisp)

zero_inf <- check_zeroinflation(mod_count)
print(zero_inf)

# DHARMa (Residuos Simulados) - Imprescindible para GLMs de conteo
res_sim <- simulateResiduals(mod_count)
cat("\n[4] DHARMA TEST (Residuos Simulados):\n")
testDispersion(res_sim)

# ============================================================================
# PARTE 5: PANEL DE GRÁFICOS PARA PAPER (Patchwork)
# ============================================================================

cat("\n--- Generando Panel Gráfico para Publicación ---\n")

# Gráfico 1: Diagnóstico de Normalidad (performance)
p1 <- plot(check_normality(mod_gauss)) + 
  labs(title = "A. Normalidad de Residuos (Trigo)") +
  theme_minimal()

# Gráfico 2: Diagnóstico de Homocedasticidad (performance)
p2 <- plot(check_heteroscedasticity(mod_gauss)) + 
  labs(title = "B. Homocedasticidad (Trigo)") +
  theme_minimal()

# Gráfico 3: Outliers
p3 <- plot(outliers) + 
  labs(title = "C. Detección de Outliers (Cook's D)") +
  theme_minimal()

# Gráfico 4: Sobredispersión en Conteo
p4 <- plot(overdisp) + 
  labs(title = "D. Check de Sobredispersión (Insectos)") +
  theme_minimal()

# Combinación con Patchwork
panel_diagnostico <- (p1 + p2) / (p3 + p4) + 
  plot_annotation(
    title = "Panel de Diagnóstico de Supuestos para Publicación",
    subtitle = "Integración de Métodos Clásicos, Easystats y Patchwork",
    caption = "Generado con el ecosistema easystats y patchwork"
  )

# Guardar el panel
ggsave("01_Exploracion_Supuestos/panel_diagnostico_experto.png", 
       panel_diagnostico, width = 12, height = 10, dpi = 300)

print(panel_diagnostico)

cat("\n--- Script de Diagnóstico Experto finalizado ---\n")
cat("Gráfico guardado en: 01_Exploracion_Supuestos/panel_diagnostico_experto.png\n")


# ============================================================================
# 5. MÓDULO AVANZADO: DISEÑO EXPERIMENTAL EN INVERNADERO (DBCA Y MIXTOS)
# ============================================================================
cat("\n[5] INICIANDO MÓDULO DE INVERNADERO (LME4 Y AGRICOLAE)...\n")

library(lme4)
library(lmerTest)
library(agricolae)
library(emmeans)

# Simulación de tu próximo ensayo de Invernadero (Abril)
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
# Útil si solo tienes 3 mesones y quieres inferir solo sobre esos 3.
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
# Útil si los bloques representan una muestra aleatoria de microclimas, o si
# tienes parcelas divididas (Split-Plot) o medidas repetidas en el tiempo.
# ----------------------------------------------------------------------------
cat("\n--- Enfoque B: Modelo Mixto (Efectos Aleatorios) ---\n")
# '(1 | bloque)' le dice al modelo: "permite que cada bloque tenga su propio intercepto (nivel base)"
mod_mixto <- lmer(rendimiento ~ genotipo * tratamiento + (1 | bloque), data = datos_invernadero)

# Resumen del modelo mixto (nota cómo la varianza se particiona)
print(summary(mod_mixto))

# Medias Marginales Estimadas (EMMeans) - El estándar moderno para post-hoc en LMM
medias_mixto <- emmeans(mod_mixto, ~ tratamiento | genotipo)
cat("\nComparación de tratamientos dentro de cada genotipo (EMMeans):\n")
print(pairs(medias_mixto))


cat("\n--- Módulo de Invernadero Finalizado ---\n")



