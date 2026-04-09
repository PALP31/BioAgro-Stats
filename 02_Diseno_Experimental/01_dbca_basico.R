# ============================================================================
# 01_dbca_basico.R (VERSIÓN EXPERTA - BATERÍA COMPLETA)
# Diseño de Bloques Completos al Azar (DBCA/RCBD)
# ============================================================================
# Paquetes: agricolae, tidyverse, performance, car, emmeans, nlme, nortest, lmtest
# ============================================================================

library(agricolae)
library(tidyverse)
library(performance)
library(car)
library(emmeans)
library(nlme)
library(nortest) # Anderson-Darling
library(lmtest)  # Breusch-Pagan

set.seed(123)

# 1. GENERACIÓN DEL DISEÑO EXPERIMENTAL (REGLA ESTRICTA)
tratamientos <- c("Control", "Bio_A", "Bio_B", "Bio_C", "Bio_D")
r <- 4
design_dbca <- design.rcbd(trt = tratamientos, r = r, serie = 2, seed = 123)
datos <- design_dbca$book

# 2. ESCENARIO AVANZADO: ESTRÉS SALINO Y PGPR
# Simulación de biomasa con heterocedasticidad inducida por el estrés
datos <- datos %>%
  mutate(
    ef_trt = case_when(
      tratamientos == "Control" ~ 12,
      tratamientos == "Bio_A" ~ 15,
      tratamientos == "Bio_B" ~ 21,
      tratamientos == "Bio_C" ~ 18,
      tratamientos == "Bio_D" ~ 14
    ),
    sd_error = ifelse(tratamientos %in% c("Control", "Bio_D"), 3.0, 0.8),
    biomasa = ef_trt + as.numeric(block)*0.4 + rnorm(n(), 0, sd_error)
  )

# 3. AJUSTE DE MODELOS
mod_linear <- lm(biomasa ~ block + tratamientos, data = datos)

# 4. BATERÍA COMPLETA DE SUPUESTOS (TOP TIER)
cat("\n", rep("=", 60), "\n")
cat("DIAGNÓSTICO EXHAUSTIVO DE SUPUESTOS\n")
cat(rep("=", 60), "\n")

# --- 4.1 NORMALIDAD (Múltiples Enfoques) ---
res <- residuals(mod_linear)
cat("\n[1] NORMALIDAD:\n")
cat("- Shapiro-Wilk (Clásico): p =", shapiro.test(res)$p.value, "\n")
cat("- Anderson-Darling (Potente): p =", ad.test(res)$p.value, "\n")
cat("- Lilliefors (Corrección KS): p =", lillie.test(res)$p.value, "\n")

# --- 4.2 HOMOCEDASTICIDAD (Varianza Constante) ---
cat("\n[2] HOMOCEDASTICIDAD:\n")
cat("- Levene Test (Robusto): p =", leveneTest(biomasa ~ tratamientos, data = datos)$`Pr(>F)`[1], "\n")
cat("- Breusch-Pagan (Linealidad): p =", bptest(mod_linear)$p.value, "\n")
cat("- Bartlett (Sensible): p =", bartlett.test(biomasa ~ tratamientos, data = datos)$p.value, "\n")

# --- 4.3 INDEPENDENCIA Y AUTOCORRELACIÓN ---
cat("\n[3] INDEPENDENCIA:\n")
cat("- Durbin-Watson: DW =", durbinWatsonTest(mod_linear)$dw, ", p =", durbinWatsonTest(mod_linear)$p, "\n")

# --- 4.4 OUTLIERS Y PUNTOS INFLUYENTES ---
cat("\n[4] DIAGNÓSTICO DE INFLUENCIA:\n")
cooks <- cooks.distance(mod_linear)
cat("- Máxima Distancia de Cook:", max(cooks), "\n")
if(max(cooks) > 1) cat("  ⚠ ALERTA: Posible punto influyente detectado.\n")

# 5. SOLUCIÓN AVANZADA CON GLS (Weights por Tratamiento)
mod_gls <- gls(biomasa ~ tratamientos + block, 
               data = datos,
               weights = varIdent(form = ~ 1 | tratamientos))

cat("\n--- Comparación AIC (LM vs GLS) ---\n")
cat("AIC LM:", AIC(mod_linear), "\n")
cat("AIC GLS (Varianza Flexible):", AIC(mod_gls), "\n")

# 6. MEDIAS ROBUSTAS
cat("\n--- Comparaciones Post-Hoc (Basadas en GLS) ---\n")
print(pairs(emmeans(mod_gls, ~ tratamientos), adjust = "tukey"))

cat("\n--- Script Finalizado (Full Assumptions) ---\n")
