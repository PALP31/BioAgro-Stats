# ============================================================================
# 02_alpha_lattice_mejoramiento.R (VERSIÓN EXPERTA - DHARMA & MIXED)
# Diseño Alpha-Lattice (GxE y Estrés Térmico)
# ============================================================================
# Paquetes: agricolae, tidyverse, lme4, lmerTest, performance, DHARMa, nortest
# ============================================================================

library(agricolae)
library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(DHARMa) # Diagnóstico avanzado con residuos simulados
library(nortest)

set.seed(123)

# 1. GENERACIÓN DEL DISEÑO EXPERIMENTAL (REGLA ESTRICTA)
genotipos <- paste0("G_", sprintf("%02d", 1:16))
k <- 4; r <- 2
design_alpha <- design.alpha(trt = genotipos, k = k, r = r, serie = 2, seed = 123)
datos_alpha <- design_alpha$book

# 2. EJERCICIO AVANZADO: MULTI-AMBIENTAL (MET) Y ESTRÉS TÉRMICO
datos_met <- bind_rows(
  datos_alpha %>% mutate(env = "Optimo", ef_env = 0),
  datos_alpha %>% mutate(env = "Calor_Extremo", ef_env = -1400)
)

# Efectos GxE Dinámicos
datos_met <- datos_met %>%
  mutate(
    yield_base = 5200 + ef_env + rnorm(16, 0, 500)[as.numeric(genotipos)],
    yield = yield_base + rnorm(n(), 0, 200)
  )

# 3. MODELO MIXTO DE ALTO NIVEL
mod_lmer <- lmer(yield ~ env + (1|genotipos) + (1|genotipos:env) + 
                   (1|env:replication) + (1|env:replication:block), data = datos_met)

# 4. DIAGNÓSTICO TOP TIER CON DHARMA (RESIDUOS SIMULADOS)
cat("\n", rep("=", 60), "\n")
cat("DIAGNÓSTICO AVANZADO: MODELO MIXTO (DHARMa)\n")
cat(rep("=", 60), "\n")

# --- 4.1 Generar Residuos Simulados ---
sim_res <- simulateResiduals(fittedModel = mod_lmer, n = 1000)

# --- 4.2 Pruebas de Diagnóstico DHARMa ---
cat("\n[1] TEST DE UNIFORMIDAD (KS Test):\n")
print(testUniformity(sim_res))

cat("\n[2] TEST DE SOBREDISPERSIÓN:\n")
print(testDispersion(sim_res))

cat("\n[3] TEST DE OUTLIERS:\n")
print(testOutliers(sim_res))

# --- 4.3 Heterocedasticidad por Variable (GxE) ---
cat("\n[4] HETEROCEDASTICIDAD POR AMBIENTE (ENV):\n")
print(testCategorical(sim_res, catPred = datos_met$env))

# 5. DIAGNÓSTICO CLÁSICO ADICIONAL
res_lmer <- residuals(mod_lmer)
cat("\n[5] NORMALIDAD (Anderson-Darling): p =", ad.test(res_lmer)$p.value, "\n")

# 6. VISUALIZACIÓN DHARMA
# plot(sim_res) # Esto abriría una ventana gráfica con 4 paneles diagnósticos

# 7. PARÁMETROS GENÉTICOS (Heredabilidad)
cat("\n--- Componentes de Varianza ---\n")
var_comp <- as.data.frame(VarCorr(mod_lmer))
print(var_comp)

cat("\n--- Script Alpha-Lattice Experto Finalizado ---\n")
