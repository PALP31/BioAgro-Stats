# ============================================================================
# 03_parcelas_divididas_splitplot.R (VERSIÓN EXPERTA - HIERARCHICAL DIAG)
# Diseño de Parcelas Sub-subdivididas (Split-Split-Plot)
# ============================================================================
# Paquetes: agricolae, tidyverse, lme4, lmerTest, performance, DHARMa, emmeans
# ============================================================================

library(agricolae)
library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(DHARMa)
library(emmeans)

set.seed(123)

# 1. GENERACIÓN DEL DISEÑO EXPERIMENTAL (REGLA ESTRICTA)
riegos <- c("Control", "Drought")
micorr <- c("No_HMA", "HMA")
r <- 3
design_sp <- design.split(trt1 = riegos, trt2 = micorr, r = r, 
                          design = "rcbd", serie = 2, seed = 123)
datos_sp <- design_sp$book

# 2. ESCALANDO A SPLIT-SPLIT-PLOT (Riego x Micorrizas x NPK)
npk_levels <- c("0%", "50%", "100%")
datos_ssp <- datos_sp %>%
  group_by(plots) %>%
  slice(rep(1, each = length(npk_levels))) %>%
  mutate(npk = rep(npk_levels, length.out = n())) %>%
  ungroup() %>%
  mutate(
    id_main = factor(paste0(block, "_", riegos)),
    id_sub = factor(paste0(id_main, "_", micorr))
  )

# Simulación de datos con variabilidad jerárquica compleja
datos_ssp <- datos_ssp %>%
  mutate(
    yield_base = case_when(
      riegos == "Control" ~ 52 + as.numeric(gsub("%", "", npk))*0.25,
      riegos == "Drought" & micorr == "HMA" ~ 48 + as.numeric(gsub("%", "", npk))*0.2,
      TRUE ~ 32 + as.numeric(gsub("%", "", npk))*0.15
    ),
    err_main = rnorm(n(), 0, 4)[as.numeric(id_main)],
    err_sub = rnorm(n(), 0, 2.5)[as.numeric(id_sub)],
    err_resid = rnorm(n(), 0, 1.2),
    biomasa = yield_base + err_main + err_sub + err_resid
  )

# 3. ANÁLISIS DE PARCELAS SUB-SUBDIVIDIDAS (MODELO JERÁRQUICO)
mod_ssp <- lmer(biomasa ~ riegos * micorr * npk + 
                  (1|block) + (1|block:riegos) + (1|block:riegos:micorr), 
                data = datos_ssp)

# 4. DIAGNÓSTICO TOP TIER: SUPUESTOS JERÁRQUICOS
cat("\n", rep("=", 60), "\n")
cat("DIAGNÓSTICO JERÁRQUICO DE SUPUESTOS\n")
cat(rep("=", 60), "\n")

# --- 4.1 Evaluación por Niveles con DHARMa ---
sim_ssp <- simulateResiduals(mod_ssp)

cat("\n[1] VALIDEZ GLOBAL DEL MODELO (DHARMa):\n")
print(testUniformity(sim_ssp))

cat("\n[2] SUPUESTO DE INDEPENDENCIA EN PARCELA PRINCIPAL (Riego):\n")
# Recalcular residuos para el factor principal
res_main <- recalculateResiduals(sim_ssp, group = datos_ssp$id_main)
print(testUniformity(res_main))

cat("\n[3] SUPUESTO DE INDEPENDENCIA EN SUBPARCELA (Micorriza):\n")
res_sub <- recalculateResiduals(sim_ssp, group = datos_ssp$id_sub)
print(testUniformity(res_sub))

# --- 4.2 Heterocedasticidad entre Factores Críticos ---
cat("\n[4] HETEROCEDASTICIDAD ENTRE NIVELES DE RIEGO:\n")
print(testCategorical(sim_ssp, catPred = datos_ssp$riegos))

# 5. RESUMEN DE COMPONENTES DE VARIANZA (ERRORES A, B, C)
cat("\n--- Componentes de Error (Varianza Aleatoria) ---\n")
var_ssp <- as.data.frame(VarCorr(mod_ssp))
print(var_ssp)

# 6. COMPARACIONES TRIPLES DE ALTO NIVEL
cat("\n--- Efecto Sinergia HMA x NPK bajo Sequía ---\n")
print(pairs(emmeans(mod_ssp, ~ micorr | riegos * npk), adjust = "tukey"))

cat("\n--- Script Split-Split-Plot Experto Finalizado ---\n")
