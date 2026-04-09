# ============================================================================
# 01_anova_trigo.R (ESTRÉS SALINO X TÉRMICO)
# ANOVA de 2 Vías con Bloques Completos al Azar (DBCA)
# ============================================================================
library(tidyverse)
library(emmeans)
library(performance)
library(car)
library(multcomp)
library(multcompView)

set.seed(789)

# 1. SIMULACIÓN DE DATOS: BIOMASA DE TRIGO (g)
# Factores: Salinidad (Control, 100mM NaCl), Calor (25°C, 38°C), Bloque (4)
datos_trigo <- expand.grid(
  sal = factor(c("Control", "Salino")),
  heat = factor(c("Optimo", "Calor")),
  block = factor(paste0("B", 1:4))
) %>%
  mutate(
    yield_base = case_when(
      sal == "Control" & heat == "Optimo" ~ 35,
      sal == "Salino" & heat == "Optimo" ~ 28,
      sal == "Control" & heat == "Calor" ~ 25,
      sal == "Salino" & heat == "Calor" ~ 15  # Sinergia negativa de estreses
    ),
    biom = yield_base + as.numeric(block)*0.5 + rnorm(n(), 0, 1.5)
  )

# 2. AJUSTE DEL MODELO (ANOVA 2 VÍAS + BLOQUES)
# Formula: Y ~ Bloque + A * B
mod_trigo <- aov(biom ~ block + sal * heat, data = datos_trigo)

cat("\n--- [1] TABLA ANOVA (INTERACCIÓN SALINIDAD X CALOR) ---\n")
print(summary(mod_trigo))

# 3. EVALUACIÓN DE SUPUESTOS (Easystats)
cat("\n--- [2] DIAGNÓSTICO DE SUPUESTOS ---\n")
print(check_model(mod_trigo))

# 4. COMPARACIONES MÚLTIPLES (EMMEANS)
cat("\n--- [3] COMPARACIONES PAIRWISE (INTERACCIÓN) ---\n")
emm_inter <- emmeans(mod_trigo, ~ sal * heat)
letras_inter <- cld(emm_inter, Letters = letters, adjust = "tukey")
print(letras_inter)

# 5. GRÁFICO PROFESIONAL: GRUPOS POR INTERACCIÓN
ggplot(datos_trigo, aes(x = heat, y = biom, fill = sal)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(data = letras_inter, aes(y = emmean + (2.5*SE), label = .group), 
            position = position_dodge(0.75), fontface = "bold") +
  scale_fill_manual(values = c("Control" = "#2E86AB", "Salino" = "#E74C3C")) +
  labs(title = "Biomasa de Trigo bajo Estrés Combinado",
       subtitle = "Interacción Salinidad x Calor (ANOVA 2 Vías + Bloques)",
       x = "Condición Térmica", y = "Biomasa (g)", fill = "Salinidad") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

cat("\n--- Script 01 finalizado correctamente ---\n")
