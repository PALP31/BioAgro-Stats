# ============================================================================
# 02_glm_edafofauna.R (GLM: CONTEOS DE ORGANISMOS)
# Modelado de Abundancia de Fauna con Binomial Negativa
# ============================================================================
library(MASS)        # Para Binomial Negativa (glm.nb)
library(tidyverse)
library(performance) # Diagnóstico avanzado
library(DHARMa)      # Validación por residuos simulados (Imprescindible en GLMs)
library(emmeans)     # Medias estimadas en escala real

set.seed(456)

# 1. SIMULACIÓN DE DATOS: ABUNDANCIA DE EDAFOFAUNA
# Escenario: 3 tratamientos (Control, Bioestimulante_A, Bioestimulante_B)
# Aplicación en huertos con 30 parcelas totales.
tratamientos <- c("Control", "Bio_A", "Bio_B")
datos_fauna <- data.frame(
  trata = factor(rep(tratamientos, each = 10)),
  abundancia = c(
    rnbinom(10, mu = 12, size = 1.2),  # Control
    rnbinom(10, mu = 45, size = 1.2),  # Bio_A (Fuerte mejora)
    rnbinom(10, mu = 28, size = 1.2)   # Bio_B
  )
)

# 2. AJUSTE DE MODELOS (POISSON VS BINOMIAL NEGATIVA)
# Probamos Poisson primero
mod_poisson <- glm(abundancia ~ trata, data = datos_fauna, family = poisson)

# Validamos sobredispersión (performance::check_overdispersion)
cat("\n--- [1] VALIDACIÓN DE SOBREDISPERSIÓN (POISSON) ---\n")
print(check_overdispersion(mod_poisson))

# Ajustamos Binomial Negativa (Maneja varianza > media)
mod_nb <- glm.nb(abundancia ~ trata, data = datos_fauna)

# 3. DIAGNÓSTICO AVANZADO CON DHARMA
cat("\n--- [2] DIAGNÓSTICO CON DHARMA (RESIDUOS SIMULADOS) ---\n")
sim_nb <- simulateResiduals(mod_nb)
# plot(sim_nb) # Visualización integral
print(testDispersion(sim_nb))

# 4. COMPARACIONES MÚLTIPLES (ESCALE="RESPONSE")
cat("\n--- [3] MEDIAS ESTIMADAS (TRANSFORMACIÓN A ESCALA NATURAL) ---\n")
emm_fauna <- emmeans(mod_nb, ~ trata, type = "response")
print(emm_fauna)
letras_fauna <- cld(emm_fauna, Letters = letters, adjust = "tukey")

# 5. GRÁFICO PROFESIONAL: VIOLIN PLOT + EMMEANS
ggplot(datos_fauna, aes(x = trata, y = abundancia, fill = trata)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.4) +
  geom_point(data = as.data.frame(emm_fauna), aes(x = trata, y = response), 
             size = 4, color = "black") +
  geom_errorbar(data = as.data.frame(emm_fauna), 
                aes(x = trata, y = response, ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, color = "black") +
  labs(title = "Abundancia de Edafofauna en Huertos",
       subtitle = "Modelado con GLM Binomial Negativa (emmeans + IC 95%)",
       x = "Tratamiento de Bioestimulante", y = "Conteo de Individuos") +
  theme_minimal() +
  scale_fill_viridis_d() +
  theme(legend.position = "none")

cat("\n--- Script 02 finalizado correctamente ---\n")
