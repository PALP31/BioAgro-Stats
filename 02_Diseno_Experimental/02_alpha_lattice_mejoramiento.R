# ============================================================================
# 02_alpha_lattice_mejoramiento.R (GxE Y ESTRÉS TÉRMICO)
# Diseño Alpha-Lattice (Bloques Incompletos)
# ============================================================================
# Este diseño es vital en mejoramiento cuando evaluamos muchos genotipos.
# Un bloque completo sería muy grande y heterogéneo; los bloques incompletos
# controlan la variabilidad espacial a menor escala.
# ============================================================================

library(agricolae)
library(tidyverse)
library(lme4)        # Para modelos mixtos
library(lmerTest)    # Para p-valores en modelos mixtos
library(performance) # Diagnóstico avanzado de modelos mixtos
library(emmeans)     # Medias ajustadas

set.seed(123)

# 1. GENERACIÓN DEL DISEÑO
# Regla: design.alpha(trt, k, r, serie=2, seed=123)
# k=4 (tamaño de bloque), r=2 (repeticiones)
genotipos <- paste0("G_", sprintf("%02d", 1:16))
design_alpha <- design.alpha(trt = genotipos, k = 4, r = 2, serie = 2, seed = 123)
datos_alpha <- design_alpha$book

# 2. EJERCICIO AVANZADO: INTERACCIÓN GENOTIPO X AMBIENTE (GxE)
# Escenario: ¿Cómo rinde el trigo en un ambiente Óptimo vs Calor Extremo?
datos_met <- bind_rows(
  datos_alpha %>% mutate(env = "Optimo", ef_env = 0),
  datos_alpha %>% mutate(env = "Calor_Extremo", ef_env = -1500)
)

# Simulamos efectos GxE: Algunos genotipos son resilientes al calor
datos_met <- datos_met %>%
  mutate(
    # agricolae genera la columna 'genotipos'
    g_id = as.numeric(as.factor(genotipos)),
    yield_base = 5000 + ef_env + rnorm(16, 0, 400)[g_id],
    # Interacción: GxE dinámico
    gxe_ef = ifelse(env == "Calor_Extremo", rnorm(16, 0, 300)[g_id], 0),
    yield = yield_base + gxe_ef + rnorm(n(), 0, 150)
  )

cat("\n--- [1] MODELO MIXTO PARA ALPHA-LATTICE (GxE) ---\n")
# Usamos LMER porque los bloques incompletos son efectos aleatorios.
# Formula: yield ~ Ambiente + (1|Genotipo) + (1|Genotipo:Ambiente) + Errores jerárquicos
modelo_mixto <- lmer(yield ~ env + (1|genotipos) + (1|genotipos:env) + 
                       (1|env:replication/block), data = datos_met)

print(summary(modelo_mixto))

# EXPLICACIÓN: El modelo mixto estima la varianza de los genotipos (V_g)
# y la varianza de la interacción (V_gxe). Si V_gxe es alta, los genotipos
# cambian de ranking según el ambiente.

cat("\n--- [2] EVALUACIÓN DE SUPUESTOS CON PERFORMANCE ---\n")
# check_model nos muestra si el modelo mixto es válido.
# Un gráfico clave aquí es la normalidad de los efectos aleatorios (QQ-plot).
print(check_model(modelo_mixto))

cat("\n--- [3] GRÁFICO DE INTERACCIÓN GxE (REACTION NORMS) ---\n")
# ¿Para qué sirve? Visualiza la estabilidad. Líneas planas = estables.
# Líneas cruzadas = Interacción significativa (Cambio de ranking).

# Medias estimadas por combinación G x E
emm_gxe <- emmeans(modelo_mixto, ~ genotipos | env) %>% as.data.frame()

ggplot(emm_gxe, aes(x = env, y = emmean, group = genotipos, color = genotipos)) +
  # Líneas que conectan las medias estimadas de cada genotipo entre ambientes
  stat_summary(fun = mean, geom = "line", alpha = 0.6, linewidth = 1) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "Interacción Genotipo x Ambiente en Trigo",
       subtitle = "Respuesta diferencial a ambientes: Óptimo vs Calor Extremo",
       x = "Ambiente (Estrés Térmico)",
       y = "Rendimiento Estimado (kg/ha)") +
  theme_minimal() +
  theme(legend.position = "none", # Muchos genotipos para leyenda
        plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold"))

cat("\n--- [4] ¿POR QUÉ USAMOS ESTE GRÁFICO? ---\n")
cat("- Identifica genotipos con alta plasticidad fenotípica.\n")
cat("- Permite seleccionar variedades que mantienen el rendimiento bajo calor.\n")
cat("- El uso de 'stat_summary' garantiza que graficamos las tendencias reales.\n")
