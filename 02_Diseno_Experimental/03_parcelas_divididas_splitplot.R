# ============================================================================
# 03_parcelas_divididas_splitplot.R (SPLIT-SPLIT-PLOT)
# Diseño de Parcelas Sub-subdivididas (Riego x HMA x NPK)
# ============================================================================
# Este diseño es el más complejo: Los factores se aplican a diferentes escalas.
# Parcela Principal: Riego (Sequía). Subparcela: HMA (Inóculo).
# Sub-subparcela: NPK (Fertilizante). Cada nivel tiene su propio error.
# ============================================================================

library(agricolae)
library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(emmeans)

set.seed(123)

# 1. GENERACIÓN DEL DISEÑO
# Regla: design.split(trt1, trt2, r, design="rcbd", serie=2, seed=123)
riegos <- c("Control", "Drought")
micorr <- c("No_HMA", "HMA")
r <- 3 # Bloques
design_sp <- design.split(trt1 = riegos, trt2 = micorr, r = r, 
                          design = "rcbd", serie = 2, seed = 123)
datos_sp <- design_sp$book

# Escalamiento manual a Sub-subparcelas (NPK)
npk_levels <- c("0%", "50%", "100%")
datos_ssp <- datos_sp %>%
  group_by(plots) %>%
  slice(rep(1, each = length(npk_levels))) %>%
  mutate(npk = rep(npk_levels, length.out = n())) %>%
  ungroup()

# 2. SIMULACIÓN DE BIOMASA (TRIPLE INTERACCIÓN)
datos_ssp <- datos_ssp %>%
  mutate(
    yield_base = case_when(
      riegos == "Control" ~ 52 + as.numeric(gsub("%", "", npk))*0.25,
      riegos == "Drought" & micorr == "HMA" ~ 48 + as.numeric(gsub("%", "", npk))*0.2,
      TRUE ~ 32 + as.numeric(gsub("%", "", npk))*0.15
    ),
    # Errores jerárquicos por nivel de parcela
    err_main = rnorm(n(), 0, 4), 
    err_sub = rnorm(n(), 0, 2),
    biomasa = yield_base + err_main + err_sub + rnorm(n(), 0, 1)
  )

cat("\n--- [1] ANÁLISIS DE VARIANZA JERÁRQUICO (Split-Split-Plot) ---\n")
# Estructura de error compleja: (1|block) + (1|block:riegos) + (1|block:riegos:micorr)
modelo_ssp <- lmer(biomasa ~ riegos * micorr * npk + 
                     (1|block) + (1|block:riegos) + (1|block:riegos:micorr), 
                   data = datos_ssp)

print(anova(modelo_ssp))

# EXPLICACIÓN: El factor Riego (Parcela Principal) es el más difícil de probar
# (menos grados de libertad). El NPK (Sub-subparcela) es el más preciso.

cat("\n--- [2] EVALUACIÓN INTEGRAL DE SUPUESTOS ---\n")
# Usar performance::check_model() para validar la estructura jerárquica.
print(check_model(modelo_ssp))

cat("\n--- [3] GRÁFICOS DE INTERACCIÓN TRIPLES ---\n")
# Medias estimadas por emmeans
emm_ssp <- emmeans(modelo_ssp, ~ npk | micorr * riegos) %>% as.data.frame()

# Gráfico visualmente potente
ggplot(emm_ssp, aes(x = npk, y = emmean, color = micorr, group = micorr)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  # El facet_wrap permite ver la interacción según Riego
  facet_wrap(~ riegos, labeller = label_both) +
  labs(title = "Sinergia Microorganismos x Fertilizante bajo Sequía",
       subtitle = "Diseño de Parcelas Sub-subdivididas (Split-Split-Plot)",
       x = "Dosis de NPK (% de la recomendación)",
       y = "Biomasa Estimada (g/planta)") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold")) +
  scale_color_manual(values = c("HMA" = "#228B22", "No_HMA" = "#B22222"))

cat("\n--- [4] ¿POR QUÉ USAMOS ESTAS PRUEBAS? ---\n")
cat("- Modelos Mixtos: Los únicos que modelan la 'dependencia' de subparcelas.\n")
cat("- Emmeans: Permiten comparar medias en diseños desbalanceados o complejos.\n")
cat("- Faceting: Permite ver si el inóculo (HMA) compensa la falta de fertilizante NPK.\n")
