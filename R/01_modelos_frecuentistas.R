# =============================================================================
# Script: 01_modelos_frecuentistas.R
# Descripción: Plantilla con modelos estadísticos frecuentistas aplicados a
#              datos experimentales biológicos / agronómicos.
#              Modelos incluidos: ANOVA, LM, LMM, GLM, GAM, GLMM
# Autor: [nombre]
# Fecha: [fecha]
# =============================================================================

# ── 0. Paquetes ──────────────────────────────────────────────────────────────
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,   # manipulación y visualización de datos
  lme4,        # modelos lineales mixtos (LMM) y generalizados mixtos (GLMM)
  lmerTest,    # p-valores para LMM (Satterthwaite)
  emmeans,     # medias marginales estimadas y comparaciones post-hoc
  mgcv,        # modelos aditivos generalizados (GAM)
  car,         # pruebas de Levene, tipo III SS
  performance, # verificación de supuestos
  see,         # gráficos diagnósticos
  ggplot2      # visualización
)

# ── 1. Datos de ejemplo ──────────────────────────────────────────────────────
# Sustituir con: datos <- read.csv("ruta/a/datos.csv")
set.seed(42)
datos <- data.frame(
  tratamiento  = rep(c("Control", "Dosis_baja", "Dosis_alta"), each = 20),
  bloque       = factor(rep(1:4, times = 15)),
  rendimiento  = c(
    rnorm(20, mean = 5.0, sd = 0.8),
    rnorm(20, mean = 6.2, sd = 0.9),
    rnorm(20, mean = 7.5, sd = 1.0)
  ),
  conteo       = rpois(60, lambda = c(3, 6, 10)),
  biomasa      = rgamma(60, shape = 2, rate = c(0.4, 0.3, 0.2))
)
datos$tratamiento <- factor(datos$tratamiento,
                            levels = c("Control", "Dosis_baja", "Dosis_alta"))

# Vista rápida
glimpse(datos)
summary(datos)

# ── 2. Estadística descriptiva ───────────────────────────────────────────────
datos %>%
  group_by(tratamiento) %>%
  summarise(
    n      = n(),
    media  = mean(rendimiento),
    DE     = sd(rendimiento),
    EE     = DE / sqrt(n),
    IC_inf = media - 1.96 * EE,
    IC_sup = media + 1.96 * EE
  )

# ── 3. ANOVA de una vía ───────────────────────────────────────────────────────
cat("\n====== ANOVA de una vía ======\n")
modelo_anova <- aov(rendimiento ~ tratamiento, data = datos)
summary(modelo_anova)

# Verificar supuestos
par(mfrow = c(2, 2))
plot(modelo_anova)
par(mfrow = c(1, 1))

# Homogeneidad de varianzas (Levene)
car::leveneTest(rendimiento ~ tratamiento, data = datos)

# Comparaciones múltiples (Tukey)
emmeans(modelo_anova, pairwise ~ tratamiento, adjust = "tukey")

# ── 4. Modelo lineal (LM) ────────────────────────────────────────────────────
cat("\n====== Modelo Lineal (LM) ======\n")
modelo_lm <- lm(rendimiento ~ tratamiento, data = datos)
summary(modelo_lm)
car::Anova(modelo_lm, type = "III")

# Diagnósticos con {performance}
performance::check_model(modelo_lm)

# ── 5. Modelo lineal mixto (LMM) ─────────────────────────────────────────────
cat("\n====== Modelo Lineal Mixto (LMM) ======\n")
# bloque como efecto aleatorio (diseño de bloques completos aleatorizados)
modelo_lmm <- lmer(rendimiento ~ tratamiento + (1 | bloque), data = datos)
summary(modelo_lmm)
anova(modelo_lmm)                        # F-test con df Satterthwaite
emmeans(modelo_lmm, pairwise ~ tratamiento, adjust = "tukey")
performance::check_model(modelo_lmm)

# ── 6. Modelo lineal generalizado (GLM) ──────────────────────────────────────
cat("\n====== GLM - Poisson (conteos) ======\n")
modelo_glm_pois <- glm(conteo ~ tratamiento, data = datos, family = poisson(link = "log"))
summary(modelo_glm_pois)
car::Anova(modelo_glm_pois, type = "II")

# Verificar sobredispersión
performance::check_overdispersion(modelo_glm_pois)

# Si hay sobredispersión, usar quasi-Poisson o binomial negativa
cat("\n====== GLM - Gamma (biomasa continua positiva) ======\n")
modelo_glm_gamma <- glm(biomasa ~ tratamiento, data = datos,
                        family = Gamma(link = "log"))
summary(modelo_glm_gamma)

# ── 7. Modelo aditivo generalizado (GAM) ─────────────────────────────────────
cat("\n====== GAM ======\n")
# Simulamos una variable continua (e.g., temperatura)
datos$temperatura <- runif(60, 15, 35)

modelo_gam <- mgcv::gam(
  rendimiento ~ tratamiento + s(temperatura, k = 5),
  data = datos,
  method = "REML"
)
summary(modelo_gam)
mgcv::gam.check(modelo_gam)

# Visualizar splines
plot(modelo_gam, pages = 1, residuals = TRUE, shade = TRUE)

# ── 8. Modelo lineal generalizado mixto (GLMM) ───────────────────────────────
cat("\n====== GLMM - Poisson con efecto de bloque ======\n")
modelo_glmm <- lme4::glmer(
  conteo ~ tratamiento + (1 | bloque),
  data   = datos,
  family = poisson(link = "log")
)
summary(modelo_glmm)
emmeans(modelo_glmm, pairwise ~ tratamiento, adjust = "tukey", type = "response")

# ── 9. Visualización de resultados ───────────────────────────────────────────
medias <- emmeans(modelo_lmm, ~ tratamiento) %>% as.data.frame()

ggplot(medias, aes(x = tratamiento, y = emmean, ymin = lower.CL, ymax = upper.CL,
                   color = tratamiento)) +
  geom_point(size = 3) +
  geom_errorbar(width = 0.2) +
  labs(
    title    = "Medias marginales estimadas ± IC 95%",
    subtitle = "Modelo lineal mixto (LMM)",
    x        = "Tratamiento",
    y        = "Rendimiento (unidades)",
    caption  = "Ajuste: Satterthwaite; comparaciones Tukey"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

# ── 10. Exportar tablas de resultados ────────────────────────────────────────
# Descomentar para guardar resultados
# write.csv(as.data.frame(summary(modelo_lmm)$coefficients),
#           "resultados/tabla_coeficientes_lmm.csv", row.names = TRUE)
