# 02_Diseno_Experimental/02_disenos_bloques.R
# RCBD y ejemplo de bloques incompletos (alpha-lattice aproximado con lme4)

suppressPackageStartupMessages({
  library(agricolae)
  library(lme4)
  library(lmerTest)
})

set.seed(42)

# -------------------------------
# 1) RCBD (Bloques Completos al Azar)
# -------------------------------
rcbd <- data.frame(
  bloque = factor(rep(1:4, each = 3)),
  tratamiento = factor(rep(c("Control", "Moderado", "Severo"), times = 4))
)

efecto_trat <- c(Control = 6.4, Moderado = 5.8, Severo = 4.9)
efecto_bloq <- c("1" = 0.20, "2" = -0.10, "3" = 0.15, "4" = -0.25)

rcbd$rendimiento <- with(rcbd,
  efecto_trat[tratamiento] + efecto_bloq[bloque] + rnorm(nrow(rcbd), 0, 0.25)
)

modelo_rcbd <- aov(rendimiento ~ bloque + tratamiento, data = rcbd)
cat("=== ANOVA RCBD ===\n")
print(summary(modelo_rcbd))

cat("\n=== Comparaciones de medias (LSD, agricolae) ===\n")
print(agricolae::LSD.test(modelo_rcbd, "tratamiento", p.adj = "none"))

# -------------------------------
# 2) Bloques incompletos (estructura tipo alpha-lattice)
# -------------------------------
incomp <- expand.grid(
  rep = factor(1:3),
  bloque_incomp = factor(1:4),
  tratamiento = factor(paste0("T", 1:4))
)

# Se simula una estructura con bloque incompleto anidado en repeticion
treatment_slope <- -1.2
incomp$y <- 50 +
  as.numeric(incomp$tratamiento) * treatment_slope +
  rnorm(nlevels(incomp$rep), 0, 1.5)[incomp$rep] +
  rnorm(nlevels(incomp$bloque_incomp), 0, 1.0)[incomp$bloque_incomp] +
  rnorm(nrow(incomp), 0, 1.8)

modelo_incomp <- lmer(y ~ tratamiento + (1 | rep/bloque_incomp), data = incomp)
cat("\n=== Modelo mixto para bloques incompletos ===\n")
print(summary(modelo_incomp))
print(anova(modelo_incomp))
