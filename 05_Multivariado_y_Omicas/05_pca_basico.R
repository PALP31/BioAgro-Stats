# 05_Multivariado_y_Omicas/05_pca_basico.R
# PCA básico con variables fisiológicas simuladas

suppressPackageStartupMessages({
  library(ggplot2)
})

set.seed(777)

n <- 80
trat <- factor(rep(c("Control", "Moderado", "Severo", "Recuperación"), each = n / 4))

datos <- data.frame(
  tratamiento = trat,
  clorofila = rnorm(n, mean = ifelse(trat == "Severo", 28, 38), sd = 3),
  conductancia = rnorm(n, mean = ifelse(trat == "Severo", 0.18, 0.35), sd = 0.05),
  biomasa = rnorm(n, mean = ifelse(trat == "Severo", 12, 18), sd = 2),
  prolina = rnorm(n, mean = ifelse(trat == "Severo", 7.5, 3.2), sd = 1.1),
  ndvi = rnorm(n, mean = ifelse(trat == "Severo", 0.54, 0.72), sd = 0.06)
)

X <- scale(datos[, -1])
pca <- prcomp(X, center = TRUE, scale. = TRUE)

cat("=== Resumen PCA ===\n")
print(summary(pca))

scores <- as.data.frame(pca$x[, 1:2])
scores$tratamiento <- datos$tratamiento

p <- ggplot(scores, aes(PC1, PC2, color = tratamiento)) +
  geom_point(size = 2, alpha = 0.85) +
  theme_minimal() +
  labs(title = "PCA: variables fisiológicas de estrés", x = "PC1", y = "PC2")
print(p)
