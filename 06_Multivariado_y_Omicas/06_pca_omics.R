# ---------------------------------------------------------------------
# 06_Multivariado_y_Omicas
# PCA para datos transcriptómicos/metabolómicos de gran dimensión
# ---------------------------------------------------------------------

set.seed(1001)

# Matriz grande simulada: 150 muestras x 2000 variables
n_muestras <- 150
n_vars <- 2000

grupo <- factor(rep(c("Control", "Estres_hidrico", "Estres_termico"), each = 50))
X <- matrix(rnorm(n_muestras * n_vars, mean = 0, sd = 1), nrow = n_muestras, ncol = n_vars)

# Insertamos señal en un subconjunto de variables para emular biomarcadores
X[grupo == "Estres_hidrico", 1:80] <- X[grupo == "Estres_hidrico", 1:80] + 1.2
X[grupo == "Estres_termico", 81:160] <- X[grupo == "Estres_termico", 81:160] + 1.0

# Escalado eficiente y PCA por SVD truncada (prcomp)
pca <- prcomp(X, center = TRUE, scale. = TRUE)

var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
cat("\n--- PCA para matrices ómicas grandes ---\n")
cat("Varianza explicada PC1-PC5:\n")
print(round(var_exp[1:5], 4))

# Visualización base (primeras dos componentes)
plot(pca$x[, 1], pca$x[, 2],
     col = as.integer(grupo), pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "PCA - perfiles ómicos simulados")
legend("topright", legend = levels(grupo), col = 1:3, pch = 19, bty = "n")
