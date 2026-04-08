# ---------------------------------------------------------------------
# 06_Multivariado_y_Omicas
# PLS-DA para clasificación de grupos en datos ómicos
# ---------------------------------------------------------------------

set.seed(1002)

n_muestras <- 180
n_vars <- 1200
clase <- factor(rep(c("Control", "Sequía", "Salinidad"), each = 60))

X <- matrix(rnorm(n_muestras * n_vars), nrow = n_muestras, ncol = n_vars)
X[clase == "Sequía", 1:60] <- X[clase == "Sequía", 1:60] + 1.5
X[clase == "Salinidad", 61:120] <- X[clase == "Salinidad", 61:120] + 1.3

if (requireNamespace("mixOmics", quietly = TRUE)) {
  plsda_fit <- mixOmics::plsda(X, clase, ncomp = 3)

  cat("\n--- PLS-DA (mixOmics) ---\n")
  print(plsda_fit)

  # Performance simple (clasificación interna)
  perf <- mixOmics::perf(plsda_fit, validation = "Mfold", folds = 5, nrepeat = 2, progressBar = FALSE)
  print(perf$error.rate)

  plotIndiv(plsda_fit, comp = c(1, 2), group = clase, legend = TRUE,
            title = "PLS-DA - Datos ómicos simulados")
} else {
  message("Paquete 'mixOmics' no disponible. Instálalo con: install.packages('mixOmics')")
}
