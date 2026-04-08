# ---------------------------------------------------------------------
# 06_Multivariado_y_Omicas
# PERMANOVA para comunidades microbianas/edafofauna
# ---------------------------------------------------------------------

set.seed(1003)

n_muestras <- 120
n_taxa <- 300

meta <- data.frame(
  tratamiento = factor(rep(c("Control", "Bioestimulante", "Deficit_hidrico"), each = 40)),
  sitio = factor(rep(paste0("Sitio", 1:4), length.out = n_muestras))
)

# Matriz de abundancia simulada (conteos)
abund <- matrix(rpois(n_muestras * n_taxa, lambda = 8), nrow = n_muestras, ncol = n_taxa)

# Señales de comunidad por tratamiento
abund[meta$tratamiento == "Bioestimulante", 1:35] <- abund[meta$tratamiento == "Bioestimulante", 1:35] + 5
abund[meta$tratamiento == "Deficit_hidrico", 36:70] <- abund[meta$tratamiento == "Deficit_hidrico", 36:70] + 6

if (requireNamespace("vegan", quietly = TRUE)) {
  dist_bc <- vegan::vegdist(abund, method = "bray")

  cat("\n--- PERMANOVA (adonis2) ---\n")
  perma <- vegan::adonis2(dist_bc ~ tratamiento + sitio, data = meta, permutations = 999)
  print(perma)

  # Ordenación auxiliar
  ord <- vegan::metaMDS(abund, distance = "bray", k = 2, trymax = 20, trace = 0)
  plot(ord$points, col = as.integer(meta$tratamiento), pch = 19,
       xlab = "MDS1", ylab = "MDS2", main = "Ordenación de comunidades simuladas")
  legend("topright", legend = levels(meta$tratamiento), col = 1:3, pch = 19, bty = "n")
} else {
  message("Paquete 'vegan' no disponible. Instálalo con: install.packages('vegan')")
}
