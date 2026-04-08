# -------------------------------------------------------------------------
# 02_Diseno_Experimental
# Bloques Completos al Azar (BCA) y Alpha-Lattice
# Contexto: ensayos simulados de sequía y salinidad en genotipos vegetales
# -------------------------------------------------------------------------

set.seed(321)

# =========================
# A) BCA: ensayo de sequía
# =========================
tratamientos <- factor(c("Riego_optimo", "Deficit_moderado", "Deficit_severo"))
bloques <- factor(paste0("Bloque_", 1:5))

bca <- expand.grid(bloque = bloques, tratamiento = tratamientos)

ef_trat <- c(Riego_optimo = 0, Deficit_moderado = -7, Deficit_severo = -13)
ef_bloque <- rnorm(length(bloques), 0, 2)
names(ef_bloque) <- levels(bloques)

bca$rendimiento <- 65 +
  ef_trat[as.character(bca$tratamiento)] +
  ef_bloque[as.character(bca$bloque)] +
  rnorm(nrow(bca), 0, 2.8)

cat("\n--- ANOVA para BCA (sequía) ---\n")
mod_bca <- aov(rendimiento ~ tratamiento + bloque, data = bca)
print(summary(mod_bca))

# =========================================
# B) Alpha-lattice: ensayo de salinidad
# =========================================
# Simulamos 24 genotipos, 3 repeticiones, bloques incompletos de tamaño 4
set.seed(322)
genotipos <- factor(paste0("G", sprintf("%02d", 1:24)))
reps <- factor(paste0("Rep_", 1:3))

# Estructura lattice simulada manualmente
alpha <- do.call(rbind, lapply(reps, function(r) {
  data.frame(
    rep = r,
    bloque_incompleto = factor(rep(paste0(r, "_B", 1:6), each = 4)),
    genotipo = sample(genotipos, 24, replace = FALSE)
  )
}))

# Efectos simulados
ef_gen <- rnorm(length(genotipos), 0, 5); names(ef_gen) <- levels(genotipos)
ef_rep <- rnorm(length(reps), 0, 1.5); names(ef_rep) <- levels(reps)
ef_bi <- rnorm(length(unique(alpha$bloque_incompleto)), 0, 2.0)
names(ef_bi) <- levels(alpha$bloque_incompleto)

alpha$indice_tolerancia_salina <- 50 +
  ef_gen[as.character(alpha$genotipo)] +
  ef_rep[as.character(alpha$rep)] +
  ef_bi[as.character(alpha$bloque_incompleto)] +
  rnorm(nrow(alpha), 0, 2.5)

# Ajuste como modelo mixto (equivalente práctico para alpha-lattice)
if (requireNamespace("lme4", quietly = TRUE)) {
  cat("\n--- Modelo mixto aproximando Alpha-lattice ---\n")
  mod_alpha <- lme4::lmer(
    indice_tolerancia_salina ~ rep + (1 | genotipo) + (1 | rep:bloque_incompleto),
    data = alpha
  )
  print(summary(mod_alpha))
} else {
  message("Paquete 'lme4' no disponible. Instálalo con: install.packages('lme4')")
}

# Nota: para diseño alpha-lattice formal se puede complementar con agricolae o sommer.
