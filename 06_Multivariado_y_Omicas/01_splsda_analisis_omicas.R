# ============================================================================
# 01_splsda_analisis_omicas.R
# Análisis sPLS-DA (Sparse Partial Least Squares Discriminant Analysis)
# Contexto: Transcriptómica de Trigo bajo Estrés Hídrico y Salino
# ============================================================================
# Este script utiliza 'mixOmics' para identificar biomarcadores (genes) 
# que discriminan entre grupos de estrés mediante reducción de dimensiones
# y selección de variables (L1 penalization).
# ============================================================================

library(mixOmics)
library(MASS)
library(ggplot2)
library(lattice)
library(tidyverse)

set.seed(123)

# 1. SIMULACIÓN DE DATOS (Omicas: 30 muestras x 500 genes)
# 3 Grupos: Control, Sequia, Salinidad (n=10 por grupo)
n_muestras <- 30
n_genes <- 500
grupos <- factor(rep(c("Control", "Sequia", "Salinidad"), each = 10))

# Simular expresión génica (log-counts)
X <- matrix(rnorm(n_muestras * n_genes, 5, 2), nrow = n_muestras, ncol = n_genes)
colnames(X) <- paste0("Gene_", 1:n_genes)
rownames(X) <- paste0("Sample_", 1:n_muestras)

# Inducir señal biológica en los primeros 50 genes para Sequia y Salinidad
X[11:20, 1:25] <- X[11:20, 1:25] + 4  # Genes sobreexpresados en Sequia
X[21:30, 26:50] <- X[21:30, 26:50] + 5 # Genes sobreexpresados en Salinidad

# 2. OPTIMIZACIÓN DEL MODELO (Tuning)
# tune.splsda evalúa el número óptimo de componentes y variables (keepX)
# usando validación cruzada (Mfold o LOO).

cat("\n--- Iniciando Tuneo del Modelo sPLS-DA (Cálculo de keepX óptimo) ---\n")
list_keepX <- c(5, 10, 20, 50) # Candidatos para selección de genes por componente

tune_splsda <- tune.splsda(
  X = X, 
  Y = grupos, 
  ncomp = 3, 
  test.keepX = list_keepX,
  validation = "Mfold", 
  folds = 5, 
  dist = "max.dist", 
  progressBar = FALSE,
  nrepeat = 10
)

# Extraer parámetros óptimos
ncomp_opt <- tune_splsda$choice.ncomp$ncomp
keepX_opt <- tune_splsda$choice.keepX[1:ncomp_opt]

cat("Número óptimo de componentes:", ncomp_opt, "\n")
cat("Variables óptimas por componente (keepX):", keepX_opt, "\n")

# 3. AJUSTE DEL MODELO FINAL
final_splsda <- splsda(X, grupos, ncomp = ncomp_opt, keepX = keepX_opt)

# 4. VISUALIZACIÓN DE RESULTADOS (Nivel Publicación)

# 4.1 Gráfico de Individuos (Indiv Plot) con Elipses de Confianza
cat("\n--- Generando Gráfico de Individuos ---\n")
plotIndiv(
  final_splsda, 
  comp = c(1, 2),
  group = grupos, 
  ind.names = FALSE,
  ellipse = TRUE, 
  legend = TRUE,
  title = "sPLS-DA: Discriminación de Estrés en Trigo",
  style = "ggplot2"
)

# 4.2 Gráfico de Loadings (Importancia de Biomarcadores)
# Identifica qué genes contribuyen más a la separación de los grupos.
cat("\n--- Generando Gráfico de Loadings (Top Genes comp 1) ---\n")
plotLoadings(
  final_splsda, 
  comp = 1, 
  method = 'median', 
  contrib = 'max', 
  size.name = 0.8, 
  legend = TRUE,
  title = "Biomarcadores Clave (Componente 1)"
)

# 5. VALIDACIÓN DEL DESEMPEÑO
perf_splsda <- perf(final_splsda, validation = "Mfold", folds = 5, nrepeat = 10)
cat("\n--- Error de Clasificación (Overall Error Rate) ---\n")
print(perf_splsda$error.rate)

cat("\n--- Script Finalizado (sPLS-DA Omics) ---\n")

