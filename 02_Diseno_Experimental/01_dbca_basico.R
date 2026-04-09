# ============================================================================
# 01_dbca_basico.R (GUÍA MAESTRA: PGPR Y ESTRÉS SALINO)
# Diseño de Bloques Completos al Azar (DBCA)
# ============================================================================
# Este script evalúa el uso de bioestimulantes (PGPR) para mitigar el estrés
# salino en plantas. El DBCA se usa cuando el área experimental no es homogénea
# (ej. gradiente de luz o humedad), agrupando parcelas similares en 'bloques'.
# ============================================================================

# Cargar librerías críticas
library(agricolae)     # Para el diseño experimental
library(tidyverse)     # Para manipulación de datos y visualización
library(performance)   # Para evaluación integral de supuestos (Top Tool)
library(car)           # Para pruebas de varianza (Levene)
library(emmeans)       # Para medias estimadas
library(multcomp)      # Para obtener letras de significancia
library(multcompView)  # Para convertir comparaciones en letras
library(nlme)          # Para modelos avanzados (GLS)

set.seed(123)

# 1. GENERACIÓN DEL DISEÑO EXPERIMENTAL
# Regla: design.rcbd(trt, r, serie=2, seed=123)
# ¿Para qué sirve? Crea una distribución aleatoria de tratamientos dentro de bloques
# para evitar sesgos por la ubicación de las plantas.
tratamientos <- c("Control", "Rhizobium", "Bacillus", "Pseudomonas", "Azospirillum")
r <- 4
design_dbca <- design.rcbd(trt = tratamientos, r = r, serie = 2, seed = 123)
datos <- design_dbca$book

# 2. SIMULACIÓN DE DATOS: BIOMASA BAJO ESTRÉS SALINO
# Escenario: El estrés salino aumenta la varianza (Heterocedasticidad).
datos <- datos %>%
  mutate(
    # agricolae usa 'tratamientos' como nombre de columna por defecto
    ef_trt = case_when(
      tratamientos == "Control" ~ 12,
      tratamientos == "Rhizobium" ~ 15,
      tratamientos == "Bacillus" ~ 21,
      tratamientos == "Pseudomonas" ~ 18,
      tratamientos == "Azospirillum" ~ 14
    ),
    # Simulamos varianza desigual: los microorganismos eficientes estabilizan la planta
    sd_error = ifelse(tratamientos == "Bacillus", 0.8, 3.0),
    biomasa = ef_trt + as.numeric(block)*0.4 + rnorm(n(), 0, sd_error)
  )

cat("\n--- [1] ANÁLISIS DE VARIANZA (ANOVA) ---\n")
# El ANOVA nos dice si AL MENOS UN tratamiento es diferente a los demás.
modelo_aov <- aov(biomasa ~ block + tratamientos, data = datos)
print(summary(modelo_aov))

# EXPLICACIÓN: Si el p-valor de 'tratamientos' es < 0.05, rechazamos la hipótesis
# nula (H0) de que todas las medias son iguales. Hay efecto de los bioestimulantes.

cat("\n--- [2] EVALUACIÓN DE SUPUESTOS CON 'PERFORMANCE' ---\n")
# ¿Por qué lo usamos? performance::check_model verifica:
# 1. Normalidad: ¿Los residuos siguen una campana de Gauss?
# 2. Homocedasticidad: ¿La varianza es igual en todos los grupos?
# 3. Outliers: ¿Hay datos extraños que arruinan el modelo?
print(check_model(modelo_aov))

cat("\n--- [3] PRUEBAS DE SIGNIFICANCIA POST-HOC (TUKEY) ---\n")
# Cuando el ANOVA es significativo, necesitamos saber qué cepa es mejor.
emm <- emmeans(modelo_aov, ~ tratamientos)
# Obtenemos las letras de significancia (CLD)
letras <- cld(emm, alpha = 0.05, Letters = letters, adjust = "tukey")
print(letras)

# 4. VISUALIZACIÓN ELEGANTE (BOXPLOT + LETRAS)
# Un gráfico dice más que mil tablas. Aquí combinamos biomasa y significancia.
ggplot(datos, aes(x = tratamientos, y = biomasa, fill = tratamientos)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  # Añadir letras de significancia sobre las cajas
  geom_text(data = letras, aes(y = emmean + (3*SE), label = .group), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Efecto de PGPR en Biomasa de Tomate bajo Salinidad",
       subtitle = "Letras diferentes indican diferencias significativas (Tukey, p < 0.05)",
       x = "Tratamiento Bioestimulante",
       y = "Biomasa de Raíz (g)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

cat("\n--- [4] ¿POR QUÉ USAMOS ESTAS PRUEBAS? ---\n")
cat("- Shapiro-Wilk: Valida que podamos usar estadística paramétrica.\n")
cat("- Levene: Clave en estrés abiótico para ver si el estrés descontrola la varianza.\n")
cat("- Tukey: Es el estándar de oro para comparar todos contra todos sin inflar el error.\n")
