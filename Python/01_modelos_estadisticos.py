"""
Script: 01_modelos_estadisticos.py
Descripción: Plantilla con modelos estadísticos frecuentistas aplicados a
             datos experimentales biológicos / agronómicos.
             Modelos incluidos: ANOVA, LM, GLM, GAM, modelos mixtos (LMM/GLMM)
Autor: [nombre]
Fecha: [fecha]
"""

# ── 0. Importaciones ──────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Modelos estadísticos
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# GAM
from pygam import LinearGAM, PoissonGAM, s, f  # pip install pygam

# Diagnósticos
from scipy import stats

# Configuración general
np.random.seed(42)
pd.set_option("display.float_format", "{:.4f}".format)
sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)

# ── 1. Datos de ejemplo ───────────────────────────────────────────────────────
n_per_group = 20
tratamientos = np.repeat(["Control", "Dosis_baja", "Dosis_alta"], n_per_group)
bloques = np.tile(np.arange(1, 5), 15)

rendimiento = np.concatenate([
    np.random.normal(5.0, 0.8, n_per_group),
    np.random.normal(6.2, 0.9, n_per_group),
    np.random.normal(7.5, 1.0, n_per_group),
])
conteo = np.random.poisson(lam=np.repeat([3, 6, 10], n_per_group))
temperatura = np.random.uniform(15, 35, 60)

df = pd.DataFrame({
    "tratamiento": pd.Categorical(
        tratamientos,
        categories=["Control", "Dosis_baja", "Dosis_alta"],
        ordered=True
    ),
    "bloque": bloques.astype(str),
    "rendimiento": rendimiento,
    "conteo": conteo,
    "temperatura": temperatura,
})

print("=== Vista general de los datos ===")
print(df.describe())
print(df.dtypes)

# ── 2. Estadística descriptiva ────────────────────────────────────────────────
desc = (df.groupby("tratamiento", observed=True)["rendimiento"]
          .agg(n="count", media="mean", DE="std")
          .assign(EE=lambda x: x["DE"] / np.sqrt(x["n"]))
          .assign(IC_inf=lambda x: x["media"] - 1.96 * x["EE"],
                  IC_sup=lambda x: x["media"] + 1.96 * x["EE"]))
print("\n=== Estadística descriptiva ===")
print(desc)

# ── 3. ANOVA de una vía ───────────────────────────────────────────────────────
print("\n====== ANOVA de una vía ======")
modelo_anova = smf.ols("rendimiento ~ C(tratamiento)", data=df).fit()
tabla_anova = anova_lm(modelo_anova, typ=1)
print(tabla_anova)

# Comparaciones múltiples (Tukey HSD)
tukey = pairwise_tukeyhsd(df["rendimiento"], df["tratamiento"], alpha=0.05)
print("\nComparaciones Tukey HSD:")
print(tukey)

# Gráfico de diagnóstico de residuos
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
residuos = modelo_anova.resid
stats.probplot(residuos, dist="norm", plot=axes[0])
axes[0].set_title("Q-Q Plot de residuos")
axes[1].scatter(modelo_anova.fittedvalues, residuos, alpha=0.6)
axes[1].axhline(0, color="red", linestyle="--")
axes[1].set_xlabel("Valores ajustados")
axes[1].set_ylabel("Residuos")
axes[1].set_title("Residuos vs. Valores ajustados")
plt.tight_layout()
plt.savefig("resultados/diagnostico_anova.png", dpi=150)
plt.show()

# ── 4. Modelo lineal (LM) ─────────────────────────────────────────────────────
print("\n====== Modelo Lineal (LM) ======")
modelo_lm = smf.ols("rendimiento ~ C(tratamiento)", data=df).fit()
print(modelo_lm.summary())

# ── 5. Modelo lineal mixto (LMM) ──────────────────────────────────────────────
print("\n====== Modelo Lineal Mixto (LMM) ======")
modelo_lmm = smf.mixedlm(
    "rendimiento ~ C(tratamiento)",
    data=df,
    groups=df["bloque"]
).fit(reml=True)
print(modelo_lmm.summary())

# ── 6. GLM – Poisson (conteos) ────────────────────────────────────────────────
print("\n====== GLM - Poisson (conteos) ======")
modelo_glm_pois = smf.glm(
    "conteo ~ C(tratamiento)",
    data=df,
    family=sm.families.Poisson(link=sm.families.links.Log())
).fit()
print(modelo_glm_pois.summary())

# Verificar sobredispersión (razón deviance/gl)
ratio_disp = modelo_glm_pois.deviance / modelo_glm_pois.df_resid
print(f"\nRazón deviance/df (>1 indica sobredispersión): {ratio_disp:.3f}")

# ── 7. GLM – Gamma (variable continua positiva) ───────────────────────────────
print("\n====== GLM - Gamma ======")
# Simulamos biomasa continua y positiva
df["biomasa"] = np.random.gamma(shape=2, scale=1.0 / np.repeat([0.4, 0.3, 0.2],
                                                                 n_per_group))
modelo_glm_gamma = smf.glm(
    "biomasa ~ C(tratamiento)",
    data=df,
    family=sm.families.Gamma(link=sm.families.links.Log())
).fit()
print(modelo_glm_gamma.summary())

# ── 8. GAM – Modelo aditivo generalizado ─────────────────────────────────────
print("\n====== GAM (LinearGAM) ======")
X = np.column_stack([
    (df["tratamiento"] == "Dosis_baja").astype(int),
    (df["tratamiento"] == "Dosis_alta").astype(int),
    df["temperatura"].values
])
y = df["rendimiento"].values

gam = LinearGAM(f(0) + f(1) + s(2, n_splines=8)).fit(X, y)
print(gam.summary())

# Visualizar efecto de temperatura
fig, ax = plt.subplots(figsize=(7, 5))
XX = gam.generate_X_grid(term=2)
pdep, confi = gam.partial_dependence(term=2, X=XX, width=0.95)
ax.plot(XX[:, 2], pdep, color="steelblue", linewidth=2)
ax.fill_between(XX[:, 2], confi[:, 0], confi[:, 1], alpha=0.3)
ax.set_xlabel("Temperatura (°C)")
ax.set_ylabel("Efecto parcial sobre rendimiento")
ax.set_title("Efecto suavizado de temperatura (GAM)")
plt.tight_layout()
plt.savefig("resultados/gam_temperatura.png", dpi=150)
plt.show()

# ── 9. GLMM – Poisson con efecto de bloque ───────────────────────────────────
print("\n====== GLMM - Poisson con efecto aleatorio de bloque ======")
modelo_glmm = smf.glm(
    "conteo ~ C(tratamiento)",
    data=df,
    family=sm.families.Poisson()
).fit()
# Nota: statsmodels no soporta efectos aleatorios en GLM directamente.
# Para GLMM completo se recomienda usar R (lme4::glmer) o pymer4 / glmmTMB.
print(modelo_glmm.summary())
print("\nNota: Para GLMM con efectos aleatorios en Python se recomienda pymer4 o R+lme4.")

# ── 10. Visualización comparativa ────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data=df, x="tratamiento", y="rendimiento", palette="Set2", ax=ax)
sns.stripplot(data=df, x="tratamiento", y="rendimiento",
              color="black", alpha=0.4, size=4, ax=ax)
ax.set_title("Rendimiento por tratamiento")
ax.set_xlabel("Tratamiento")
ax.set_ylabel("Rendimiento (unidades)")
plt.tight_layout()
plt.savefig("resultados/boxplot_tratamientos.png", dpi=150)
plt.show()
