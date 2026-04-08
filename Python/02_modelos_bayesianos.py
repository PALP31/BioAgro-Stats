"""
Script: 02_modelos_bayesianos.py
Descripción: Plantilla para análisis bayesiano de datos experimentales
             biológicos / agronómicos con PyMC.
             Modelos: regresión lineal bayesiana, ANOVA bayesiano, GLMM bayesiano.
Autor: [nombre]
Fecha: [fecha]
"""

# ── 0. Importaciones ──────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import arviz as az          # pip install arviz
import pymc as pm           # pip install pymc

# Configuración
np.random.seed(42)
az.style.use("arviz-darkgrid")

# ── 1. Datos de ejemplo ───────────────────────────────────────────────────────
n_per_group = 20
n_bloques = 4
tratamientos_str = np.repeat(["Control", "Dosis_baja", "Dosis_alta"], n_per_group)
bloques_num = np.tile(np.arange(n_bloques), 15)

rendimiento = np.concatenate([
    np.random.normal(5.0, 0.8, n_per_group),
    np.random.normal(6.2, 0.9, n_per_group),
    np.random.normal(7.5, 1.0, n_per_group),
])
conteo = np.random.poisson(lam=np.repeat([3, 6, 10], n_per_group))

# Codificación numérica de tratamiento
trat_map = {"Control": 0, "Dosis_baja": 1, "Dosis_alta": 2}
tratamiento_idx = np.array([trat_map[t] for t in tratamientos_str])
n_tratamientos = len(trat_map)

# ── 2. Modelo lineal bayesiano (ANOVA bayesiano) ─────────────────────────────
print("\n====== Modelo Lineal Bayesiano (ANOVA) ======")

with pm.Model(coords={"tratamiento": list(trat_map.keys())}) as modelo_bayes_lm:

    # Priors
    mu_tratamiento = pm.Normal("mu_tratamiento", mu=5.0, sigma=3.0,
                               dims="tratamiento")
    sigma = pm.HalfNormal("sigma", sigma=2.0)

    # Esperanza
    mu_obs = mu_tratamiento[tratamiento_idx]

    # Verosimilitud
    y_obs = pm.Normal("y_obs", mu=mu_obs, sigma=sigma, observed=rendimiento)

    # Muestreo
    traza_lm = pm.sample(
        draws=2000,
        tune=1000,
        chains=4,
        cores=4,
        random_seed=42,
        target_accept=0.9,
        return_inferencedata=True
    )

# ── 3. Diagnósticos MCMC ─────────────────────────────────────────────────────
print("\nDiagnósticos MCMC:")
print(az.summary(traza_lm, var_names=["mu_tratamiento", "sigma"]))

# Trazas
az.plot_trace(traza_lm, var_names=["mu_tratamiento", "sigma"])
plt.tight_layout()
plt.savefig("resultados/trazas_modelo_lm_bayes.png", dpi=150)
plt.show()

# Distribuciones posteriores
az.plot_posterior(traza_lm, var_names=["mu_tratamiento"])
plt.savefig("resultados/posterior_mu_tratamientos.png", dpi=150)
plt.show()

# Posterior predictive check
with modelo_bayes_lm:
    ppc = pm.sample_posterior_predictive(traza_lm, random_seed=42)

az.plot_ppc(az.from_pymc(
    posterior_predictive=ppc,
    model=modelo_bayes_lm
), num_pp_samples=100)
plt.savefig("resultados/ppc_modelo_lm_bayes.png", dpi=150)
plt.show()

# ── 4. Contrastes entre tratamientos ─────────────────────────────────────────
post_mu = traza_lm.posterior["mu_tratamiento"]

diferencia_alta_control = (
    post_mu.sel(tratamiento="Dosis_alta") -
    post_mu.sel(tratamiento="Control")
)
prob_pos = float((diferencia_alta_control > 0).mean())
print(f"\nP(Dosis_alta > Control) = {prob_pos:.3f}")

# HDI al 94%
hdi_diff = az.hdi(diferencia_alta_control.values, hdi_prob=0.94)
print(f"HDI 94% de la diferencia: [{hdi_diff[0]:.3f}, {hdi_diff[1]:.3f}]")

# ── 5. Modelo Lineal Mixto Bayesiano (LMM) ───────────────────────────────────
print("\n====== LMM Bayesiano (efecto aleatorio de bloque) ======")

with pm.Model(coords={
    "tratamiento": list(trat_map.keys()),
    "bloque": [str(i) for i in range(n_bloques)]
}) as modelo_lmm_bayes:

    # Priors para efectos fijos
    alpha = pm.Normal("alpha", mu=5.0, sigma=3.0)
    beta  = pm.Normal("beta", mu=0.0, sigma=2.0, dims="tratamiento")

    # Prior para efectos aleatorios de bloque
    sigma_bloque = pm.HalfNormal("sigma_bloque", sigma=1.0)
    u_bloque = pm.Normal("u_bloque", mu=0.0, sigma=sigma_bloque, dims="bloque")

    sigma_resid = pm.HalfNormal("sigma_resid", sigma=2.0)

    mu_obs = alpha + beta[tratamiento_idx] + u_bloque[bloques_num]

    y_obs = pm.Normal("y_obs", mu=mu_obs, sigma=sigma_resid, observed=rendimiento)

    traza_lmm = pm.sample(
        draws=2000,
        tune=1000,
        chains=4,
        cores=4,
        random_seed=42,
        target_accept=0.9,
        return_inferencedata=True
    )

print(az.summary(traza_lmm, var_names=["alpha", "beta", "sigma_bloque", "sigma_resid"]))

az.plot_forest(traza_lmm, var_names=["beta"], combined=True,
               hdi_prob=0.95, r_hat=True)
plt.title("Efectos de tratamiento – LMM Bayesiano")
plt.savefig("resultados/forest_beta_lmm_bayes.png", dpi=150)
plt.show()

# ── 6. GLMM Bayesiano – Poisson ───────────────────────────────────────────────
print("\n====== GLMM Bayesiano - Poisson ======")

with pm.Model(coords={
    "tratamiento": list(trat_map.keys()),
    "bloque": [str(i) for i in range(n_bloques)]
}) as modelo_glmm_pois:

    alpha_p = pm.Normal("alpha_p", mu=1.0, sigma=1.5)
    beta_p  = pm.Normal("beta_p",  mu=0.0, sigma=1.0, dims="tratamiento")

    sigma_b = pm.HalfNormal("sigma_b", sigma=0.5)
    u_b = pm.Normal("u_b", mu=0.0, sigma=sigma_b, dims="bloque")

    log_mu = alpha_p + beta_p[tratamiento_idx] + u_b[bloques_num]
    mu_p = pm.math.exp(log_mu)

    y_cnt = pm.Poisson("y_cnt", mu=mu_p, observed=conteo)

    traza_glmm_pois = pm.sample(
        draws=2000,
        tune=1000,
        chains=4,
        cores=4,
        random_seed=42,
        target_accept=0.9,
        return_inferencedata=True
    )

print(az.summary(traza_glmm_pois, var_names=["alpha_p", "beta_p", "sigma_b"]))

# ── 7. Comparación de modelos (LOO-CV) ────────────────────────────────────────
print("\n====== Comparación de modelos (LOO) ======")
loo_lm  = az.loo(traza_lm,  pointwise=True)
loo_lmm = az.loo(traza_lmm, pointwise=True)
comp = az.compare({"LM_bayes": traza_lm, "LMM_bayes": traza_lmm})
print(comp)

az.plot_compare(comp)
plt.savefig("resultados/comparacion_modelos_loo.png", dpi=150)
plt.show()

# ── 8. Exportar resultados ────────────────────────────────────────────────────
# Guardar traza como NetCDF (formato estándar ArviZ)
# traza_lmm.to_netcdf("resultados/traza_lmm_bayes.nc")

# Exportar resumen a CSV
# summary_df = az.summary(traza_lmm).reset_index()
# summary_df.to_csv("resultados/resumen_lmm_bayes.csv", index=False)
