"""Cinética enzimática: leyes de velocidad, linealizaciones, ajustes y simulación mecanística.

Módulo pedagógico (numpy/scipy, sin gráficos) que acompaña al cuaderno sobre la
cinética de la glucoquinasa humana. Todas las funciones que reciben una
concentración de sustrato ``s`` aceptan escalares o arreglos de numpy.
"""

from __future__ import annotations

import re

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

# ---------------------------------------------------------------------------
# Constantes físicas
# ---------------------------------------------------------------------------
R_KCAL = 1.987204259e-3  # kcal / (mol K)
R_J = 8.314462618  # J / (mol K)
KB_J = 1.380649e-23  # J / K
H_J = 6.62607015e-34  # J s
N_A = 6.02214076e23  # 1 / mol

_EPS = 1e-12


def _as_array(x):
    return np.asarray(x, dtype=float)


# ---------------------------------------------------------------------------
# 1. Leyes de velocidad
# ---------------------------------------------------------------------------
def michaelis_menten(s, vmax, km):
    """Ecuación de Michaelis–Menten: ``v = vmax·s / (km + s)``."""
    s = _as_array(s)
    return vmax * s / (km + s)


def hill(s, vmax, s_half, n):
    """Ecuación de Hill: ``v = vmax·s^n / (s_half^n + s^n)``.

    ``s_half`` (S0.5) es la concentración a la que ``v = vmax/2``; ``n`` es el
    coeficiente de Hill (n > 1: cooperatividad positiva, n = 1: Michaelis–Menten).
    """
    s = _as_array(s)
    sn = np.power(s, n)
    return vmax * sn / (np.power(s_half, n) + sn)


def apparent_parameters(kind, vmax, km, i, ki, ki_prime=None):
    """Parámetros aparentes en presencia de inhibidor a concentración ``i``.

    Con ``alpha = 1 + i/ki`` y ``alpha' = 1 + i/ki_prime``:

    - competitive:    ``vmax_app = vmax``,          ``km_app = km·alpha``
    - uncompetitive:  ``vmax_app = vmax/alpha'``,   ``km_app = km/alpha'``  (alpha' usa ki)
    - noncompetitive: ``vmax_app = vmax/alpha``,    ``km_app = km``          (ki = ki')
    - mixed:          ``vmax_app = vmax/alpha'``,   ``km_app = km·alpha/alpha'``

    Devuelve ``dict(vmax_app, km_app)``.
    """
    kind = str(kind).lower()
    alpha = 1.0 + i / ki
    if kind == "competitive":
        return {"vmax_app": vmax, "km_app": km * alpha}
    if kind == "uncompetitive":
        return {"vmax_app": vmax / alpha, "km_app": km / alpha}
    if kind == "noncompetitive":
        return {"vmax_app": vmax / alpha, "km_app": km}
    if kind == "mixed":
        if ki_prime is None:
            raise ValueError("mixed inhibition requires ki_prime")
        alpha_p = 1.0 + i / ki_prime
        return {"vmax_app": vmax / alpha_p, "km_app": km * alpha / alpha_p}
    raise ValueError(f"unknown inhibition kind: {kind!r}")


def competitive_inhibition(s, vmax, km, i, ki):
    """Inhibición competitiva: ``v = vmax·s / (km·(1 + i/ki) + s)``."""
    p = apparent_parameters("competitive", vmax, km, i, ki)
    return michaelis_menten(s, p["vmax_app"], p["km_app"])


def uncompetitive_inhibition(s, vmax, km, i, ki):
    """Inhibición acompetitiva: ``v = vmax·s / (km + s·(1 + i/ki))``."""
    p = apparent_parameters("uncompetitive", vmax, km, i, ki)
    return michaelis_menten(s, p["vmax_app"], p["km_app"])


def noncompetitive_inhibition(s, vmax, km, i, ki):
    """Inhibición no competitiva pura: ``v = vmax·s / ((km + s)·(1 + i/ki))``."""
    p = apparent_parameters("noncompetitive", vmax, km, i, ki)
    return michaelis_menten(s, p["vmax_app"], p["km_app"])


def mixed_inhibition(s, vmax, km, i, ki, ki_prime):
    """Inhibición mixta: ``v = vmax·s / (km·(1 + i/ki) + s·(1 + i/ki'))``."""
    p = apparent_parameters("mixed", vmax, km, i, ki, ki_prime)
    return michaelis_menten(s, p["vmax_app"], p["km_app"])


# ---------------------------------------------------------------------------
# 2. Linealizaciones
# ---------------------------------------------------------------------------
def _positive_mask(s, v):
    s = _as_array(s)
    v = _as_array(v)
    mask = (s > 0) & (v > 0) & np.isfinite(s) & np.isfinite(v)
    return s[mask], v[mask]


def lineweaver_burk(s, v):
    """Lineweaver–Burk: ``1/v = (km/vmax)·(1/s) + 1/vmax``.

    Devuelve ``(x, y) = (1/s, 1/v)``; pendiente = km/vmax, ordenada = 1/vmax.
    Se descartan puntos con s ≤ 0 o v ≤ 0.
    """
    s, v = _positive_mask(s, v)
    return 1.0 / s, 1.0 / v


def eadie_hofstee(s, v):
    """Eadie–Hofstee: ``v = -km·(v/s) + vmax``.

    Devuelve ``(x, y) = (v/s, v)``; pendiente = -km, ordenada = vmax.
    """
    s, v = _positive_mask(s, v)
    return v / s, v


def hanes_woolf(s, v):
    """Hanes–Woolf: ``s/v = (1/vmax)·s + km/vmax``.

    Devuelve ``(x, y) = (s, s/v)``; pendiente = 1/vmax, ordenada = km/vmax.
    """
    s, v = _positive_mask(s, v)
    return s, s / v


def linear_fit(x, y):
    """Regresión lineal por mínimos cuadrados: ``y = slope·x + intercept``.

    Devuelve ``dict(slope, intercept, r2)`` con ``r2 = 1 - SS_res/SS_tot``.
    """
    x = _as_array(x)
    y = _as_array(y)
    slope, intercept = np.polyfit(x, y, 1)
    return {"slope": float(slope), "intercept": float(intercept), "r2": _r2(y, slope * x + intercept)}


def _r2(y, y_pred):
    y = _as_array(y)
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 0.0:
        return 1.0 if ss_res <= _EPS else 0.0
    return 1.0 - ss_res / ss_tot


def hill_plot(s, v, vmax):
    """Gráfico de Hill: ``log10(v/(vmax - v)) = n·log10(s) - n·log10(s_half)``.

    Devuelve ``(log10(s), log10(v/(vmax - v)))`` restringido a ``s > 0`` y
    ``0 < v < vmax``; la pendiente es el coeficiente de Hill ``n``.
    """
    s = _as_array(s)
    v = _as_array(v)
    mask = (s > 0) & (v > 0) & (v < vmax)
    s, v = s[mask], v[mask]
    return np.log10(s), np.log10(v / (vmax - v))


# ---------------------------------------------------------------------------
# 3. Ajustes no lineales
# ---------------------------------------------------------------------------
def _guess_vmax_km(s, v):
    vmax0 = float(np.max(v)) if np.size(v) else 1.0
    if vmax0 <= 0:
        vmax0 = 1.0
    idx = int(np.argmin(np.abs(v - 0.5 * vmax0)))
    km0 = float(s[idx]) if np.size(s) else 1.0
    if km0 <= 0:
        km0 = float(np.median(s[s > 0])) if np.any(s > 0) else 1.0
    return vmax0, km0


def fit_michaelis_menten(s, v):
    """Ajuste no lineal de ``v = vmax·s/(km + s)`` (scipy ``curve_fit``).

    Devuelve ``dict(vmax, km, vmax_err, km_err, r2, residuals)``; los errores
    son desviaciones estándar de la matriz de covarianza.
    """
    s = _as_array(s)
    v = _as_array(v)
    vmax0, km0 = _guess_vmax_km(s, v)
    popt, pcov = curve_fit(michaelis_menten, s, v, p0=[vmax0, km0], bounds=(_EPS, np.inf))
    err = np.sqrt(np.diag(pcov))
    pred = michaelis_menten(s, *popt)
    return {
        "vmax": float(popt[0]),
        "km": float(popt[1]),
        "vmax_err": float(err[0]),
        "km_err": float(err[1]),
        "r2": _r2(v, pred),
        "residuals": v - pred,
    }


def fit_hill(s, v):
    """Ajuste no lineal de ``v = vmax·s^n/(s_half^n + s^n)``.

    Devuelve ``dict(vmax, s_half, n, vmax_err, s_half_err, n_err, r2, residuals)``.
    """
    s = _as_array(s)
    v = _as_array(v)
    vmax0, s_half0 = _guess_vmax_km(s, v)
    popt, pcov = curve_fit(
        hill, s, v, p0=[vmax0, s_half0, 1.0], bounds=([_EPS, _EPS, 0.1], [np.inf, np.inf, 10.0])
    )
    err = np.sqrt(np.diag(pcov))
    pred = hill(s, *popt)
    return {
        "vmax": float(popt[0]),
        "s_half": float(popt[1]),
        "n": float(popt[2]),
        "vmax_err": float(err[0]),
        "s_half_err": float(err[1]),
        "n_err": float(err[2]),
        "r2": _r2(v, pred),
        "residuals": v - pred,
    }


_INHIBITION_MODELS = {
    "competitive": competitive_inhibition,
    "uncompetitive": uncompetitive_inhibition,
    "noncompetitive": noncompetitive_inhibition,
    "mixed": mixed_inhibition,
}


def fit_inhibition(s, v, i, kind):
    """Ajuste global de un modelo de inhibición a datos con varias [I].

    ``s``, ``v`` e ``i`` son arreglos planos de igual longitud. Se ajusta
    ``v = f(s, i; vmax, km, ki[, ki'])`` con ``f`` la ley de velocidad de
    ``kind`` ∈ {competitive, uncompetitive, noncompetitive, mixed}.
    Devuelve ``dict(vmax, km, ki, [ki_prime], *_err, r2, residuals)``.
    """
    kind = str(kind).lower()
    if kind not in _INHIBITION_MODELS:
        raise ValueError(f"unknown inhibition kind: {kind!r}")
    s = _as_array(s)
    v = _as_array(v)
    i = _as_array(i)
    model = _INHIBITION_MODELS[kind]

    control = i <= 0
    vmax0, km0 = _guess_vmax_km(s[control], v[control]) if np.any(control) else _guess_vmax_km(s, v)
    ki0 = float(np.median(i[i > 0])) if np.any(i > 0) else 1.0

    xdata = np.vstack([s, i])
    if kind == "mixed":
        def f(x, vmax, km, ki, ki_prime):
            return model(x[0], vmax, km, x[1], ki, ki_prime)
        p0 = [vmax0, km0, ki0, ki0]
        names = ["vmax", "km", "ki", "ki_prime"]
    else:
        def f(x, vmax, km, ki):
            return model(x[0], vmax, km, x[1], ki)
        p0 = [vmax0, km0, ki0]
        names = ["vmax", "km", "ki"]

    popt, pcov = curve_fit(f, xdata, v, p0=p0, bounds=(_EPS, np.inf), maxfev=20000)
    err = np.sqrt(np.diag(pcov))
    pred = f(xdata, *popt)
    out = {name: float(val) for name, val in zip(names, popt)}
    out.update({f"{name}_err": float(e) for name, e in zip(names, err)})
    out["r2"] = _r2(v, pred)
    out["residuals"] = v - pred
    out["kind"] = kind
    return out


# ---------------------------------------------------------------------------
# 4. Simulación mecanística
# ---------------------------------------------------------------------------
def steady_state_parameters(k1, k_minus1, k2):
    """Parámetros de estado estacionario del mecanismo E + S ⇌ ES → E + P.

    ``km = (k_minus1 + k2)/k1``, ``kcat = k2``, ``kd = k_minus1/k1``.
    """
    return {"km": (k_minus1 + k2) / k1, "kcat": float(k2), "kd": k_minus1 / k1}


def simulate_mechanism(e0, s0, k1, k_minus1, k2, t_end, n_points=400, k_minus2=0.0):
    """Integra las EDO del mecanismo E + S ⇌(k1, k-1) ES ⇌(k2, k-2) E + P.

    ``dE/dt = -k1·E·S + (k_minus1 + k2)·ES - k_minus2·E·P``
    ``dS/dt = -k1·E·S + k_minus1·ES``
    ``dES/dt = k1·E·S - (k_minus1 + k2)·ES + k_minus2·E·P``
    ``dP/dt = k2·ES - k_minus2·E·P``

    Devuelve ``dict(t, E, S, ES, P, v0)`` donde ``v0`` es la velocidad máxima
    de formación de producto ``max(dP/dt)``, alcanzada justo después del
    estado pre-estacionario (en t = 0, dP/dt = 0 porque ES(0) = 0).
    """

    def rhs(_t, y):
        e, s, es, p = y
        binding = k1 * e * s
        rev = k_minus2 * e * p
        return [
            -binding + (k_minus1 + k2) * es - rev,
            -binding + k_minus1 * es,
            binding - (k_minus1 + k2) * es + rev,
            k2 * es - rev,
        ]

    t_eval = np.linspace(0.0, t_end, int(n_points))
    sol = solve_ivp(rhs, (0.0, t_end), [e0, s0, 0.0, 0.0], t_eval=t_eval, method="LSODA", rtol=1e-8, atol=1e-12)
    e, s, es, p = sol.y
    rate = k2 * es - k_minus2 * e * p
    return {"t": sol.t, "E": e, "S": s, "ES": es, "P": p, "v0": float(np.max(rate))}


def initial_rates_from_simulation(e0, s_values, k1, k_minus1, k2, t_window, k_minus2=0.0, n_points=200):
    """Velocidad inicial ``v0`` para cada ``s0`` como pendiente de ``P(t)`` en ``t_window``.

    ``v0 = dP/dt ≈ ΔP/Δt`` (pendiente de la regresión lineal de ``P`` frente a
    ``t``). ``t_window = (t_start, t_end)`` debe empezar tras el estado
    pre-estacionario y terminar antes de que el sustrato se agote. Devuelve un
    arreglo de ``v0``, uno por cada ``s0``.
    """
    t_start, t_end = float(t_window[0]), float(t_window[1])
    rates = []
    for s0 in np.atleast_1d(_as_array(s_values)):
        sim = simulate_mechanism(e0, s0, k1, k_minus1, k2, t_end, n_points=n_points, k_minus2=k_minus2)
        mask = sim["t"] >= t_start
        rates.append(linear_fit(sim["t"][mask], sim["P"][mask])["slope"])
    return np.asarray(rates)


# ---------------------------------------------------------------------------
# 5. Cooperatividad
# ---------------------------------------------------------------------------
def hill_coefficient_from_fit(s, v):
    """Coeficiente de Hill ``n`` por ajuste no lineal de ``v = vmax·s^n/(s_half^n + s^n)``."""
    return fit_hill(s, v)["n"]


def hill_coefficient_from_ratio(s90, s10):
    """Coeficiente de Hill a partir del cociente de saturación: ``n = log10(81) / log10(s90/s10)``."""
    return float(np.log10(81.0) / np.log10(s90 / s10))


def substrate_at_fraction(vmax, s_half, n, fraction):
    """Inversa de Hill: sustrato para ``v = fraction·vmax``: ``s = s_half·(f/(1 - f))^(1/n)``."""
    f = float(fraction)
    if not 0.0 < f < 1.0:
        raise ValueError("fraction must be in (0, 1)")
    return float(s_half * (f / (1.0 - f)) ** (1.0 / n))


# ---------------------------------------------------------------------------
# 6. Teoría del estado de transición
# ---------------------------------------------------------------------------
def eyring_rate(delta_g_kcal, temperature_k=298.15, kappa=1.0):
    """Ecuación de Eyring: ``k = kappa·(kB·T/h)·exp(-ΔG‡/(R·T))`` (ΔG‡ en kcal/mol, k en 1/s)."""
    t = _as_array(temperature_k)
    return kappa * (KB_J * t / H_J) * np.exp(-_as_array(delta_g_kcal) / (R_KCAL * t))


def barrier_from_rate(k, temperature_k=298.15):
    """Inversa de Eyring (kappa = 1): ``ΔG‡ = -R·T·ln(k·h/(kB·T))`` en kcal/mol."""
    t = _as_array(temperature_k)
    return -R_KCAL * t * np.log(_as_array(k) * H_J / (KB_J * t))


def arrhenius(temperature_k, a, ea_kcal):
    """Ecuación de Arrhenius: ``k = A·exp(-Ea/(R·T))`` (Ea en kcal/mol)."""
    t = _as_array(temperature_k)
    return a * np.exp(-ea_kcal / (R_KCAL * t))


def fit_arrhenius(temperatures_k, rates):
    """Ajuste lineal de ``ln k = ln A - (Ea/R)·(1/T)``. Devuelve ``dict(a, ea_kcal, r2)``."""
    t = _as_array(temperatures_k)
    k = _as_array(rates)
    fit = linear_fit(1.0 / t, np.log(k))
    return {"a": float(np.exp(fit["intercept"])), "ea_kcal": -fit["slope"] * R_KCAL, "r2": fit["r2"]}


def fit_eyring(temperatures_k, rates):
    """Ajuste lineal de Eyring: ``ln(k/T) = ln(kB/h) + ΔS‡/R - (ΔH‡/R)·(1/T)``.

    ``ΔH‡ = -pendiente·R``, ``ΔS‡ = R·(ordenada - ln(kB/h))``,
    ``ΔG‡(298.15 K) = ΔH‡ - T·ΔS‡``. Devuelve
    ``dict(delta_h_kcal, delta_s_cal, delta_g_kcal, r2)`` (ΔS‡ en cal/mol/K).
    """
    t = _as_array(temperatures_k)
    k = _as_array(rates)
    fit = linear_fit(1.0 / t, np.log(k / t))
    delta_h = -fit["slope"] * R_KCAL
    delta_s_cal = R_KCAL * 1000.0 * (fit["intercept"] - np.log(KB_J / H_J))
    delta_g = delta_h - 298.15 * delta_s_cal / 1000.0
    return {
        "delta_h_kcal": float(delta_h),
        "delta_s_cal": float(delta_s_cal),
        "delta_g_kcal": float(delta_g),
        "r2": fit["r2"],
    }


def catalytic_proficiency(kcat_over_km, k_uncat):
    """Proficiencia catalítica: ``(kcat/km) / k_uncat`` (unidades 1/M)."""
    return kcat_over_km / k_uncat


def rate_enhancement(kcat, k_uncat):
    """Aceleración de la reacción: ``kcat / k_uncat`` (adimensional)."""
    return kcat / k_uncat


# ---------------------------------------------------------------------------
# 7. Perfil de pH
# ---------------------------------------------------------------------------
def bell_shaped_ph_profile(ph, vmax, pka1, pka2):
    """Perfil de pH en campana: ``v = vmax / (1 + 10^(pka1 - pH) + 10^(pH - pka2))``."""
    ph = _as_array(ph)
    return vmax / (1.0 + 10.0 ** (pka1 - ph) + 10.0 ** (ph - pka2))


# ---------------------------------------------------------------------------
# 8. Valores de referencia de la glucoquinasa humana
# ---------------------------------------------------------------------------
def glucokinase_reference():
    """Valores experimentales de la glucoquinasa humana (hexoquinasa IV, GCK) con su fuente primaria.

    Cada parámetro es un diccionario ``{label, value, unit, conditions, source}``.  Los valores
    se citan tal como los publican los autores; las conversiones marcadas como "derivado" son
    aritmética nuestra (p. ej. ΔG‡ a partir de k_cat con la ecuación de Eyring).  Los números
    varían entre laboratorios (k_cat entre 38 y 66 s⁻¹ para la misma enzima), así que deben
    leerse como valores representativos y contrastarse con las fuentes primarias antes de usarse
    cuantitativamente.
    """
    V12 = "Valentínová et al., PLoS One 2012, 7:e34541 (GST-GCK, 5 mM ATP)"
    S09 = "Sayed et al., Diabetes 2009, 58:1419 (GST-GCK purificada)"
    return _traducir_referencia({
        "enzyme": "glucoquinasa humana (hexoquinasa IV, GCK)",
        "s_half_mm": {"label": "S₀.₅ (glucosa)", "value": 7.69, "unit": "mM", "conditions": "silvestre, 0-100 mM glucosa, 5 mM ATP",
                      "source": V12, "approx": True, "alt": {"value": 7.6, "source": S09}},
        "hill_n": {"label": "coeficiente de Hill n (glucosa)", "value": 1.67, "unit": "adimensional", "conditions": "silvestre",
                   "source": V12, "approx": True, "alt": {"value": 1.7, "source": S09}},
        "kcat_s": {"label": "k_cat", "value": 65.8, "unit": "s⁻¹", "conditions": "silvestre; 62.3 s⁻¹ en Sayed 2009; 38 s⁻¹ en Heredia 2006 (stopped-flow)",
                   "source": V12, "approx": True},
        "km_atp_mm": {"label": "K_M (ATP)", "value": 0.44, "unit": "mM", "conditions": "silvestre", "source": V12, "approx": True},
        "dG_ddagger_from_kcat_kcal": {"label": "ΔG‡ que implica k_cat (Eyring, κ = 1)", "value": 15.0, "unit": "kcal/mol",
                                      "conditions": "derivado de k_cat = 62-66 s⁻¹ a 25 °C; 15.6 a 37 °C", "source": "aritmética nuestra"},
        "ph_optimum": {"label": "pH óptimo", "value": 8.6, "unit": "pH", "conditions": "8.5-8.7 tras corregir la acidificación por ATP",
                       "source": "Šimčíková & Heneberg, Sci Rep 2019, 9:11422"},
        "gkrp_ic50_glucose_mm": {"label": "IC₅₀ de GKRP (frente a glucosa)", "value": 13.3, "unit": "mM glucosa", "conditions": "silvestre",
                                 "source": V12},
        "kex_conformational_s": {"label": "k de intercambio conformacional (abierta ⇌ cerrada)", "value": 50.0, "unit": "s⁻¹",
                                 "conditions": "entre 5 y 100 s⁻¹ por RMN, escala de ms; comparable a k_cat", "source": "Larion et al., PLoS Biol 2012, 10:e1001452"},
        "qmmm_barrier_literature_kcal": {"label": "ΔE‡ QM/MM publicada (silvestre)", "value": 18.3, "unit": "kcal/mol",
                                         "conditions": "Asp205 base general; K169A: 32.1 kcal/mol", "source": "Zhang et al., PLoS One 2009, 4:e6304"},
        "uncatalyzed_halflife_years": {"label": "vida media de la hidrólisis no catalizada de un fosfato dianión", "value": 1.1e12, "unit": "años",
                                       "conditions": "k = 2×10⁻²⁰ s⁻¹ a 25 °C (ΔG‡ ≈ 44 kcal/mol); hidrólisis, no transferencia a glucosa",
                                       "source": "Lad, Williams & Wolfenden, PNAS 2003, 100:5607"},
        "inhibitors": [
            {"name": "N-acetilglucosamina (GlcNAc)", "kind": "competitive", "versus": "glucosa", "ki": "no localizada en fuente primaria", "unit": "",
             "note": "inhibidor de alta afinidad; competitivo y baja n hacia 1.0", "source": "Xu et al. 1995; Cárdenas et al. 1984"},
            {"name": "manoheptulosa", "kind": "mixed", "versus": "glucosa", "ki": "discutida (0.25-20 mM según fuentes secundarias)", "unit": "mM",
             "note": "tipo mixto, predominantemente no competitivo/competitivo según Scruel 1998", "source": "Xu et al. 1995; Scruel et al. 1998"},
            {"name": "glucosamina", "kind": "competitive", "versus": "glucosa", "ki": "baja afinidad; también sustrato", "unit": "", "source": "Xu et al. 1995"},
            {"name": "glucosa-6-fosfato (producto)", "kind": "none", "versus": "-", "ki": "sin inhibición fisiológica (a diferencia de HK I-III)", "unit": "",
             "source": "Viñuela, Salas & Sols 1963; Storer & Cornish-Bowden 1977"},
            {"name": "palmitoil-CoA", "kind": "slow-binding", "versus": "-", "ki": "15 µM → >90 % inactivación (horas)", "unit": "",
             "source": "Tippett & Neet 1982; Lin et al. 1989"},
            {"name": "proteína reguladora GKRP", "kind": "competitive", "versus": "glucosa", "ki": "IC₅₀ 13.3 mM glucosa; K_i no localizada", "unit": "",
             "note": "sube S₀.₅ sin cambiar V; F6P refuerza (k_off ÷60), F1P antagoniza", "source": "Vandercammen & Van Schaftingen 1991; Casey & Miller 2016"},
        ],
        "activators": [
            {"name": "RO-28-1675 (RO0281675, GKA)", "fold": 15.8, "ec50": "6.9 µM", "s_half_factor": "baja", "vmax_factor": "sube",
             "note": "activador mixto no esencial; sitio alostérico a ~20 Å del sitio de glucosa", "source": "Grimsby et al., Science 2003; Sayed 2009; Kamata 2004"},
        ],
        "mutants": [
            {"name": "L315H", "kind": "GCK-MODY", "s_half_mm": 8.09, "hill_n": 1.67, "kcat_s": 69.2, "kcat_rel": 69.2 / 65.8},
            {"name": "V244G", "kind": "GCK-MODY", "s_half_mm": 12.77, "hill_n": 1.52, "kcat_s": 75.7, "kcat_rel": 75.7 / 65.8},
            {"name": "G223S", "kind": "GCK-MODY", "s_half_mm": 16.5, "hill_n": 1.46, "kcat_s": 74.5, "kcat_rel": 74.5 / 65.8},
            {"name": "I110N", "kind": "GCK-MODY", "s_half_mm": 14.87, "hill_n": 1.37, "kcat_s": 3.07, "kcat_rel": 3.07 / 65.8},
            {"name": "V200A", "kind": "GCK-MODY", "s_half_mm": 78.3, "hill_n": 1.28, "kcat_s": 56.9, "kcat_rel": 56.9 / 65.8},
            {"name": "W99L", "kind": "activadora (hiperinsulinismo)", "s_half_mm": 2.9, "hill_n": 1.6, "kcat_s": 85.6, "kcat_rel": 85.6 / 62.3},
            {"name": "M197I", "kind": "activadora (hiperinsulinismo)", "s_half_mm": 2.6, "hill_n": 1.6, "kcat_s": 38.1, "kcat_rel": 38.1 / 62.3},
        ],
        "mutants_source": "GCK-MODY: Valentínová et al. 2012 (tabla 3); activadoras: Sayed et al. 2009 (tabla 2). k_cat relativo a la silvestre del mismo estudio.",
        "references": [
            "Valentínová L. et al. (2012) PLoS One 7:e34541. doi:10.1371/journal.pone.0034541",
            "Sayed S. et al. (2009) Diabetes 58:1419-1427. doi:10.2337/db08-1792",
            "Heredia V. V. et al. (2006) Biochemistry 45:7553-7562. doi:10.1021/bi060253q",
            "Larion M. & Miller B. G. (2012) Arch. Biochem. Biophys. 519:103-111. doi:10.1016/j.abb.2011.11.007",
            "Larion M. et al. (2012) PLoS Biol. 10:e1001452. doi:10.1371/journal.pbio.1001452",
            "Xu L. Z. et al. (1995) Biochemistry 34:6083-6092. doi:10.1021/bi00018a011",
            "Kamata K. et al. (2004) Structure 12:429-438. doi:10.1016/j.str.2004.02.005",
            "Petit P. et al. (2011) Acta Cryst. D 67:929-935 (PDB 3FGU)",
            "Grimsby J. et al. (2003) Science 301:370-373. doi:10.1126/science.1084073",
            "Gloyn A. L. et al. (2003) Diabetes 52:2433-2440. doi:10.2337/diabetes.52.9.2433",
            "Vandercammen A. & Van Schaftingen E. (1991) Eur. J. Biochem. 200:545-551. doi:10.1111/j.1432-1033.1991.tb16217.x",
            "Casey A. K. & Miller B. G. (2016) Biochemistry 55:2899-2902. doi:10.1021/acs.biochem.6b00349",
            "Cárdenas M. L. et al. (1984) Eur. J. Biochem. 145:163-171. doi:10.1111/j.1432-1033.1984.tb08536.x",
            "Storer A. C. & Cornish-Bowden A. (1976) Biochem. J. 159:7-14; (1977) Biochem. J. 165:61-69",
            "Scruel O. et al. (1998) Mol. Cell. Biochem. 187:113-120. doi:10.1023/a:1006812300200",
            "Viñuela E., Salas M. & Sols A. (1963) J. Biol. Chem. 238:1175-1177",
            "Tippett P. S. & Neet K. E. (1982) J. Biol. Chem. 257:12839-12845",
            "Šimčíková D. & Heneberg P. (2019) Sci. Rep. 9:11422. doi:10.1038/s41598-019-47883-1",
            "Zhang J. et al. (2009) PLoS One 4:e6304. doi:10.1371/journal.pone.0006304",
            "Lad C., Williams N. H. & Wolfenden R. (2003) PNAS 100:5607-5610. doi:10.1073/pnas.0631607100",
            "Wolfenden R. & Snider M. J. (2001) Acc. Chem. Res. 34:938-945. doi:10.1021/ar000058i",
            "Cornish-Bowden A. (2012) Fundamentals of Enzyme Kinetics, 4.ª ed., Wiley-Blackwell",
            "Senn H. M. & Thiel W. (2009) Angew. Chem. Int. Ed. 48:1198-1229 (métodos QM/MM)",
            "Stewart J. J. P. (2013) J. Mol. Model. 19:1-32 (PM7); Larsen A. H. et al. (2017) J. Phys.: Condens. Matter 29:273002 (ASE)",
        ],
        "notes": (
            "Valores de fuentes primarias (ver 'source' y 'references'); varían entre laboratorios y ensayos "
            "(acoplados a G6PDH, proteínas de fusión GST, temperatura no siempre indicada). Todos los números son "
            "aproximados (approx) y deben verificarse contra las fuentes primarias antes de usarse cuantitativamente."
        ),
    })


# Campos de glucokinase_reference() que el estudiante lee (se traducen si el curso está en inglés); las fuentes
# bibliográficas, los nombres de mutantes y los tipos de inhibición en inglés ("competitive"…) se dejan tal cual.
_CAMPOS_TEXTO = {"enzyme", "label", "unit", "conditions", "name", "versus", "ki", "note", "s_half_factor",
                 "vmax_factor", "mutants_source", "notes"}
_NO_TRADUCIR = {"competitive", "uncompetitive", "noncompetitive", "mixed", "none", "slow-binding", "GCK-MODY", "mM", "s⁻¹",
                "kcal/mol", "pH", ""}


def _traducir_referencia(ref):
    from .textos import idioma, t

    if idioma() == "es":
        return ref

    def tr(valor, clave=None):
        if isinstance(valor, dict):
            return {k: tr(v, k) for k, v in valor.items()}
        if isinstance(valor, list):
            return [tr(v, clave) for v in valor]
        if isinstance(valor, str) and valor not in _NO_TRADUCIR and not re.fullmatch(r"[A-Z]\d+[A-Z]", valor):
            if clave in _CAMPOS_TEXTO or (clave == "kind" and valor.startswith("activadora")) or \
                    (clave == "source" and valor == "aritmética nuestra") or (clave == "references" and "métodos" in valor):
                return t(valor)
        return valor

    return tr(ref)

# ---------------------------------------------------------------------------
# 9. Inhibidores: Dixon, Cheng–Prusoff, IC50 y gráficos secundarios
# ---------------------------------------------------------------------------
def dixon_plot(s, v, i):
    """Gráfico de Dixon: ``1/v`` frente a ``[I]``, una recta por cada ``[S]``.

    Para inhibición competitiva
    ``1/v = (Km/(Vmax·Ki·[S]))·[I] + (1/Vmax)·(1 + Km/[S])``,
    de modo que las rectas de distintas ``[S]`` se cruzan en ``[I] = −Ki`` y
    ``1/v = 1/Vmax``. ``s``, ``v`` e ``i`` son arreglos planos de igual
    longitud (datos globales con varias series de ``[S]``).

    Devuelve ``dict(lines, fits)`` donde ``lines[s] = (i_array, inv_v_array)``
    y ``fits[s]`` es :func:`linear_fit` de esa recta (claves = valores de ``[S]``
    como ``float``). Se descartan los puntos con ``v ≤ 0``.
    """
    s = _as_array(s)
    v = _as_array(v)
    i = _as_array(i)
    lines, fits = {}, {}
    for s_value in np.unique(s):
        mask = (s == s_value) & (v > 0) & np.isfinite(v)
        if np.count_nonzero(mask) < 2:
            continue
        i_s, inv_v = i[mask], 1.0 / v[mask]
        key = float(s_value)
        lines[key] = (i_s, inv_v)
        fits[key] = linear_fit(i_s, inv_v)
    return {"lines": lines, "fits": fits}


def ki_from_dixon(s_values, slopes_or_lines, km, vmax):
    """``Ki`` (competitiva) a partir de las rectas de Dixon.

    De la pendiente de cada recta, ``pendiente_s = Km/(Vmax·Ki·[S])`` →
    ``Ki = Km/(Vmax·pendiente_s·[S])``; y del cruce de cada par de rectas,
    ``[I]_cruce = −Ki``. ``slopes_or_lines`` puede ser el diccionario devuelto
    por :func:`dixon_plot`, un diccionario ``{s: fit}`` o una secuencia de
    pendientes alineada con ``s_values``.

    Devuelve ``dict(ki, ki_from_slopes, ki_from_intersection, intersection)``
    con ``intersection = (i, inv_v)`` promedio de los cruces por pares
    (``None`` si solo hay una recta); ``ki`` es la media de las estimaciones
    por pendiente.
    """
    if isinstance(slopes_or_lines, dict) and "fits" in slopes_or_lines:
        fits = slopes_or_lines["fits"]
        s_values = list(fits.keys())
        fits = [fits[key] for key in s_values]
    elif isinstance(slopes_or_lines, dict):
        s_values = list(slopes_or_lines.keys())
        fits = [slopes_or_lines[key] for key in s_values]
    else:
        fits = [{"slope": float(sl), "intercept": None} for sl in slopes_or_lines]
    s_arr = _as_array(list(s_values))
    slopes = np.array([float(f["slope"]) for f in fits])
    ki_slopes = km / (vmax * slopes * s_arr)

    crossings = []
    for a in range(len(fits)):
        for b in range(a + 1, len(fits)):
            fa, fb = fits[a], fits[b]
            if fa["intercept"] is None or fb["intercept"] is None:
                continue
            d_slope = fa["slope"] - fb["slope"]
            if abs(d_slope) < _EPS:
                continue
            x = (fb["intercept"] - fa["intercept"]) / d_slope
            crossings.append((x, fa["slope"] * x + fa["intercept"]))
    intersection = tuple(map(float, np.mean(crossings, axis=0))) if crossings else None
    ki_cross = -intersection[0] if intersection is not None else None
    return {
        "ki": float(np.mean(ki_slopes)),
        "ki_from_slopes": ki_slopes,
        "ki_from_intersection": ki_cross,
        "intersection": intersection,
    }


def cheng_prusoff(ic50, s, km, kind="competitive"):
    """Ecuación de Cheng–Prusoff: convierte ``IC50`` en ``Ki`` según el tipo de inhibición.

    - competitive:    ``Ki = IC50 / (1 + [S]/Km)``
    - uncompetitive:  ``Ki = IC50 / (1 + Km/[S])``
    - noncompetitive: ``Ki = IC50``
    """
    kind = str(kind).lower()
    if kind == "competitive":
        return ic50 / (1.0 + s / km)
    if kind == "uncompetitive":
        return ic50 / (1.0 + km / s)
    if kind == "noncompetitive":
        return float(ic50)
    raise ValueError(f"unknown inhibition kind: {kind!r}")


def ic50_curve(i, v0, ic50, hill=1.0):
    """Curva dosis–respuesta: ``v = v0 / (1 + ([I]/IC50)^hill)``."""
    i = _as_array(i)
    return v0 / (1.0 + np.power(i / ic50, hill))


def fit_ic50(i, v):
    """Ajuste no lineal de ``v = v0/(1 + ([I]/IC50)^hill)``.

    Devuelve ``dict(v0, ic50, hill, v0_err, ic50_err, hill_err, r2, residuals)``.
    """
    i = _as_array(i)
    v = _as_array(v)
    v0_0 = float(np.max(v)) if np.size(v) else 1.0
    positive = i > 0
    if np.any(positive):
        idx = int(np.argmin(np.abs(v[positive] - 0.5 * v0_0)))
        ic50_0 = float(i[positive][idx])
    else:
        ic50_0 = 1.0
    popt, pcov = curve_fit(
        ic50_curve, i, v, p0=[v0_0, ic50_0, 1.0], bounds=([_EPS, _EPS, 0.1], [np.inf, np.inf, 10.0]), maxfev=20000
    )
    err = np.sqrt(np.diag(pcov))
    pred = ic50_curve(i, *popt)
    return {
        "v0": float(popt[0]),
        "ic50": float(popt[1]),
        "hill": float(popt[2]),
        "v0_err": float(err[0]),
        "ic50_err": float(err[1]),
        "hill_err": float(err[2]),
        "r2": _r2(v, pred),
        "residuals": v - pred,
    }


def secondary_plot_competitive(i_values, km_app_values, km):
    """Gráfico secundario competitivo: ``Km_app = Km·(1 + [I]/Ki)``.

    Regresión de ``Km_app`` frente a ``[I]``: ``pendiente = Km/Ki`` →
    ``Ki = Km/pendiente``. Devuelve ``dict(ki, slope, intercept, r2)``.
    """
    fit = linear_fit(i_values, km_app_values)
    return {"ki": float(km / fit["slope"]), "slope": fit["slope"], "intercept": fit["intercept"], "r2": fit["r2"]}


def secondary_plot_uncompetitive(i_values, vmax_app_values, vmax):
    """Gráfico secundario acompetitivo: ``1/Vmax_app = (1/Vmax)·(1 + [I]/Ki')``.

    Regresión de ``1/Vmax_app`` frente a ``[I]``: ``pendiente = 1/(Vmax·Ki')``
    → ``Ki' = 1/(Vmax·pendiente)``. Devuelve ``dict(ki, slope, intercept, r2)``
    (``ki`` es ``Ki'``).
    """
    fit = linear_fit(i_values, 1.0 / _as_array(vmax_app_values))
    return {"ki": float(1.0 / (vmax * fit["slope"])), "slope": fit["slope"], "intercept": fit["intercept"], "r2": fit["r2"]}


# ---------------------------------------------------------------------------
# 10. Constantes catalíticas y activadores
# ---------------------------------------------------------------------------
def kcat_km_from_fit(fit, e0, km_to_molar=1e-3):
    """``kcat`` y ``kcat/Km`` a partir de un ajuste (``dict`` con ``vmax`` y ``km``) y ``[E]₀``.

    ``kcat = Vmax/[E]₀`` y ``kcat/Km``. Si ``Vmax`` está en µM/s y ``[E]₀`` en
    µM, ``kcat`` queda en 1/s; ``kcat_over_km`` se da en las unidades de ``Km``
    del ajuste y ``kcat_over_km_molar`` en M⁻¹·s⁻¹ usando ``km_to_molar``
    (factor que convierte ``Km`` a mol/L; 1e-3 si ``Km`` está en mM).
    Devuelve ``dict(kcat, kcat_over_km, kcat_over_km_molar, units)``.
    """
    vmax = float(fit["vmax"])
    km = float(fit["km"] if "km" in fit else fit["s_half"])
    kcat = vmax / e0
    return {
        "kcat": kcat,
        "kcat_over_km": kcat / km,
        "kcat_over_km_molar": kcat / (km * km_to_molar),
        "units": (
            "kcat en 1/s si Vmax y [E]₀ comparten unidad de concentración; "
            "kcat_over_km en (unidad de Km)⁻¹·s⁻¹; kcat_over_km_molar en M⁻¹·s⁻¹ "
            f"(Km × {km_to_molar:g} = mol/L)"
        ),
    }


def turnover_time(kcat):
    """Tiempo de recambio: ``τ = 1/kcat`` (segundos por ciclo catalítico si ``kcat`` está en 1/s)."""
    return 1.0 / _as_array(kcat)


def activator_effect(s, vmax, s_half, n, vmax_factor=1.0, s_half_factor=1.0):
    """Curva de Hill en presencia de un activador alostérico (p. ej. un GKA).

    ``v = (vmax·vmax_factor)·s^n / ((s_half·s_half_factor)^n + s^n)``: el
    activador multiplica ``Vmax`` por ``vmax_factor`` (≥ 1) y ``S0.5`` por
    ``s_half_factor`` (≤ 1, mayor afinidad aparente). Con ambos factores = 1 se
    recupera :func:`hill`.
    """
    return hill(s, vmax * vmax_factor, s_half * s_half_factor, n)
