"""Cinética enzimática: leyes de velocidad, linealizaciones, ajustes y simulación mecanística.

Módulo pedagógico (numpy/scipy, sin gráficos) que acompaña al cuaderno sobre la
cinética de la glucoquinasa humana. Todas las funciones que reciben una
concentración de sustrato ``s`` aceptan escalares o arreglos de numpy.
"""

from __future__ import annotations

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
    """Valores de referencia (orden de magnitud) para la glucoquinasa humana (hexoquinasa IV).

    Devuelve un diccionario con parámetros cinéticos aproximados y una lista de
    inhibidores; ``notes`` advierte que deben contrastarse con fuentes primarias.
    """
    return {
        "enzyme": "glucoquinasa humana (hexoquinasa IV, GCK)",
        "s_half_mm": {"value": 7.5, "unit": "mM glucosa", "approx": True},
        "hill_n": {"value": 1.7, "unit": "adimensional", "approx": True},
        "kcat_s": {"value": 60.0, "unit": "1/s", "approx": True},
        "km_atp_mm": {"value": 0.4, "unit": "mM ATP", "approx": True},
        "inhibitors": [
            {"name": "manoheptulosa", "kind": "competitive", "versus": "glucosa", "approx": True},
            {"name": "N-acetilglucosamina", "kind": "competitive", "versus": "glucosa", "approx": True},
            {
                "name": "proteína reguladora de la glucoquinasa (GKRP)",
                "kind": "competitive",
                "versus": "glucosa",
                "context": "hígado",
                "approx": True,
            },
        ],
        "notes": (
            "Magnitudes típicas de libros de texto y literatura para la glucoquinasa "
            "humana a 25-30 °C; todos los números son aproximados (approx) y deben "
            "verificarse contra fuentes primarias antes de usarse cuantitativamente."
        ),
    }
