import numpy as np
import pytest

from enzimas import kinetics as kin


def test_michaelis_menten_limits():
    assert kin.michaelis_menten(1e9, 10.0, 2.0) == pytest.approx(10.0, rel=1e-6)
    assert kin.michaelis_menten(2.0, 10.0, 2.0) == pytest.approx(5.0)
    assert kin.michaelis_menten(0.0, 10.0, 2.0) == 0.0


def test_hill_reduces_to_michaelis_menten():
    s = np.linspace(0.0, 20.0, 50)
    np.testing.assert_allclose(kin.hill(s, 3.0, 4.0, 1.0), kin.michaelis_menten(s, 3.0, 4.0))
    assert kin.hill(4.0, 3.0, 4.0, 2.5) == pytest.approx(1.5)


def test_apparent_parameters_and_inhibition_laws():
    vmax, km, i, ki, kip = 10.0, 2.0, 3.0, 1.5, 4.0
    alpha, alpha_p = 1 + i / ki, 1 + i / kip
    c = kin.apparent_parameters("competitive", vmax, km, i, ki)
    assert c == {"vmax_app": vmax, "km_app": pytest.approx(km * alpha)}
    u = kin.apparent_parameters("uncompetitive", vmax, km, i, ki)
    assert u["vmax_app"] == pytest.approx(vmax / alpha) and u["km_app"] == pytest.approx(km / alpha)
    n = kin.apparent_parameters("noncompetitive", vmax, km, i, ki)
    assert n["vmax_app"] == pytest.approx(vmax / alpha) and n["km_app"] == pytest.approx(km)
    m = kin.apparent_parameters("mixed", vmax, km, i, ki, kip)
    assert m["vmax_app"] == pytest.approx(vmax / alpha_p)
    assert m["km_app"] == pytest.approx(km * alpha / alpha_p)

    s = np.linspace(0.1, 10, 20)
    np.testing.assert_allclose(kin.competitive_inhibition(s, vmax, km, i, ki), vmax * s / (km * alpha + s))
    np.testing.assert_allclose(kin.uncompetitive_inhibition(s, vmax, km, i, ki), vmax * s / (km + s * alpha))
    np.testing.assert_allclose(kin.noncompetitive_inhibition(s, vmax, km, i, ki), vmax * s / ((km + s) * alpha))
    np.testing.assert_allclose(kin.mixed_inhibition(s, vmax, km, i, ki, kip), vmax * s / (km * alpha + s * alpha_p))
    np.testing.assert_allclose(kin.mixed_inhibition(s, vmax, km, i, ki, ki), kin.noncompetitive_inhibition(s, vmax, km, i, ki))
    with pytest.raises(ValueError):
        kin.apparent_parameters("mixed", vmax, km, i, ki)


def test_linearizations_recover_parameters():
    vmax, km = 8.0, 3.0
    s = np.array([0.0, 0.5, 1, 2, 4, 8, 16])
    v = kin.michaelis_menten(s, vmax, km)
    x, y = kin.lineweaver_burk(s, v)
    assert len(x) == 6
    fit = kin.linear_fit(x, y)
    assert 1 / fit["intercept"] == pytest.approx(vmax)
    assert fit["slope"] / fit["intercept"] == pytest.approx(km)
    assert fit["r2"] == pytest.approx(1.0)
    fit = kin.linear_fit(*kin.eadie_hofstee(s, v))
    assert -fit["slope"] == pytest.approx(km) and fit["intercept"] == pytest.approx(vmax)
    fit = kin.linear_fit(*kin.hanes_woolf(s, v))
    assert 1 / fit["slope"] == pytest.approx(vmax) and fit["intercept"] * vmax == pytest.approx(km)


def test_hill_plot_slope_is_n():
    s = np.logspace(-1, 2, 30)
    v = kin.hill(s, 5.0, 7.5, 1.7)
    x, y = kin.hill_plot(np.concatenate([[0.0], s]), np.concatenate([[0.0], v]), 5.0)
    assert len(x) == 30
    assert kin.linear_fit(x, y)["slope"] == pytest.approx(1.7)


def test_fit_michaelis_menten_and_hill_recover_parameters():
    rng = np.random.default_rng(0)
    s = np.array([0.25, 0.5, 1, 2, 4, 8, 16, 32])
    v = kin.michaelis_menten(s, 12.0, 3.0) * (1 + 0.01 * rng.standard_normal(s.size))
    fit = kin.fit_michaelis_menten(s, v)
    assert fit["vmax"] == pytest.approx(12.0, rel=0.05)
    assert fit["km"] == pytest.approx(3.0, rel=0.05)
    assert fit["r2"] > 0.99 and fit["residuals"].shape == s.shape
    assert fit["km_err"] > 0

    s = np.logspace(-0.5, 1.8, 15)
    v = kin.hill(s, 60.0, 7.5, 1.7) * (1 + 0.01 * rng.standard_normal(s.size))
    fit = kin.fit_hill(s, v)
    assert fit["vmax"] == pytest.approx(60.0, rel=0.05)
    assert fit["s_half"] == pytest.approx(7.5, rel=0.05)
    assert fit["n"] == pytest.approx(1.7, rel=0.05)
    assert kin.hill_coefficient_from_fit(s, v) == pytest.approx(fit["n"])


@pytest.mark.parametrize("kind,tol", [("competitive", 0.05), ("uncompetitive", 0.05), ("noncompetitive", 0.05), ("mixed", 0.15)])
def test_fit_inhibition_recovers_parameters(kind, tol):
    rng = np.random.default_rng(1)
    vmax, km, ki, kip = 10.0, 2.0, 1.5, 4.0
    s_grid = np.array([0.5, 1, 2, 4, 8, 16])
    i_grid = np.array([0.0, 0.75, 1.5, 3.0])
    s, i = np.meshgrid(s_grid, i_grid)
    s, i = s.ravel(), i.ravel()
    if kind == "mixed":
        v = kin.mixed_inhibition(s, vmax, km, i, ki, kip)
    else:
        v = getattr(kin, f"{kind}_inhibition")(s, vmax, km, i, ki)
    v = v * (1 + 0.01 * rng.standard_normal(v.size))
    fit = kin.fit_inhibition(s, v, i, kind)
    assert fit["vmax"] == pytest.approx(vmax, rel=tol)
    assert fit["km"] == pytest.approx(km, rel=tol)
    assert fit["ki"] == pytest.approx(ki, rel=tol)
    if kind == "mixed":
        assert fit["ki_prime"] == pytest.approx(kip, rel=tol)
    assert fit["r2"] > 0.99


def test_simulation_conserves_mass_and_reproduces_michaelis_menten():
    k1, km1, k2, e0 = 100.0, 50.0, 10.0, 1e-3
    ss = kin.steady_state_parameters(k1, km1, k2)
    assert ss == {"km": pytest.approx(0.6), "kcat": 10.0, "kd": pytest.approx(0.5)}

    sim = kin.simulate_mechanism(e0, 2.0, k1, km1, k2, 0.5)
    np.testing.assert_allclose(sim["E"] + sim["ES"], e0, rtol=1e-6)
    np.testing.assert_allclose(sim["S"] + sim["ES"] + sim["P"], 2.0, rtol=1e-6)
    assert sim["v0"] == pytest.approx(kin.michaelis_menten(2.0, k2 * e0, ss["km"]), rel=0.02)

    s_values = np.array([0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 5.0])
    v0 = kin.initial_rates_from_simulation(e0, s_values, k1, km1, k2, (0.05, 0.3))
    fit = kin.fit_michaelis_menten(s_values, v0)
    assert fit["km"] == pytest.approx(ss["km"], rel=0.10)
    assert fit["vmax"] == pytest.approx(k2 * e0, rel=0.10)


def test_hill_coefficient_from_ratio_and_substrate_at_fraction():
    vmax, s_half, n = 60.0, 7.5, 1.7
    s90 = kin.substrate_at_fraction(vmax, s_half, n, 0.9)
    s10 = kin.substrate_at_fraction(vmax, s_half, n, 0.1)
    assert kin.hill(s90, vmax, s_half, n) == pytest.approx(0.9 * vmax)
    assert kin.hill_coefficient_from_ratio(s90, s10) == pytest.approx(n)
    assert kin.substrate_at_fraction(vmax, s_half, n, 0.5) == pytest.approx(s_half)


def test_eyring_round_trip_and_fits():
    k = kin.eyring_rate(15.0, 298.15)
    assert kin.barrier_from_rate(k, 298.15) == pytest.approx(15.0)
    assert kin.eyring_rate(0.0) == pytest.approx(kin.KB_J * 298.15 / kin.H_J)

    temps = np.linspace(280.0, 320.0, 9)
    dh, ds = 12.0, -8.0  # kcal/mol, cal/mol/K
    rates = kin.eyring_rate(dh - temps * ds / 1000.0, temps)
    fit = kin.fit_eyring(temps, rates)
    assert fit["delta_h_kcal"] == pytest.approx(dh, rel=1e-6)
    assert fit["delta_s_cal"] == pytest.approx(ds, rel=1e-5)
    assert fit["delta_g_kcal"] == pytest.approx(dh - 298.15 * ds / 1000.0, rel=1e-6)

    rates = kin.arrhenius(temps, 1e10, 14.0)
    fit = kin.fit_arrhenius(temps, rates)
    assert fit["a"] == pytest.approx(1e10, rel=1e-6) and fit["ea_kcal"] == pytest.approx(14.0)
    assert kin.rate_enhancement(60.0, 1e-6) == pytest.approx(6e7)
    assert kin.catalytic_proficiency(60.0 / 7.5e-3, 1e-6) == pytest.approx(8e9)


def test_ph_profile_maximum_near_mean_pka():
    ph = np.linspace(3.0, 11.0, 8001)
    v = kin.bell_shaped_ph_profile(ph, 1.0, 6.0, 9.0)
    assert ph[np.argmax(v)] == pytest.approx(7.5, abs=0.01)
    assert v.max() < 1.0


def test_glucokinase_reference():
    ref = kin.glucokinase_reference()
    assert ref["s_half_mm"]["approx"] and ref["hill_n"]["value"] == pytest.approx(1.7, abs=0.1)
    assert "fuentes primarias" in ref["notes"]
    assert all("source" in inh for inh in ref["inhibitors"]) and ref["references"]
    assert all(v["source"] for k, v in ref.items() if isinstance(v, dict) and "value" in v)


# ---------------------------------------------------------------------------
# Inhibidores: Dixon, Cheng–Prusoff, IC50, gráficos secundarios
# ---------------------------------------------------------------------------
def test_cheng_prusoff_identities():
    ic50, s, km = 3.0, 2.0, 2.0
    assert kin.cheng_prusoff(ic50, s, km, "competitive") == pytest.approx(ic50 / 2)
    assert kin.cheng_prusoff(ic50, s, km, "uncompetitive") == pytest.approx(ic50 / 2)
    assert kin.cheng_prusoff(ic50, s, km, "noncompetitive") == pytest.approx(ic50)
    assert kin.cheng_prusoff(ic50, 1e-9, km, "competitive") == pytest.approx(ic50, rel=1e-6)
    # IC50 medida a [S] = 3·Km para un competitivo con Ki dado → Cheng–Prusoff recupera Ki
    vmax, km, ki, s = 10.0, 2.0, 1.5, 6.0
    ic50_true = ki * (1 + s / km)
    assert kin.competitive_inhibition(s, vmax, km, ic50_true, ki) == pytest.approx(kin.michaelis_menten(s, vmax, km) / 2)
    assert kin.cheng_prusoff(ic50_true, s, km) == pytest.approx(ki)
    with pytest.raises(ValueError):
        kin.cheng_prusoff(ic50, s, km, "mixed")


def test_dixon_intersection_recovers_ki():
    vmax, km, ki = 10.0, 2.0, 1.5
    s_grid = np.array([1.0, 2.0, 4.0, 8.0])
    i_grid = np.array([0.0, 1.0, 2.0, 4.0])
    s, i = np.meshgrid(s_grid, i_grid)
    s, i = s.ravel(), i.ravel()
    v = kin.competitive_inhibition(s, vmax, km, i, ki)
    dixon = kin.dixon_plot(s, v, i)
    assert set(dixon["lines"]) == set(map(float, s_grid))
    i_arr, inv_v = dixon["lines"][2.0]
    np.testing.assert_allclose(inv_v, 1.0 / kin.competitive_inhibition(2.0, vmax, km, i_arr, ki))
    assert all(fit["r2"] == pytest.approx(1.0) for fit in dixon["fits"].values())
    res = kin.ki_from_dixon(None, dixon, km, vmax)
    assert res["ki"] == pytest.approx(ki)
    np.testing.assert_allclose(res["ki_from_slopes"], ki)
    assert res["ki_from_intersection"] == pytest.approx(ki)
    assert res["intersection"][0] == pytest.approx(-ki) and res["intersection"][1] == pytest.approx(1 / vmax)
    slopes = [dixon["fits"][float(sv)]["slope"] for sv in s_grid]
    assert kin.ki_from_dixon(s_grid, slopes, km, vmax)["ki"] == pytest.approx(ki)


def test_fit_ic50_round_trip():
    rng = np.random.default_rng(3)
    i = np.array([0.0, 0.1, 0.3, 1.0, 2.0, 3.0, 10.0, 30.0, 100.0])
    v = kin.ic50_curve(i, 5.0, 2.0, 1.3) * (1 + 0.01 * rng.standard_normal(i.size))
    fit = kin.fit_ic50(i, v)
    assert fit["v0"] == pytest.approx(5.0, rel=0.03)
    assert fit["ic50"] == pytest.approx(2.0, rel=0.05)
    assert fit["hill"] == pytest.approx(1.3, rel=0.05)
    assert fit["r2"] > 0.99 and fit["residuals"].shape == i.shape
    assert kin.ic50_curve(2.0, 5.0, 2.0, 1.3) == pytest.approx(2.5)


def test_secondary_plots_recover_ki():
    vmax, km, ki = 10.0, 2.0, 1.5
    i = np.array([0.0, 0.5, 1.0, 2.0, 4.0])
    km_app = np.array([kin.apparent_parameters("competitive", vmax, km, x, ki)["km_app"] for x in i])
    res = kin.secondary_plot_competitive(i, km_app, km)
    assert res["ki"] == pytest.approx(ki) and res["r2"] == pytest.approx(1.0)
    assert res["intercept"] == pytest.approx(km)
    vmax_app = np.array([kin.apparent_parameters("uncompetitive", vmax, km, x, ki)["vmax_app"] for x in i])
    res = kin.secondary_plot_uncompetitive(i, vmax_app, vmax)
    assert res["ki"] == pytest.approx(ki) and res["r2"] == pytest.approx(1.0)
    assert res["intercept"] == pytest.approx(1 / vmax)


def test_kcat_turnover_and_activator():
    res = kin.kcat_km_from_fit({"vmax": 6.0, "km": 7.5}, 0.1)
    assert res["kcat"] == pytest.approx(60.0)
    assert res["kcat_over_km"] == pytest.approx(8.0)
    assert res["kcat_over_km_molar"] == pytest.approx(8000.0)
    assert "M⁻¹" in res["units"]
    assert kin.turnover_time(60.0) == pytest.approx(1 / 60)
    s = np.linspace(0, 20, 30)
    np.testing.assert_allclose(kin.activator_effect(s, 60.0, 7.5, 1.7), kin.hill(s, 60.0, 7.5, 1.7))
    np.testing.assert_allclose(kin.activator_effect(s, 60.0, 7.5, 1.7, 1.5, 0.4), kin.hill(s, 90.0, 3.0, 1.7))
