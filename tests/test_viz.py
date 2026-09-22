import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from enzimas import kinetics as kin  # noqa: E402
from enzimas import viz  # noqa: E402

VMAX, KM = 10.0, 7.5


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def mm_data():
    rng = np.random.default_rng(0)
    s = np.array([0.5, 1, 2, 4, 6, 8, 12, 16, 24, 32])
    v = kin.michaelis_menten(s, VMAX, KM) * (1 + rng.normal(0, 0.02, s.size))
    return s, v, kin.fit_michaelis_menten(s, v)


@pytest.fixture
def md_df():
    rng = np.random.default_rng(1)
    t = np.linspace(0, 100, 200)
    return pd.DataFrame({
        "tiempo_ps": t,
        "temperatura_K": 300 + rng.normal(0, 2, t.size),
        "rmsd_CA_A": 1.2 * (1 - np.exp(-t / 30)) + rng.normal(0, 0.05, t.size),
        "d_PG_O6_A": 3.4 + rng.normal(0, 0.2, t.size),
        "d_PG_O3B_A": 1.6 + rng.normal(0, 0.02, t.size),
        "d_O6_OD1asp205_A": 2.7 + rng.normal(0, 0.15, t.size),
        "d_Mg_PG_A": 3.0 + rng.normal(0, 0.1, t.size),
        "d_NZlys169_O6_A": 3.2 + rng.normal(0, 0.2, t.size),
    })


def _check(fig, n_axes=1):
    assert isinstance(fig, Figure)
    assert len(fig.axes) >= n_axes
    return fig


# ---------------------------------------------------------------------------
# Sistema de diseño
# ---------------------------------------------------------------------------
def test_style_and_palette():
    path = viz.apply_style()
    assert path.endswith("enzimas.mplstyle")
    assert matplotlib.rcParams["axes.facecolor"] == viz.SURFACE
    assert not matplotlib.rcParams["axes.spines.top"]
    assert not matplotlib.rcParams["legend.frameon"]
    assert matplotlib.rcParams["lines.linewidth"] == 2.0
    cycle = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
    assert [c.lower() for c in cycle] == viz.PALETTE
    assert len(viz.PALETTE) == 8 and viz.PALETTE[0] == "#2a78d6" and viz.PALETTE[-1] == "#e34948"
    assert viz.COLORS["reactivo"] == viz.PALETTE[0] and viz.COLORS["ts"] == viz.PALETTE[1]
    assert set(viz.STATUS) == {"good", "warning", "critical"}
    assert not set(viz.STATUS.values()) & set(viz.PALETTE)
    ramp = viz.sequential_blue(4)
    assert len(ramp) == 4 and len(set(ramp)) == 4 and ramp[-1] == "#0d366b"
    assert viz.sequential_blue(1) == ["#0d366b"]


def test_figure_and_save(tmp_path):
    fig, ax = viz.figure(5, 3)
    assert isinstance(fig, Figure) and fig.get_size_inches().tolist() == [5, 3]
    ax.plot([0, 1], [0, 1])
    out = viz.save_figure(fig, tmp_path / "sub" / "fig.png", dpi=50)
    assert (tmp_path / "sub" / "fig.png").exists() and str(out).endswith("fig.png")


def test_plot_on_external_axes_returns_parent_figure():
    fig, ax = plt.subplots()
    out = viz.plot_energy_profile(np.linspace(0, 1, 5), [0, 4, 9, 3, -2], ax=ax)
    assert out is fig
    assert not ax.spines["top"].get_visible()


# ---------------------------------------------------------------------------
# Energía
# ---------------------------------------------------------------------------
def test_plot_energy_profile_variants():
    x = np.linspace(-1.5, 1.5, 13)
    e = 18 * np.exp(-x ** 2 / 0.35) - 2 * x
    fig = _check(viz.plot_energy_profile(x, e, xlabel="ξ (Å)", title="Perfil", subtitle="PM7",
                                         labels=["R"] + [None] * 11 + ["P"]))
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert any("ΔE‡" in t and "kcal/mol" in t for t in texts)
    _check(viz.plot_energy_profile([0, 1, 2], [0, 10, -5], smooth=True, annotate_barrier=False, ts_index=1))
    _check(viz.plot_energy_profile(x, e, relative=False, smooth=False, label="serie"))
    # entrada desordenada: etiquetas y ts_index siguen a sus puntos
    fig = viz.plot_energy_profile([2, 0, 1], [-5, 0, 10], ts_index=2, labels=["P", "R", "TS"], smooth=False)
    ax = fig.axes[0]
    by_text = {t.get_text(): t.xy for t in ax.texts if t.get_text() in ("P", "R", "TS")}
    assert by_text["R"][0] == 0 and by_text["TS"][0] == 1 and by_text["P"][0] == 2
    ts_marker = [ln for ln in ax.lines if ln.get_markerfacecolor() == viz.COLORS["ts"]][0]
    assert ts_marker.get_xdata()[0] == 1 and ts_marker.get_ydata()[0] == 10


def test_plot_energy_profiles():
    x = np.linspace(0, 1, 9)
    fig = _check(viz.plot_energy_profiles([(x, np.sin(np.pi * x) * 10, "enzima"), (x, np.sin(np.pi * x) * 20, "agua")]))
    assert fig.axes[0].get_legend() is not None
    with pytest.raises(ValueError):
        viz.plot_energy_profiles([(x, x, str(i)) for i in range(9)])


def test_plot_energy_levels():
    levels = [("E + S", 0), ("ES", -3), ("TS", 15), ("EP", -5), ("E + P", -8)]
    fig = _check(viz.plot_energy_levels(levels, title="Niveles", subtitle="con y sin enzima",
                                        compare=[("S", 0), ("S", 0), ("TS‡", 24), ("P", -8), ("P", -8)]))
    ax = fig.axes[0]
    assert [t.get_text() for t in ax.get_xticklabels()] == ["E + S", "ES", "TS", "EP", "E + P"]
    assert any("ΔE‡ = 18.0" in t.get_text() for t in ax.texts)
    assert ax.get_legend() is not None
    _check(viz.plot_energy_levels(levels[:2], connect=False, show_values=False))


# ---------------------------------------------------------------------------
# Dinámica molecular
# ---------------------------------------------------------------------------
def test_plot_md_timeseries(md_df):
    fig = _check(viz.plot_md_timeseries(md_df, ["d_PG_O6_A", "d_O6_OD1asp205_A"], ylabel="d (Å)", title="MD"))
    assert fig.axes[0].get_legend() is not None
    fig = _check(viz.plot_md_timeseries(md_df, "rmsd_CA_A", ylabel="RMSD (Å)"))
    assert fig.axes[0].get_legend() is None


def test_plot_md_summary(md_df):
    fig = _check(viz.plot_md_summary(md_df, crystal_value=3.1, subtitle="prueba"), n_axes=4)
    assert len(fig.axes) == 4
    _check(viz.plot_md_summary(md_df.drop(columns=["rmsd_CA_A"])), n_axes=4)


# ---------------------------------------------------------------------------
# Cinética
# ---------------------------------------------------------------------------
def test_plot_michaelis_menten(mm_data):
    s, v, fit = mm_data
    fig = _check(viz.plot_michaelis_menten(s, v, fit, title="MM", subtitle="datos"))
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert "Vmax" in texts and "Km" in texts and "Vmax/2" in texts
    assert any("R²" in t for t in texts)
    _check(viz.plot_michaelis_menten(s, v))  # solo puntos
    _check(viz.plot_michaelis_menten(s, v, fit, show_km=False, show_params=False, point_label="datos"))


def test_plot_hill_vs_mm():
    s = np.linspace(0, 40, 100)
    fig = _check(viz.plot_hill_vs_mm(s, kin.michaelis_menten(s, VMAX, KM), kin.hill(s, VMAX, KM, 1.7), n_hill=1.7, s_half=KM))
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert labels == ["Michaelis–Menten (n = 1)", "Hill (n = 1.7)"]


def test_plot_linearizations(mm_data):
    s, v, _ = mm_data
    fig = _check(viz.plot_linearizations(s, v), n_axes=3)
    assert fig.axes[0].get_xlabel() == "1/[S]" and fig.axes[0].get_ylabel() == "1/v₀"
    assert fig.axes[1].get_xlabel() == "v₀/[S]" and fig.axes[2].get_ylabel() == "[S]/v₀"
    fits = {"lineweaver_burk": kin.linear_fit(*kin.lineweaver_burk(s, v))}
    _check(viz.plot_linearizations(s, v, fits=fits, subtitle="x"), n_axes=3)


@pytest.mark.parametrize("kind", ["competitive", "uncompetitive", "noncompetitive", "mixed"])
def test_plot_inhibition(kind, mm_data):
    s, _, _ = mm_data
    fn = {
        "competitive": lambda i: kin.competitive_inhibition(s, VMAX, KM, i, 2.0),
        "uncompetitive": lambda i: kin.uncompetitive_inhibition(s, VMAX, KM, i, 2.0),
        "noncompetitive": lambda i: kin.noncompetitive_inhibition(s, VMAX, KM, i, 2.0),
        "mixed": lambda i: kin.mixed_inhibition(s, VMAX, KM, i, 2.0, 5.0),
    }[kind]
    curves = [(i, fn(i)) for i in (0, 1, 2, 5)]
    fig = _check(viz.plot_inhibition_family(s, curves, kind))
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert labels[0].startswith("[I] = 0 mM") and labels[-1] == "[I] = 5 mM"
    assert fig.axes[0].get_title(loc="left").startswith("Inhibición")
    fig = _check(viz.plot_inhibition_lineweaver(s, curves, kind))
    assert fig.axes[0].get_xlabel() == "1/[S]"


def test_plot_ode_and_initial_rates():
    sim = kin.simulate_mechanism(0.01, 10.0, 5.0, 2.0, 1.0, 40.0, n_points=200)
    fig = _check(viz.plot_ode_simulation(sim, t_zoom=2.0), n_axes=2)
    assert fig.axes[1].get_xlim()[1] == pytest.approx(2.0)
    _check(viz.plot_ode_simulation(sim, subtitle="sin zoom"), n_axes=2)
    s_values = np.array([0.5, 1, 2, 4, 8, 16])
    v0 = kin.initial_rates_from_simulation(0.01, s_values, 5.0, 2.0, 1.0, (1.0, 3.0), n_points=100)
    fit = kin.fit_michaelis_menten(s_values, v0)
    fig = _check(viz.plot_initial_rates_from_ode(s_values, v0, fit))
    assert fig.axes[0].get_legend() is not None


def test_plot_arrhenius_and_eyring():
    t = np.array([283.15, 293.15, 298.15, 303.15, 313.15])
    k = kin.arrhenius(t, 1e10, 12.0)
    fig = _check(viz.plot_arrhenius(t, k, kin.fit_arrhenius(t, k)))
    assert fig.axes[0].get_xlabel().startswith("1000/T")
    assert any("Ea = 12.0" in tx.get_text() for tx in fig.axes[0].texts)
    fig = _check(viz.plot_eyring(t, k, kin.fit_eyring(t, k)))
    assert fig.axes[0].get_ylabel() == "ln(k/T)"
    assert any("ΔG‡" in tx.get_text() for tx in fig.axes[0].texts)
    _check(viz.plot_arrhenius(t, k))
    _check(viz.plot_eyring(t, k))


def test_plot_ph_profile_and_hill_plot():
    ph = np.linspace(4, 10, 13)
    v = kin.bell_shaped_ph_profile(ph, VMAX, 6.0, 8.5)
    fig = _check(viz.plot_ph_profile(ph, v, pkas=(6.0, 8.5)))
    assert any("pKa₁" in t.get_text() for t in fig.axes[0].texts)
    dense = np.linspace(4, 10, 100)
    _check(viz.plot_ph_profile(dense, kin.bell_shaped_ph_profile(dense, VMAX, 6.0, 8.5)))
    s = np.logspace(-0.5, 1.8, 12)
    fig = _check(viz.plot_hill_plot(s, kin.hill(s, VMAX, KM, 1.7), VMAX))
    assert any("n = 1.70" in t.get_text() for t in fig.axes[0].texts)


# ---------------------------------------------------------------------------
# 3D y utilidades de geometría (sin py3Dmol)
# ---------------------------------------------------------------------------
def test_mode_animation_frames_and_xyz_block():
    rng = np.random.default_rng(2)
    xyz = rng.normal(size=(5, 3))
    disp = rng.normal(size=(5, 3))
    frames = viz.mode_animation_frames(xyz, disp, n_frames=24, amplitude=0.5)
    assert len(frames) == 24 and all(f.shape == (5, 3) for f in frames)
    np.testing.assert_allclose(frames[0], xyz)
    excursions = [np.max(np.linalg.norm(f - xyz, axis=1)) for f in frames]
    assert max(excursions) == pytest.approx(0.5, rel=1e-6)
    block = viz.xyz_block(["C", "H", "H", "H", "O"], xyz, comment="prueba")
    lines = block.strip().split("\n")
    assert lines[0] == "5" and lines[1] == "prueba" and len(lines) == 7
    assert lines[2].split()[0] == "C" and len(lines[2].split()) == 4
    with pytest.raises(ValueError):
        viz.xyz_block(["C"], xyz)


def test_view_functions_return_py3dmol_views():
    pytest.importorskip("py3Dmol")
    rng = np.random.default_rng(3)
    symbols = ["C", "H", "H", "H", "O", "H"]
    xyz = np.array([[0, 0, 0], [1.09, 0, 0], [-0.36, 1.03, 0], [-0.36, -0.51, 0.89], [-0.5, -0.7, -1.2], [-1.4, -0.9, -1.3]])
    frames = viz.mode_animation_frames(xyz, rng.normal(size=(6, 3)), n_frames=6, amplitude=0.3)
    view = viz.view_frames(frames, symbols, interval_ms=50)
    assert hasattr(view, "addModelsAsFrames") and hasattr(view, "addModel")
    view = viz.view_frames(frames, symbols, bonds_from=0, loop=False)
    assert hasattr(view, "addModelsAsFrames")
    view = viz.view_frames(frames, symbols, bonds_from=[(0, 1), (0, 4)])
    assert hasattr(view, "addModelsAsFrames")
    pdb = (
        "ATOM      1  N   ASP A 201      10.000  10.000  10.000  1.00  0.00           N\n"
        "ATOM      2  CA  ASP A 201      11.400  10.000  10.000  1.00  0.00           C\n"
        "HETATM    3  C1  GLC A 501      12.000  12.000  12.000  1.00  0.00           C\n"
        "HETATM    4 MG    MG A 601      14.000  12.000  12.000  1.00  0.00          MG\n"
        "END\n"
    )
    view = viz.view_complex(pdb, highlight_residues=(("ASP", np.int64(201)),), label_map={201: "Asp205"})
    assert hasattr(view, "addModel") and hasattr(view, "addLabel")
    view = viz.view_qm_region(pdb, pdb, link_atom_indices=[np.int64(1)])
    assert hasattr(view, "addModel")
    view = viz.view_complex(pdb, ligand_resnames=(), show_labels=False)
    assert hasattr(view, "addModel")
