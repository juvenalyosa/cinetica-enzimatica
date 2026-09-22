import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

ipywidgets = pytest.importorskip("ipywidgets")

from enzimas import interactivo  # noqa: E402

EXPLORERS = [
    "explorar_michaelis_menten",
    "explorar_mecanismo",
    "explorar_hill",
    "explorar_inhibicion",
    "explorar_eyring",
    "explorar_temperatura",
    "explorar_ph",
    "explorar_activador",
    "explorar_perfil_energia",
]


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _first_slider(box):
    return next(w for w in box.controles.values() if isinstance(w, ipywidgets.FloatSlider))


def _nudge(slider):
    """Cambia el valor del deslizador (siempre a uno distinto) para disparar el redibujado."""
    new = slider.min if slider.value != slider.min else slider.max
    slider.value = new
    return new


@pytest.mark.parametrize("name", EXPLORERS)
def test_explorer_returns_widget_and_updates(name):
    box = getattr(interactivo, name)()
    assert isinstance(box, ipywidgets.Widget)
    assert isinstance(box.salida, ipywidgets.Output)
    assert box.estado["error"] is None
    assert isinstance(box.estado["resumen"], str) and box.estado["resumen"]
    first_summary = box.estado["resumen"]
    _nudge(_first_slider(box))
    assert box.estado["error"] is None
    assert box.estado["resumen"] != first_summary
    assert not plt.get_fignums()  # cada figura se cierra tras mostrarse


@pytest.mark.parametrize("name", EXPLORERS)
def test_draw_function_returns_figure_and_summary(name):
    box = getattr(interactivo, name)()
    draw = getattr(interactivo, name.replace("explorar_", "_dibujar_"))
    fig, summary = draw(**{k: w.value for k, w in box.controles.items()})
    assert isinstance(fig, Figure)
    assert isinstance(summary, str) and len(summary) > 20
    assert all(ax.get_xlabel() for ax in fig.axes)


def test_inhibition_dropdown_and_all_kinds():
    box = interactivo.explorar_inhibicion()
    kind = box.controles["kind"]
    assert box.controles["ki_prime"].disabled
    for value in ("uncompetitive", "noncompetitive", "mixed"):
        kind.value = value
        assert box.estado["error"] is None
        assert box.controles["ki_prime"].disabled == (value != "mixed")
    box.controles["i"].value = 0.0  # sin inhibidor: una sola curva, sin excepción
    assert box.estado["error"] is None


def test_energy_profile_clamps_barrier_below_reaction_energy():
    fig, summary = interactivo._dibujar_perfil_energia(5.0, 8.0)
    assert isinstance(fig, Figure)
    assert "se elevó" in summary


def test_explorer_helper_reports_errors():
    def bad(x):
        raise RuntimeError("boom")

    slider = ipywidgets.FloatSlider(value=1.0, min=0.0, max=2.0)
    with pytest.raises(RuntimeError):
        interactivo._explorer(bad, {"x": slider})


def test_tiempo_formatting():
    assert interactivo._tiempo(2.5e-9) == "2.5 ns"
    assert interactivo._tiempo(0.0159) == "15.9 ms"
    assert interactivo._tiempo(90.0) == "1.5 min"
    assert interactivo._tiempo(3.15576e8) == "10 años"


def test_todo_builds_tab_with_all_explorers():
    tab = interactivo.todo()
    assert isinstance(tab, ipywidgets.Tab)
    assert len(tab.children) == len(EXPLORERS)
    assert tab.get_title(0) == "Michaelis–Menten"
    assert all(child.estado["error"] is None for child in tab.children)
