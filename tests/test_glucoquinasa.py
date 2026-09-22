"""Pruebas del modelo QM/MM que no necesitan MOPAC (partición, términos MM y geometría)."""
import numpy as np
import pytest

pytest.importorskip("openmm")

from enzimas import datos  # noqa: E402
from enzimas.glucoquinasa import QM_FORMAL_CHARGE, ModeloGlucoquinasa  # noqa: E402


@pytest.fixture(scope="module")
def modelo(tmp_path_factory):
    try:
        pdb = datos.ruta("md/minimizado_completo.pdb")
        prmtop = datos.ruta("sistema/complejo.prmtop")
    except FileNotFoundError:
        pytest.skip("faltan los datos precalculados")
    return ModeloGlucoquinasa(pdb, prmtop, workdir=str(tmp_path_factory.mktemp("qmmm")))


def _fd(fn, xyz, k, c, h=1e-4):
    xp, xm = xyz.copy(), xyz.copy()
    xp[k, c] += h
    xm[k, c] -= h
    return (fn(xp) - fn(xm)) / (2 * h)


def test_particion(modelo):
    d = modelo.describe()
    assert d["qm_charge"] == QM_FORMAL_CHARGE
    assert d["n_link"] == 4 and d["n_qm"] > 60
    assert len(d["fixed"]) == 2 * d["n_link"]  # H de enlace + ancla por cada corte
    resn = {a["resname"] for a in d["qm_atoms"]}
    assert {"GLC", "ATP", "ASP", "LYS", "THR", "MG"} <= resn
    # la carga total del clúster se conserva
    assert abs(d["mm_charge"] + QM_FORMAL_CHARGE - modelo.q_all[sorted(set(modelo.qm_global) | set(modelo.excluded) | set(modelo.mm_global))].sum()) < 1e-6


def test_gradiente_potencial_y_lj(modelo):
    xyz = modelo.qm_xyz.copy()
    v, dv = modelo.potential(xyz)
    e_lj, g_lj = modelo.lennard_jones(xyz)
    assert np.isfinite(e_lj)
    for k in (modelo.O6, modelo.PG, 0):
        for c in range(3):
            assert abs(_fd(lambda x: modelo.potential(x)[0][k], xyz, k, c) - dv[k, c]) < 1e-3 * max(1.0, abs(dv[k, c]))
            assert abs(_fd(lambda x: modelo.lennard_jones(x)[0], xyz, k, c) - g_lj[k, c]) < 1e-3 * max(1.0, abs(g_lj[k, c]))


def test_coordenada_de_reaccion(modelo):
    xyz = modelo.qm_xyz.copy()
    xi, dxi = modelo.reaction_coordinate(xyz)
    assert xi < 0  # reactivo: PG unido a O3B
    for k in (modelo.PG, modelo.O6, modelo.O3B):
        for c in range(3):
            assert abs(_fd(lambda x: modelo.reaction_coordinate(x)[0], xyz, k, c) - dxi[k, c]) < 1e-5


def test_construccion_producto(modelo):
    p = modelo.build_product(modelo.qm_xyz)
    rc = modelo.reaction_coordinates(p)
    assert abs(rc["d_PG_O6"] - 1.62) < 1e-6
    assert rc["d_PG_O3B"] > 2.5 and rc["xi"] > 0
    assert abs(rc["d_OD1_H"] - 0.98) < 1e-6
    # los enlaces P-O del grupo gamma se conservan tras la reflexión
    for k in (modelo.O1G, modelo.O2G, modelo.O3G):
        d0 = np.linalg.norm(modelo.qm_xyz[k] - modelo.qm_xyz[modelo.PG])
        d1 = np.linalg.norm(p[k] - p[modelo.PG])
        assert abs(d0 - d1) < 1e-6


def test_salida_pdb_xyz(modelo, tmp_path):
    modelo.write_pdb(tmp_path / "qm.pdb", modelo.qm_xyz)
    txt = (tmp_path / "qm.pdb").read_text()
    assert txt.count("HETATM") == modelo.n_qm
    xyz_txt = modelo.xyz_text(modelo.qm_xyz, "prueba")
    assert xyz_txt.splitlines()[0] == str(modelo.n_qm)
