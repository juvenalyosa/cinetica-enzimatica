"""Instalación del entorno en Google Colab (o en cualquier Linux sin conda).

Uso en la primera celda del notebook::

    !git clone -q https://github.com/juvenalyosa/cinetica-enzimatica.git
    %cd cinetica-enzimatica
    from enzimas import colab_setup
    colab_setup.instalar(modo="rapido")      # solo lectura de datos precalculados (~1 min)
    colab_setup.instalar(modo="completo")    # + OpenMM, ASE, MOPAC (~3 min)

MOPAC se descarga como binario oficial de GitHub (openmopac/mopac) y se deja en
``work/mopac/.../bin/mopac``; :func:`mopac_ejecutable` devuelve la ruta.
"""
from __future__ import annotations

import glob
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

from .textos import t

RAIZ = Path(__file__).resolve().parents[1]
MOPAC_VERSION = "23.2.5"
MOPAC_URL = f"https://github.com/openmopac/mopac/releases/download/v{MOPAC_VERSION}/mopac-{MOPAC_VERSION}-linux.tar.gz"
PAQUETES_RAPIDO = ["numpy", "scipy", "pandas", "matplotlib", "py3Dmol", "ipywidgets"]
PAQUETES_COMPLETO = PAQUETES_RAPIDO + ["openmm", "ase", "mdtraj"]


def _pip(paquetes):
    faltan = []
    for p in paquetes:
        modulo = {"py3Dmol": "py3Dmol", "openmm": "openmm", "ase": "ase", "mdtraj": "mdtraj"}.get(p, p)
        try:
            __import__(modulo)
        except ImportError:
            faltan.append(p)
    if faltan:
        print(t("instalando:"), " ".join(faltan), flush=True)
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", *faltan], check=True)


def en_colab():
    try:
        import google.colab  # noqa: F401
        return True
    except ImportError:
        return False


def mopac_ejecutable():
    """Ruta a MOPAC: variable MOPAC_EXE, luego el PATH, luego la descarga local en work/mopac."""
    env = os.environ.get("MOPAC_EXE")
    if env and Path(env).exists():
        return env
    en_path = shutil.which("mopac")
    if en_path:
        return en_path
    local = sorted(glob.glob(str(RAIZ / "work" / "mopac" / "mopac-*-linux" / "bin" / "mopac")))
    return local[-1] if local else None


def instalar_mopac():
    exe = mopac_ejecutable()
    if exe:
        return exe
    if platform.system() != "Linux":
        raise RuntimeError("Descarga automática de MOPAC solo para Linux; en macOS/Windows instálalo desde "
                           "https://github.com/openmopac/mopac/releases y define MOPAC_EXE.")
    destino = RAIZ / "work" / "mopac"
    destino.mkdir(parents=True, exist_ok=True)
    tgz = destino / "mopac-linux.tar.gz"
    print(t("descargando MOPAC"), MOPAC_VERSION, flush=True)
    urllib.request.urlretrieve(MOPAC_URL, tgz)
    with tarfile.open(tgz) as tf:
        tf.extractall(destino)
    exe = mopac_ejecutable()
    if exe is None:
        raise RuntimeError("no se encontró el binario de MOPAC tras la descarga")
    os.chmod(exe, 0o755)
    # la biblioteca compartida vive junto al binario
    lib = str(Path(exe).parents[1] / "lib")
    os.environ["LD_LIBRARY_PATH"] = lib + ":" + os.environ.get("LD_LIBRARY_PATH", "")
    return exe


def instalar(modo="rapido", idioma=None):
    """``rapido``: solo gráficos y datos precalculados.  ``completo``: además OpenMM, ASE y MOPAC.

    ``idioma`` ("es"/"en", opcional) elige el idioma de los mensajes y del resto del curso (``textos.usar``).
    """
    if idioma is not None:
        from .textos import usar
        usar(idioma)
    modo = str(modo).lower()
    _pip(PAQUETES_COMPLETO if modo == "completo" else PAQUETES_RAPIDO)
    if str(RAIZ) not in sys.path:
        sys.path.insert(0, str(RAIZ))
    info = dict(modo=modo, colab=en_colab(), python=platform.python_version(), mopac=None)
    if modo == "completo":
        try:
            info["mopac"] = instalar_mopac()
            out = subprocess.run([info["mopac"], "--version"], capture_output=True, text=True, timeout=30)
            info["mopac_version"] = (out.stdout or out.stderr).strip().splitlines()[0]
        except Exception as exc:  # el notebook sigue funcionando en modo rápido
            info["mopac_error"] = str(exc)
    print(t("entorno listo:"), info, flush=True)
    return info
