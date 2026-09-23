"""Visores moleculares 3D del curso: la enzima, la reacción y la dinámica, manipulables con el ratón.

Cada función devuelve un ``IPython.display.HTML`` autónomo (3Dmol.js desde CDN + un reproductor propio):

* :func:`complejo_cristal`  la glucoquinasa del cristal 3FGU: dominios, sitio activo, superficie.
* :func:`region_qm`         la partición QM/MM: región cuántica y cargas de la enzima.
* :func:`pelicula_reaccion` la transferencia del fosfato **sincronizada con el perfil de energía**.
* :func:`modo_imaginario`   la vibración del estado de transición (la frecuencia imaginaria).
* :func:`pelicula_md`       la dinámica molecular del sitio activo sincronizada con d(Pγ–O6).

Los enlaces NO se adivinan por distancia: salen de la topología de Amber (``qmmm/topologia_qm.json`` y
``md/pelicula_md.json.gz``, generados por ``scripts/05_datos_visor.py``). Los cuatro enlaces que cambian en la
reacción (Pγ–O3β, Pγ–O6, O6–H y H–O de Asp205) se dibujan aparte, con un grosor y un brillo que siguen su
orden de enlace de Pauling, n = exp((r₀ − r)/0.6): enteros cuando están formados, discontinuos y tenues cuando
están a medio formar o romper, e invisibles cuando ya no existen. Así, en el estado de transición se ve al fósforo
«entre» los dos oxígenos, que es exactamente lo que ocurre. La coordinación del Mg²⁺ (no covalente) va en verde
discontinuo fino.
"""
from __future__ import annotations

import base64
import functools
import gzip
import itertools
import json
from pathlib import Path

import numpy as np

from . import datos

URL_3DMOL = "https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js"
DESPLAZAMIENTO = 4                     # numeración de tleap = cristal − 4

# paleta (misma identidad que viz.py y las ilustraciones)
AZUL, NARANJA, AGUA, AMARILLO, MAGENTA, VERDE, VIOLETA, ROJO = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
ELEMENTOS = {"O": "#ff5a5f", "N": "#6b9cff", "P": "#ff9f1c", "S": "#ffd166", "H": "#f4f4f6", "MG": "#3ddc84",
             "Mg": "#3ddc84", "K": "#b69cff", "NA": "#b69cff"}
CARBONO = {"glc": "#ffb38a", "atp": "#ffd27a", "prot": "#c9d4e6", "asp": "#9ec5ff", "lys": "#c7b8ff",
           "thr": "#a8e6cf", "otro": "#c9d4e6"}
ROMPE, FORMA, COORD = "#ff4d6d", "#2ee6a6", "#3ddc84"   # enlaces que se rompen / se forman; coordinación del Mg
_contador = itertools.count()


# ---------------------------------------------------------------------------------------------- utilidades
def _b64_int16(coords):
    """Coordenadas (Å) → base64 de int16 en centésimas de Å (≈ 4 veces más ligero que texto)."""
    arr = np.round(np.asarray(coords, dtype=float) * 100.0).astype("<i2")
    return base64.b64encode(arr.tobytes()).decode("ascii")


def _leer_pdb(texto, solo=None):
    atomos = []
    for l in texto.splitlines():
        if l.startswith(("ATOM", "HETATM")):
            elem = (l[76:78].strip() or l[12:16].strip()[0]).capitalize()
            a = dict(n=l[12:16].strip(), r=l[17:20].strip(), c=l[21], s=int(l[22:26]), e=elem, het=l.startswith("HETATM"),
                     x=float(l[30:38]), y=float(l[38:46]), z=float(l[46:54]),
                     q=float(l[60:66]) if l[60:66].strip() else 0.0)
            if solo is None or solo(a):
                atomos.append(a)
    return atomos


def _rol(resname, resseq_cristal):
    r = resname.upper()
    if r in ("GLC", "BGC"):
        return "glc"
    if r in ("ATP", "ANP", "ADP"):
        return "atp"
    if r == "ASP" and resseq_cristal == 205:
        return "asp"
    if r == "LYS" and resseq_cristal == 169:
        return "lys"
    if r == "THR" and resseq_cristal == 228:
        return "thr"
    return "prot"


@functools.lru_cache(maxsize=1)
def _orientacion_sitio():
    """(R, origen) común a todos los visores del sitio activo (cristal minimizado).

    La cámara mira desde la boca de la hendidura, pero perpendicular al salto del fosfato: el eje
    O3β (ATP) → O6 (glucosa) queda en el plano de la pantalla, de izquierda a derecha, como en el perfil
    de energía. Así los enlaces que se rompen y se forman se ven de lado y no de frente.
    """
    return _orientacion(_leer_pdb(Path(datos.ruta("md/minimizado_soluto.pdb")).read_text()), ("GLC", "O6"), ("ATP", "O3B"))


def _orientacion(todos, destino, origen_fosfato, ligandos=("GLC", "ATP", "MG", "BGC", "ANP"), de_frente=False):
    """Rotación (R, origen) que pone el sitio activo de frente a la cámara (+z).

    Con ``de_frente`` (vista general) se prueban ~600 direcciones y se elige la que tiene menos átomos de
    proteína en el «tubo» (radio 7 Å) entre el sitio activo y la cámara: la vista más despejada. Sin él (vistas
    del sitio activo) se mira perpendicular al salto del fosfato, para ver de lado los enlaces que cambian.
    El eje O3β → O6 queda horizontal en la pantalla (de izquierda a derecha, como en el perfil de energía).
    """
    xyz = np.array([[a["x"], a["y"], a["z"]] for a in todos if a["r"] not in ligandos and a["e"] != "H"])
    lig = np.array([[a["x"], a["y"], a["z"]] for a in todos if a["r"] in ligandos])
    pos = {(a["r"], a["n"]): np.array([a["x"], a["y"], a["z"]]) for a in todos if a["r"] in ligandos}
    centro = lig.mean(0)
    u = pos[destino] - pos[origen_fosfato]
    u /= np.linalg.norm(u)
    if not de_frente:
        # vistas del sitio activo: desde la boca de la hendidura, pero perpendicular al salto del fosfato (el plano
        # de corte del visor quita la proteína que queda delante)
        prot = np.array([[a["x"], a["y"], a["z"]] for a in todos if a["n"] == "CA"])
        fuera = centro - prot.mean(0)
        vista = fuera - np.dot(fuera, u) * u
        vista /= np.linalg.norm(vista)
        return np.vstack([u, np.cross(vista, u), vista]), centro
    n = 600
    k = np.arange(n) + 0.5
    phi, theta = np.arccos(1 - 2 * k / n), np.pi * (1 + 5 ** 0.5) * k
    dirs = np.c_[np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]
    rel = xyz - centro
    t = rel @ dirs.T                                               # avance hacia la cámara
    perp2 = (rel ** 2).sum(1)[:, None] - t ** 2
    tapan = ((t > 2.0) & (perp2 < 49.0)).sum(0)
    vista = dirs[int(np.argmin(tapan))]
    x = u - np.dot(u, vista) * vista
    x /= np.linalg.norm(x)
    R = np.vstack([x, np.cross(vista, x), vista])                 # filas: ejes x, y, z de la pantalla
    return R, centro


def _girar(xyz, R, origen):
    return (np.asarray(xyz, float) - origen) @ R.T


def _pdb_girado(lineas, R, origen):
    out = []
    for l in lineas:
        x = _girar([[float(l[30:38]), float(l[38:46]), float(l[46:54])]], R, origen)[0]
        out.append(f"{l[:30]}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}{l[54:]}")
    return out


def _orden_pauling(r, r0, b=0.6):
    return float(np.exp((r0 - r) / b))


def _html(config, alto=560):
    """Envuelve la configuración en el reproductor HTML/JS. Devuelve IPython.display.HTML."""
    from IPython.display import HTML

    uid = f"enz3d{next(_contador)}"
    cfg = json.dumps(config, separators=(",", ":"), ensure_ascii=False)
    return HTML(_PLANTILLA.replace("__UID__", uid).replace("__ALTO__", str(alto)).replace("__URL__", URL_3DMOL)
                .replace("__CFG__", cfg.replace("</", "<\\/")))


# ---------------------------------------------------------------------------------------------- datos de la región QM
def _topologia():
    return datos.json_("qmmm/topologia_qm.json")


def _atomos_qm(topo, coords):
    """Átomos de la región QM para ``model.addAtoms`` con enlaces fijos explícitos (sin H de enlace)."""
    visibles = [k for k, a in enumerate(topo["atomos"]) if a["nombre"] != "HL"]
    nuevo = {k: i for i, k in enumerate(visibles)}
    vecinos = {k: [] for k in visibles}
    for i, j in topo["enlaces"]:
        if i in nuevo and j in nuevo:
            vecinos[i].append(nuevo[j])
            vecinos[j].append(nuevo[i])
    atomos = []
    for k in visibles:
        a = topo["atomos"][k]
        cristal = a["resseq"] + DESPLAZAMIENTO if a["residuo"] not in ("GLC", "ATP", "MG") else a["resseq"]
        x, y, z = (float(v) for v in coords[k])
        atomos.append(dict(elem=a["elemento"].capitalize(), x=x, y=y, z=z, bonds=vecinos[k], bondOrder=[1] * len(vecinos[k]),
                           atom=a["nombre"], resn=a["residuo"], resi=cristal, hetflag=True, serial=len(atomos),
                           rol=_rol(a["residuo"], cristal)))
    return atomos, visibles, nuevo


def _enlaces_dinamicos(topo, frames, nuevo):
    """Por fotograma: [[i, j, orden, color], ...] de los enlaces que cambian y de la coordinación del Mg."""
    salida = []
    colores = {"se rompe": ROMPE, "se forma": FORMA}
    mg = topo["mg"]
    for xyz in frames:
        fila = []
        for d in topo["dinamicos"]:
            r = float(np.linalg.norm(xyz[d["i"]] - xyz[d["j"]]))
            n = _orden_pauling(r, d["r0"])
            if n > 0.08:
                fila.append([nuevo[d["i"]], nuevo[d["j"]], round(min(n, 1.0), 3), colores[d["papel"]]])
        for o in topo["coordinacion_mg"]:
            if o in nuevo and np.linalg.norm(xyz[o] - xyz[mg]) < 2.6:
                fila.append([nuevo[mg], nuevo[o], -1, COORD])        # −1 = coordinación (no covalente)
        salida.append(fila)
    return salida


def _etiquetas_qm(topo, nuevo):
    """[índice, texto, desplazamiento x, desplazamiento y en píxeles] de los átomos protagonistas."""
    d = topo["destacados"]
    return [[nuevo[d["PG"]], "Pγ", 0, -26], [nuevo[d["O6"]], "O6 (glucosa)", 34, 20], [nuevo[d["O3B"]], "O3β (ATP)", -34, 20],
            [nuevo[d["OD"]], "Asp205", 30, -18], [nuevo[d["MG"]], "Mg²⁺", 0, 28]]


def _entorno_proteina(radio=5.0, excluir_qm=True):
    """Proteína (caricatura completa) + residuos cercanos a la región QM (varillas), del cristal minimizado."""
    part = datos.json_("qmmm/particion.json")
    qm_glob = {a["global_index"] for a in part["qm_atoms"] if a.get("global_index") is not None}
    texto = Path(datos.ruta("md/minimizado_soluto.pdb")).read_text()
    todos = _leer_pdb(texto)
    xyz = np.array([[a["x"], a["y"], a["z"]] for a in todos])
    qm_xyz = xyz[sorted(i for i in qm_glob if i < len(todos))]   # el agua del Mg no está en el soluto
    cerca = np.min(np.linalg.norm(xyz[:, None, :] - qm_xyz[None, :, :], axis=2), axis=1) < radio
    res_cerca = {(todos[i]["s"]) for i in np.where(cerca)[0] if todos[i]["r"] not in ("WAT", "HOH", "Na+", "K+")}
    lineas = []
    for i, (a, l) in enumerate(zip(todos, [l for l in texto.splitlines() if l.startswith(("ATOM", "HETATM"))])):
        if a["r"] in ("GLC", "ATP", "MG", "WAT", "HOH", "Na+", "K+"):
            continue
        if excluir_qm and i in qm_glob and a["n"] not in ("N", "CA", "C", "O"):
            continue
        es_bb = a["n"] in ("N", "CA", "C", "O")
        if not es_bb and (a["s"] not in res_cerca or a["e"] == "H"):
            continue
        lineas.append(l[:22] + f"{a['s'] + DESPLAZAMIENTO:4d}" + l[26:])
    # enlaces frontera: el átomo MM real que el H de enlace sustituye, unido al átomo QM
    frontera = [[l["mm"], l["qm"]] for l in part["links"]]
    return "\n".join(lineas), frontera, todos


# ---------------------------------------------------------------------------------------------- visores públicos
def pelicula_reaccion(camino="enzima", proteina=True, titulo=None):
    """La transferencia del fosfato, fotograma a fotograma, con el perfil de energía sincronizado.

    ``camino``: "enzima" (camino de mínima energía descendiendo desde el TS, QM/MM) o "agua" (el mismo sitio
    activo en disolvente implícito COSMO).
    """
    topo = _topologia()
    nombre_xyz, nombre_csv, col_e = {
        "enzima": ("qmmm/camino_descenso.xyz", "qmmm/camino_descenso.csv", "energia_rel_kcal"),
        "agua": ("qmmm/agua_camino_descenso.xyz", "qmmm/agua_camino_descenso.csv", "energia_kcal"),
    }[camino]
    _, frames, _ = datos.leer_xyz_multiple(nombre_xyz)
    tabla = datos.csv(nombre_csv)
    energia = tabla[col_e].to_numpy(dtype=float)
    xi = tabla["xi"].to_numpy(dtype=float)
    # quitar fotogramas repetidos (el descenso guarda dos veces el punto de partida)
    keep = [0] + [k for k in range(1, len(frames)) if np.abs(frames[k] - frames[k - 1]).max() > 1e-4]
    frames, energia, xi = [frames[k] for k in keep], energia[keep], xi[keep]
    energia = energia - energia[0]
    d = topo["destacados"]
    lecturas, fases = [], []
    k_ts = int(np.argmax(energia))
    for k, f in enumerate(frames):
        d1 = np.linalg.norm(f[d["PG"]] - f[d["O3B"]])
        d2 = np.linalg.norm(f[d["PG"]] - f[d["O6"]])
        dh = np.linalg.norm(f[d["H6"]] - f[d["OD"]])
        lecturas.append([f"{energia[k]:.1f} kcal/mol", f"{xi[k]:+.2f} Å", f"{d1:.2f} Å", f"{d2:.2f} Å", f"{dh:.2f} Å"])
        if abs(k - k_ts) <= 1:
            fases.append("Estado de transición: el fósforo está a medio camino entre los dos oxígenos")
        elif k < k_ts:
            fases.append("Reactivo → cima: el fosfato se separa del ATP y se acerca a la glucosa")
        elif dh > 1.3:
            fases.append("Después de la cima: el fosfato ya está en la glucosa; el protón de O6 aún no se ha movido")
        else:
            fases.append("Hacia el producto: el fosfato está en la glucosa y Asp205 empieza a tomar el protón de O6")
    R_, o_ = _orientacion_sitio()
    frames = [_girar(f, R_, o_) for f in frames]
    atomos, visibles, nuevo = _atomos_qm(topo, frames[0])
    config = dict(
        tipo="reaccion",
        titulo=titulo or ("La reacción dentro de la enzima" if camino == "enzima" else "La reacción en agua (sin el resto de la enzima)"),
        subtitulo="Arrastra para girar · rueda para acercar · ▶ para reproducir",
        modelos=[dict(atomos=atomos, frames=_b64_int16(np.array([f[visibles] for f in frames]).ravel()),
                      n=len(visibles), estilo="bola")],
        dinamicos=_enlaces_dinamicos(topo, frames, nuevo),
        etiquetas=_etiquetas_qm(topo, nuevo),
        grafica=dict(x=xi.round(3).tolist(), y=energia.round(3).tolist(), xlabel="ξ = d(Pγ–O3β) − d(Pγ–O6)  (Å)",
                     ylabel="energía (kcal/mol)", titulo="Perfil de energía", ts=k_ts,
                     puntos=[[0, "R", AZUL], [k_ts, "TS", NARANJA], [len(xi) - 1, "P", AGUA]]),
        lecturas=dict(nombres=["Energía", "ξ", "Pγ–O3β", "Pγ–O6", "H···O Asp205"], valores=lecturas),
        fases=fases,
        leyenda=[[ROMPE, "enlace que se rompe"], [FORMA, "enlace que se forma"], [COORD, "Mg²⁺ coordinado (no covalente)"]],
    )
    if camino == "enzima" and proteina:
        pdb_env, frontera, todos = _entorno_proteina()
        part = datos.json_("qmmm/particion.json")
        glob2k = {a["global_index"]: a["k"] for a in part["qm_atoms"] if a.get("global_index") is not None}
        config["entorno"] = "\n".join(_pdb_girado(pdb_env.splitlines(), R_, o_))
        # plano de corte justo delante del átomo QM más cercano a la cámara: nada de la proteína tapa la química
        zmax = max(float(np.max(f[visibles][:, 2])) for f in frames)
        config["corte"] = [zmax + 1.2, 16.0]
        config["entorno"] = "\n".join(l for l in config["entorno"].splitlines() if float(l[46:54]) < zmax + 1.2
                                      or l[12:16].strip() in ("N", "CA", "C", "O"))
        config["frontera"] = [[_girar([[todos[mm]["x"], todos[mm]["y"], todos[mm]["z"]]], R_, o_)[0].round(3).tolist(),
                               nuevo[glob2k[qm]]] for mm, qm in frontera if glob2k.get(qm) in nuevo]
    return _html(config, alto=600)


def modo_imaginario(n_cuadros=30, amplitud=0.35):
    """La vibración del estado de transición a lo largo de su modo imaginario (−149 cm⁻¹)."""
    topo = _topologia()
    vib = np.load(datos.ruta("qmmm/frecuencias_ts.npz"))
    ts, modo = vib["ts_xyz"], vib["modo_imaginario"]
    modo = modo / np.abs(modo).max()
    fases_t = np.sin(np.linspace(0, 2 * np.pi, n_cuadros, endpoint=False))
    d = topo["destacados"]
    crudos = [ts + amplitud * s * modo for s in fases_t]
    xi = np.array([np.linalg.norm(f[d["PG"]] - f[d["O3B"]]) - np.linalg.norm(f[d["PG"]] - f[d["O6"]]) for f in crudos])
    R_, o_ = _orientacion_sitio()
    frames = [_girar(f, R_, o_) for f in crudos]
    atomos, visibles, nuevo = _atomos_qm(topo, frames[0])
    freq = float(vib["freq_cm_signed"].min())
    rel = xi - xi.mean()
    config = dict(
        tipo="reaccion", titulo=f"El modo imaginario del estado de transición ({freq:.0f} cm⁻¹)",
        subtitulo="No es una vibración normal: hacia un lado cae al reactivo, hacia el otro al producto",
        modelos=[dict(atomos=atomos, frames=_b64_int16(np.array([f[visibles] for f in frames]).ravel()), n=len(visibles), estilo="bola")],
        dinamicos=_enlaces_dinamicos(topo, frames, nuevo), etiquetas=_etiquetas_qm(topo, nuevo),
        grafica=dict(x=list(range(n_cuadros)), y=rel.round(3).tolist(), xlabel="un ciclo de la vibración",
                     ylabel="desplazamiento en ξ (Å)", titulo="¿Hacia dónde empuja el modo?", ts=None, puntos=[],
                     bandas_y=[[0, 10, "hacia el producto", AGUA], [-10, 0, "hacia el reactivo", AZUL]]),
        lecturas=dict(nombres=["ξ relativo"], valores=[[f"{v:+.2f} Å"] for v in rel]),
        fases=["El fósforo se acerca al O6 de la glucosa: hacia el producto" if v > 0 else
               "El fósforo vuelve hacia el O3β del ATP: hacia el reactivo" for v in rel],
        leyenda=[[ROMPE, "enlace que se rompe"], [FORMA, "enlace que se forma"], [COORD, "Mg²⁺ coordinado (no covalente)"]],
        auto=True,
    )
    return _html(config, alto=560)


def pelicula_md(umbral=3.5, cristal=2.68):
    """La dinámica molecular: la enzima se mueve y los reactivos «vibran» sin dejar de apuntarse."""
    with gzip.open(datos.PRECALCULADO / "md" / "pelicula_md.json.gz", "rt", encoding="utf-8") as f:
        peli = json.load(f)
    vecinos = [[] for _ in peli["atomos"]]
    for i, j in peli["enlaces"]:
        vecinos[i].append(j)
        vecinos[j].append(i)
    R_, o_ = _orientacion_sitio()
    frames = np.array(peli["frames"], dtype=float)
    frames = _girar(frames.reshape(len(frames), -1, 3), R_, o_).reshape(len(frames), -1)
    atomos = []
    for k, a in enumerate(peli["atomos"]):
        atomos.append(dict(elem=a["e"].capitalize(), x=frames[0][3 * k], y=frames[0][3 * k + 1], z=frames[0][3 * k + 2],
                           bonds=vecinos[k], bondOrder=[1] * len(vecinos[k]), atom=a["n"], resn=a["r"], resi=a["s"],
                           hetflag=bool(a["h"]), serial=k, sitio=a["sitio"], rol=_rol(a["r"], a["s"] if not a["h"] else 0)))
    t, dd = peli["tiempo_ps"], peli["d_PG_O6"]
    # núcleo químico (glucosa, fosfatos, Mg): el plano de corte va justo delante, la adenina puede quedar cortada
    nucleo = [k for k, a in enumerate(peli["atomos"]) if a["r"] in ("GLC", "MG") or (a["r"] == "ATP" and a["n"][:1] in ("P", "O"))]
    zmax = float(np.percentile(frames.reshape(len(frames), -1, 3)[:, nucleo, 2].max(axis=1), 90))
    config = dict(
        tipo="md", titulo="La película de la enzima: 1 ns de dinámica molecular",
        subtitulo="Arrastra para girar · rueda para acercar · ▶ para reproducir",
        modelos=[dict(atomos=atomos, frames=_b64_int16(frames.ravel()), n=len(atomos), estilo="md")],
        dinamicos=[[] for _ in t], etiquetas=[],
        grafica=dict(x=t, y=dd, xlabel="tiempo (ps)", ylabel="d(Pγ–O6) (Å)", titulo="¿Se apuntan los reactivos?",
                     ts=None, puntos=[], bandas_y=[[0, umbral, f"ataque cercano (< {umbral} Å)", AGUA]],
                     linea_y=[cristal, f"cristal {cristal} Å"]),
        lecturas=dict(nombres=["tiempo", "d(Pγ–O6)"], valores=[[f"{a:.0f} ps", f"{b:.2f} Å"] for a, b in zip(t, dd)]),
        fases=["Lista para reaccionar: el fósforo del ATP apunta al O6 de la glucosa" if b < umbral else
               "Un momento de separación: los reactivos se alejan un poco" for b in dd],
        leyenda=[[CARBONO["glc"], "glucosa"], [CARBONO["atp"], "ATP"], [ELEMENTOS["Mg"], "Mg²⁺"],
                 ["#7f9cc9", "dominio grande"], ["#e2a2c3", "dominio pequeño"]],
        corte=[zmax + 1.5, 40.0],
    )
    return _html(config, alto=600)


def complejo_cristal(pdb_texto=None):
    """La glucoquinasa del cristal 3FGU: dos dominios que abrazan a la glucosa y al ATP (AMP‑PNP)."""
    if pdb_texto is None:
        pdb_texto = (datos.RAIZ / "data" / "raw" / "3FGU.pdb").read_text()
    lineas = [l for l in pdb_texto.splitlines() if l.startswith(("ATOM", "HETATM")) and l[21] == "A"
              and l[17:20].strip() not in ("HOH",) and l[16] in " A"]
    # misma puesta en escena que los demás visores: la hendidura mira a la cámara
    R_, o_ = _orientacion(_leer_pdb("\n".join(lineas)), ("BGC", "O6"), ("ANP", "N3B"), de_frente=True)
    lineas = _pdb_girado(lineas, R_, o_)
    config = dict(tipo="complejo", titulo="La glucoquinasa humana (cristal 3FGU)",
                  subtitulo="Arrastra para girar · rueda para acercar · botones para cambiar la vista",
                  pdb="\n".join(lineas),
                  leyenda=[["#7f9cc9", "dominio grande"], ["#e2a2c3", "dominio pequeño"], [CARBONO["glc"], "glucosa"],
                           [CARBONO["atp"], "AMP‑PNP (análogo del ATP)"], [ELEMENTOS["Mg"], "Mg²⁺"]],
                  residuos=[[205, "Asp205 · base"], [169, "Lys169"], [228, "Thr228"], [151, "Ser151"]])
    return _html(config, alto=560)


def region_qm():
    """La partición QM/MM: 78 átomos cuánticos dentro del campo eléctrico de la enzima."""
    topo = _topologia()
    xyz = np.array([[a["x"], a["y"], a["z"]] for a in _leer_pdb(Path(datos.ruta("qmmm/region_qm_inicial.pdb")).read_text())])
    R_, o_ = _orientacion_sitio()
    xyz = _girar(xyz, R_, o_)
    atomos, visibles, nuevo = _atomos_qm(topo, xyz)
    env = _leer_pdb(Path(datos.ruta("qmmm/entorno_mm_8A.pdb")).read_text())
    cargas = [[*_girar([[a["x"], a["y"], a["z"]]], R_, o_)[0].round(2).tolist(), round(a["q"], 2)] for a in env if abs(a["q"]) > 0.3]
    config = dict(tipo="region", titulo="QM/MM: la región cuántica dentro de la enzima",
                  subtitulo="Bolas y varillas = región QM (PM7) · esferas rojas/azules = cargas de la enzima (−/+)",
                  modelos=[dict(atomos=atomos, frames=None, n=len(visibles), estilo="bola")],
                  cargas=cargas, etiquetas=_etiquetas_qm(topo, nuevo), dinamicos=[_enlaces_dinamicos(topo, [xyz], nuevo)[0]],
                  leyenda=[[ROJO, "carga negativa de la enzima"], [AZUL, "carga positiva de la enzima"],
                           [CARBONO["glc"], "glucosa"], [CARBONO["atp"], "trifosfato"]])
    return _html(config, alto=540)


# ---------------------------------------------------------------------------------------------- plantilla HTML/JS
_PLANTILLA = r"""
<div id="__UID__" class="enz3d" style="font-family:Figtree,'Avenir Next','Segoe UI',Roboto,Helvetica,Arial,sans-serif;
  background:radial-gradient(120% 140% at 20% 0%,#1b2436 0%,#0c111b 55%,#080b12 100%);border-radius:18px;padding:16px 16px 12px;
  color:#e8ecf3;max-width:1100px;box-shadow:0 10px 30px rgba(8,11,18,.25)">
 <div style="display:flex;justify-content:space-between;align-items:baseline;gap:12px;flex-wrap:wrap">
  <div><div class="t" style="font-size:19px;font-weight:650;letter-spacing:.01em"></div>
       <div class="st" style="font-size:13px;color:#9aa6ba;margin-top:2px"></div></div>
  <div class="bt" style="display:flex;gap:6px;flex-wrap:wrap"></div>
 </div>
 <div style="display:flex;gap:14px;margin-top:12px;flex-wrap:wrap">
  <div class="v" style="position:relative;flex:1.35 1 420px;height:__ALTO__px;border-radius:14px;overflow:hidden;
       background:radial-gradient(90% 90% at 50% 40%,#1a2334 0%,#0b0f18 100%)"></div>
  <div class="panel" style="flex:1 1 300px;display:flex;flex-direction:column;gap:10px;min-width:0"></div>
 </div>
 <div class="ctl" style="display:flex;align-items:center;gap:12px;margin-top:12px"></div>
 <div class="ley" style="display:flex;gap:16px;flex-wrap:wrap;margin-top:10px;font-size:12.5px;color:#aab4c5"></div>
 <div class="err" style="display:none;margin-top:8px;font-size:13px;color:#ffb4a8"></div>
</div>
<script>
(function(){
const root=document.getElementById("__UID__"); const cfg=__CFG__;
const URL="__URL__";
function load(cb){
  if(window.$3Dmol){cb(window.$3Dmol);return;}
  if(typeof window.define==="function"&&window.define.amd&&window.require){
    window.require.config({paths:{"3Dmol":URL.replace(/\.js$/,"")}});
    window.require(["3Dmol"],m=>cb(m||window.$3Dmol));return;}
  const s=document.createElement("script");s.src=URL;s.onload=()=>cb(window.$3Dmol);
  s.onerror=()=>{const e=root.querySelector(".err");e.style.display="block";
    e.textContent="No se pudo cargar 3Dmol.js (¿sin conexión a internet?). Vuelve a ejecutar la celda.";};
  document.head.appendChild(s);}
const EL={O:"#ff5a5f",N:"#6b9cff",P:"#ff9f1c",S:"#ffd166",H:"#f4f4f6",Mg:"#3ddc84",MG:"#3ddc84",K:"#b69cff",Na:"#b69cff"};
const CARB={glc:"#ffb38a",atp:"#ffd27a",prot:"#c9d4e6",asp:"#9ec5ff",lys:"#c7b8ff",thr:"#a8e6cf",otro:"#c9d4e6"};
function esquema(rol){const m=Object.assign({},EL);m.C=CARB[rol]||CARB.prot;return {prop:"elem",map:m};}
function el(tag,st,txt){const e=document.createElement(tag);if(st)e.style.cssText=st;if(txt!==undefined)e.textContent=txt;return e;}
function boton(txt,on){const b=el("button","background:rgba(255,255,255,.08);color:#e8ecf3;border:1px solid rgba(255,255,255,.14);"+
  "border-radius:999px;padding:6px 12px;font-weight:600;font-size:12.5px;font-family:inherit;cursor:pointer",txt);b.onclick=on;return b;}
function activo(b,on){b.style.background=on?"rgba(42,120,214,.55)":"rgba(255,255,255,.08)";b.style.borderColor=on?"rgba(120,170,255,.8)":"rgba(255,255,255,.14)";}
function tono(hex,n){const k=.45+.55*Math.min(1,Math.max(0,n));const c=parseInt(hex.slice(1),16);
  const r=(c>>16)&255,g=(c>>8)&255,b=c&255,f=v=>Math.round(26+(v-26)*k).toString(16).padStart(2,"0");return "#"+f(r)+f(g)+f(b);}
function decode(b64){const bin=atob(b64);const n=bin.length/2;const out=new Array(n);
  for(let i=0;i<n;i++){let v=bin.charCodeAt(2*i)|(bin.charCodeAt(2*i+1)<<8);if(v>32767)v-=65536;out[i]=v/100;}return out;}
root.querySelector(".t").textContent=cfg.titulo||"";root.querySelector(".st").textContent=cfg.subtitulo||"";
const ley=root.querySelector(".ley");
(cfg.leyenda||[]).forEach(([c,t])=>{const s=el("span","display:inline-flex;align-items:center;gap:6px");
  s.appendChild(el("span","width:10px;height:10px;border-radius:50%;background:"+c+";display:inline-block"));s.appendChild(el("span","",t));ley.appendChild(s);});

// ------------------------------------------------------------------ gráfica sincronizada (SVG)
function grafica(g){
  const W=420,H=280,L=52,R=14,T=34,B=44; const NS="http://www.w3.org/2000/svg";
  const svg=document.createElementNS(NS,"svg");svg.setAttribute("viewBox",`0 0 ${W} ${H}`);svg.style.cssText="width:100%;height:auto;display:block";
  const xs=g.x,ys=g.y; let x0=Math.min(...xs),x1=Math.max(...xs),y0=Math.min(...ys),y1=Math.max(...ys);
  if(g.linea_y){y0=Math.min(y0,g.linea_y[0]);}
  const py=(y1-y0)*0.12||1; y0-=py; y1+=py*1.6;
  const X=v=>L+(v-x0)/(x1-x0||1)*(W-L-R), Y=v=>T+(1-(v-y0)/(y1-y0))*(H-T-B);
  const mk=(t,a,txt)=>{const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);if(txt!==undefined)e.textContent=txt;svg.appendChild(e);return e;};
  const defs=mk("defs",{}); defs.innerHTML='<linearGradient id="gr__UID__" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#eb6834" stop-opacity=".35"/><stop offset="1" stop-color="#2a78d6" stop-opacity=".02"/></linearGradient>'+
    '<filter id="gl__UID__" x="-50%" y="-50%" width="200%" height="200%"><feGaussianBlur stdDeviation="3" result="b"/><feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge></filter>';
  (g.bandas_y||[]).forEach(b=>{const a=Math.max(b[0],y0),z=Math.min(b[1],y1);if(z<=a)return;
    mk("rect",{x:L,y:Y(z),width:W-L-R,height:Y(a)-Y(z),fill:b[3],opacity:.14});
    mk("text",{x:W-R-6,y:Y(a)-7,"text-anchor":"end",fill:"#b9c3d3","font-size":11},b[2]);});
  for(let k=0;k<=4;k++){const v=y0+(y1-y0)*k/4;mk("line",{x1:L,x2:W-R,y1:Y(v),y2:Y(v),stroke:"rgba(255,255,255,.07)"});
    mk("text",{x:L-8,y:Y(v)+4,"text-anchor":"end",fill:"#8793a8","font-size":11},v.toFixed(Math.abs(y1-y0)<6?1:0));}
  mk("text",{x:L,y:18,fill:"#e8ecf3","font-size":13.5,"font-weight":650},g.titulo||"");
  mk("text",{x:(L+W-R)/2,y:H-8,"text-anchor":"middle",fill:"#9aa6ba","font-size":11.5},g.xlabel||"");
  mk("text",{x:14,y:(T+H-B)/2,"text-anchor":"middle",fill:"#9aa6ba","font-size":11.5,transform:`rotate(-90 14 ${(T+H-B)/2})`},g.ylabel||"");
  if(g.linea_y){mk("line",{x1:L,x2:W-R,y1:Y(g.linea_y[0]),y2:Y(g.linea_y[0]),stroke:"#e8ecf3","stroke-dasharray":"4 4",opacity:.55});
    mk("text",{x:L+6,y:Y(g.linea_y[0])+14,fill:"#c9d2df","font-size":11},g.linea_y[1]);}
  let d="";xs.forEach((v,i)=>{d+=(i?"L":"M")+X(v).toFixed(1)+" "+Y(ys[i]).toFixed(1);});
  if(g.ts!==null&&g.ts!==undefined){mk("path",{d:d+`L${X(xs[xs.length-1])} ${Y(y0)}L${X(xs[0])} ${Y(y0)}Z`,fill:`url(#gr__UID__)`});}
  mk("path",{d:d,fill:"none",stroke:"#8fb6ff","stroke-width":2.2,"stroke-linejoin":"round","stroke-linecap":"round"});
  (g.puntos||[]).forEach(([i,t,c])=>{mk("circle",{cx:X(xs[i]),cy:Y(ys[i]),r:5,fill:c,stroke:"#0c111b","stroke-width":1.5});
    mk("text",{x:X(xs[i]),y:Y(ys[i])-11,"text-anchor":"middle",fill:"#e8ecf3","font-size":12,"font-weight":650},t);});
  if(g.ts!==null&&g.ts!==undefined){const i=g.ts;const yb=ys[0];
    mk("line",{x1:X(xs[i]),x2:X(xs[i]),y1:Y(yb),y2:Y(ys[i])+7,stroke:"#eb6834","stroke-width":1.4,"stroke-dasharray":"3 3"});
    mk("line",{x1:X(xs[0]),x2:X(xs[i]),y1:Y(yb),y2:Y(yb),stroke:"rgba(255,255,255,.25)","stroke-dasharray":"2 4"});
    mk("text",{x:X(xs[i])-8,y:(Y(yb)+Y(ys[i]))/2,"text-anchor":"end",fill:"#ffc4a3","font-size":12,"font-weight":650},"ΔE‡ = "+(ys[i]-yb).toFixed(1)+" kcal/mol");}
  const vl=mk("line",{x1:0,x2:0,y1:T,y2:H-B,stroke:"rgba(255,255,255,.35)","stroke-dasharray":"2 4"});
  const halo=mk("circle",{r:11,fill:"#ffd166",opacity:.25,filter:`url(#gl__UID__)`});
  const dot=mk("circle",{r:6,fill:"#ffd166",stroke:"#0c111b","stroke-width":1.5});
  return {svg,set(i){const x=X(xs[i]),y=Y(ys[i]);vl.setAttribute("x1",x);vl.setAttribute("x2",x);
    halo.setAttribute("cx",x);halo.setAttribute("cy",y);dot.setAttribute("cx",x);dot.setAttribute("cy",y);}};
}

load(function($3Dmol){
 try{
  const vdiv=root.querySelector(".v");
  const viewer=$3Dmol.createViewer(vdiv,{backgroundColor:"black",backgroundAlpha:0,antialias:true,cartoonQuality:10});
  viewer.enableFog(true);
  const panel=root.querySelector(".panel"), ctl=root.querySelector(".ctl"), bt=root.querySelector(".bt");
  let modelos=[], nframes=1, girando=false, mostrarEtiquetas=true, mostrarProt=true, env=null;
  root._viewer=viewer;

  if(cfg.tipo==="complejo"){
    const m=viewer.addModel(cfg.pdb,"pdb");
    m.setStyle({hetflag:false},{cartoon:{color:"#7f9cc9"}}); m.setStyle({resi:["65-204","441-465"],hetflag:false},{cartoon:{color:"#e2a2c3"}});
    m.setStyle({resn:"BGC"},{stick:{radius:.28,colorscheme:esquema("glc")},sphere:{scale:.38,colorscheme:esquema("glc")}});
    m.setStyle({resn:"ANP"},{stick:{radius:.28,colorscheme:esquema("atp")},sphere:{scale:.38,colorscheme:esquema("atp")}});
    m.setStyle({resn:"MG"},{sphere:{scale:.7,color:"#3ddc84"}}); m.setStyle({resn:"K"},{sphere:{scale:.6,color:"#b69cff"}});
    const sel=cfg.residuos.map(r=>r[0]);
    m.setStyle({resi:sel,hetflag:false,not:{atom:["N","C","O"]}},{stick:{radius:.18,colorscheme:esquema("asp")}},true);
    let sup=null;
    const vista=(modo)=>{viewer.removeAllLabels();
      if(modo==="sitio"){viewer.zoomTo({resn:["BGC","ANP","MG"]},800);
        cfg.residuos.forEach(([r,t])=>{viewer.addLabel(t,{fontSize:12,fontColor:"#0c111b",backgroundColor:"#e8ecf3",backgroundOpacity:.85,
          borderRadius:6,inFront:true},{resi:r,atom:"CB"});});
        viewer.addLabel("glucosa",{fontSize:12,fontColor:"#0c111b",backgroundColor:"#ffb38a",backgroundOpacity:.9,inFront:true},{resn:"BGC",atom:"C4"});
        viewer.addLabel("ATP",{fontSize:12,fontColor:"#0c111b",backgroundColor:"#ffd27a",backgroundOpacity:.9,inFront:true},{resn:"ANP",atom:"PB"});}
      else{viewer.zoomTo();viewer.zoom(1.3);}
      viewer.render();};
    const bG=boton("Toda la enzima",()=>{vista("todo");activo(bG,true);activo(bS,false);});
    const bS=boton("Sitio activo",()=>{vista("sitio");activo(bS,true);activo(bG,false);});
    const bSup=boton("Superficie",()=>{if(sup===null){sup=viewer.addSurface($3Dmol.SurfaceType.MS,{opacity:.6,color:"#9fb3d6"},{hetflag:false});activo(bSup,true);}
      else{viewer.removeSurface(sup);sup=null;activo(bSup,false);}});
    const bR=boton("Girar",()=>{girando=!girando;viewer.spin(girando?"y":false,.4);activo(bR,girando);});
    [bG,bS,bSup,bR].forEach(b=>bt.appendChild(b)); activo(bG,true);
    const info=el("div","background:rgba(255,255,255,.05);border:1px solid rgba(255,255,255,.1);border-radius:14px;padding:14px 16px;font-size:13.5px;line-height:1.55;color:#c9d2df");
    info.innerHTML="<div style='font-weight:650;color:#e8ecf3;font-size:15px;margin-bottom:6px'>Qué estás viendo</div>"+
      "Las cintas son la cadena de la proteína: <b style='color:#9fb6de'>dominio grande</b> y <b style='color:#efb3d1'>dominio pequeño</b>. "+
      "En la hendidura entre ambos están la <b style='color:#ffb38a'>glucosa</b> y el <b style='color:#ffd27a'>ATP</b> (aquí AMP‑PNP), con el <b style='color:#3ddc84'>Mg²⁺</b>.<br><br>"+
      "Pulsa <b>Sitio activo</b> para acercarte: verás Asp205, la base que tomará el protón de la glucosa. "+
      "<b>Superficie</b> muestra la forma de la proteína: la glucosa queda casi enterrada.";
    panel.appendChild(info);
    // un giro lento de bienvenida que se detiene solo a los 10 s (o al pulsar «Girar»)
    vista("todo"); viewer.spin("y",.25); girando=true; activo(bR,true);
    setTimeout(()=>{if(girando){girando=false;viewer.spin(false);activo(bR,false);}},10000);
    return;
  }

  // ---------------------------------------------------------------- modelos con fotogramas (reacción, MD, región)
  cfg.modelos.forEach(md=>{
    const m=viewer.addModel(); m.addAtoms(md.atomos);
    if(md.frames){const flat=decode(md.frames); nframes=flat.length/(3*md.n);
      const arr=[];for(let f=0;f<nframes;f++){const fr=[];for(let i=0;i<md.n;i++){const o=3*(f*md.n+i);fr.push([flat[o],flat[o+1],flat[o+2]]);}arr.push(fr);}
      m.setCoordinates(arr,"array");}
    if(md.estilo==="bola"){
      ["glc","atp","asp","lys","thr","prot"].forEach(r=>{m.setStyle({rol:r},{stick:{radius:.16,colorscheme:esquema(r)},sphere:{scale:.27,colorscheme:esquema(r)}});});
      m.setStyle({elem:"Mg"},{sphere:{scale:.45,color:"#3ddc84"}});
      m.setStyle({elem:"H"},{stick:{radius:.12,colorscheme:esquema("prot")},sphere:{scale:.18,color:"#f4f4f6"}});
    }else{
      m.setStyle({hetflag:false},{cartoon:{color:"#7f9cc9"}});
      m.setStyle({hetflag:false,resi:["65-204","441-465"]},{cartoon:{color:"#e2a2c3"}});
      m.setStyle({hetflag:false,sitio:1},{stick:{radius:.12,colorscheme:esquema("prot")}},true);
      m.setStyle({resn:"GLC"},{stick:{radius:.2,colorscheme:esquema("glc")},sphere:{scale:.26,colorscheme:esquema("glc")}});
      m.setStyle({resn:"ATP"},{stick:{radius:.2,colorscheme:esquema("atp")},sphere:{scale:.26,colorscheme:esquema("atp")}});
      m.setStyle({resn:"MG"},{sphere:{scale:.5,color:"#3ddc84"}});
    }
    modelos.push(m);});
  if(cfg.entorno){env=viewer.addModel(cfg.entorno,"pdb",{keepH:false});
    const envC=Object.assign({},EL);envC.C="#7d8aa3";
    env.setStyle({},{cartoon:{color:"#33456a"},stick:{radius:.06,colorscheme:{prop:"elem",map:envC}}});
    env.setStyle({atom:["N","CA","C","O"]},{cartoon:{color:"#33456a"}});}
  // cargas puntuales de la enzima: pequeñas y apagadas para no confundirlas con átomos
  if(cfg.cargas){cfg.cargas.forEach(([x,y,z,q])=>{viewer.addSphere({center:{x,y,z},radius:.12+Math.min(Math.abs(q),1)*.16,
      color:q<0?"#9c3346":"#3160b8",wireframe:true,linewidth:1});});}

  // cilindros de los enlaces que cambian, fotograma a fotograma
  const main=modelos[0]; let shapes=[], labels=[];
  function pos(f,i){const a=(main.frames&&main.frames.length?main.frames[f]:main.selectedAtoms({}))[i];return {x:a.x,y:a.y,z:a.z};}
  function dibujar(f){
    shapes.forEach(s=>viewer.removeShape(s)); shapes=[]; labels.forEach(l=>viewer.removeLabel(l)); labels=[];
    const din=cfg.dinamicos[Math.min(f,cfg.dinamicos.length-1)]||[];
    din.forEach(([i,j,n,c])=>{const a=pos(f,i),b=pos(f,j);
      if(n<0){shapes.push(viewer.addCylinder({start:a,end:b,radius:.045,color:"#2a8f5c",dashed:true,dashLength:.14,gapLength:.14}));return;}
      const full=n>0.72;
      shapes.push(viewer.addCylinder({start:a,end:b,radius:.06+.12*n,color:tono(c,n),dashed:!full,dashLength:.2,gapLength:.14,
        fromCap:2,toCap:2}));});
    (cfg.frontera||[]).forEach(([p,i])=>{shapes.push(viewer.addCylinder({start:{x:p[0],y:p[1],z:p[2]},end:pos(f,i),radius:.1,color:"#c9d4e6"}));});
    if(mostrarEtiquetas){(cfg.etiquetas||[]).forEach(([i,t,dx,dy])=>{labels.push(viewer.addLabel(t,{position:pos(f,i),fontSize:11,fontColor:"#0c111b",
      backgroundColor:"#e8ecf3",backgroundOpacity:.82,borderRadius:5,inFront:true,alignment:"center",screenOffset:{x:dx||0,y:dy||0}}));});}
  }

  // panel lateral: gráfica + fase + lecturas
  let G=null; if(cfg.grafica){const box=el("div","background:rgba(255,255,255,.04);border:1px solid rgba(255,255,255,.09);border-radius:14px;padding:8px 8px 4px");
    G=grafica(cfg.grafica); box.appendChild(G.svg); panel.appendChild(box);}
  const fase=el("div","background:rgba(255,209,102,.08);border:1px solid rgba(255,209,102,.25);border-radius:12px;padding:10px 12px;font-size:13.5px;line-height:1.45;color:#ffe3a3;min-height:40px");
  const lect=el("div","display:grid;grid-template-columns:repeat(auto-fit,minmax(104px,1fr));gap:8px");
  let celdas=[]; if(cfg.lecturas){cfg.lecturas.nombres.forEach(nm=>{const c=el("div","background:rgba(255,255,255,.05);border-radius:10px;padding:7px 9px");
    c.appendChild(el("div","font-size:11px;letter-spacing:.02em;color:#8793a8;white-space:nowrap",nm));
    const v=el("div","font-size:15px;font-weight:650;color:#f1f4f9;margin-top:2px;white-space:nowrap","");
    c.appendChild(v);lect.appendChild(c);celdas.push(v);});}
  if(cfg.fases)panel.appendChild(fase); if(cfg.lecturas)panel.appendChild(lect);
  if(cfg.tipo==="region"){const info=el("div","background:rgba(255,255,255,.05);border:1px solid rgba(255,255,255,.1);border-radius:14px;padding:14px 16px;font-size:13.5px;line-height:1.55;color:#c9d2df");
    info.innerHTML="<div style='font-weight:650;color:#e8ecf3;font-size:15px;margin-bottom:6px'>Qué estás viendo</div>"+
      "En <b>bolas y varillas</b>, los 78 átomos que se tratan con mecánica cuántica (PM7): glucosa, trifosfato, Mg²⁺ con su agua "+
      "y las cadenas laterales de Asp205, Lys169 y Thr228.<br><br>Las <b style='color:#ff8a8a'>esferas rojas</b> y <b style='color:#8ab4ff'>azules</b> son "+
      "cargas negativas y positivas de la enzima (región MM): sus electrones no se calculan, pero su campo eléctrico sí entra en el cálculo cuántico.";
    panel.appendChild(info);}

  // controles
  let f=0, timer=null, vel=1;
  const play=boton("▶  Reproducir",()=>{if(timer){parar();}else{iniciar();}});
  const slider=el("input","flex:1;accent-color:#ffd166");slider.type="range";slider.min=0;slider.max=Math.max(0,nframes-1);slider.value=0;
  const cont=el("span","font-size:12.5px;color:#9aa6ba;min-width:70px;text-align:right","");
  slider.oninput=()=>{parar();ir(+slider.value);};
  const bV=boton("×1",()=>{vel=vel===1?2:(vel===2?0.5:1);bV.textContent="×"+vel;if(timer){parar();iniciar();}});
  if(nframes>1){ctl.appendChild(play);ctl.appendChild(slider);ctl.appendChild(cont);ctl.appendChild(bV);}
  const bE=boton("Etiquetas",()=>{mostrarEtiquetas=!mostrarEtiquetas;activo(bE,mostrarEtiquetas);dibujar(f);viewer.render();});
  if(cfg.etiquetas&&cfg.etiquetas.length){bt.appendChild(bE);activo(bE,true);}
  if(env){const bP=boton("Proteína",()=>{mostrarProt=!mostrarProt;if(mostrarProt){env.show();}else{env.hide();}activo(bP,mostrarProt);viewer.render();});
    bt.appendChild(bP);activo(bP,true);}
  const bR=boton("Girar",()=>{girando=!girando;viewer.spin(girando?"y":false,.5);activo(bR,girando);}); bt.appendChild(bR);
  const bC=boton("Centrar",()=>{centrar(600);}); bt.appendChild(bC);
  // «corte» cinematográfico: solo se dibuja una rebanada alrededor del sitio activo, así la proteína que está
  // delante no tapa la química
  function centrar(ms){if(cfg.tipo==="md"){viewer.zoomTo({resn:["GLC","ATP","MG"]},ms);viewer.zoom(.8,ms);}
    else{viewer.zoomTo({model:main},ms);viewer.zoom(1.2,ms);}
    if(cfg.corte){setTimeout(()=>{viewer.setSlab(-cfg.corte[0],cfg.corte[1]);viewer.render();},ms+20);}}
  function ir(k){f=k;slider.value=k;cont.textContent=(k+1)+" / "+nframes;
    // sin fotogramas (vista estática) no se llama a setFrame: 3Dmol dejaría el modelo sin átomos
    if(nframes>1){viewer.setFrame(k).then(()=>{dibujar(k);viewer.render();});}else{dibujar(k);viewer.render();}
    if(G)G.set(Math.min(k,cfg.grafica.x.length-1));
    if(cfg.lecturas){const vals=cfg.lecturas.valores[Math.min(k,cfg.lecturas.valores.length-1)];celdas.forEach((c,i)=>c.textContent=vals[i]);}
    if(cfg.fases)fase.textContent=cfg.fases[Math.min(k,cfg.fases.length-1)];}
  function iniciar(){play.textContent="⏸  Pausa";activo(play,true);
    timer=setInterval(()=>{ir((f+1)%nframes);},(cfg.tipo==="md"?140:220)/vel);}
  function parar(){if(timer){clearInterval(timer);timer=null;}play.textContent="▶  Reproducir";activo(play,false);}
  centrar(0); ir(0); root._ir=ir;
  if(cfg.auto&&nframes>1)iniciar();
  // pausar cuando el visor sale de la pantalla (ahorra batería)
  if("IntersectionObserver" in window){new IntersectionObserver(es=>{es.forEach(e=>{if(!e.isIntersecting)parar();});}).observe(root);}
 }catch(err){const e=root.querySelector(".err");e.style.display="block";e.textContent="Error en el visor 3D: "+err;console.error(err);}
});
})();
</script>
"""
