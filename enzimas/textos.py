"""Idioma del curso: el mismo código sirve al cuaderno en español y al cuaderno en inglés.

El texto fuente está en español y es la propia clave de traducción (al estilo de gettext)::

    from enzimas.textos import t
    ax.set_title(t("Perfil de energía"))
    t("media {x:.2f} Å").format(x=x)          # con números: traducir la plantilla, formatear después

``usar("en")`` activa el inglés (lo hace la segunda celda del cuaderno en inglés, antes de dibujar nada).
En español ``t`` devuelve el texto tal cual, así el cuaderno en español no cambia. En inglés, un texto sin
traducción se devuelve en español y se anota en ``faltan`` (los tests exigen que quede vacío).
"""
from __future__ import annotations

_IDIOMA = "es"
faltan: set[str] = set()
_EN: dict[str, str] | None = None


def usar(idioma: str) -> str:
    """Elige el idioma de gráficas, tarjetas, visores y exploradores: "es" o "en"."""
    global _IDIOMA
    if idioma not in ("es", "en"):
        raise ValueError(f"idioma no soportado: {idioma!r} (use 'es' o 'en')")
    _IDIOMA = idioma
    return _IDIOMA


def idioma() -> str:
    return _IDIOMA


def _diccionario():
    global _EN
    if _EN is None:
        from .traducciones_en import EN
        _EN = EN
    return _EN


def t(texto: str) -> str:
    """Traduce ``texto`` (escrito en español) al idioma activo."""
    if _IDIOMA == "es" or not isinstance(texto, str) or not texto.strip():
        return texto
    trad = _diccionario().get(texto)
    if trad is None:
        faltan.add(texto)
        return texto
    return trad
