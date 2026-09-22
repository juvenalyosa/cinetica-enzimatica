# Cinética enzimática, paso a paso, con simulaciones QM/MM

[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/juvenalyosa/cinetica-enzimatica/blob/main/notebooks/Cinetica_Enzimatica_Glucoquinasa.ipynb)
[![Licencia MIT](https://img.shields.io/badge/licencia-MIT-blue.svg)](LICENSE)

Curso abierto, en español, que explica la cinética enzimática **desde cero** usando una
enzima real, la **glucoquinasa humana** (hexoquinasa IV, EC 2.7.1.2), y simulaciones
moleculares que el estudiante ejecuta celda a celda en Google Colab:

* estructura cristalina del complejo catalítico (PDB 3FGU: glucosa + ATP + Mg²⁺),
* dinámica molecular clásica (OpenMM, Amber ff14SB),
* reacción química por QM/MM (MOPAC PM7 + campo de la enzima): perfil de energía,
  estado de transición, frecuencia imaginaria, camino de reacción,
* de la barrera de energía a *k*<sub>cat</sub> (teoría del estado de transición, Eyring),
* Michaelis–Menten, linealizaciones, inhibición, cooperatividad (Hill), temperatura y pH.

Cada sección empieza con una explicación muy simple y sube de nivel poco a poco.

## Cómo usarlo

**En Colab (recomendado):** pulsa el botón *Abrir en Colab*. El notebook tiene dos modos:

| Modo | Qué hace | Tiempo |
|---|---|---|
| `rapido` | Lee los resultados precalculados que viajan en `data/precalculado/` y dibuja todo | ~10 min |
| `completo` | Reejecuta MD y QM/MM en la máquina de Colab (CPU) | horas |

**En tu computadora:**

```bash
git clone https://github.com/juvenalyosa/cinetica-enzimatica.git
cd cinetica-enzimatica
pip install -e ".[sim,dev]"
jupyter lab notebooks/Cinetica_Enzimatica_Glucoquinasa.ipynb
```

Para el modo completo hace falta además MOPAC (https://github.com/openmopac/mopac) y,
para regenerar la preparación del sistema, AmberTools (antechamber, tleap).

## Estructura del repositorio

```
notebooks/   Cinetica_Enzimatica_Glucoquinasa.ipynb   el curso, celda a celda
enzimas/     paquete Python
   kinetics.py      leyes de velocidad, ajustes, simulación del mecanismo, Eyring
   viz.py           gráficos y visores 3D con un estilo único
   glucoquinasa.py  modelo QM/MM de la reacción (PM7 + Amber, NEB, dímero, frecuencias)
   qmmm_mopac.py    interfaz MOPAC (SADDLE = QST2, TS, FORCETS, IRC), de la suite Leonardo
   openmm_qmmm.py   utilidades OpenMM/QM-MM, de la suite Leonardo
   cluster.py       modelos de clúster con átomos de enlace, de la suite Leonardo
   datos.py         acceso a los datos precalculados
   colab_setup.py   instalación en Colab
scripts/     01_preparar_sistema.py, 02_dinamica_molecular.py, 03_qmmm_reaccion.py
data/        raw/ (PDB) y precalculado/ (sistema, MD, QM/MM)
tests/       pytest
```

## Cómo se generaron los datos precalculados

1. `scripts/01_preparar_sistema.py`: 3FGU, cadena A; bucles faltantes con PDBFixer;
   AMP-PNP convertido en ATP; glucosa y ATP con GAFF2/AM1-BCC; proteína ff14SB; Mg²⁺ y K⁺
   (Li–Merz); caja de agua TIP3P y Na⁺ neutralizantes (tleap).
2. `scripts/02_dinamica_molecular.py`: minimización, calentamiento, equilibración NPT y
   1 ns de producción con OpenMM.
3. `scripts/03_qmmm_reaccion.py`: región QM (glucosa, trifosfato, Asp205, Lys169, Thr228,
   Mg²⁺ y su agua) con PM7 en el potencial electrostático de la enzima más
   Lennard-Jones QM–MM (optimizador externo con ASE); escaneo de la coordenada de reacción,
   NEB con imagen trepadora, método del dímero, frecuencias numéricas y descenso desde el TS.
4. `scripts/03b_qmmm_agua.py` y `scripts/03c_qmmm_agua_validacion.py`: el mismo clúster en
   agua implícita (COSMO) con los métodos nativos de MOPAC: `SADDLE` (QST2) como estimación,
   refinamiento con el método del dímero y la palabra clave `TS`, `FORCETS` y frecuencias
   numéricas; el `IRC` de MOPAC abortó en esa superficie y el camino se obtuvo por descenso
   desde el TS.

Resultados principales (PM7, entorno fijo, energías electrónicas): barrera en la enzima
19.1 kcal/mol con una sola frecuencia imaginaria (−149 cm⁻¹); sitio activo en agua
21.6 kcal/mol; ΔE de reacción +13.7 kcal/mol; distancia Pγ–O6 media en la MD 3.44 Å.

Las barreras PM7 son **semicuantitativas**: el objetivo es entender los conceptos, no
reproducir el valor experimental con precisión química.

## Créditos y licencia

Los módulos QM/MM provienen de la suite [Leonardo](https://github.com/juvenalyosa/Leonardo)
(Juvenal Yosa, MIT). Estructura: Petit *et al.*, PDB 3FGU. Código con licencia MIT.
