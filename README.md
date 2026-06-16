# osa-ha-martini3

A small Python tool that generates primitive starting geometries and Martini 3
coarse-grained force-field files (GROMACS `.gro` + `.itp`) for **hyaluronic acid
(HA)** and **OSA-grafted hyaluronic acid (OSA-HA)** chains of arbitrary length.

It extends the published `carbo2martini3_2.0.py` generator with an octenyl
succinic anhydride (OSA) graft, letting you build OSA-HA at a chosen degree of
substitution (DS) with a reproducible random seed.

## Provenance and license

This repository is a **derivative work**. The core machinery and the HA /
glucopyranose Martini 3 parameters come from `carbo2martini3_2.0.py` by
**Valery Lutsyk and Wojciech Plazinski**, distributed under **CC BY 4.0** as
Supporting Information of:

> V. Lutsyk, W. Plazinski. *Extending the Martini 3 Coarse-Grained Force Field
> to Hyaluronic Acid.* J. Phys. Chem. B **2025**, 129 (9), 2408–2425.
> DOI: [10.1021/acs.jpcb.4c08043](https://doi.org/10.1021/acs.jpcb.4c08043)

The underlying carbohydrate model was first described in:

> V. Lutsyk, P. Wolski, W. Plazinski. *Extending the Martini 3 Coarse-Grained
> Force Field to Carbohydrates.* J. Chem. Theory Comput. **2022**, 18 (8),
> 5089–5107. DOI: [10.1021/acs.jctc.2c00553](https://doi.org/10.1021/acs.jctc.2c00553)

In keeping with CC BY 4.0, this repo is also released under **CC BY 4.0** (see
[`LICENSE`](LICENSE)). You may reuse and modify it, including commercially,
provided you keep attribution to the original authors and indicate changes.

**This derivative is not endorsed by the original authors.**

## What is original vs. derived

**Derived (unchanged from `carbo2martini3_2.0.py`):**
- The `Polysaccharide` class and the `.gro` / `.itp` writers.
- `rotation_matrix`, `pickAngle`, and the geometry placement scheme.
- The HA and glucopyranose generators (`HA`, `GLCA14`, `GLCB12`, `GLCB13`,
  `GLCB14`, `GLCA16`) and their force-field parameters.

**Original additions in this repo:**
- `OSA_HA(units_num, ds_percent, seed)` — builds an OSA-grafted HA chain.
- Supporting functions: `_ha_osa_templates`, `_ha_absolute_coord`,
  `_place_osa_fragment`, `_resolve_repeat_local`.
- OSA bead definitions and the ester-linkage bonded terms attaching the OSA
  fragment to the GlcNAc residue.
- Degree-of-substitution control with a reproducible random seed, and menu
  option 7 in `main()`.

## ⚠️ Scientific caveat (please read)

The published HA parameters were validated against atomistic MD and
octanol–water transfer free energies. **The OSA graft bead types and the
ester-linkage bonded parameters added here have NOT been independently
validated.** Treat them as a preliminary starting point and validate against an
atomistic reference before any production simulation. Outputs are *primitive*
starting geometries intended to be minimised, solvated, and equilibrated before
use.

## Requirements

- Python 3.8+
- NumPy

```bash
pip install numpy
```

## Usage

Run interactively:

```bash
python osa_ha_martini3.py
```

You will be prompted to:
1. Select the linkage / molecule type (choose **7** for OSA_HA).
2. Enter the number of HA disaccharide repeat units (> 1).
3. For OSA_HA: enter the **DS percentage** (fraction of repeats that receive an
   OSA graft) and a **random seed** (for reproducible placement of grafts).

This writes `OSA_HA.gro` and `OSA_HA.itp` (or `HA.gro` / `HA.itp` for plain HA)
to the working directory.

### Downstream

The generated `.gro`/`.itp` are *starting* files. A typical next step is to
solvate and add ions (e.g. with `insane.py`), then energy-minimise and
equilibrate (NVT → NPT) in GROMACS with the standard Martini 3 settings
(reaction-field electrostatics, `epsilon_r = 15`) before production MD.

## Citing

If you use this tool, please cite **both** original papers above. A
machine-readable citation is in [`CITATION.cff`](CITATION.cff).
