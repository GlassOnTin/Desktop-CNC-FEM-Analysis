# TTC450 FEM Baseline (Parametric Geometry + FEniCSx)

This repo builds a **null baseline FEM model** of a TTC450 Pro‑style CNC, with a focus on:
- Parametric, mesh‑stable C‑beam geometry (no STL dependency)
- Separate parts with tagged coupling surfaces
- FEniCSx workflow, including MPC‑based tie coupling

The geometry is intentionally simplified to support robust meshing and repeatable optimization.

---

## Quick start

Create the **combined mesh** (single `.msh`, multiple parts, tagged interfaces):
```bash
python fem/analysis/run_ttc450_mpc.py --mesh-only
```

Create the mesh **and** write XDMF + tag counts:
```bash
python fem/analysis/run_ttc450_mpc.py --mesh-only --self-check --write-xdmf --write-tags-json
```

Run the **MPC tie solve** (requires `dolfinx_mpc`):
```bash
python fem/analysis/run_ttc450_mpc.py
```

---

## What’s in the model

### Geometry (parametric)
- **Base frame**: 600×600 of 2020 extrusion (optional in mesh)
- **Y‑axis beams**: 2× C‑beam 4080, aligned along Y
- **X‑gantry beam**: 1× C‑beam 4080, aligned along X
- **Riser plates**: 6 mm thick, attached to gantry ends
- **Z‑axis beam**: 150 mm C‑beam section (simplified as block in combined mesh)

### Meshing strategy
- **C‑beam profile** is **parametric** and mesh‑stable (open “C” shape)
- **Swept volumes** for beams (prismatic, robust)
- **Riser plates** as flat blocks
- **Interface facets tagged** for tie constraints

---

## Key scripts

### 1) Mesh generation

Generate per‑part meshes (for separate coupling in FEniCSx):
```bash
python - <<'PY'
from fem.geometry.mesh_generator import generate_ttc450_hybrid_parts
from fem.config import OUTPUT_DIR

generate_ttc450_hybrid_parts(
    OUTPUT_DIR / "parts_cshape",
    mesh_size_min=8.0,
    mesh_size_max=25.0,
    profile="cshape",
    include_base=False,
    beam_shape="profile",
    interface_gap=0.0,
)
PY
```

Generate a **combined mesh** (multiple volumes + facet tags):
```bash
python - <<'PY'
from fem.geometry.mesh_generator import generate_ttc450_combined_mesh
from fem.config import OUTPUT_DIR

generate_ttc450_combined_mesh(
    OUTPUT_DIR / "ttc450_combined.msh",
    mesh_size_min=8.0,
    mesh_size_max=25.0,
    profile="cshape",
    interface_gap=0.2,
    include_base=False,
)
PY
```

---

### 2) Facet tag inspection + XDMF export

```bash
python fem/analysis/fenicsx_tie_parts.py
```

This prints tag IDs and writes XDMF files with facet tags for each part.

---

### 3) MPC tie solve (FEniCSx)

```bash
python fem/analysis/run_ttc450_mpc.py
```

Defaults:
- `--gap 0.2`
- `--mesh-min 8`
- `--mesh-max 25`
- `--profile cshape`

Optional flags:
- `--mesh-only` generate mesh only
- `--self-check` print facet tag counts
- `--write-xdmf` export XDMF for ParaView
- `--write-tags-json` write tag IDs + counts JSON

---

## Tag names used for coupling

These are emitted as **Physical Groups** in the mesh:
- `x_beam_end_left`
- `x_beam_end_right`
- `riser_left_inner`
- `riser_right_inner`

Volume names:
- `x_beam`, `y_beam_left`, `y_beam_right`
- `riser_left`, `riser_right`
- `z_block`, `z_interface`
- `base_frame` (optional)

---

## Dependencies

Required:
- `gmsh`
- `meshio`
- `numpy`, `scipy`

For MPC ties:
- `dolfinx_mpc` (see `requirements.txt`)

FEniCSx should already be installed in your environment.

---

## Notes / gotchas

- **Parametric C‑shape** is used for robustness (no fragile STL or boolean cuts).
- The combined mesh uses a small **interface gap** (default 0.2 mm) to avoid coincident‑facet issues.
- If you change the gap, update:
  - `--gap` in `run_ttc450_mpc.py`
  - `INTERFACE_GAP` in `fenicsx_mpc_tie_combined.py`

---

## Next steps

Common extensions:
- Replace `z_block` with a swept C‑beam when needed
- Add refined mesh sizing near interfaces
- Add real bearing blocks or rails as separate parts
- Integrate loads/boundaries in `fenicsx_static.py`
