# Hexapod Stewart Platform FEM Analysis

Finite Element Analysis of a Stewart platform hexapod using the same 4080 C-beam actuators as the TTC450 gantry.

## Overview

This analysis explores an alternative CNC kinematic configuration - a 6-axis parallel manipulator (Stewart platform) - using the same structural components as the TTC450 Cartesian gantry. The goal is to compare stiffness characteristics between serial (gantry) and parallel (hexapod) kinematic configurations.

### Key Advantage: Symmetric Stiffness

Unlike the TTC450 gantry which has a 3.4x weaker Y-axis, the hexapod provides identical stiffness in X and Y directions due to its symmetric parallel kinematic architecture.

## Comparison: Hexapod vs TTC450 Gantry

Both designs use identical **600mm 4080 C-beam actuators** with Aluminum 6061-T6.

| Load Case | TTC450 Gantry | Hexapod | Winner |
|-----------|---------------|---------|--------|
| **Z-axis weight** (50N down) | 10.3 µm | 1.2 µm | Hexapod (8.6x) |
| **Heavy cut X** (100N + weight) | 15.7 µm | 15.2 µm | Tie |
| **Heavy cut Y** (100N + weight) | 52.8 µm | 15.3 µm | Hexapod (3.5x) |
| **Plunge cut** (100-150N Z) | 19.7 µm | 3.0 µm | Hexapod (6.6x) |

### Key Findings

1. **Hexapod has perfect X/Y symmetry** - 15.2 µm in X, 15.3 µm in Y (essentially identical)

2. **Hexapod is much stiffer in Z** - 8.6x better for static weight, 6.6x better for plunge cuts due to axial strut loading

3. **X-direction is now equal** - Both designs achieve ~15-16 µm for 100N lateral load

4. **Hexapod eliminates the weak axis** - TTC450's Y-direction (53 µm) is 3.5x worse than hexapod (15 µm)

---

## Results Summary

### Static Analysis (Cutting Loads)

Tool loads applied at platform center with 5kg Z-axis + spindle weight (50N down).

| Load Case | Cutting Force | Tool Deflection |
|-----------|---------------|-----------------|
| Z-axis weight only | 0 N | 1.2 µm |
| **Heavy cut (X)** | 100 N | **15.2 µm** |
| **Heavy cut (Y)** | 100 N | **15.3 µm** |
| **Heavy plunge (Z)** | 100 N | **3.0 µm** |

### Undeformed Mesh

![Hexapod Mesh](docs/images/hexapod_mesh.png?v=3)

### Heavy Cut - X Direction (100N)

![Heavy Cut X](docs/images/hexapod_heavy_cut_x.png?v=3)

### Heavy Cut - Y Direction (100N)

![Heavy Cut Y](docs/images/hexapod_heavy_cut_y.png?v=3)

### Heavy Plunge Cut - Z Direction (150N)

![Heavy Cut Z](docs/images/hexapod_heavy_cut_z.png?v=3)

---

## Geometry

### Stewart Platform Configuration

The hexapod uses a symmetric 6-6 Stewart platform configuration with paired joints:

| Parameter | Value |
|-----------|-------|
| Base radius | 300 mm |
| Platform radius | 120 mm |
| Strut length | 600 mm (all identical) |
| Platform height | 564 mm |
| Strut angle | 20° from vertical (all identical) |

### Joint Layout

- **Base joints**: 3 pairs at 120° intervals, each pair separated by 40° (±20°)
- **Platform joints**: 3 pairs at 120° intervals, rotated 30° from base, same 40° spread
- **Strut pairing**: Direct connection (base[i] to platform[i]) creates 3 pairs of parallel struts

### Strut Cross-Section

Each strut uses a simplified 40x40mm hollow rectangular section:
- Wall thickness: 1.5mm
- Equivalent to the strong axis of the 4080 C-beam

### FEM Mesh

- **Nodes:** 7,616
- **Elements:** 22,662 (tetrahedra)
- **Boundary conditions:** Base plate bottom face fixed

---

## Design Analysis

### Why Hexapod is Stiffer

1. **Axial vs Bending Loading**
   - Gantry: Cutting forces cause beam *bending* (I/c matters)
   - Hexapod: Cutting forces cause strut *tension/compression* (A matters)
   - Axial stiffness scales with area; bending stiffness scales with I/L³

2. **Load Distribution**
   - Gantry: Load path through sequential joints creates compliance accumulation
   - Hexapod: Parallel struts share load simultaneously

3. **Symmetric Architecture**
   - Gantry: X and Y axes have fundamentally different structural paths
   - Hexapod: 6 struts arranged symmetrically about vertical axis

### Trade-offs

| Aspect | TTC450 Gantry | Hexapod |
|--------|---------------|---------|
| Stiffness symmetry | Poor (Y 3.4x weaker) | Excellent |
| Workspace volume | Large rectangular | Smaller, complex shape |
| Control complexity | Simple Cartesian | 6-DOF inverse kinematics |
| Singularities | None | Exist at workspace limits |
| Self-interference | None | Struts may collide |
| Cable management | Simple | Complex (moving platform) |

### Practical Considerations

For a hobby CNC, the hexapod's advantages in stiffness may be offset by:
- Complex motion control (6 coordinated axes)
- Reduced workspace efficiency
- Difficult cable routing to spindle
- Higher actuator count (6 vs 3 linear axes)

---

## Quick Start

### Generate Hexapod Mesh
```bash
python fem/geometry/generate_hexapod.py
```

### Run Load Case Analysis
```bash
python fem/analysis/run_hexapod_analysis.py
```

### Generate Visualizations
```bash
python fem/visualization/render_hexapod.py
```

### View Results in ParaView
```bash
paraview fem/results/hexapod_heavy_cut_x.xdmf
```

---

## File Structure

```
fem/
├── geometry/
│   └── generate_hexapod.py      # Stewart platform geometry generator
├── analysis/
│   └── run_hexapod_analysis.py  # Static load case analysis
├── visualization/
│   └── render_hexapod.py        # PyVista batch rendering
└── results/
    ├── hexapod.msh              # Generated mesh
    ├── hexapod_geometry.json    # Joint positions
    ├── hexapod_*.xdmf           # Load case results
    └── hexapod_load_case_summary.txt
```

---

## Conclusions

The hexapod Stewart platform provides:
- **2.5x better stiffness** in the weak (Y) direction compared to TTC450
- **Symmetric stiffness** in X and Y (~21 µm for 100N)
- **2x better Z-axis stiffness** due to axial strut loading

However, for a practical hobby CNC, the Cartesian gantry remains more practical due to simpler control, larger workspace, and easier construction.

The hexapod configuration would be more appropriate for:
- High-precision finishing operations
- Applications requiring isotropic stiffness
- Situations where workspace size is less critical than rigidity
