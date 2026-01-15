# Hexapod Stewart Platform FEM Analysis

Finite Element Analysis of a Stewart platform hexapod using the same 4080 C-beam actuators as the TTC450 gantry.

## Overview

This analysis explores an alternative CNC kinematic configuration - a 6-axis parallel manipulator (Stewart platform) - using the same structural components as the TTC450 Cartesian gantry. The goal is to compare stiffness characteristics between serial (gantry) and parallel (hexapod) kinematic configurations.

### Key Advantage: Symmetric Stiffness

Unlike the TTC450 gantry which has a 3.4x weaker Y-axis, the hexapod provides identical stiffness in X and Y directions due to its symmetric parallel kinematic architecture.

## Comparison: Hexapod vs TTC450 Gantry

Both designs use identical **600mm 4080 C-beam actuators** with Aluminum 6061-T6.

| Load Case | TTC450 Gantry | Hexapod (50% stroke) | Winner |
|-----------|---------------|----------------------|--------|
| **Z-axis weight** (50N down) | 10.3 µm | 1.3 µm | Hexapod (8x) |
| **Heavy cut X** (100N + weight) | 15.7 µm | 3.4 µm | Hexapod (5x) |
| **Heavy cut Y** (100N + weight) | 52.8 µm | 3.3 µm | Hexapod (16x) |
| **Plunge cut** (100-150N Z) | 19.7 µm | 3.1 µm | Hexapod (6x) |

### Key Findings

1. **Hexapod has near-perfect X/Y symmetry** - 3.3-3.4 µm in X and Y

2. **Hexapod is dramatically stiffer** - 5-17x better than TTC450 in all directions

3. **Mid-stroke is the realistic operating point** - Beams are fixed-length with sliding carriages; 50% stroke gives maximum workspace while maintaining stiffness

4. **Steeper strut angles improve lateral stiffness** - At 43° from vertical, struts provide excellent resistance to cutting forces

---

## Results Summary

### Static Analysis (Cutting Loads)

Tool loads applied at platform center with 5kg Z-axis + spindle weight (50N down).

| Load Case | Cutting Force | Tool Deflection |
|-----------|---------------|-----------------|
| Z-axis weight only | 0 N | 1.3 µm |
| **Heavy cut (X)** | 100 N | **3.4 µm** |
| **Heavy cut (Y)** | 100 N | **3.3 µm** |
| **Heavy plunge (Z)** | 100 N | **3.1 µm** |

### Undeformed Mesh

![Hexapod Mesh](docs/images/hexapod_mesh.png?v=9)

### Heavy Cut - X Direction (100N)

![Heavy Cut X](docs/images/hexapod_heavy_cut_x.png?v=9)

### Heavy Cut - Y Direction (100N)

![Heavy Cut Y](docs/images/hexapod_heavy_cut_y.png?v=9)

### Heavy Plunge Cut - Z Direction (150N)

![Heavy Cut Z](docs/images/hexapod_heavy_cut_z.png?v=9)

### Helical Bore Animation

Animation showing the hexapod retracting from an 80mm radius × 180mm deep helical bore (4 revolutions). Platform starts low (struts extended ~2/3) and lifts up (struts retract to ~1/3), demonstrating the parallel kinematic motion. Ball joints at platform, universal joints at base allow strut reorientation without axial rotation.

![Hexapod Helical Bore](docs/images/hexapod_helix.gif?v=10)

---

## Geometry

### Stewart Platform Configuration

The hexapod uses a symmetric 6-6 Stewart platform configuration with paired joints.
Beams are fixed-length 4080 C-beams with sliding carriages for actuation.
The base is a 600mm tall hollow hexagonal prism that accommodates the beam extensions,
with edges aligned to allow struts to pass without intersection.

| Parameter | Value |
|-----------|-------|
| Beam length | 600 mm (4080 C-beam) |
| Stroke position | 50% (mid-stroke) |
| Effective strut length | 300 mm |
| Base pillar | Hexagonal, 700mm across × 600mm tall, 6mm wall |
| Base joint radius | 300 mm |
| Platform radius | 120 mm |
| Platform height | 219 mm |
| Strut angle | 43° from vertical |

### Joint Layout

- **Base joints**: 3 pairs at 120° intervals, each pair separated by 40° (±20°)
- **Platform joints**: 3 pairs at 120° intervals, rotated 30° from base, same 40° spread
- **Strut pairing**: Direct connection (base[i] to platform[i]) creates 3 pairs of parallel struts
- **Beam extension**: 300mm below base joints (inside pillar)

### Strut Cross-Section

Each strut uses a 40×80mm hollow rectangular section (4080 C-beam profile):
- Wall thickness: 1.5mm
- Full 600mm beam shown with base attachment at midpoint

### FEM Mesh

- **Nodes:** 17,541
- **Elements:** 53,013 (tetrahedra)
- **Boundary conditions:** Pillar bottom face fixed

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
- **5-16x better stiffness** than TTC450 gantry across all load cases
- **Near-perfect X/Y symmetry** - 3.3-3.4 µm in both directions
- **Excellent Z-axis stiffness** - 3.1 µm under 150N plunge load

However, for a practical hobby CNC, the Cartesian gantry remains more practical due to simpler control, larger workspace, and easier construction.

The hexapod configuration would be more appropriate for:
- High-precision finishing operations
- Applications requiring isotropic stiffness
- Situations where workspace size is less critical than rigidity
