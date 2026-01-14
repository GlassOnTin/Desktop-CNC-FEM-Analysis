#!/usr/bin/env python3
"""Run static load case variations on TTC450 gantry.

Applies realistic loading at the tool position:
- Z-axis + spindle weight (~5 kg = 50N)
- Cutting forces in X, Y, Z directions

Based on TTC450 Pro specifications:
- 500W spindle: ~1 kg
- Z-axis carriage + rails: ~2-3 kg
- Total moving mass: ~4-5 kg
"""

import numpy as np
from pathlib import Path
from dataclasses import dataclass
from mpi4py import MPI

import dolfinx
from dolfinx import fem, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
import ufl

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "fem" / "results"
MESH_FILE = RESULTS_DIR / "ttc450_hollow.msh"

# Material: Aluminum 6061-T6
E_MPA = 69000.0     # MPa (N/mm²)
NU = 0.33           # Poisson's ratio
RHO_MM = 2.7e-9     # kg/mm³ (2700 kg/m³)
GRAVITY_MM = 9810.0 # mm/s²


@dataclass
class LoadCase:
    """Definition of a load case."""
    name: str
    description: str
    # Point load at tool position [Fx, Fy, Fz] in Newtons
    # Positive Z is up, so gravity/weight is negative Z
    point_load_N: tuple[float, float, float]
    include_gravity: bool = True


# Define load cases for TTC450 Pro
# Tool position is approximately at gantry midpoint: (0, 20, 200) mm
LOAD_CASES = [
    LoadCase(
        name="z_axis_weight",
        description="Z-axis + 500W spindle weight (5kg)",
        point_load_N=(0, 0, -50),  # 5 kg * 9.81 ≈ 50N down
        include_gravity=True,
    ),
    LoadCase(
        name="light_cut_x",
        description="Light cut: 5kg weight + 20N in X",
        point_load_N=(20, 0, -50),
        include_gravity=True,
    ),
    LoadCase(
        name="light_cut_y",
        description="Light cut: 5kg weight + 20N in Y",
        point_load_N=(0, 20, -50),
        include_gravity=True,
    ),
    LoadCase(
        name="moderate_cut_x",
        description="Moderate cut: 5kg weight + 50N in X",
        point_load_N=(50, 0, -50),
        include_gravity=True,
    ),
    LoadCase(
        name="moderate_cut_y",
        description="Moderate cut: 5kg weight + 50N in Y",
        point_load_N=(0, 50, -50),
        include_gravity=True,
    ),
    LoadCase(
        name="heavy_cut_x",
        description="Heavy cut: 5kg weight + 100N in X",
        point_load_N=(100, 0, -50),
        include_gravity=True,
    ),
    LoadCase(
        name="heavy_cut_y",
        description="Heavy cut: 5kg weight + 100N in Y",
        point_load_N=(0, 100, -50),
        include_gravity=True,
    ),
    LoadCase(
        name="plunge_cut",
        description="Plunge: 5kg weight + 50N downward cutting force",
        point_load_N=(0, 0, -100),  # 50N weight + 50N cutting
        include_gravity=True,
    ),
]


def find_closest_node(coords: np.ndarray, target: np.ndarray) -> int:
    """Find the index of the node closest to the target point."""
    distances = np.linalg.norm(coords - target, axis=1)
    return int(np.argmin(distances))


def run_load_case(
    domain,
    V,
    bc,
    load_case: LoadCase,
    tool_position: np.ndarray,
    coords: np.ndarray,
) -> dict:
    """Run a single load case and return results."""
    from petsc4py import PETSc

    print(f"\n{'─' * 50}")
    print(f"Load Case: {load_case.name}")
    print(f"  {load_case.description}")
    print(f"  Point load: {load_case.point_load_N} N")

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Lamé parameters
    lmbda = E_MPA * NU / ((1 + NU) * (1 - 2 * NU))
    mu = E_MPA / (2 * (1 + NU))

    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(3) + 2 * mu * epsilon(u)

    # Bilinear form
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx

    # Start with gravity or zero load
    if load_case.include_gravity:
        f_body = fem.Constant(domain, default_scalar_type((0.0, 0.0, -RHO_MM * GRAVITY_MM)))
        L = ufl.dot(f_body, v) * ufl.dx
    else:
        f_body = fem.Constant(domain, default_scalar_type((0.0, 0.0, 0.0)))
        L = ufl.dot(f_body, v) * ufl.dx

    # Find node for point load
    Fx, Fy, Fz = load_case.point_load_N
    tool_node_idx = find_closest_node(coords, tool_position)

    if any(f != 0 for f in [Fx, Fy, Fz]):
        actual_pos = coords[tool_node_idx]
        print(f"  Load applied at node {tool_node_idx}: {actual_pos}")

    # Solve using LinearProblem first (handles BCs properly)
    problem = fem.petsc.LinearProblem(
        a, L, bcs=[bc],
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"}
    )
    uh = problem.solve()

    # If we have point loads, we need to solve again with modified RHS
    if any(f != 0 for f in [Fx, Fy, Fz]):
        # Get the assembled system
        A = fem.petsc.assemble_matrix(fem.form(a), bcs=[bc])
        A.assemble()

        b = fem.petsc.assemble_vector(fem.form(L))

        # Add point loads before applying BCs
        # DOFs are ordered as [u0_x, u0_y, u0_z, u1_x, u1_y, u1_z, ...]
        dof_x = 3 * tool_node_idx
        dof_y = 3 * tool_node_idx + 1
        dof_z = 3 * tool_node_idx + 2

        b.array[dof_x] += Fx
        b.array[dof_y] += Fy
        b.array[dof_z] += Fz

        # Apply BCs
        fem.petsc.apply_lifting(b, [fem.form(a)], [[bc]])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem.petsc.set_bc(b, [bc])

        # Solve
        solver = PETSc.KSP().create(domain.comm)
        solver.setOperators(A)
        solver.setType(PETSc.KSP.Type.PREONLY)
        solver.getPC().setType(PETSc.PC.Type.LU)

        uh = fem.Function(V)
        solver.solve(b, uh.x.petsc_vec)
        uh.x.scatter_forward()

        solver.destroy()
        A.destroy()
        b.destroy()

    # Extract results
    u_array = uh.x.array.reshape((-1, 3))
    disp_mag = np.linalg.norm(u_array, axis=1)
    max_disp = disp_mag.max()
    max_idx = disp_mag.argmax()
    max_pos = coords[max_idx]

    # Displacement at tool position
    tool_disp = u_array[tool_node_idx]
    tool_disp_mag = np.linalg.norm(tool_disp)

    print(f"  Results:")
    print(f"    Max displacement: {max_disp:.4f} mm at {max_pos}")
    print(f"    Tool displacement: [{tool_disp[0]:.4f}, {tool_disp[1]:.4f}, {tool_disp[2]:.4f}] mm")
    print(f"    Tool displacement magnitude: {tool_disp_mag:.4f} mm")

    return {
        'name': load_case.name,
        'description': load_case.description,
        'point_load_N': load_case.point_load_N,
        'max_displacement_mm': max_disp,
        'max_position': max_pos,
        'tool_displacement_mm': tool_disp,
        'tool_displacement_mag_mm': tool_disp_mag,
        'u_array': u_array,
        'uh': uh,
    }


def main():
    print("=" * 60)
    print("TTC450 LOAD CASE ANALYSIS")
    print("=" * 60)

    if not MESH_FILE.exists():
        print(f"Mesh file not found: {MESH_FILE}")
        print("Run generate_ttc450_simple.py first")
        return

    # Load mesh
    domain, cell_tags, facet_tags = io.gmshio.read_from_msh(
        str(MESH_FILE), MPI.COMM_WORLD, 0
    )

    coords = domain.geometry.x
    n_nodes = domain.topology.index_map(0).size_local
    n_cells = domain.topology.index_map(3).size_local
    bbox_min = coords.min(axis=0)
    bbox_max = coords.max(axis=0)

    print(f"Mesh: {n_nodes} nodes, {n_cells} elements")
    print(f"Bounding box: [{bbox_min}] to [{bbox_max}]")

    # Tool position: gantry midpoint (center of X-beam)
    # From previous analysis, max displacement was at [0, 20, 200]
    tool_position = np.array([0.0, 20.0, 200.0])
    print(f"Tool position: {tool_position}")

    # Function space
    V = fem.functionspace(domain, ("Lagrange", 1, (3,)))

    # Boundary condition: fix Z=0 face
    z_min = bbox_min[2]
    tol = 1.0

    def base_boundary(x):
        return x[2] < z_min + tol

    bc_dofs = fem.locate_dofs_geometrical(V, base_boundary)
    u_zero = fem.Function(V)
    u_zero.x.array[:] = 0.0
    bc = fem.dirichletbc(u_zero, bc_dofs)

    print(f"Fixed {len(bc_dofs)} DOFs at Z < {z_min + tol:.1f} mm")

    # Run all load cases
    results = []
    for load_case in LOAD_CASES:
        result = run_load_case(domain, V, bc, load_case, tool_position, coords)
        results.append(result)

    # Summary table
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"{'Load Case':<20} {'Load (N)':<20} {'Tool Disp (mm)':<15} {'Max Disp (mm)':<15}")
    print("-" * 70)

    for r in results:
        load_str = f"({r['point_load_N'][0]:.0f}, {r['point_load_N'][1]:.0f}, {r['point_load_N'][2]:.0f})"
        print(f"{r['name']:<20} {load_str:<20} {r['tool_displacement_mag_mm']:<15.4f} {r['max_displacement_mm']:<15.4f}")

    # Save summary to file
    summary_file = RESULTS_DIR / "load_case_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("TTC450 Load Case Analysis Results\n")
        f.write("=" * 70 + "\n\n")
        f.write("Tool position: (0, 20, 200) mm (gantry midpoint)\n")
        f.write("Material: Aluminum 6061-T6 (E=69 GPa, nu=0.33)\n\n")

        f.write(f"{'Load Case':<20} {'Fx':<8} {'Fy':<8} {'Fz':<8} {'Tool Disp':<12} {'Max Disp':<12}\n")
        f.write(f"{'':<20} {'(N)':<8} {'(N)':<8} {'(N)':<8} {'(mm)':<12} {'(mm)':<12}\n")
        f.write("-" * 70 + "\n")

        for r in results:
            Fx, Fy, Fz = r['point_load_N']
            f.write(f"{r['name']:<20} {Fx:<8.0f} {Fy:<8.0f} {Fz:<8.0f} "
                   f"{r['tool_displacement_mag_mm']:<12.4f} {r['max_displacement_mm']:<12.4f}\n")

        f.write("\n\nDetailed tool displacements (X, Y, Z components):\n")
        f.write("-" * 70 + "\n")
        for r in results:
            td = r['tool_displacement_mm']
            f.write(f"{r['name']:<20}: X={td[0]:+.4f}, Y={td[1]:+.4f}, Z={td[2]:+.4f} mm\n")

    print(f"\nSummary saved to: {summary_file}")

    # Save the heavy cut case for visualization
    for r in results:
        if r['name'] == 'heavy_cut_x':
            out_path = RESULTS_DIR / "ttc450_heavy_cut_x.xdmf"
            with io.XDMFFile(MPI.COMM_WORLD, str(out_path), "w") as xdmf:
                xdmf.write_mesh(domain)
                r['uh'].name = "displacement"
                xdmf.write_function(r['uh'])
            print(f"Saved heavy_cut_x results to: {out_path}")
            break


if __name__ == "__main__":
    main()
