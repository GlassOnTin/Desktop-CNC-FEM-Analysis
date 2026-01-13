#!/usr/bin/env python3
"""Run static and modal FEM analysis on TTC450 gantry mesh.

Boundary conditions: Base frame bottom face fixed (Z=0)
Load case: Cutting forces at spindle mount point
"""

import numpy as np
from pathlib import Path
from mpi4py import MPI
from petsc4py import PETSc

import dolfinx
from dolfinx import fem, mesh, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem, assemble_matrix
import ufl

try:
    from slepc4py import SLEPc
    HAS_SLEPC = True
except ImportError:
    HAS_SLEPC = False
    print("Warning: SLEPc not available for modal analysis")

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "fem" / "results"
MESH_FILE = RESULTS_DIR / "ttc450_hollow.msh"

# Material: Aluminum 6061-T6
E_MPA = 69000.0     # MPa (N/mm²)
NU = 0.33           # Poisson's ratio
RHO_MM = 2.7e-9     # kg/mm³ (2700 kg/m³)
GRAVITY_MM = 9810.0 # mm/s²


def run_static_analysis(mesh_path: Path, output_path: Path):
    """Run static linear elasticity analysis.

    - Fixed BC at Z=0 (base frame bottom)
    - Gravity load on entire structure
    - Optional point load at tool position
    """
    print("=" * 60)
    print("STATIC ANALYSIS - TTC450 Gantry")
    print("=" * 60)

    # Load mesh
    domain, cell_tags, facet_tags = io.gmshio.read_from_msh(
        str(mesh_path), MPI.COMM_WORLD, 0
    )

    coords = domain.geometry.x
    n_nodes = domain.topology.index_map(0).size_local
    n_cells = domain.topology.index_map(3).size_local

    bbox_min = coords.min(axis=0)
    bbox_max = coords.max(axis=0)

    print(f"Mesh: {n_nodes} nodes, {n_cells} elements")
    print(f"Bounding box: [{bbox_min}] to [{bbox_max}]")

    # Function space (P1 vector)
    V = fem.functionspace(domain, ("Lagrange", 1, (3,)))

    # Trial and test functions
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Lamé parameters
    lmbda = E_MPA * NU / ((1 + NU) * (1 - 2 * NU))
    mu = E_MPA / (2 * (1 + NU))

    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(3) + 2 * mu * epsilon(u)

    # Bilinear form (stiffness)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx

    # Body force (gravity in -Z)
    f = fem.Constant(domain, default_scalar_type((0.0, 0.0, -RHO_MM * GRAVITY_MM)))
    L = ufl.dot(f, v) * ufl.dx

    # Boundary condition: fix Z=0 face (base frame bottom)
    z_min = bbox_min[2]
    tol = 1.0  # mm tolerance

    def base_boundary(x):
        return x[2] < z_min + tol

    bc_dofs = fem.locate_dofs_geometrical(V, base_boundary)
    u_zero = fem.Function(V)
    u_zero.x.array[:] = 0.0
    bc = fem.dirichletbc(u_zero, bc_dofs)

    print(f"Fixed {len(bc_dofs)} DOFs at Z < {z_min + tol:.1f} mm")

    # Solve
    problem = LinearProblem(
        a, L, bcs=[bc],
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"}
    )
    uh = problem.solve()

    # Results
    u_array = uh.x.array.reshape((-1, 3))
    disp_mag = np.linalg.norm(u_array, axis=1)
    max_disp = disp_mag.max()
    max_idx = disp_mag.argmax()
    max_pos = coords[max_idx]

    print(f"\nResults (gravity only):")
    print(f"  Max displacement: {max_disp:.4f} mm at {max_pos}")
    print(f"  Max Z-displacement: {np.abs(u_array[:, 2]).max():.4f} mm")

    # Save results
    with io.XDMFFile(MPI.COMM_WORLD, str(output_path), "w") as xdmf:
        xdmf.write_mesh(domain)
        uh.name = "displacement"
        xdmf.write_function(uh)

    print(f"\nSaved to: {output_path}")

    return {
        'max_displacement_mm': max_disp,
        'max_position': max_pos,
        'u_array': u_array,
        'domain': domain,
        'V': V,
    }


def run_modal_analysis(mesh_path: Path, output_dir: Path, n_modes: int = 6):
    """Run modal analysis to find natural frequencies.

    - Fixed BC at Z=0 (base frame bottom)
    - Uses shift-invert spectral transformation for robust low-frequency extraction
    """
    if not HAS_SLEPC:
        print("SLEPc not available - skipping modal analysis")
        return None

    print("\n" + "=" * 60)
    print("MODAL ANALYSIS - TTC450 Gantry")
    print("=" * 60)

    # Load mesh
    domain, _, _ = io.gmshio.read_from_msh(str(mesh_path), MPI.COMM_WORLD, 0)

    coords = domain.geometry.x
    n_nodes = domain.topology.index_map(0).size_local
    bbox_min = coords.min(axis=0)

    print(f"Mesh: {n_nodes} nodes")

    # Function space
    V = fem.functionspace(domain, ("Lagrange", 1, (3,)))

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Lamé parameters
    lmbda = E_MPA * NU / ((1 + NU) * (1 - 2 * NU))
    mu = E_MPA / (2 * (1 + NU))

    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(3) + 2 * mu * epsilon(u)

    # Stiffness and mass forms
    k_form = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    m_form = RHO_MM * ufl.dot(u, v) * ufl.dx

    # Get boundary DOFs for penalty
    z_min = bbox_min[2]
    tol = 1.0

    def base_boundary(x):
        return x[2] < z_min + tol

    bc_dofs = fem.locate_dofs_geometrical(V, base_boundary)
    bc_dofs_set = set(bc_dofs)

    print(f"Constrained {len(bc_dofs_set)} DOFs at base")

    # Assemble matrices
    K = assemble_matrix(fem.form(k_form))
    K.assemble()

    M = assemble_matrix(fem.form(m_form))
    M.assemble()

    # Apply BCs via penalty - use smaller penalty to avoid numerical issues
    # The penalty should be large relative to max stiffness but not extreme
    diag_vals = K.getDiagonal()
    k_max = diag_vals.array.max()
    penalty = k_max * 1e6  # Penalty relative to max stiffness
    print(f"Using penalty: {penalty:.2e} (k_max={k_max:.2e})")

    for dof in bc_dofs_set:
        K.setValueLocal(dof, dof, penalty, addv=PETSc.InsertMode.ADD_VALUES)
    K.assemble()

    # Eigenvalue solver with shift-invert for low frequencies
    eps = SLEPc.EPS().create(MPI.COMM_WORLD)
    eps.setOperators(K, M)
    eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)

    # Use shift-invert spectral transformation
    # Target eigenvalues near sigma (shift) - helps find lowest modes
    sigma_shift = 1.0  # Small positive shift
    st = eps.getST()
    st.setType(SLEPc.ST.Type.SINVERT)
    st.setShift(sigma_shift)

    # Configure Krylov solver for shift-invert
    ksp = st.getKSP()
    ksp.setType(PETSc.KSP.Type.PREONLY)
    pc = ksp.getPC()
    pc.setType(PETSc.PC.Type.LU)

    eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
    eps.setTarget(sigma_shift)
    eps.setType(SLEPc.EPS.Type.KRYLOVSCHUR)

    nev = n_modes + 10
    ncv = min(nev * 4, K.getSize()[0])
    eps.setDimensions(nev=nev, ncv=ncv)
    eps.setTolerances(tol=1e-6, max_it=2000)
    eps.setFromOptions()

    print(f"Solving for {n_modes} modes (shift-invert, sigma={sigma_shift})...")
    eps.solve()

    n_conv = eps.getConverged()
    print(f"Converged: {n_conv} eigenvalues")

    # Extract results
    vr = K.createVecRight()
    vi = K.createVecRight()

    eigenvalues = []
    frequencies = []
    mode_shapes = []

    # Expected frequency range for gantry: 50-500 Hz
    freq_min = 10.0   # Hz - filter out rigid body modes
    freq_max = 2000.0  # Hz - filter out penalty-induced modes

    for i in range(n_conv):
        eigval = eps.getEigenpair(i, vr, vi)

        # Skip very small or negative eigenvalues (rigid body modes)
        if eigval.real < 1e-3:
            continue

        omega = np.sqrt(eigval.real)
        freq = omega / (2 * np.pi)

        # Filter to expected frequency range
        if freq < freq_min or freq > freq_max:
            continue

        eigenvalues.append(eigval.real)
        frequencies.append(freq)

        mode = vr.array.copy().reshape((-1, 3))
        mode = mode / np.linalg.norm(mode)
        mode_shapes.append(mode)

        if len(frequencies) >= n_modes:
            break

    eps.destroy()
    vr.destroy()
    vi.destroy()

    print(f"\nNatural Frequencies:")
    for i, f in enumerate(frequencies):
        print(f"  Mode {i+1}: {f:.1f} Hz")

    # Save mode shapes
    output_dir.mkdir(parents=True, exist_ok=True)

    for i, (mode, freq) in enumerate(zip(mode_shapes, frequencies)):
        u_mode = fem.Function(V)
        u_mode.x.array[:] = mode.flatten()
        u_mode.name = f"mode_{i+1}"

        out_path = output_dir / f"mode_{i+1:02d}_{freq:.1f}Hz.xdmf"
        with io.XDMFFile(MPI.COMM_WORLD, str(out_path), "w") as xdmf:
            xdmf.write_mesh(domain)
            xdmf.write_function(u_mode)
        print(f"  Saved: {out_path.name}")

    return {
        'frequencies': frequencies,
        'mode_shapes': mode_shapes,
        'n_modes': len(frequencies),
    }


def main():
    if not MESH_FILE.exists():
        print(f"Mesh file not found: {MESH_FILE}")
        print("Run generate_ttc450_simple.py first")
        return

    # Static analysis
    static_out = RESULTS_DIR / "ttc450_static_gravity.xdmf"
    static_results = run_static_analysis(MESH_FILE, static_out)

    # Modal analysis
    modes_dir = RESULTS_DIR / "ttc450_modes"
    modal_results = run_modal_analysis(MESH_FILE, modes_dir, n_modes=6)

    # Summary
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"Static: Max deflection = {static_results['max_displacement_mm']:.4f} mm")

    if modal_results and modal_results['n_modes'] > 0:
        print(f"Modal: {modal_results['n_modes']} modes found")
        print(f"  Mode 1: {modal_results['frequencies'][0]:.1f} Hz")
        if len(modal_results['frequencies']) > 1:
            print(f"  Mode 2: {modal_results['frequencies'][1]:.1f} Hz")
    else:
        print("Modal: No modes found in expected frequency range")


if __name__ == "__main__":
    main()
