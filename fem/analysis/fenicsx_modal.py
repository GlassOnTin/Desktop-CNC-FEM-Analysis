"""Modal analysis (natural frequencies) using FEniCSx and SLEPc.

Solves the generalized eigenvalue problem: K*u = omega^2 * M * u
"""

import numpy as np
from pathlib import Path
from typing import Optional, Dict, Tuple, List
from mpi4py import MPI
from petsc4py import PETSc

import dolfinx
from dolfinx import fem, mesh, io, default_scalar_type
from dolfinx.fem.petsc import assemble_matrix
import ufl

try:
    from slepc4py import SLEPc
    HAS_SLEPC = True
except ImportError:
    HAS_SLEPC = False
    print("Warning: SLEPc not available. Modal analysis will not work.")

from ..config import (
    GRAVITY, MATERIALS, CBEAM_40X80, DEFAULT_GANTRY,
    BOUNDARY_TOLERANCE_MM, N_MODES_DEFAULT, SPINDLE_MASS_KG,
    analytical_beam_frequency
)


def create_stiffness_mass_forms(
    domain: mesh.Mesh,
    E: float,
    nu: float,
    rho: float
) -> Tuple[ufl.Form, ufl.Form, fem.FunctionSpace]:
    """Create stiffness and mass matrix forms for modal analysis.

    Args:
        domain: DOLFINx mesh
        E: Young's modulus (N/mm² = MPa)
        nu: Poisson's ratio
        rho: Density (kg/mm³)

    Returns:
        (stiffness_form, mass_form, function_space)
    """
    # Vector function space for displacement (P1 elements)
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    # Trial and test functions
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Lamé parameters
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))

    # Strain tensor
    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    # Stress tensor
    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

    # Stiffness form
    k = ufl.inner(sigma(u), epsilon(v)) * ufl.dx

    # Mass form
    m = rho * ufl.dot(u, v) * ufl.dx

    return k, m, V


def locate_boundary_dofs(
    V: fem.FunctionSpace,
    domain: mesh.Mesh,
    boundary_type: str = 'fixed'
) -> List[np.ndarray]:
    """Locate DOFs on boundary regions for different support conditions.

    Args:
        V: Function space
        domain: Mesh
        boundary_type: 'cantilever', 'fixed', or 'simply_supported'

    Returns:
        List of dof index arrays for boundary conditions
    """
    coords = domain.geometry.x
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    tol = BOUNDARY_TOLERANCE_MM

    bc_dofs = []

    def left_boundary(x):
        return x[0] < x_min + tol

    def right_boundary(x):
        return x[0] > x_max - tol

    if boundary_type == 'cantilever':
        left_dofs = fem.locate_dofs_geometrical(V, left_boundary)
        bc_dofs.append(left_dofs)

    elif boundary_type == 'fixed':
        left_dofs = fem.locate_dofs_geometrical(V, left_boundary)
        right_dofs = fem.locate_dofs_geometrical(V, right_boundary)
        bc_dofs.append(left_dofs)
        bc_dofs.append(right_dofs)

    elif boundary_type == 'simply_supported':
        # Approximate as fixed-fixed for modal analysis
        left_dofs = fem.locate_dofs_geometrical(V, left_boundary)
        right_dofs = fem.locate_dofs_geometrical(V, right_boundary)
        bc_dofs.append(left_dofs)
        bc_dofs.append(right_dofs)

    return bc_dofs


def solve_modal(
    mesh_path: Path,
    n_modes: int = N_MODES_DEFAULT,
    material_name: str = 'aluminum_6061_t6',
    boundary_type: str = 'fixed',
    point_masses: Optional[List[Tuple[float, float, float, float, str]]] = None
) -> Dict:
    """Solve modal analysis problem.

    Args:
        mesh_path: Path to mesh file (XDMF or Gmsh .msh)
        n_modes: Number of modes to compute
        material_name: Key from MATERIALS dict
        boundary_type: 'cantilever', 'fixed', or 'simply_supported'
        point_masses: List of (mass_kg, x, y, z, name) tuples

    Returns:
        Dict with eigenvalues, frequencies, and mode shapes
    """
    if not HAS_SLEPC:
        raise RuntimeError("SLEPc required for modal analysis. Install with: apt install python3-slepc4py")

    # Get material properties
    mat = MATERIALS[material_name]
    E_mpa = mat['E'] / 1e6  # Pa -> MPa (N/mm²)
    nu = mat['nu']
    rho_mm = mat['rho'] / 1e9  # kg/m³ -> kg/mm³

    # Load mesh
    mesh_path = Path(mesh_path)
    if mesh_path.suffix == '.xdmf':
        with io.XDMFFile(MPI.COMM_WORLD, str(mesh_path), "r") as xdmf:
            domain = xdmf.read_mesh(name="Grid")
    elif mesh_path.suffix == '.msh':
        domain, _, _ = io.gmshio.read_from_msh(str(mesh_path), MPI.COMM_WORLD, 0)
    else:
        raise ValueError(f"Unsupported mesh format: {mesh_path.suffix}")

    n_nodes = domain.topology.index_map(0).size_local
    n_cells = domain.topology.index_map(3).size_local
    print(f"Loaded mesh: {n_nodes} nodes, {n_cells} cells")

    # Create forms
    k_form, m_form, V = create_stiffness_mass_forms(domain, E_mpa, nu, rho_mm)

    # Get boundary DOFs
    bc_dof_lists = locate_boundary_dofs(V, domain, boundary_type)

    # Create boundary conditions - collect BC DOFs
    bc_dofs_set = set()
    for dofs in bc_dof_lists:
        bc_dofs_set.update(dofs)

    print(f"Applied {len(bc_dof_lists)} boundary condition(s) ({boundary_type}), {len(bc_dofs_set)} constrained DOFs")

    # Assemble matrices without BCs
    k_compiled = fem.form(k_form)
    m_compiled = fem.form(m_form)

    K = assemble_matrix(k_compiled)
    K.assemble()

    M = assemble_matrix(m_compiled)
    M.assemble()

    # Apply boundary conditions by setting diagonal to large values
    # This effectively fixes those DOFs (large stiffness, normal mass -> high frequency)
    large_k = 1e20  # Very large stiffness
    for dof in bc_dofs_set:
        K.setValueLocal(dof, dof, large_k, addv=PETSc.InsertMode.ADD_VALUES)
    K.assemble()

    # Set up eigenvalue solver
    eps = SLEPc.EPS().create(MPI.COMM_WORLD)
    eps.setOperators(K, M)
    eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)  # Generalized Hermitian

    # Request smallest magnitude eigenvalues (lowest frequencies)
    eps.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)

    # Use Krylov-Schur - the default and most robust
    eps.setType(SLEPc.EPS.Type.KRYLOVSCHUR)

    # Set dimensions: nev=requested, ncv=working subspace size (should be much larger)
    nev = n_modes + 20  # Extra to filter
    ncv = min(nev * 3, K.getSize()[0])  # Working space
    eps.setDimensions(nev=nev, ncv=ncv)

    # Set tolerances
    eps.setTolerances(tol=1e-8, max_it=1000)

    # Configure solver
    eps.setFromOptions()

    # Solve
    print(f"Solving for {n_modes} modes (dim: {K.getSize()})...")
    eps.solve()

    n_conv = eps.getConverged()
    print(f"Converged eigenvalues: {n_conv}")

    # Extract results
    eigenvalues = []
    frequencies = []
    mode_shapes = []

    # Create vectors for eigenvector extraction
    vr = K.createVecRight()
    vi = K.createVecRight()

    for i in range(min(n_conv, n_modes + 6)):
        eigval = eps.getEigenpair(i, vr, vi)

        # Skip rigid body modes (near-zero eigenvalues)
        if abs(eigval.real) < 1e-6:
            continue

        # Natural frequency: omega^2 = eigenvalue, f = omega / (2*pi)
        if eigval.real > 0:
            omega = np.sqrt(eigval.real)
            freq = omega / (2 * np.pi)

            eigenvalues.append(eigval.real)
            frequencies.append(freq)

            # Extract mode shape
            mode = vr.array.copy().reshape((-1, 3))
            # Normalize
            mode = mode / np.linalg.norm(mode)
            mode_shapes.append(mode)

        if len(frequencies) >= n_modes:
            break

    # Clean up
    eps.destroy()
    vr.destroy()
    vi.destroy()

    print(f"\nNatural frequencies:")
    for i, f in enumerate(frequencies[:6]):
        print(f"  Mode {i+1}: {f:.2f} Hz")

    # Compare with analytical (beam theory)
    coords = domain.geometry.x
    L = coords[:, 0].max() - coords[:, 0].min()
    f_analytical = analytical_beam_frequency(
        L=L,
        E=mat['E'],
        I=CBEAM_40X80['Iy'],
        rho=mat['rho'],
        A=CBEAM_40X80['area'],
        mode=1,
        support='fixed_fixed' if boundary_type != 'cantilever' else 'simply_supported'
    )
    print(f"\nAnalytical Mode 1 (beam theory): {f_analytical:.2f} Hz")

    return {
        'eigenvalues': eigenvalues,
        'frequencies': frequencies,
        'mode_shapes': mode_shapes,
        'n_modes': len(frequencies),
        'mesh': domain,
        'function_space': V,
        'boundary_type': boundary_type,
        'material': material_name,
        'analytical_mode1': f_analytical,
    }


def save_mode_shapes(results: Dict, output_dir: Path) -> List[Path]:
    """Save mode shapes to VTK files.

    Args:
        results: Results from solve_modal()
        output_dir: Directory for output files

    Returns:
        List of paths to saved files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    domain = results['mesh']
    V = results['function_space']
    mode_shapes = results['mode_shapes']
    frequencies = results['frequencies']

    saved_files = []

    for i, (mode, freq) in enumerate(zip(mode_shapes, frequencies)):
        u = fem.Function(V)
        u.x.array[:] = mode.flatten()
        u.name = f"mode_{i+1}"

        output_path = output_dir / f"mode_{i+1:02d}_{freq:.1f}Hz.xdmf"

        with io.XDMFFile(MPI.COMM_WORLD, str(output_path), "w") as xdmf:
            xdmf.write_mesh(domain)
            xdmf.write_function(u)

        saved_files.append(output_path)
        print(f"Saved mode {i+1} ({freq:.1f} Hz) to {output_path}")

    return saved_files


def chatter_analysis(
    frequencies: List[float],
    spindle_rpm_range: Tuple[int, int] = (8000, 24000),
    flute_counts: List[int] = [1, 2, 3, 4]
) -> Dict:
    """Analyze potential chatter conditions.

    Args:
        frequencies: Natural frequencies from modal analysis
        spindle_rpm_range: (min_rpm, max_rpm) for spindle
        flute_counts: List of endmill flute counts to check

    Returns:
        Dict with chatter analysis results
    """
    min_rpm, max_rpm = spindle_rpm_range
    warnings = []
    safe_zones = []

    for n_flutes in flute_counts:
        for freq in frequencies[:3]:  # Check first 3 modes
            # Tooth passing frequency = (RPM * n_flutes) / 60
            # Chatter occurs when TPF ≈ natural frequency

            # Critical RPM where TPF = natural frequency
            critical_rpm = (freq * 60) / n_flutes

            if min_rpm <= critical_rpm <= max_rpm:
                warnings.append({
                    'frequency': freq,
                    'flutes': n_flutes,
                    'critical_rpm': critical_rpm,
                    'warning': f"Mode at {freq:.1f} Hz resonates with {n_flutes}-flute at {critical_rpm:.0f} RPM"
                })

    # Find safe RPM ranges (away from any critical RPM)
    critical_rpms = sorted([w['critical_rpm'] for w in warnings])

    if not critical_rpms:
        safe_zones.append((min_rpm, max_rpm))
    else:
        # Add margins around critical RPMs
        margin = 500  # RPM margin

        prev = min_rpm
        for rpm in critical_rpms:
            if rpm - margin > prev:
                safe_zones.append((prev, rpm - margin))
            prev = rpm + margin

        if prev < max_rpm:
            safe_zones.append((prev, max_rpm))

    return {
        'warnings': warnings,
        'safe_zones': safe_zones,
        'spindle_range': spindle_rpm_range,
        'modes_analyzed': len(frequencies),
    }


if __name__ == "__main__":
    from ..config import OUTPUT_DIR

    msh_path = OUTPUT_DIR / "x_gantry_hybrid.msh"
    xdmf_path = OUTPUT_DIR / "x_gantry_hybrid.xdmf"

    mesh_file = xdmf_path if xdmf_path.exists() else msh_path

    if mesh_file.exists():
        results = solve_modal(mesh_file, n_modes=6, boundary_type='fixed')
        save_mode_shapes(results, OUTPUT_DIR / "modes")

        chatter = chatter_analysis(results['frequencies'])
        print(f"\nChatter warnings: {len(chatter['warnings'])}")
        for w in chatter['warnings']:
            print(f"  {w['warning']}")
    else:
        print(f"No mesh found at {mesh_file}")
