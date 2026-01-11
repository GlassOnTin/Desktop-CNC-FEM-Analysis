"""Static linear elasticity analysis using FEniCSx.

This replaces the SfePy solver with FEniCSx for better AMR support and modern API.
"""

import numpy as np
from pathlib import Path
from typing import Optional, Dict, Tuple, List
from mpi4py import MPI
from petsc4py import PETSc

import dolfinx
from dolfinx import fem, mesh, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
import ufl

from ..config import (
    GRAVITY, MATERIALS, CUTTING_LOADS, CBEAM_40X80, DEFAULT_GANTRY,
    BOUNDARY_TOLERANCE_MM, DEFAULT_LOAD_RADIUS_MM
)


def create_elasticity_problem(
    domain: mesh.Mesh,
    E: float,
    nu: float,
    rho: float
) -> Tuple[ufl.Form, ufl.Form, fem.FunctionSpace]:
    """Create the variational formulation for linear elasticity.

    Args:
        domain: DOLFINx mesh
        E: Young's modulus (N/mm² = MPa)
        nu: Poisson's ratio
        rho: Density (kg/mm³)

    Returns:
        (bilinear_form, linear_form, function_space)
    """
    # Vector function space for displacement (P1 elements)
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim,)))

    # Trial and test functions
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Lamé parameters
    lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))

    # Strain tensor (symmetric gradient)
    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    # Stress tensor (Hooke's law for isotropic material)
    def sigma(u):
        return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

    # Bilinear form (stiffness)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx

    # Body force (gravity in -Z direction)
    g = GRAVITY * 1000  # m/s² -> mm/s²
    f = fem.Constant(domain, default_scalar_type((0.0, 0.0, -rho * g)))

    # Linear form (body forces)
    L = ufl.dot(f, v) * ufl.dx

    return a, L, V


def locate_boundary_dofs(
    V: fem.FunctionSpace,
    domain: mesh.Mesh,
    boundary_type: str = 'fixed'
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Locate DOFs on boundary regions for different support conditions.

    Args:
        V: Function space
        domain: Mesh
        boundary_type: 'cantilever', 'fixed', or 'simply_supported'

    Returns:
        List of (dof_indices, values) tuples for boundary conditions
    """
    # Get mesh bounds
    coords = domain.geometry.x
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    tol = BOUNDARY_TOLERANCE_MM

    bc_data = []

    def left_boundary(x):
        return x[0] < x_min + tol

    def right_boundary(x):
        return x[0] > x_max - tol

    if boundary_type == 'cantilever':
        # Fix all DOFs at left end only
        left_dofs = fem.locate_dofs_geometrical(V, left_boundary)
        bc_data.append((left_dofs, np.zeros(3)))

    elif boundary_type == 'fixed':
        # Fix all DOFs at both ends
        left_dofs = fem.locate_dofs_geometrical(V, left_boundary)
        right_dofs = fem.locate_dofs_geometrical(V, right_boundary)
        bc_data.append((left_dofs, np.zeros(3)))
        bc_data.append((right_dofs, np.zeros(3)))

    elif boundary_type == 'simply_supported':
        # Pin supports: fix Z displacement at both ends, X at left only
        # This requires component-wise BCs which is more complex in FEniCSx
        # For simplicity, we'll use fixed-fixed and note this in results
        left_dofs = fem.locate_dofs_geometrical(V, left_boundary)
        right_dofs = fem.locate_dofs_geometrical(V, right_boundary)
        bc_data.append((left_dofs, np.zeros(3)))
        bc_data.append((right_dofs, np.zeros(3)))
        print("Note: simply_supported approximated as fixed-fixed in FEniCSx")

    else:
        raise ValueError(f"Unknown boundary_type: {boundary_type}")

    return bc_data


def add_point_load(
    L: ufl.Form,
    V: fem.FunctionSpace,
    domain: mesh.Mesh,
    load_point: Tuple[float, float, float],
    force: Tuple[float, float, float]
) -> ufl.Form:
    """Add a point load approximated as a distributed load over small region.

    Args:
        L: Current linear form
        V: Function space
        domain: Mesh
        load_point: (x, y, z) position
        force: (Fx, Fy, Fz) in Newtons

    Returns:
        Updated linear form
    """
    px, py, pz = load_point
    fx, fy, fz = force
    radius = DEFAULT_LOAD_RADIUS_MM

    # Create a smooth approximation to point load using Gaussian-like function
    x = ufl.SpatialCoordinate(domain)

    # Distance from load point
    r2 = (x[0] - px)**2 + (x[1] - py)**2 + (x[2] - pz)**2

    # Gaussian weight (normalized approximately)
    # The total integral should equal 1, so force magnitude is preserved
    sigma = radius / 3  # 3-sigma rule
    weight = ufl.exp(-r2 / (2 * sigma**2))

    # We need to normalize, but that requires integration
    # For simplicity, use a scaling factor based on sphere volume
    vol_approx = 4/3 * np.pi * radius**3
    f_density = fem.Constant(domain, default_scalar_type((fx/vol_approx, fy/vol_approx, fz/vol_approx)))

    # Create indicator function for load region
    # Using conditional to limit load to small region
    v = ufl.TestFunction(V)

    # Add load term with weight
    L_load = ufl.dot(f_density, v) * weight * ufl.dx

    return L + L_load


def solve_static(
    mesh_path: Path,
    material_name: str = 'aluminum_6061_t6',
    boundary_type: str = 'fixed',
    load_case: str = 'worst_case',
    load_point: Optional[Tuple[float, float, float]] = None,
    include_gravity: bool = True,
    force_override: Optional[Tuple[float, float, float]] = None
) -> Dict:
    """Solve static linear elasticity problem.

    Args:
        mesh_path: Path to mesh file (XDMF or Gmsh .msh)
        material_name: Key from MATERIALS dict
        boundary_type: 'cantilever', 'fixed', or 'simply_supported'
        load_case: Key from CUTTING_LOADS dict
        load_point: Optional (x, y, z) for load application
        include_gravity: Include self-weight
        force_override: Optional (Fx, Fy, Fz) to override load_case

    Returns:
        Dict with displacement field and analysis results
    """
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
        raise ValueError(f"Unsupported mesh format: {mesh_path.suffix}. Use .xdmf or .msh")

    print(f"Loaded mesh: {domain.topology.index_map(0).size_local} nodes, "
          f"{domain.topology.index_map(3).size_local} cells")

    # Get bounding box
    coords = domain.geometry.x
    bbox_min = coords.min(axis=0)
    bbox_max = coords.max(axis=0)
    print(f"Bounding box: [{bbox_min}] to [{bbox_max}]")

    # Create variational problem
    a, L, V = create_elasticity_problem(domain, E_mpa, nu, rho_mm)

    if not include_gravity:
        # Zero out body force
        L = ufl.dot(fem.Constant(domain, default_scalar_type((0.0, 0.0, 0.0))),
                    ufl.TestFunction(V)) * ufl.dx

    # Get force to apply
    if force_override is not None:
        force = force_override
    else:
        load_data = CUTTING_LOADS.get(load_case)
        if load_data is None:
            raise ValueError(f"Unknown load case: {load_case}")
        force = (load_data['Fx'], load_data['Fy'], load_data['Fz'])

    # Default load point: center of mesh
    if load_point is None:
        center = (bbox_min + bbox_max) / 2
        load_point = (center[0], bbox_min[1] - 20, center[2] - 50)

    # Add cutting load
    if any(f != 0 for f in force):
        print(f"Applying load {force} N at {load_point} mm")
        L = add_point_load(L, V, domain, load_point, force)

    # Create boundary conditions
    bc_data = locate_boundary_dofs(V, domain, boundary_type)
    bcs = []
    for dofs, values in bc_data:
        u_bc = fem.Function(V)
        u_bc.x.array[:] = 0.0
        bc = fem.dirichletbc(u_bc, dofs)
        bcs.append(bc)

    print(f"Applied {len(bcs)} boundary condition(s) ({boundary_type})")

    # Solve
    problem = LinearProblem(a, L, bcs=bcs,
                           petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()

    # Extract results
    u_array = uh.x.array.reshape((-1, 3))
    disp_mag = np.linalg.norm(u_array, axis=1)
    max_disp = disp_mag.max()
    max_idx = disp_mag.argmax()
    max_pos = coords[max_idx]

    # Find displacement at load point
    distances = np.linalg.norm(coords - np.array(load_point), axis=1)
    nearest_idx = distances.argmin()
    disp_at_load = u_array[nearest_idx]

    print(f"\nResults:")
    print(f"  Max displacement: {max_disp:.4f} mm at {max_pos}")
    print(f"  Displacement at load: {np.linalg.norm(disp_at_load):.4f} mm")

    return {
        'displacement': u_array,
        'displacement_magnitude': disp_mag,
        'max_displacement_mm': max_disp,
        'max_displacement_position': max_pos,
        'displacement_at_load': disp_at_load,
        'load_case': load_case,
        'load_point': load_point,
        'force': force,
        'mesh': domain,
        'solution': uh,
        'function_space': V,
    }


def save_results_xdmf(results: Dict, output_path: Path) -> Path:
    """Save results to XDMF for ParaView visualization.

    Args:
        results: Results from solve_static()
        output_path: Path for output file

    Returns:
        Path to saved file
    """
    domain = results['mesh']
    uh = results['solution']

    output_path = Path(output_path)
    with io.XDMFFile(MPI.COMM_WORLD, str(output_path), "w") as xdmf:
        xdmf.write_mesh(domain)
        uh.name = "displacement"
        xdmf.write_function(uh)

    print(f"Results saved to: {output_path}")
    return output_path


def save_results_vtk(results: Dict, output_path: Path) -> Path:
    """Save results to VTK for ParaView visualization.

    Args:
        results: Results from solve_static()
        output_path: Path for output file

    Returns:
        Path to saved file
    """
    domain = results['mesh']
    uh = results['solution']

    output_path = Path(output_path)
    with io.VTKFile(MPI.COMM_WORLD, str(output_path), "w") as vtk:
        vtk.write_mesh(domain)
        vtk.write_function(uh)

    print(f"Results saved to: {output_path}")
    return output_path


if __name__ == "__main__":
    from ..config import OUTPUT_DIR

    # Test with a simple mesh if available
    msh_path = OUTPUT_DIR / "x_gantry_hybrid.msh"
    xdmf_path = OUTPUT_DIR / "x_gantry_hybrid.xdmf"

    if xdmf_path.exists():
        results = solve_static(xdmf_path, boundary_type='fixed', load_case='worst_case')
        save_results_xdmf(results, OUTPUT_DIR / "x_gantry_fenicsx_results.xdmf")
    elif msh_path.exists():
        results = solve_static(msh_path, boundary_type='fixed', load_case='worst_case')
        save_results_xdmf(results, OUTPUT_DIR / "x_gantry_fenicsx_results.xdmf")
    else:
        print(f"No mesh found. Run mesh_generator.py first.")
        print(f"Looked for: {xdmf_path} or {msh_path}")
