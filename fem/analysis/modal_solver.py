"""Modal analysis (eigenvalue problem) for CNC gantry using SfePy."""

import numpy as np
from pathlib import Path
from typing import Optional, List, Dict, Tuple

from scipy.sparse.linalg import eigsh
from sfepy.discrete import Problem, FieldVariable, Material, Integral, Equation, Equations
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.discrete.conditions import Conditions, EssentialBC
from sfepy.terms import Term
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson

from ..config import (
    MATERIALS, N_MODES_DEFAULT, SPINDLE_RPM_RANGE, FLUTE_COUNTS,
    CBEAM_40X80, analytical_beam_frequency, BOUNDARY_TOLERANCE_MM
)


def solve_modal(
    mesh_path: Path,
    material_name: str = 'aluminum_6061_t6',
    n_modes: int = N_MODES_DEFAULT,
    boundary_type: str = 'simply_supported',
    point_masses: Optional[List[Tuple[float, float, float, float, str]]] = None
) -> Dict:
    """Solve modal analysis eigenvalue problem for CNC gantry.

    Solves: K * phi = omega^2 * M * phi

    Args:
        mesh_path: Path to VTK mesh
        material_name: Key from MATERIALS dict
        n_modes: Number of modes to compute
        boundary_type: 'free', 'cantilever', 'simply_supported', or 'fixed'
        point_masses: List of (mass_kg, x_mm, y_mm, z_mm, description) tuples
                     for lumped masses at specific locations (e.g., spindle)

    Returns:
        Dict with frequencies, mode_shapes
    """
    material_props = MATERIALS[material_name]

    # Load mesh
    mesh = Mesh.from_file(str(mesh_path))

    # Create domain
    domain = FEDomain('domain', mesh)

    # Get bounding box
    bbox = domain.get_mesh_bounding_box()
    min_coords = bbox[0]
    max_coords = bbox[1]
    eps = 0.1

    # Regions
    omega_region = domain.create_region('Omega', 'all')

    x_min, x_max = min_coords[0], max_coords[0]
    z_min, z_max = min_coords[2], max_coords[2]
    bc_tol = BOUNDARY_TOLERANCE_MM  # mm tolerance for boundary selection

    # Build boundary conditions based on support type
    # SfePy doesn't support OR in region expressions, so we create multiple regions
    fixed_regions = []
    bc_dofs = {}  # Maps region name to DOF constraints

    if boundary_type == 'cantilever':
        # Fix all DOFs at left end only
        fixed_regions.append(('FixedLeft', f'vertices in (x < {x_min + bc_tol})'))
        bc_dofs['FixedLeft'] = {'u.all': 0.0}

    elif boundary_type == 'fixed':
        # Fix all DOFs at both ends (fixed-fixed beam)
        fixed_regions.append(('FixedLeft', f'vertices in (x < {x_min + bc_tol})'))
        fixed_regions.append(('FixedRight', f'vertices in (x > {x_max - bc_tol})'))
        bc_dofs['FixedLeft'] = {'u.all': 0.0}
        bc_dofs['FixedRight'] = {'u.all': 0.0}

    elif boundary_type == 'simply_supported':
        # Pin both ends: fix vertical (Z) displacement, allow rotation
        # This approximates simple supports for a beam along X
        fixed_regions.append(('PinLeft', f'vertices in (x < {x_min + bc_tol})'))
        fixed_regions.append(('PinRight', f'vertices in (x > {x_max - bc_tol})'))
        # Fix Z displacement at both ends, X at left to prevent rigid body motion
        bc_dofs['PinLeft'] = {'u.0': 0.0, 'u.2': 0.0}  # X and Z fixed
        bc_dofs['PinRight'] = {'u.2': 0.0}  # Only Z fixed

    # 'free' case: no fixed_regions added

    # Create region objects
    created_regions = {}
    for region_name, region_expr in fixed_regions:
        try:
            created_regions[region_name] = domain.create_region(region_name, region_expr, 'facet')
        except Exception:
            created_regions[region_name] = domain.create_region(region_name, region_expr, 'vertex')

    # Field and variables
    field = Field.from_args('displacement', np.float64, 'vector',
                            omega_region, approx_order=1)

    u = FieldVariable('u', 'unknown', field)
    v = FieldVariable('v', 'test', field, primary_var_name='u')

    # Material - convert to mm-based units for consistency with mesh
    E_mpa = material_props['E'] / 1e6  # Pa -> MPa
    nu = material_props['nu']
    rho_mm = material_props['rho'] / 1e9  # kg/m³ -> kg/mm³

    D = stiffness_from_youngpoisson(dim=3, young=E_mpa, poisson=nu)
    m = Material('m', D=D, rho=rho_mm)

    integral = Integral('i', order=2)

    # Stiffness term
    t_stiff = Term.new('dw_lin_elastic(m.D, v, u)',
                       integral, omega_region, m=m, v=v, u=u)

    # Mass term
    t_mass = Term.new('dw_volume_dot(m.rho, v, u)',
                      integral, omega_region, m=m, v=v, u=u)

    # Create problem with stiffness equation
    eq_stiff = Equation('stiff', t_stiff)
    eqs = Equations([eq_stiff])

    pb = Problem('modal', equations=eqs)

    # Apply boundary conditions
    if created_regions:
        ebcs = []
        for region_name, region_obj in created_regions.items():
            dof_spec = bc_dofs[region_name]
            bc = EssentialBC(f'bc_{region_name}', region_obj, dof_spec)
            ebcs.append(bc)
        pb.time_update(ebcs=Conditions(ebcs))
    else:
        pb.time_update()

    pb.update_materials()

    # Create integrals container for evaluate()
    from sfepy.discrete import Integrals
    integrals = Integrals([integral])

    # Use evaluate() to get stiffness matrix
    K_sparse = pb.evaluate(
        'dw_lin_elastic.i.Omega(m.D, v, u)',
        mode='weak', dw_mode='matrix',
        copy_materials=False,
        integrals=integrals
    )

    # Use evaluate() to get mass matrix
    M_sparse = pb.evaluate(
        'dw_volume_dot.i.Omega(m.rho, v, u)',
        mode='weak', dw_mode='matrix',
        copy_materials=False,
        integrals=integrals
    )

    # Apply point masses if provided (e.g., spindle motor)
    if point_masses:
        from scipy.sparse import lil_matrix
        M_lil = lil_matrix(M_sparse)
        mesh_coords = mesh.coors

        variables = pb.get_variables()
        u_var_local = variables['u']
        eq_map_arr = u_var_local.eq_map.eq

        for mass_kg, x_mm, y_mm, z_mm, desc in point_masses:
            target = np.array([x_mm, y_mm, z_mm])
            distances = np.linalg.norm(mesh_coords - target, axis=1)
            nearest_node = np.argmin(distances)
            min_dist = distances[nearest_node]

            mass_value = mass_kg

            for dof_comp in range(3):
                global_dof = nearest_node * 3 + dof_comp
                if global_dof < len(eq_map_arr):
                    eq_num = eq_map_arr[global_dof]
                    if eq_num >= 0:
                        M_lil[eq_num, eq_num] += mass_value

            print(f"  Point mass '{desc}': {mass_kg*1000:.0f}g at node {nearest_node} "
                  f"(dist={min_dist:.1f}mm from target)")

        M_sparse = M_lil.tocsr()

    # Solve eigenvalue problem using shift-invert
    try:
        eigenvalues, eigenvectors = eigsh(
            K_sparse, k=n_modes, M=M_sparse,
            sigma=1.0, which='LM'
        )
    except Exception as e:
        print(f"Shift-invert failed, trying standard mode: {e}")
        eigenvalues, eigenvectors = eigsh(
            K_sparse, k=n_modes, M=M_sparse,
            which='SM'
        )

    # Convert eigenvalues to frequencies
    omega = np.sqrt(np.abs(eigenvalues))
    frequencies_hz = omega / (2 * np.pi)

    # Sort by frequency
    idx = np.argsort(frequencies_hz)
    frequencies_hz = frequencies_hz[idx]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Filter out near-zero eigenvalues (numerical artifacts)
    valid_mask = np.abs(eigenvalues) > 1.0
    frequencies_hz = frequencies_hz[valid_mask]
    eigenvalues = eigenvalues[valid_mask]
    eigenvectors = eigenvectors[:, valid_mask]
    omega = omega[idx][valid_mask]

    # Get equation-to-DOF mapping
    variables = pb.get_variables()
    u_var = variables['u']
    eq_map = u_var.eq_map.eq

    # Mass-normalize eigenvectors
    eigenvectors_normalized = np.zeros_like(eigenvectors)
    for i in range(eigenvectors.shape[1]):
        phi = eigenvectors[:, i]
        modal_mass = phi @ M_sparse @ phi
        if modal_mass > 0:
            eigenvectors_normalized[:, i] = phi / np.sqrt(modal_mass)
        else:
            eigenvectors_normalized[:, i] = phi

    return {
        'frequencies_hz': frequencies_hz,
        'eigenvalues': eigenvalues,
        'mode_shapes': eigenvectors_normalized,
        'omega_rad_s': omega,
        'mesh': mesh,
        'field': field,
        'n_dofs': K_sparse.shape[0],
        'eq_map': eq_map,
        'boundary_type': boundary_type,
    }


def check_chatter_risk(
    frequencies_hz: np.ndarray,
    spindle_rpm_range: Tuple[float, float] = SPINDLE_RPM_RANGE,
    flute_counts: List[int] = FLUTE_COUNTS,
    margin: float = 0.2
) -> Dict:
    """Check for potential chatter risk based on modal frequencies.

    Chatter occurs when tooth passing frequency (TPF) coincides with
    a natural frequency of the structure. TPF = (RPM / 60) * flutes

    Args:
        frequencies_hz: Array of natural frequencies
        spindle_rpm_range: (min_rpm, max_rpm) operating range
        flute_counts: List of flute counts to check (common endmills)
        margin: Frequency ratio margin (0.2 = within 20% is risky)

    Returns:
        Dict with chatter analysis results
    """
    min_rpm, max_rpm = spindle_rpm_range
    results = {
        'safe': True,
        'warnings': [],
        'analysis': []
    }

    for n_flutes in flute_counts:
        # Calculate tooth passing frequency range
        tpf_min = (min_rpm / 60) * n_flutes
        tpf_max = (max_rpm / 60) * n_flutes

        for i, freq in enumerate(frequencies_hz):
            # Check if any TPF is within margin of this mode
            ratio_min = freq / tpf_max if tpf_max > 0 else float('inf')
            ratio_max = freq / tpf_min if tpf_min > 0 else float('inf')

            # If mode frequency is within TPF range (or close to it)
            if (1 - margin) <= ratio_min <= (1 + margin) or \
               (1 - margin) <= ratio_max <= (1 + margin):
                warning = (
                    f"Mode {i+1} ({freq:.1f} Hz) may resonate with "
                    f"{n_flutes}-flute endmill at {int(freq*60/n_flutes)} RPM"
                )
                results['warnings'].append(warning)
                results['safe'] = False

            results['analysis'].append({
                'mode': i + 1,
                'frequency_hz': freq,
                'flutes': n_flutes,
                'critical_rpm': int(freq * 60 / n_flutes),
            })

    # Minimum mode 1 frequency recommendation
    mode_1_freq = frequencies_hz[0] if len(frequencies_hz) > 0 else 0
    if mode_1_freq < 50:
        results['warnings'].append(
            f"Mode 1 frequency ({mode_1_freq:.1f} Hz) is below 50 Hz - "
            f"consider stiffening the structure"
        )
        results['safe'] = False
    elif mode_1_freq < 80:
        results['warnings'].append(
            f"Mode 1 frequency ({mode_1_freq:.1f} Hz) is marginal - "
            f"aim for >80 Hz for reliable operation"
        )

    results['mode_1_frequency'] = mode_1_freq
    results['stiffness_margin'] = mode_1_freq / 50.0  # Ratio to minimum target

    return results


def save_mode_shapes_vtk(
    modal_results: Dict,
    output_dir: Path,
    scale_factor: float = 1.0
) -> List[Path]:
    """Save mode shapes as VTK files for visualization.

    Args:
        modal_results: Results from solve_modal()
        output_dir: Directory for output files
        scale_factor: Displacement scale for visualization

    Returns:
        List of output file paths
    """
    import meshio
    import shutil

    output_dir = Path(output_dir)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mesh = modal_results['mesh']
    mode_shapes = modal_results['mode_shapes']
    frequencies = modal_results['frequencies_hz']
    eq_map = modal_results.get('eq_map')
    n_modes = len(frequencies)

    points = mesh.coors.copy()
    cells = []
    for group in mesh.descs:
        if group == '3_4':
            conn = mesh.get_conn(group)
            cells.append(meshio.CellBlock("tetra", conn))

    n_nodes = len(points)
    n_total_dofs = n_nodes * 3

    output_files = []

    for i in range(n_modes):
        mode = mode_shapes[:, i]
        n_free_dofs = len(mode)

        if n_free_dofs == n_total_dofs:
            mode_3d = mode.reshape((n_nodes, 3))
        elif eq_map is not None:
            mode_3d = np.zeros((n_nodes, 3))
            for node in range(n_nodes):
                for comp in range(3):
                    dof_idx = node * 3 + comp
                    if dof_idx < len(eq_map):
                        eq_num = eq_map[dof_idx]
                        if eq_num >= 0 and eq_num < n_free_dofs:
                            mode_3d[node, comp] = mode[eq_num]
        else:
            mode_3d = np.zeros((n_nodes, 3))

        mode_3d = mode_3d * scale_factor
        mode_mag = np.linalg.norm(mode_3d, axis=1)

        meshio_mesh = meshio.Mesh(
            points=points,
            cells=cells,
            point_data={
                'mode_shape': mode_3d,
                'mode_magnitude': mode_mag,
            }
        )

        output_path = output_dir / f'mode_{i+1:02d}_f{frequencies[i]:.1f}Hz.vtk'
        meshio.write(str(output_path), meshio_mesh)
        output_files.append(output_path)

        print(f"Mode {i+1}: {frequencies[i]:.2f} Hz saved to {output_path.name}")

    return output_files


def compare_with_analytical(
    fem_frequencies: np.ndarray,
    beam_length: float = 600.0,
    beam_material: str = 'aluminum_6061_t6',
    support: str = 'simply_supported'
) -> Dict:
    """Compare FEM results with analytical beam theory.

    Args:
        fem_frequencies: Array of FEM frequencies in Hz
        beam_length: Beam span in mm
        beam_material: Material key
        support: Support type ('simply_supported' or 'fixed_fixed')

    Returns:
        Dict with comparison results
    """
    mat = MATERIALS[beam_material]
    cb = CBEAM_40X80

    comparisons = []
    for i, fem_freq in enumerate(fem_frequencies[:3], start=1):
        analytical_freq = analytical_beam_frequency(
            L=beam_length,
            E=mat['E'],
            I=cb['Iy'],  # Strong axis
            rho=mat['rho'],
            A=cb['area'],
            mode=i,
            support=support
        )

        error_pct = (fem_freq - analytical_freq) / analytical_freq * 100

        comparisons.append({
            'mode': i,
            'fem_hz': fem_freq,
            'analytical_hz': analytical_freq,
            'error_pct': error_pct,
        })

        print(f"Mode {i}: FEM={fem_freq:.1f} Hz, Analytical={analytical_freq:.1f} Hz, "
              f"Error={error_pct:+.1f}%")

    return {'comparisons': comparisons}


if __name__ == "__main__":
    from ..config import OUTPUT_DIR

    vtk_path = OUTPUT_DIR / "x_gantry_hybrid.vtk"
    if vtk_path.exists():
        print("Solving modal analysis...")
        results = solve_modal(vtk_path, n_modes=6, boundary_type='simply_supported')

        print("\nNatural Frequencies:")
        for i, freq in enumerate(results['frequencies_hz']):
            print(f"  Mode {i+1}: {freq:.2f} Hz")

        # Check chatter risk
        print("\nChatter Analysis:")
        chatter = check_chatter_risk(results['frequencies_hz'])
        if chatter['safe']:
            print("  No chatter risk detected in normal operating range")
        else:
            for warning in chatter['warnings']:
                print(f"  WARNING: {warning}")

        # Compare with analytical
        print("\nAnalytical Comparison:")
        compare_with_analytical(results['frequencies_hz'])

        # Save mode shapes
        mode_dir = OUTPUT_DIR / "modes"
        save_mode_shapes_vtk(results, mode_dir)
    else:
        print(f"Mesh not found: {vtk_path}")
        print("Run export and meshing first")
