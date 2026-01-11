#!/usr/bin/env python3
"""Main entry point for CNC gantry FEM analysis."""

import argparse
from pathlib import Path
from typing import Optional

from .config import (
    OUTPUT_DIR, MATERIALS, CUTTING_LOADS, DEFAULT_GANTRY,
    analytical_beam_deflection, analytical_beam_frequency,
    CBEAM_40X80, get_x_gantry_refinement, SPINDLE_MASS_KG
)


def run_full_analysis(
    component: str = 'x_gantry',
    include_rails: bool = True,
    load_case: str = 'worst_case',
    mesh_size: float = 2.0,
    n_modes: int = 6,
    output_dir: Optional[Path] = None,
    solver: str = 'fenicsx'
):
    """Run complete FEM analysis workflow.

    Args:
        component: Component to analyze ('x_gantry', 'y_extension', 'x_gantry_system')
        include_rails: Include HGR20 rails in X-gantry model
        load_case: Load case key from CUTTING_LOADS
        mesh_size: Target mesh element size in mm
        n_modes: Number of modal frequencies to compute
        output_dir: Output directory for results
        solver: 'fenicsx' (default) or 'sfepy'
    """
    from .geometry.export_stl import export_x_gantry, export_y_extension, export_x_gantry_with_extensions, export_ttc450_pro
    from .geometry.mesh_generator import generate_mesh_from_stl, convert_vtk_to_xdmf

    if output_dir is None:
        output_dir = OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"CNC Gantry FEM Analysis (solver: {solver})")
    print("=" * 60)

    # Step 1: Export STL
    print(f"\n[1/4] Exporting {component} to STL...")
    if component == 'x_gantry':
        stl_path = export_x_gantry(output_dir, include_rails=include_rails)
    elif component == 'y_extension':
        stl_path = export_y_extension(output_dir)
    elif component == 'x_gantry_system':
        stl_path = export_x_gantry_with_extensions(output_dir)
    elif component == 'ttc450_pro':
        stl_path = export_ttc450_pro(output_dir, simplified=True)
    else:
        raise ValueError(f"Unknown component: {component}")

    # Step 2: Generate mesh
    print(f"\n[2/4] Generating tetrahedral mesh (size={mesh_size}mm)...")
    refinement = get_x_gantry_refinement() if 'x_gantry' in component else None

    # Generate both VTK (for legacy SfePy) and XDMF (for FEniCSx)
    mesh_result = generate_mesh_from_stl(
        stl_path, output_dir,
        mesh_size=mesh_size,
        refinement_regions=refinement,
        output_formats=['vtk', 'xdmf', 'msh']
    )
    print(f"  Mesh: {mesh_result['n_nodes']} nodes, {mesh_result['n_elements']} elements")

    # Select mesh format based on solver
    if solver == 'fenicsx':
        mesh_path = mesh_result.get('xdmf_path') or mesh_result.get('msh_path')
    else:
        mesh_path = mesh_result['vtk_path']

    # Step 3: Static analysis
    print(f"\n[3/4] Running static analysis (load case: {load_case})...")

    import numpy as np

    if solver == 'fenicsx':
        from .analysis.fenicsx_static import solve_static, save_results_xdmf, save_results_vtk
        static_results = solve_static(
            mesh_path,
            material_name='aluminum_6061_t6',
            boundary_type='fixed',
            load_case=load_case,
            include_gravity=True
        )
        # Save in both formats for compatibility
        static_xdmf = output_dir / f"{component}_static_results.xdmf"
        save_results_xdmf(static_results, static_xdmf)
        static_vtk = output_dir / f"{component}_static_results.vtk"
        try:
            save_results_vtk(static_results, static_vtk)
        except Exception as e:
            print(f"  Note: VTK export skipped ({e})")
    else:
        from .analysis.static_solver import setup_gantry_problem, solve_with_cutting_load, save_results_vtk
        setup = setup_gantry_problem(mesh_path, material_name='aluminum_6061_t6', boundary_type='fixed')
        static_results = solve_with_cutting_load(setup, load_case=load_case, include_gravity=True)
        static_vtk = output_dir / f"{component}_static_results.vtk"
        save_results_vtk(static_results, static_vtk)

    print(f"\n  Static Analysis Results:")
    print(f"    Max displacement: {static_results['max_displacement_mm']:.4f} mm")
    print(f"    At position: ({static_results['max_displacement_position'][0]:.1f}, "
          f"{static_results['max_displacement_position'][1]:.1f}, "
          f"{static_results['max_displacement_position'][2]:.1f}) mm")

    disp_at_load = static_results['displacement_at_load']
    print(f"    Displacement at load point: {np.linalg.norm(disp_at_load):.4f} mm")

    # Compare with analytical
    load_data = CUTTING_LOADS[load_case]
    total_force = np.sqrt(load_data['Fx']**2 + load_data['Fy']**2 + load_data['Fz']**2)
    mat = MATERIALS['aluminum_6061_t6']
    cb = CBEAM_40X80

    analytical_defl = analytical_beam_deflection(
        P=total_force,
        L=600.0,
        E=mat['E'],
        I=cb['Iy'],
        support='fixed_fixed'
    )
    print(f"\n  Analytical beam deflection (C-Beam only, fixed-fixed): {analytical_defl:.4f} mm")
    if static_results['max_displacement_mm'] > 0:
        print(f"  Stiffness improvement from hybrid beam: "
              f"{analytical_defl / static_results['max_displacement_mm']:.1f}x")

    # Step 4: Modal analysis
    # Note: Using SfePy for modal analysis as FEniCSx/SLEPc eigenvalue solver
    # needs additional configuration for this problem type
    print(f"\n[4/4] Running modal analysis ({n_modes} modes)...")
    print("  (Using SfePy for modal - FEniCSx eigenvalue solver WIP)")

    # Always use SfePy modal solver (works reliably)
    from .analysis.modal_solver import solve_modal, check_chatter_risk, save_mode_shapes_vtk, compare_with_analytical
    # Include spindle as point mass
    spindle_position = (
        0.0,
        -CBEAM_40X80['width'] / 2,
        DEFAULT_GANTRY.spindle_z_offset
    )
    point_masses = [(SPINDLE_MASS_KG, *spindle_position, 'Spindle/Router')]

    # Use VTK mesh for SfePy modal analysis
    vtk_path = mesh_result['vtk_path']
    modal_results = solve_modal(
        vtk_path, n_modes=n_modes, boundary_type='fixed',
        point_masses=point_masses
    )
    frequencies = modal_results['frequencies_hz']

    print(f"\n  Natural Frequencies:")
    for i, freq in enumerate(frequencies[:6]):
        print(f"    Mode {i+1}: {freq:.2f} Hz")

    chatter = check_chatter_risk(frequencies)
    print(f"\n  Chatter Analysis:")
    print(f"    Mode 1 frequency: {chatter['mode_1_frequency']:.1f} Hz")
    print(f"    Stiffness margin: {chatter['stiffness_margin']:.1f}x minimum target (50 Hz)")

    compare_with_analytical(frequencies)
    mode_dir = output_dir / "modes"
    save_mode_shapes_vtk(modal_results, mode_dir)

    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    print(f"\nOutput files in: {output_dir}")
    print(f"  - {component}.stl (geometry)")
    print(f"  - {component}.vtk / .xdmf (mesh)")
    if solver == 'fenicsx':
        print(f"  - {component}_static_results.xdmf (deflection)")
    else:
        print(f"  - {component}_static_results.vtk (deflection)")
    print(f"  - modes/ (modal analysis)")

    return {
        'static': static_results,
        'modal': modal_results,
        'mesh': mesh_result,
    }


def run_analytical_validation():
    """Run analytical validation of beam theory vs expected values."""
    print("\n" + "=" * 60)
    print("Analytical Validation")
    print("=" * 60)

    mat = MATERIALS['aluminum_6061_t6']
    cb = CBEAM_40X80

    print(f"\nC-Beam 40x80 Properties:")
    print(f"  Iy = {cb['Iy']:.2e} mm^4")
    print(f"  Area = {cb['area']:.0f} mm^2")
    print(f"  E = {mat['E']/1e9:.0f} GPa")

    print(f"\nDeflection under 200N center load (600mm span):")
    for support in ['simply_supported', 'fixed_fixed']:
        defl = analytical_beam_deflection(
            P=200.0, L=600.0, E=mat['E'], I=cb['Iy'], support=support
        )
        print(f"  {support}: {defl:.4f} mm")

    print(f"\nFirst 3 natural frequencies (600mm span):")
    for mode in [1, 2, 3]:
        freq = analytical_beam_frequency(
            L=600.0, E=mat['E'], I=cb['Iy'],
            rho=mat['rho'], A=cb['area'],
            mode=mode, support='simply_supported'
        )
        print(f"  Mode {mode}: {freq:.1f} Hz")


def main():
    parser = argparse.ArgumentParser(
        description='CNC Gantry FEM Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python -m fem.run_analysis --component x_gantry --load-case worst_case
  python -m fem.run_analysis --component x_gantry --solver sfepy  # Use legacy solver
  python -m fem.run_analysis --component x_gantry --no-rails  # Baseline comparison
  python -m fem.run_analysis --validate  # Analytical validation only
        '''
    )

    parser.add_argument(
        '--component', '-c',
        choices=['x_gantry', 'y_extension', 'x_gantry_system', 'ttc450_pro'],
        default='x_gantry',
        help='Component to analyze (default: x_gantry)'
    )

    parser.add_argument(
        '--no-rails',
        action='store_true',
        help='Exclude HGR20 rails from X-gantry (baseline analysis)'
    )

    parser.add_argument(
        '--load-case', '-l',
        choices=list(CUTTING_LOADS.keys()),
        default='worst_case',
        help='Cutting load case (default: worst_case)'
    )

    parser.add_argument(
        '--mesh-size', '-m',
        type=float,
        default=5.0,
        help='Target mesh element size in mm (default: 5.0)'
    )

    parser.add_argument(
        '--n-modes', '-n',
        type=int,
        default=6,
        help='Number of modal frequencies to compute (default: 6)'
    )

    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        default=None,
        help='Output directory (default: fem/results)'
    )

    parser.add_argument(
        '--solver', '-s',
        choices=['fenicsx', 'sfepy'],
        default='fenicsx',
        help='FEM solver to use (default: fenicsx)'
    )

    parser.add_argument(
        '--validate',
        action='store_true',
        help='Run analytical validation only'
    )

    args = parser.parse_args()

    if args.validate:
        run_analytical_validation()
    else:
        run_full_analysis(
            component=args.component,
            include_rails=not args.no_rails,
            load_case=args.load_case,
            mesh_size=args.mesh_size,
            n_modes=args.n_modes,
            output_dir=args.output_dir,
            solver=args.solver
        )


if __name__ == "__main__":
    main()
