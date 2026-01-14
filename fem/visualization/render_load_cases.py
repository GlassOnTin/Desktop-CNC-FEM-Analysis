#!/usr/bin/env python3
"""Generate load case visualizations using PyVista.

Creates isometric screenshots with WarpByVector deformation colored by displacement.
"""

import numpy as np
import pyvista as pv
from pathlib import Path
import h5py

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "fem" / "results"
IMAGES_DIR = PROJECT_ROOT / "docs" / "images"

# Ensure output directory exists
IMAGES_DIR.mkdir(parents=True, exist_ok=True)

# Use off-screen rendering
pv.OFF_SCREEN = True


def load_xdmf_result(h5_path: Path, array_name: str = "displacement") -> pv.UnstructuredGrid:
    """Load FEM result from HDF5 file and convert to PyVista mesh."""
    with h5py.File(h5_path, 'r') as f:
        # Read mesh topology and geometry
        points = f['/Mesh/mesh/geometry'][:]
        cells_data = f['/Mesh/mesh/topology'][:]

        # Read function data
        func_group = f'/Function/{array_name}/0'
        if func_group in f:
            array_data = f[func_group][:]
        else:
            array_data = None

    # Create cell array for tetrahedra
    n_cells = cells_data.shape[0]
    cells = np.hstack([np.full((n_cells, 1), 4, dtype=np.int64), cells_data]).ravel()
    cell_types = np.full(n_cells, pv.CellType.TETRA, dtype=np.uint8)

    # Create PyVista mesh
    grid = pv.UnstructuredGrid(cells, cell_types, points)

    # Add point data
    if array_data is not None:
        grid.point_data[array_name] = array_data

    return grid


def setup_isometric_camera(plotter, mesh):
    """Set up isometric camera view."""
    bounds = mesh.bounds
    cx = (bounds[0] + bounds[1]) / 2
    cy = (bounds[2] + bounds[3]) / 2
    cz = (bounds[4] + bounds[5]) / 2

    size = max(bounds[1] - bounds[0],
               bounds[3] - bounds[2],
               bounds[5] - bounds[4])

    # Isometric view angles
    dist = size * 2.5
    angle_z = np.radians(35.264)  # arctan(1/sqrt(2))
    angle_xy = np.radians(45)

    cam_x = cx + dist * np.cos(angle_xy) * np.cos(angle_z)
    cam_y = cy + dist * np.sin(angle_xy) * np.cos(angle_z)
    cam_z = cz + dist * np.sin(angle_z)

    plotter.camera_position = [
        (cam_x, cam_y, cam_z),  # Position
        (cx, cy, cz),           # Focal point
        (0, 0, 1)               # View up
    ]


def render_load_case(h5_path: Path, title: str, output_name: str, scale_mm: float = 50.0):
    """Render a load case result with warped displacement."""
    print(f"Rendering {output_name}...")

    mesh = load_xdmf_result(h5_path, "displacement")

    # Get displacement array
    disp = mesh.point_data["displacement"]

    # Calculate magnitude for coloring
    disp_mag = np.linalg.norm(disp, axis=1)
    mesh.point_data["displacement_mag"] = disp_mag

    max_disp = disp_mag.max()

    # Apply warp - scale so max displacement = scale_mm for visibility
    scale_factor = scale_mm / max(max_disp, 1e-10)
    warped = mesh.warp_by_vector("displacement", factor=scale_factor)

    # Create plotter
    plotter = pv.Plotter(off_screen=True, window_size=[1920, 1080])
    plotter.background_color = 'white'

    # Add warped mesh colored by displacement magnitude
    plotter.add_mesh(warped, scalars="displacement_mag",
                     cmap='coolwarm', show_edges=False,
                     scalar_bar_args={
                         'title': 'Displacement (mm)',
                         'title_font_size': 16,
                         'label_font_size': 14,
                         'position_x': 0.85,
                         'position_y': 0.3,
                         'width': 0.1,
                         'height': 0.4,
                         'color': 'black',
                     })

    # Add title
    plotter.add_text(f"{title}\n"
                     f"Max displacement: {max_disp*1000:.1f} \u03bcm (scaled {scale_factor:.0f}x)",
                     position='upper_left', font_size=14, color='black')

    setup_isometric_camera(plotter, warped)

    output_path = IMAGES_DIR / f"{output_name}.png"
    plotter.screenshot(str(output_path))
    print(f"  Saved: {output_path}")
    plotter.close()

    return max_disp


def main():
    print("=" * 60)
    print("Generating Load Case Visualizations")
    print("=" * 60)

    # First, we need to run the load cases to generate the XDMF files
    # Check if they exist, if not, run the analysis

    # Define load cases to visualize
    load_cases = [
        ("z_axis_weight", "Z-Axis Weight (5kg)", (0, 0, -50)),
        ("light_cut_y", "Light Cut Y (20N)", (0, 20, -50)),
        ("moderate_cut_y", "Moderate Cut Y (50N)", (0, 50, -50)),
        ("heavy_cut_x", "Heavy Cut X (100N)", (100, 0, -50)),
        ("heavy_cut_y", "Heavy Cut Y (100N)", (0, 100, -50)),
    ]

    # We need to generate XDMF files for each load case
    # Let's modify the load case runner to save all cases
    print("\nRunning load case analysis to generate result files...")

    import sys
    sys.path.insert(0, str(PROJECT_ROOT / "fem" / "analysis"))

    from run_load_cases import LOAD_CASES, run_load_case, find_closest_node
    from run_load_cases import MESH_FILE, E_MPA, NU, RHO_MM
    from mpi4py import MPI
    from dolfinx import fem, io
    import ufl

    # Load mesh once
    domain, _, _ = io.gmshio.read_from_msh(str(MESH_FILE), MPI.COMM_WORLD, 0)
    coords = domain.geometry.x
    bbox_min = coords.min(axis=0)

    V = fem.functionspace(domain, ("Lagrange", 1, (3,)))

    # Boundary condition
    z_min = bbox_min[2]
    tol = 1.0

    def base_boundary(x):
        return x[2] < z_min + tol

    bc_dofs = fem.locate_dofs_geometrical(V, base_boundary)
    u_zero = fem.Function(V)
    u_zero.x.array[:] = 0.0
    bc = fem.dirichletbc(u_zero, bc_dofs)

    tool_position = np.array([0.0, 20.0, 200.0])

    # Run each load case and save XDMF
    for lc in LOAD_CASES:
        result = run_load_case(domain, V, bc, lc, tool_position, coords)

        # Save to XDMF
        out_path = RESULTS_DIR / f"ttc450_{lc.name}.xdmf"
        with io.XDMFFile(MPI.COMM_WORLD, str(out_path), "w") as xdmf:
            xdmf.write_mesh(domain)
            result['uh'].name = "displacement"
            xdmf.write_function(result['uh'])
        print(f"  Saved: {out_path.name}")

    # Now render visualizations
    print("\n" + "=" * 60)
    print("Rendering visualizations...")
    print("=" * 60)

    cases_to_render = [
        ("ttc450_z_axis_weight", "Z-Axis + Spindle Weight (5kg down)"),
        ("ttc450_heavy_cut_x", "Heavy Cut - X Direction (100N)"),
        ("ttc450_heavy_cut_y", "Heavy Cut - Y Direction (100N)"),
    ]

    for filename, title in cases_to_render:
        h5_path = RESULTS_DIR / f"{filename}.h5"
        if h5_path.exists():
            render_load_case(h5_path, title, filename, scale_mm=50.0)
        else:
            print(f"  Warning: {h5_path} not found")

    print("\n" + "=" * 60)
    print("All visualizations complete!")
    print(f"Output directory: {IMAGES_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
