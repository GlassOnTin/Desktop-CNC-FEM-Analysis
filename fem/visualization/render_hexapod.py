#!/usr/bin/env python3
"""Generate hexapod Stewart platform visualizations using PyVista.

Creates isometric screenshots with WarpByVector deformation colored by displacement.
"""

import numpy as np
import pyvista as pv
from pathlib import Path
import h5py
import json

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


def setup_isometric_camera(plotter, mesh, elevation=35, azimuth=45):
    """Set up isometric camera view."""
    bounds = mesh.bounds
    cx = (bounds[0] + bounds[1]) / 2
    cy = (bounds[2] + bounds[3]) / 2
    cz = (bounds[4] + bounds[5]) / 2

    size = max(bounds[1] - bounds[0],
               bounds[3] - bounds[2],
               bounds[5] - bounds[4])

    # Camera position
    dist = size * 2.0
    angle_z = np.radians(elevation)
    angle_xy = np.radians(azimuth)

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

    setup_isometric_camera(plotter, warped, elevation=25, azimuth=45)

    output_path = IMAGES_DIR / f"{output_name}.png"
    plotter.screenshot(str(output_path))
    print(f"  Saved: {output_path}")
    plotter.close()

    return max_disp


def render_undeformed_mesh(msh_path: Path, output_name: str):
    """Render the undeformed hexapod mesh."""
    print(f"Rendering undeformed mesh...")

    # Load mesh from gmsh format
    mesh = pv.read(str(msh_path))

    # Create plotter
    plotter = pv.Plotter(off_screen=True, window_size=[1920, 1080])
    plotter.background_color = 'white'

    # Add mesh with edges
    plotter.add_mesh(mesh, color='lightblue', show_edges=True, edge_color='darkblue',
                     opacity=0.9)

    # Add title
    n_points = mesh.n_points
    n_cells = mesh.n_cells
    plotter.add_text(f"Hexapod Stewart Platform\n"
                     f"Nodes: {n_points:,}  Elements: {n_cells:,}",
                     position='upper_left', font_size=14, color='black')

    setup_isometric_camera(plotter, mesh, elevation=25, azimuth=45)

    output_path = IMAGES_DIR / f"{output_name}.png"
    plotter.screenshot(str(output_path))
    print(f"  Saved: {output_path}")
    plotter.close()


def main():
    print("=" * 60)
    print("Generating Hexapod Visualizations")
    print("=" * 60)

    # First render undeformed mesh
    msh_path = RESULTS_DIR / "hexapod.msh"
    if msh_path.exists():
        render_undeformed_mesh(msh_path, "hexapod_mesh")
    else:
        print(f"  Warning: {msh_path} not found - run generate_hexapod.py first")

    # Render load cases
    cases_to_render = [
        ("hexapod_z_axis_weight", "Z-Axis + Spindle Weight (5kg down)"),
        ("hexapod_heavy_cut_x", "Heavy Cut - X Direction (100N)"),
        ("hexapod_heavy_cut_y", "Heavy Cut - Y Direction (100N)"),
        ("hexapod_heavy_cut_z", "Heavy Plunge Cut - Z Direction (150N)"),
    ]

    for filename, title in cases_to_render:
        h5_path = RESULTS_DIR / f"{filename}.h5"
        if h5_path.exists():
            render_load_case(h5_path, title, filename, scale_mm=50.0)
        else:
            print(f"  Warning: {h5_path} not found - run run_hexapod_analysis.py first")

    print("\n" + "=" * 60)
    print("All visualizations complete!")
    print(f"Output directory: {IMAGES_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
