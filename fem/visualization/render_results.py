#!/usr/bin/env python3
"""Generate FEM result visualizations using PyVista.

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


def load_fem_result(h5_path: Path, array_name: str = "displacement") -> pv.UnstructuredGrid:
    """Load FEM result from HDF5 file and convert to PyVista mesh."""
    with h5py.File(h5_path, 'r') as f:
        # Read mesh topology and geometry
        points = f['/Mesh/mesh/geometry'][:]
        cells_data = f['/Mesh/mesh/topology'][:]

        # Read function data
        if 'Function' in f:
            func_group = f'/Function/{array_name}/0'
            if func_group in f:
                array_data = f[func_group][:]
            else:
                # Try direct path
                array_data = None
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


def load_mode_result(h5_path: Path, mode_num: int) -> pv.UnstructuredGrid:
    """Load mode shape from HDF5 file."""
    with h5py.File(h5_path, 'r') as f:
        # Read mesh topology and geometry
        points = f['/Mesh/mesh/geometry'][:]
        cells_data = f['/Mesh/mesh/topology'][:]

        # Read mode shape data
        array_name = f'mode_{mode_num}'
        func_path = f'/Function/{array_name}/0'
        array_data = f[func_path][:]

    # Create cell array for tetrahedra
    n_cells = cells_data.shape[0]
    cells = np.hstack([np.full((n_cells, 1), 4, dtype=np.int64), cells_data]).ravel()
    cell_types = np.full(n_cells, pv.CellType.TETRA, dtype=np.uint8)

    # Create PyVista mesh
    grid = pv.UnstructuredGrid(cells, cell_types, points)
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


def render_mesh_only():
    """Render undeformed mesh for reference."""
    print("Rendering undeformed mesh...")

    h5_path = RESULTS_DIR / "ttc450_static_gravity.h5"
    mesh = load_fem_result(h5_path, "displacement")

    # Create plotter
    plotter = pv.Plotter(off_screen=True, window_size=[1920, 1080])
    plotter.background_color = 'white'

    # Add mesh with edges
    plotter.add_mesh(mesh, color='lightsteelblue', show_edges=True,
                     edge_color='darkgray', opacity=1.0)

    # Add title
    plotter.add_text("TTC450 Gantry - FEM Mesh\n15,810 nodes, 54,426 elements",
                     position='upper_left', font_size=14, color='black')

    setup_isometric_camera(plotter, mesh)

    output_path = IMAGES_DIR / "mesh_undeformed.png"
    plotter.screenshot(str(output_path))
    print(f"  Saved: {output_path}")
    plotter.close()


def render_static_result():
    """Render static analysis with warped displacement."""
    print("Rendering static analysis...")

    h5_path = RESULTS_DIR / "ttc450_static_gravity.h5"
    mesh = load_fem_result(h5_path, "displacement")

    # Get displacement array
    disp = mesh.point_data["displacement"]

    # Calculate magnitude for coloring
    disp_mag = np.linalg.norm(disp, axis=1)
    mesh.point_data["displacement_mag"] = disp_mag

    # Apply warp - scale so max displacement = 100mm
    # Actual max is ~0.001mm, so scale factor = 100000
    scale_factor = 100.0 / max(disp_mag.max(), 1e-10)
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
    max_disp = disp_mag.max()
    plotter.add_text(f"Static Analysis - Gravity Load\n"
                     f"Max displacement: {max_disp:.4f} mm (scaled {scale_factor:.0f}x)",
                     position='upper_left', font_size=14, color='black')

    setup_isometric_camera(plotter, warped)

    output_path = IMAGES_DIR / "static_gravity_warped.png"
    plotter.screenshot(str(output_path))
    print(f"  Saved: {output_path}")
    plotter.close()


def render_mode_shape(mode_num: int, frequency: float):
    """Render modal analysis result with warped mode shape."""
    print(f"Rendering mode {mode_num} ({frequency:.1f} Hz)...")

    h5_path = RESULTS_DIR / "ttc450_modes" / f"mode_{mode_num:02d}_{frequency:.1f}Hz.h5"
    array_name = f"mode_{mode_num}"

    mesh = load_mode_result(h5_path, mode_num)

    # Get mode shape
    mode = mesh.point_data[array_name]

    # Calculate magnitude for coloring
    mode_mag = np.linalg.norm(mode, axis=1)
    mesh.point_data["mode_mag"] = mode_mag

    # Apply warp - scale for 100mm max displacement
    scale_factor = 100.0 / max(mode_mag.max(), 1e-10)
    warped = mesh.warp_by_vector(array_name, factor=scale_factor)

    # Create plotter
    plotter = pv.Plotter(off_screen=True, window_size=[1920, 1080])
    plotter.background_color = 'white'

    # Add warped mesh colored by mode magnitude
    plotter.add_mesh(warped, scalars="mode_mag",
                     cmap='coolwarm', show_edges=False,
                     scalar_bar_args={
                         'title': 'Displacement',
                         'title_font_size': 16,
                         'label_font_size': 14,
                         'position_x': 0.85,
                         'position_y': 0.3,
                         'width': 0.1,
                         'height': 0.4,
                         'color': 'black',
                     })

    # Add title
    plotter.add_text(f"Mode {mode_num} - {frequency:.1f} Hz\n"
                     f"(Normalized mode shape, scaled for visibility)",
                     position='upper_left', font_size=14, color='black')

    setup_isometric_camera(plotter, warped)

    output_path = IMAGES_DIR / f"mode_{mode_num:02d}_{frequency:.1f}Hz.png"
    plotter.screenshot(str(output_path))
    print(f"  Saved: {output_path}")
    plotter.close()


def main():
    print("=" * 60)
    print("Generating FEM Result Visualizations")
    print("=" * 60)

    # Render undeformed mesh
    render_mesh_only()

    # Render static result
    render_static_result()

    # Render mode shapes
    modes = [
        (1, 114.7),
        (2, 120.8),
        (3, 165.0),
        (4, 359.5),
        (5, 373.9),
        (6, 443.2),
    ]

    for mode_num, freq in modes:
        render_mode_shape(mode_num, freq)

    print("\n" + "=" * 60)
    print("All visualizations complete!")
    print(f"Output directory: {IMAGES_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
