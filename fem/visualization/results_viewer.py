"""Visualization of FEM results using PyVista."""

import numpy as np
from pathlib import Path
from typing import List, Optional

import pyvista as pv


def plot_displacement(vtk_path: Path, show_deformed: bool = True,
                      scale_factor: float = 100.0, save_path: Optional[Path] = None):
    """Plot displacement field on mesh.

    Args:
        vtk_path: Path to results VTK file
        show_deformed: Show deformed mesh
        scale_factor: Deformation scale factor
        save_path: Optional path to save screenshot
    """
    mesh = pv.read(str(vtk_path))

    plotter = pv.Plotter(off_screen=save_path is not None)

    if show_deformed and 'displacement' in mesh.point_data:
        # Apply deformation
        displacement = mesh.point_data['displacement']
        mesh.points = mesh.points + displacement * scale_factor

    # Plot displacement magnitude
    if 'displacement_magnitude' in mesh.point_data:
        scalars = 'displacement_magnitude'
        title = 'Displacement Magnitude (mm)'
    elif 'displacement' in mesh.point_data:
        disp = mesh.point_data['displacement']
        mesh.point_data['disp_mag'] = np.linalg.norm(disp, axis=1)
        scalars = 'disp_mag'
        title = 'Displacement Magnitude (mm)'
    else:
        scalars = None
        title = None

    plotter.add_mesh(
        mesh,
        scalars=scalars,
        cmap='viridis',
        scalar_bar_args={'title': title} if title else None,
        show_edges=False
    )

    plotter.add_axes()
    plotter.add_title(f"Static Analysis - Scale: {scale_factor}x")

    if save_path:
        plotter.screenshot(str(save_path))
        print(f"Screenshot saved: {save_path}")
    else:
        plotter.show()

    return plotter


def plot_mode_shape(vtk_path: Path, scale_factor: float = 1.0,
                    save_path: Optional[Path] = None):
    """Plot a single mode shape.

    Args:
        vtk_path: Path to mode shape VTK file
        scale_factor: Additional scale factor
        save_path: Optional path to save screenshot
    """
    mesh = pv.read(str(vtk_path))

    plotter = pv.Plotter(off_screen=save_path is not None)

    # Apply mode shape deformation
    if 'mode_shape' in mesh.point_data:
        mode = mesh.point_data['mode_shape']
        mesh.points = mesh.points + mode * scale_factor

    # Color by mode magnitude
    if 'mode_magnitude' in mesh.point_data:
        scalars = 'mode_magnitude'
    else:
        scalars = None

    plotter.add_mesh(
        mesh,
        scalars=scalars,
        cmap='coolwarm',
        show_edges=False
    )

    plotter.add_axes()

    # Extract frequency from filename
    filename = Path(vtk_path).stem
    plotter.add_title(filename.replace('_', ' ').title())

    if save_path:
        plotter.screenshot(str(save_path))
    else:
        plotter.show()

    return plotter


def plot_mode_gallery(mode_dir: Path, n_modes: int = 6,
                      save_path: Optional[Path] = None):
    """Plot gallery of mode shapes.

    Args:
        mode_dir: Directory containing mode VTK files
        n_modes: Number of modes to display
        save_path: Optional path to save screenshot
    """
    mode_files = sorted(Path(mode_dir).glob('mode_*.vtk'))[:n_modes]

    if not mode_files:
        print(f"No mode files found in {mode_dir}")
        return None

    # Calculate grid layout
    n_cols = min(3, len(mode_files))
    n_rows = (len(mode_files) + n_cols - 1) // n_cols

    plotter = pv.Plotter(
        shape=(n_rows, n_cols),
        off_screen=save_path is not None
    )

    for i, vtk_path in enumerate(mode_files):
        row = i // n_cols
        col = i % n_cols

        plotter.subplot(row, col)

        mesh = pv.read(str(vtk_path))

        # Apply mode shape deformation
        if 'mode_shape' in mesh.point_data:
            mode = mesh.point_data['mode_shape']
            mesh.points = mesh.points + mode

        # Color by magnitude
        scalars = 'mode_magnitude' if 'mode_magnitude' in mesh.point_data else None

        plotter.add_mesh(mesh, scalars=scalars, cmap='coolwarm', show_edges=False)

        # Extract frequency from filename
        filename = vtk_path.stem
        plotter.add_title(filename.replace('_', ' ').replace('mode ', 'Mode '), font_size=8)

    plotter.link_views()

    if save_path:
        plotter.screenshot(str(save_path))
        print(f"Gallery saved: {save_path}")
    else:
        plotter.show()

    return plotter


def create_mode_animation(vtk_path: Path, output_gif: Path,
                          n_frames: int = 30, scale: float = 1.0):
    """Create animated GIF of mode shape oscillation.

    Args:
        vtk_path: Path to mode shape VTK file
        output_gif: Output GIF path
        n_frames: Number of frames in animation
        scale: Base scale for mode shape
    """
    mesh = pv.read(str(vtk_path))
    original_points = mesh.points.copy()

    if 'mode_shape' not in mesh.point_data:
        print("No mode_shape data found")
        return

    mode = mesh.point_data['mode_shape']

    plotter = pv.Plotter(off_screen=True)
    plotter.open_gif(str(output_gif))

    for i in range(n_frames):
        # Sinusoidal oscillation
        phase = 2 * np.pi * i / n_frames
        amplitude = np.sin(phase) * scale

        mesh.points = original_points + mode * amplitude

        plotter.clear()
        plotter.add_mesh(
            mesh,
            scalars='mode_magnitude',
            cmap='coolwarm',
            clim=[0, np.max(mesh.point_data['mode_magnitude'])],
            show_edges=False
        )
        plotter.add_axes()
        plotter.write_frame()

    plotter.close()
    print(f"Animation saved: {output_gif}")


def plot_gantry_overview(
    static_vtk: Path,
    mode_dir: Optional[Path] = None,
    scale_factor: float = 100.0,
    save_path: Optional[Path] = None
):
    """Plot overview of gantry analysis results.

    Shows deflection and first few modes in a single figure.

    Args:
        static_vtk: Path to static analysis results
        mode_dir: Directory with mode shape files
        scale_factor: Displacement scale factor
        save_path: Optional path to save screenshot
    """
    # Determine layout
    n_cols = 3
    n_rows = 1
    if mode_dir and mode_dir.exists():
        mode_files = sorted(mode_dir.glob('mode_*.vtk'))[:3]
        if mode_files:
            n_rows = 2

    plotter = pv.Plotter(
        shape=(n_rows, n_cols),
        off_screen=save_path is not None
    )

    # Row 1: Static results (3 views)
    mesh = pv.read(str(static_vtk))

    views = [
        ('Isometric View', (1, 1, 1)),
        ('Side View (Y)', (0, 1, 0)),
        ('Top View (Z)', (0, 0, 1)),
    ]

    for col, (title, view_up) in enumerate(views):
        plotter.subplot(0, col)

        mesh_copy = mesh.copy()
        if 'displacement' in mesh_copy.point_data:
            disp = mesh_copy.point_data['displacement']
            mesh_copy.points = mesh_copy.points + disp * scale_factor

        scalars = 'displacement_magnitude' if 'displacement_magnitude' in mesh_copy.point_data else None
        plotter.add_mesh(mesh_copy, scalars=scalars, cmap='viridis', show_edges=False)
        plotter.add_title(title, font_size=8)
        plotter.add_axes()

    # Row 2: Mode shapes
    if n_rows > 1 and mode_dir:
        mode_files = sorted(mode_dir.glob('mode_*.vtk'))[:3]
        for col, mode_file in enumerate(mode_files):
            plotter.subplot(1, col)

            mode_mesh = pv.read(str(mode_file))
            if 'mode_shape' in mode_mesh.point_data:
                mode = mode_mesh.point_data['mode_shape']
                mode_mesh.points = mode_mesh.points + mode

            scalars = 'mode_magnitude' if 'mode_magnitude' in mode_mesh.point_data else None
            plotter.add_mesh(mode_mesh, scalars=scalars, cmap='coolwarm', show_edges=False)
            plotter.add_title(mode_file.stem.replace('_', ' ').title(), font_size=8)

    if save_path:
        plotter.screenshot(str(save_path))
        print(f"Overview saved: {save_path}")
    else:
        plotter.show()

    return plotter


if __name__ == "__main__":
    from ..config import OUTPUT_DIR

    # Test visualization
    results_vtk = OUTPUT_DIR / "x_gantry_static_results.vtk"
    if results_vtk.exists():
        plot_displacement(results_vtk, save_path=OUTPUT_DIR / "displacement.png")

    mode_dir = OUTPUT_DIR / "modes"
    if mode_dir.exists():
        plot_mode_gallery(mode_dir, save_path=OUTPUT_DIR / "mode_gallery.png")
