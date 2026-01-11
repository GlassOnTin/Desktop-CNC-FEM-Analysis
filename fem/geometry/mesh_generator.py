"""Generate tetrahedral mesh from STL using gmsh."""

import math
from pathlib import Path
from typing import List, Tuple, Optional

import gmsh
import meshio
import numpy as np

from ..config import DEFAULT_MESH_SIZE_MIN, DEFAULT_MESH_SIZE_MAX

# Mesh quality thresholds
MIN_QUALITY_THRESHOLD = 0.1  # Warn if any element below this
TARGET_QUALITY = 0.3  # Target minimum quality for good mesh


def compute_tet_quality(points: np.ndarray, cells: np.ndarray) -> np.ndarray:
    """Compute quality metric for tetrahedral elements.

    Uses the aspect ratio metric: 3 * r_in / r_out
    where r_in is inscribed sphere radius and r_out is circumscribed.
    Quality = 1.0 for regular tetrahedron, approaches 0 for degenerate.

    Args:
        points: (N, 3) array of vertex coordinates
        cells: (M, 4) array of vertex indices for each tetrahedron

    Returns:
        (M,) array of quality values in [0, 1]
    """
    n_cells = len(cells)
    qualities = np.zeros(n_cells)

    for i, cell in enumerate(cells):
        # Get vertices
        v0, v1, v2, v3 = points[cell]

        # Edge vectors
        e01 = v1 - v0
        e02 = v2 - v0
        e03 = v3 - v0
        e12 = v2 - v1
        e13 = v3 - v1

        # Volume (6x actual volume)
        vol6 = abs(np.dot(e01, np.cross(e02, e03)))

        if vol6 < 1e-15:
            qualities[i] = 0.0
            continue

        # Edge lengths
        edges = [e01, e02, e03, e12, e13, v3 - v2]
        edge_lengths = [np.linalg.norm(e) for e in edges]
        sum_edges = sum(edge_lengths)

        # Face areas (via cross products)
        face_normals = [
            np.cross(e01, e02),  # face 0-1-2
            np.cross(e01, e03),  # face 0-1-3
            np.cross(e02, e03),  # face 0-2-3
            np.cross(e12, e13),  # face 1-2-3
        ]
        face_areas = [0.5 * np.linalg.norm(n) for n in face_normals]
        sum_areas = sum(face_areas)

        # Inscribed sphere radius: r_in = 3V / A (where A is total surface area)
        volume = vol6 / 6.0
        r_in = 3.0 * volume / sum_areas if sum_areas > 0 else 0.0

        # Circumscribed sphere approximation using edge length
        # For a regular tet: r_out = edge * sqrt(6) / 4
        avg_edge = sum_edges / 6.0
        r_out_approx = avg_edge * 0.612  # sqrt(6)/4 ≈ 0.612

        # Quality metric: normalized ratio
        if r_out_approx > 0:
            qualities[i] = min(1.0, 3.0 * r_in / r_out_approx)
        else:
            qualities[i] = 0.0

    return qualities


def validate_mesh_quality(
    vtk_path: Path,
    min_threshold: float = MIN_QUALITY_THRESHOLD,
    target: float = TARGET_QUALITY
) -> dict:
    """Validate mesh quality and report statistics.

    Args:
        vtk_path: Path to VTK mesh file
        min_threshold: Warn if any element below this quality
        target: Target minimum quality

    Returns:
        Dict with quality statistics and pass/fail status
    """
    mesh = meshio.read(str(vtk_path))

    # Find tetrahedral cells
    tet_cells = None
    for cell_block in mesh.cells:
        if cell_block.type == "tetra":
            tet_cells = cell_block.data
            break

    if tet_cells is None:
        return {
            'valid': False,
            'error': 'No tetrahedral elements found',
            'n_elements': 0,
        }

    # Compute quality
    qualities = compute_tet_quality(mesh.points, tet_cells)

    # Statistics
    min_q = np.min(qualities)
    max_q = np.max(qualities)
    mean_q = np.mean(qualities)
    n_below_threshold = np.sum(qualities < min_threshold)
    n_below_target = np.sum(qualities < target)
    n_inverted = np.sum(qualities <= 0)

    result = {
        'valid': True,
        'n_elements': len(tet_cells),
        'min_quality': float(min_q),
        'max_quality': float(max_q),
        'mean_quality': float(mean_q),
        'n_inverted': int(n_inverted),
        'n_below_threshold': int(n_below_threshold),
        'n_below_target': int(n_below_target),
        'pass': n_inverted == 0 and min_q >= min_threshold,
    }

    # Report
    status = "PASS" if result['pass'] else "WARNING"
    print(f"  Mesh quality: {status}")
    print(f"    Elements: {result['n_elements']}")
    print(f"    Quality: min={min_q:.3f}, mean={mean_q:.3f}, max={max_q:.3f}")

    if n_inverted > 0:
        print(f"    ERROR: {n_inverted} inverted elements (quality <= 0)")
    if n_below_threshold > 0:
        print(f"    WARNING: {n_below_threshold} elements below quality threshold ({min_threshold})")
    if n_below_target > 0 and n_below_target != n_below_threshold:
        print(f"    INFO: {n_below_target} elements below target quality ({target})")

    return result


def stl_to_tet_mesh(
    stl_path: Path,
    output_msh: Path,
    mesh_size_min: float = DEFAULT_MESH_SIZE_MIN,
    mesh_size_max: float = DEFAULT_MESH_SIZE_MAX,
    refinement_regions: Optional[List[Tuple[float, float, float, float, float]]] = None
) -> Path:
    """Convert closed STL surface to tetrahedral volume mesh.

    Args:
        stl_path: Path to input STL file
        output_msh: Path for output .msh file
        mesh_size_min: Minimum element size (mm)
        mesh_size_max: Maximum element size (mm)
        refinement_regions: List of (x, y, z, radius, size) for local refinement

    Returns:
        Path to generated mesh file
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)  # Suppress output

    try:
        # Merge STL file
        gmsh.merge(str(stl_path))

        # Classify surfaces to identify geometric features
        angle = 40 * math.pi / 180  # Feature detection angle
        includeBoundary = True
        forceParametrizablePatches = False
        curveAngle = 180 * math.pi / 180

        gmsh.model.mesh.classifySurfaces(
            angle, includeBoundary, forceParametrizablePatches, curveAngle
        )

        # Create geometry from discrete surfaces
        gmsh.model.mesh.createGeometry()

        # Get all surfaces and create volume
        surfaces = gmsh.model.getEntities(2)
        if len(surfaces) == 0:
            raise ValueError("No surfaces found in STL")

        surface_loop = gmsh.model.geo.addSurfaceLoop([e[1] for e in surfaces])
        gmsh.model.geo.addVolume([surface_loop])
        gmsh.model.geo.synchronize()

        # Set mesh size
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size_min)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size_max)

        # Optional: Add refinement near stress concentrations
        if refinement_regions:
            field_ids = []
            for i, (x, y, z, radius, size) in enumerate(refinement_regions):
                field_id = i + 1
                gmsh.model.mesh.field.add("Ball", field_id)
                gmsh.model.mesh.field.setNumber(field_id, "Radius", radius)
                gmsh.model.mesh.field.setNumber(field_id, "VIn", size)
                gmsh.model.mesh.field.setNumber(field_id, "VOut", mesh_size_max)
                gmsh.model.mesh.field.setNumber(field_id, "XCenter", x)
                gmsh.model.mesh.field.setNumber(field_id, "YCenter", y)
                gmsh.model.mesh.field.setNumber(field_id, "ZCenter", z)
                field_ids.append(field_id)

            # Combine fields using minimum
            min_field = len(refinement_regions) + 1
            gmsh.model.mesh.field.add("Min", min_field)
            gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", field_ids)
            gmsh.model.mesh.field.setAsBackgroundMesh(min_field)

        # Generate 3D tetrahedral mesh
        gmsh.model.mesh.generate(3)

        # Write mesh
        gmsh.write(str(output_msh))

        # Get mesh statistics
        nodes = gmsh.model.mesh.getNodes()
        n_nodes = len(nodes[0])
        elements = gmsh.model.mesh.getElements(3)  # 3D elements
        n_elements = sum(len(e) for e in elements[1]) if elements[1] else 0

        print(f"Mesh generated: {n_nodes} nodes, {n_elements} tetrahedra")

    finally:
        gmsh.finalize()

    return output_msh


def convert_msh_to_vtk(msh_path: Path, vtk_path: Path) -> Path:
    """Convert gmsh .msh to VTK format for SfePy.

    Args:
        msh_path: Path to gmsh mesh file
        vtk_path: Path for output VTK file

    Returns:
        Path to VTK file
    """
    mesh = meshio.read(str(msh_path))

    # Filter to only keep 3D elements (tetrahedra)
    cells_3d = []
    for cell_block in mesh.cells:
        if cell_block.type == "tetra":
            cells_3d.append(cell_block)

    if not cells_3d:
        raise ValueError("No tetrahedral elements found in mesh")

    mesh_3d = meshio.Mesh(
        points=mesh.points,
        cells=cells_3d
    )

    meshio.write(str(vtk_path), mesh_3d)
    return vtk_path


def convert_msh_to_xdmf(msh_path: Path, xdmf_path: Path) -> Path:
    """Convert gmsh .msh to XDMF format for FEniCSx.

    Args:
        msh_path: Path to gmsh mesh file
        xdmf_path: Path for output XDMF file

    Returns:
        Path to XDMF file
    """
    mesh = meshio.read(str(msh_path))

    # Filter to only keep 3D elements (tetrahedra)
    cells_3d = []
    for cell_block in mesh.cells:
        if cell_block.type == "tetra":
            cells_3d.append(cell_block)

    if not cells_3d:
        raise ValueError("No tetrahedral elements found in mesh")

    mesh_3d = meshio.Mesh(
        points=mesh.points,
        cells=cells_3d
    )

    meshio.write(str(xdmf_path), mesh_3d)
    print(f"XDMF mesh written to: {xdmf_path}")
    return xdmf_path


def convert_vtk_to_xdmf(vtk_path: Path, xdmf_path: Path) -> Path:
    """Convert VTK to XDMF format for FEniCSx.

    Args:
        vtk_path: Path to VTK mesh file
        xdmf_path: Path for output XDMF file

    Returns:
        Path to XDMF file
    """
    mesh = meshio.read(str(vtk_path))

    # Filter to only keep 3D elements (tetrahedra)
    cells_3d = []
    for cell_block in mesh.cells:
        if cell_block.type == "tetra":
            cells_3d.append(cell_block)

    if not cells_3d:
        raise ValueError("No tetrahedral elements found in mesh")

    mesh_3d = meshio.Mesh(
        points=mesh.points,
        cells=cells_3d
    )

    meshio.write(str(xdmf_path), mesh_3d)
    print(f"XDMF mesh written to: {xdmf_path}")
    return xdmf_path


def generate_mesh_from_stl(
    stl_path: Path,
    output_dir: Path,
    mesh_size: float = 2.0,
    refinement_regions: Optional[List] = None,
    output_formats: Optional[List[str]] = None
) -> dict:
    """Complete workflow: STL -> MSH -> VTK/XDMF.

    Args:
        stl_path: Path to input STL
        output_dir: Output directory
        mesh_size: Target mesh size in mm
        refinement_regions: Optional refinement zones
        output_formats: List of formats to output ('vtk', 'xdmf', 'msh')
                       Default: ['vtk', 'xdmf', 'msh']

    Returns:
        Dict with mesh paths and statistics
    """
    if output_formats is None:
        output_formats = ['vtk', 'xdmf', 'msh']

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = stl_path.stem
    msh_path = output_dir / f"{stem}.msh"
    vtk_path = output_dir / f"{stem}.vtk"
    xdmf_path = output_dir / f"{stem}.xdmf"

    # Generate tetrahedral mesh
    stl_to_tet_mesh(
        stl_path,
        msh_path,
        mesh_size_min=mesh_size * 0.5,
        mesh_size_max=mesh_size,
        refinement_regions=refinement_regions
    )

    result = {
        'msh_path': msh_path,
    }

    # Convert to VTK (for SfePy and visualization)
    if 'vtk' in output_formats:
        convert_msh_to_vtk(msh_path, vtk_path)
        result['vtk_path'] = vtk_path

    # Convert to XDMF (for FEniCSx)
    if 'xdmf' in output_formats:
        convert_msh_to_xdmf(msh_path, xdmf_path)
        result['xdmf_path'] = xdmf_path

    # Load mesh to get stats
    mesh = meshio.read(str(msh_path))
    n_nodes = len(mesh.points)
    n_elements = sum(len(c.data) for c in mesh.cells if c.type == 'tetra')

    result['n_nodes'] = n_nodes
    result['n_elements'] = n_elements

    # Validate mesh quality (use VTK if available, otherwise MSH)
    check_path = vtk_path if 'vtk' in output_formats else msh_path
    quality_result = validate_mesh_quality(check_path)
    result['quality'] = quality_result

    return result


if __name__ == "__main__":
    # Test with X-gantry
    from .export_stl import export_component
    from ..config import get_x_gantry_refinement, OUTPUT_DIR

    stl_path = export_component('x_gantry', OUTPUT_DIR)
    result = generate_mesh_from_stl(
        stl_path,
        OUTPUT_DIR,
        mesh_size=2.0,
        refinement_regions=get_x_gantry_refinement()
    )
    print(f"Generated mesh: {result}")
