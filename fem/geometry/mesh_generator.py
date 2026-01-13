"""Generate tetrahedral mesh from STL using gmsh."""

import math
from pathlib import Path
from typing import List, Tuple, Optional

import gmsh
import meshio
import numpy as np

from ..config import DEFAULT_MESH_SIZE_MIN, DEFAULT_MESH_SIZE_MAX, CBEAM_40X80

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
    refinement_regions: Optional[List[Tuple[float, float, float, float, float]]] = None,
    optimize: bool = True,
    use_hxt: bool = True,
    use_curvature: bool = True,
    classify_angle_deg: float = 40.0,
    classify_force: bool = False
) -> Path:
    """Convert closed STL surface to tetrahedral volume mesh.

    Args:
        stl_path: Path to input STL file
        output_msh: Path for output .msh file
        mesh_size_min: Minimum element size (mm)
        mesh_size_max: Maximum element size (mm)
        refinement_regions: List of (x, y, z, radius, size) for local refinement
        optimize: Enable mesh optimization passes
        use_hxt: Use HXT 3D meshing algorithm (robust for complex STL)
        use_curvature: Enable curvature-based sizing

    Returns:
        Path to generated mesh file
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)  # Suppress output

    def _try_set_option(name: str, value: float) -> None:
        try:
            gmsh.option.setNumber(name, value)
        except Exception:
            pass

    try:
        # STL cleanup options (helps with non-manifold or duplicated facets)
        _try_set_option("Mesh.StlRemoveDuplicateNodes", 1)
        _try_set_option("Mesh.StlRemoveDuplicateTriangles", 1)
        _try_set_option("Mesh.StlOneSurfacePerCell", 1)

        # Merge STL file
        gmsh.merge(str(stl_path))

        # Remove duplicate nodes/elements after merge
        try:
            gmsh.model.mesh.removeDuplicateNodes()
            gmsh.model.mesh.removeDuplicateElements()
        except Exception:
            pass

        # Classify surfaces to identify geometric features
        angle = classify_angle_deg * math.pi / 180  # Feature detection angle
        includeBoundary = True
        forceParametrizablePatches = classify_force
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
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1 if use_curvature else 0)
        gmsh.option.setNumber("Mesh.Smoothing", 10)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10 if use_hxt else 1)
        gmsh.option.setNumber("Mesh.Optimize", 1 if optimize else 0)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1 if optimize else 0)

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


def _add_cbeam_profile_surface(
    width: float = 40.0,
    height: float = 80.0,
    outer_wall: float = 3.0,
    inner_wall: float = 2.0,
    slot_depth: float = 6.0,
    slot_width: float = 11.0,
    profile: str = "cshape",
    slot_dir: int = -1,
) -> int:
    """Create a C-beam cross-section in the XZ plane.

    profile:
      - "solid": solid rectangle (fastest/most robust)
      - "hollow": rectangular tube (hole via inner rectangle)
      - "cshape": single-wire C-shaped polygon (robust, no booleans)
      - "cslot": C-slot + internal voids (most detailed; uses booleans)
    Returns the surface tag for the resulting profile.
    """
    def _rect_wire(xmin: float, xmax: float, zmin: float, zmax: float) -> int:
        pts = [
            (xmin, 0.0, zmin),
            (xmax, 0.0, zmin),
            (xmax, 0.0, zmax),
            (xmin, 0.0, zmax),
        ]
        p_tags = [gmsh.model.occ.addPoint(x, y, z) for x, y, z in pts]
        lines = [gmsh.model.occ.addLine(p_tags[i], p_tags[(i + 1) % 4]) for i in range(4)]
        return gmsh.model.occ.addWire(lines)

    def _rect_surface(xmin: float, xmax: float, zmin: float, zmax: float) -> int:
        wire = _rect_wire(xmin, xmax, zmin, zmax)
        return gmsh.model.occ.addPlaneSurface([wire])

    if profile == "solid":
        wire = _rect_wire(-width / 2, width / 2, -height / 2, height / 2)
        surface = gmsh.model.occ.addPlaneSurface([wire])
        gmsh.model.occ.synchronize()
        return surface

    if profile == "hollow":
        outer_wire = _rect_wire(-width / 2, width / 2, -height / 2, height / 2)
        inner_wire = _rect_wire(
            -width / 2 + outer_wall,
            width / 2 - outer_wall,
            -height / 2 + outer_wall,
            height / 2 - outer_wall,
        )
        surface = gmsh.model.occ.addPlaneSurface([outer_wire, inner_wire])
        gmsh.model.occ.synchronize()
        return surface

    if profile == "cshape":
        t = outer_wall
        s = slot_depth
        # C-shaped polygon (single wire, no holes/booleans)
        pts = [
            (-width / 2, -height / 2),
            (width / 2, -height / 2),
            (width / 2, height / 2),
            (-width / 2, height / 2),
            (-width / 2, height / 2 - t),
            (-width / 2 + s, height / 2 - t),
            (-width / 2 + s, -height / 2 + t),
            (-width / 2, -height / 2 + t),
        ]
        if slot_dir > 0:
            pts = [(-x, y) for x, y in pts]
        p_tags = [gmsh.model.occ.addPoint(x, 0.0, z) for x, z in pts]
        lines = []
        for i in range(len(p_tags)):
            p1 = p_tags[i]
            p2 = p_tags[(i + 1) % len(p_tags)]
            lines.append(gmsh.model.occ.addLine(p1, p2))
        wire = gmsh.model.occ.addWire(lines)
        surface = gmsh.model.occ.addPlaneSurface([wire])
        gmsh.model.occ.synchronize()
        return surface

    # Upper and lower internal voids
    void_width = width - 2 * outer_wall
    void_height = height / 2 - outer_wall - inner_wall / 2

    upper_center_z = height / 4
    lower_center_z = -height / 4

    upper = _rect_surface(
        -void_width / 2,
        void_width / 2,
        upper_center_z - void_height / 2,
        upper_center_z + void_height / 2,
    )
    lower = _rect_surface(
        -void_width / 2,
        void_width / 2,
        lower_center_z - void_height / 2,
        lower_center_z + void_height / 2,
    )

    # C-slot opening on -X face
    slot_height = height - 4 * outer_wall
    slot_center_x = (-width / 2 + slot_depth / 2) * (-1 if slot_dir > 0 else 1)
    slot = _rect_surface(
        slot_center_x - slot_depth / 2,
        slot_center_x + slot_depth / 2,
        -slot_height / 2,
        slot_height / 2,
    )

    # Subtract voids and slot from outer profile
    outer = _rect_surface(-width / 2, width / 2, -height / 2, height / 2)
    gmsh.model.occ.synchronize()
    cut = gmsh.model.occ.cut(
        [(2, outer)],
        [(2, upper), (2, lower), (2, slot)],
        removeObject=True,
        removeTool=True,
    )
    gmsh.model.occ.synchronize()
    if not cut or not cut[0]:
        raise RuntimeError("Failed to construct C-beam profile")

    return cut[0][0][1]


def _extrude_profile(
    surface_tag: int,
    dx: float,
    dy: float,
    dz: float,
    recombine: bool = True,
    layers: int = 24,
    center: bool = True,
) -> int:
    """Extrude a 2D surface into a 3D volume and return the volume tag."""
    gmsh.model.occ.synchronize()
    if center:
        gmsh.model.occ.translate([(2, surface_tag)], -dx / 2, -dy / 2, -dz / 2)
    entities = gmsh.model.occ.extrude(
        [(2, surface_tag)],
        dx, dy, dz,
        numElements=[layers],
        heights=[],
        recombine=recombine,
    )
    # The volume is the first 3D entity in the result list
    for dim, tag in entities:
        if dim == 3:
            return tag
    raise RuntimeError("Extrusion did not create a volume")


def generate_ttc450_hybrid_mesh(
    output_msh: Path,
    mesh_size_min: float = 2.5,
    mesh_size_max: float = 12.0,
    beam_length: float = 600.0,
    base_size: float = 600.0,
    base_ext: float = 20.0,
    riser_thk: float = 6.0,
    z_beam_length: float = 150.0,
    recombine: bool = False,
    union_all: bool = False,
    profile: str = "hollow",
    interface_depth: float = 20.0,
    interface_gap: float = 0.2,
    include_base: bool = False,
) -> Path:
    """Generate a hybrid mesh using swept C-beams and Delaunay interface blocks.

    - C-beams are swept from a simplified C-profile.
    - Interface blocks (flat-faced cuboids) are meshed with standard Delaunay.
    - Volumes touch at interfaces (no overlap); suitable for tied contact in solver.
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    try:
        # Geometry parameters
        width = CBEAM_40X80["width"]
        height = CBEAM_40X80["height"]
        outer_wall = CBEAM_40X80["outer_wall"]
        inner_wall = CBEAM_40X80["inner_wall"]
        slot_depth = CBEAM_40X80["slot_depth"]
        slot_width = CBEAM_40X80["slot_width"]

        y_beam_x = base_size / 2 - width / 2
        y_beam_z = base_ext + height / 2
        x_beam_z = y_beam_z + 100

        # Y-axis rails: length along Y, slot faces ±X (outward)
        profile_tag = _add_cbeam_profile_surface(
            width=width,
            height=height,
            outer_wall=outer_wall,
            inner_wall=inner_wall,
            slot_depth=slot_depth,
            slot_width=slot_width,
            profile=profile,
        )
        y_beam = _extrude_profile(profile_tag, 0, beam_length, 0)
        gmsh.model.occ.translate([(3, y_beam)], -y_beam_x, 0, y_beam_z)

        profile_tag = _add_cbeam_profile_surface(
            width=width,
            height=height,
            outer_wall=outer_wall,
            inner_wall=inner_wall,
            slot_depth=slot_depth,
            slot_width=slot_width,
            profile=profile,
        )
        y_beam_2 = _extrude_profile(profile_tag, 0, beam_length, 0)
        gmsh.model.occ.rotate([(3, y_beam_2)], 0, 0, 0, 0, 0, 1, math.pi)
        gmsh.model.occ.translate([(3, y_beam_2)], y_beam_x, 0, y_beam_z)

        # X-gantry beam: length along X, slot faces -Y
        profile_tag = _add_cbeam_profile_surface(
            width=width,
            height=height,
            outer_wall=outer_wall,
            inner_wall=inner_wall,
            slot_depth=slot_depth,
            slot_width=slot_width,
            profile=profile,
        )
        x_beam = _extrude_profile(profile_tag, 0, beam_length, 0)
        gmsh.model.occ.rotate([(3, x_beam)], 0, 0, 0, 0, 0, 1, math.pi / 2)
        gmsh.model.occ.translate([(3, x_beam)], 0, 0, x_beam_z)

        # Z-axis beam: simplified solid block for stable meshing
        z_beam = gmsh.model.occ.addBox(
            -width / 2,
            -(width + interface_depth + 2 * interface_gap) - width / 2,
            x_beam_z - height / 2,
            width,
            width,
            z_beam_length,
        )
        z_center_z = (x_beam_z - height / 2) + z_beam_length / 2
        gmsh.model.occ.translate([(3, z_beam)], 0, 0, 0)

        base_vols = []
        if include_base:
            # Base frame (2020 square)
            base_half = base_size / 2
            base_z = base_ext / 2
            base_parts = []
            # Rails along X (front/back)
            for sy in (-1, 1):
                box = gmsh.model.occ.addBox(
                    -base_half + base_ext, sy * (base_half - base_ext / 2) - base_ext / 2, base_z - base_ext / 2,
                    base_size - 2 * base_ext, base_ext, base_ext
                )
                base_parts.append(box)
            # Rails along Y (left/right)
            for sx in (-1, 1):
                box = gmsh.model.occ.addBox(
                    sx * (base_half - base_ext / 2) - base_ext / 2, -base_half + base_ext, base_z - base_ext / 2,
                    base_ext, base_size - 2 * base_ext, base_ext
                )
                base_parts.append(box)
            # Corner blocks
            for sx in (-1, 1):
                for sy in (-1, 1):
                    box = gmsh.model.occ.addBox(
                        sx * (base_half - base_ext / 2) - base_ext / 2,
                        sy * (base_half - base_ext / 2) - base_ext / 2,
                        base_z - base_ext / 2,
                        base_ext, base_ext, base_ext
                    )
                    base_parts.append(box)

            # Fuse base parts into a single volume to avoid overlapping facets
            if base_parts:
                base_fuse, _ = gmsh.model.occ.fuse(
                    [(3, base_parts[0])],
                    [(3, tag) for tag in base_parts[1:]],
                    removeObject=True,
                    removeTool=True,
                )
                for dim, tag in base_fuse:
                    if dim == 3:
                        base_vols.append(tag)

        # Interface blocks for Delaunay bridging (touching, no overlap)
        interface_blocks = []
        for sx in (-1, 1):
            # Riser plates outside Y-beams: touch beam ends at x=±beam_length/2
            block = gmsh.model.occ.addBox(
                sx * (beam_length / 2 + interface_gap + riser_thk / 2) - riser_thk / 2,
                -width / 2,
                x_beam_z - height / 2,
                riser_thk,
                width,
                height,
            )
            interface_blocks.append(block)

        # Block between X-beam back face and Z-beam front face
        z_interface = gmsh.model.occ.addBox(
            -width / 2,
            -width / 2 - interface_gap - interface_depth,
            x_beam_z - height / 2,
            width,
            interface_depth,
            height,
        )
        interface_blocks.append(z_interface)

        all_vols = [
            (3, y_beam),
            (3, y_beam_2),
            (3, x_beam),
            (3, z_beam),
        ]
        all_vols.extend([(3, v) for v in base_vols])
        all_vols.extend([(3, v) for v in interface_blocks])
        if union_all:
            fused, _ = gmsh.model.occ.fuse(all_vols[:1], all_vols[1:], removeObject=True, removeTool=True)
            # Keep fused volumes (single in most cases)
            for dim, tag in fused:
                if dim == 3:
                    all_vols = [(3, tag)]
                    break
        gmsh.model.occ.synchronize()

        # Mesh settings: structured for swept beams, Delaunay for interface blocks
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size_min)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size_max)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay
        if recombine:
            gmsh.option.setNumber("Mesh.RecombineAll", 1)
            for dim, tag in gmsh.model.getEntities(2):
                gmsh.model.mesh.setRecombine(dim, tag)

        gmsh.model.mesh.generate(3)
        gmsh.write(str(output_msh))

        return output_msh
    finally:
        gmsh.finalize()


def _mesh_single_volume(
    create_volume_fn,
    output_msh: Path,
    mesh_size_min: float,
    mesh_size_max: float,
    tag_surfaces_fn=None,
) -> Path:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        vol = create_volume_fn()
        gmsh.model.occ.synchronize()
        if tag_surfaces_fn is not None:
            tag_surfaces_fn(vol)
        gmsh.option.setNumber("Mesh.SaveAll", 1)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size_min)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size_max)
        gmsh.model.mesh.generate(3)
        gmsh.write(str(output_msh))
        return output_msh
    finally:
        gmsh.finalize()


def generate_ttc450_hybrid_parts(
    output_dir: Path,
    mesh_size_min: float = 8.0,
    mesh_size_max: float = 25.0,
    beam_length: float = 600.0,
    base_size: float = 600.0,
    base_ext: float = 20.0,
    riser_thk: float = 6.0,
    z_beam_length: float = 150.0,
    profile: str = "solid",
    interface_depth: float = 20.0,
    interface_gap: float = 0.2,
    include_base: bool = False,
    beam_shape: str = "profile",
    couple_x_beam_ends: bool = True,
    couple_riser_inner_faces: bool = True,
    surface_tol: float = 1e-3,
) -> dict:
    """Mesh each major part separately for later tie constraints in the solver."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    width = CBEAM_40X80["width"]
    height = CBEAM_40X80["height"]
    outer_wall = CBEAM_40X80["outer_wall"]
    inner_wall = CBEAM_40X80["inner_wall"]
    slot_depth = CBEAM_40X80["slot_depth"]
    slot_width = CBEAM_40X80["slot_width"]

    y_beam_x = base_size / 2 - width / 2
    y_beam_z = base_ext + height / 2
    x_beam_z = y_beam_z + 100

    def _make_y_beam(x_sign: int):
        if beam_shape == "box":
            return gmsh.model.occ.addBox(
                -width / 2 + x_sign * y_beam_x,
                -beam_length / 2,
                y_beam_z - height / 2,
                width,
                beam_length,
                height,
            )
        profile_tag = _add_cbeam_profile_surface(
            width=width,
            height=height,
            outer_wall=outer_wall,
            inner_wall=inner_wall,
            slot_depth=slot_depth,
            slot_width=slot_width,
            profile=profile,
            slot_dir=1 if x_sign > 0 else -1,
        )
        gmsh.model.occ.translate([(2, profile_tag)], x_sign * y_beam_x, 0, y_beam_z)
        return _extrude_profile(profile_tag, 0, beam_length, 0)

    def _make_x_beam():
        if beam_shape == "box":
            return gmsh.model.occ.addBox(
                -beam_length / 2,
                -width / 2,
                x_beam_z - height / 2,
                beam_length,
                width,
                height,
            )
        profile_tag = _add_cbeam_profile_surface(
            width=width,
            height=height,
            outer_wall=outer_wall,
            inner_wall=inner_wall,
            slot_depth=slot_depth,
            slot_width=slot_width,
            profile=profile,
            slot_dir=-1,
        )
        # Rotate profile so slot faces -Y, then extrude along X
        gmsh.model.occ.rotate([(2, profile_tag)], 0, 0, 0, 0, 0, 1, math.pi / 2)
        gmsh.model.occ.translate([(2, profile_tag)], 0, 0, x_beam_z)
        return _extrude_profile(profile_tag, beam_length, 0, 0)

    def _make_z_block():
        return gmsh.model.occ.addBox(
            -width / 2,
            -(width + interface_depth + 2 * interface_gap) - width / 2,
            x_beam_z - height / 2,
            width,
            width,
            z_beam_length,
        )

    def _make_riser(x_sign: int):
        return gmsh.model.occ.addBox(
            x_sign * (beam_length / 2 + interface_gap + riser_thk / 2) - riser_thk / 2,
            -width / 2,
            x_beam_z - height / 2,
            riser_thk,
            width,
            height,
        )

    def _make_z_interface():
        return gmsh.model.occ.addBox(
            -width / 2,
            -width / 2 - interface_gap - interface_depth,
            x_beam_z - height / 2,
            width,
            interface_depth,
            height,
        )

    parts = {
        "y_beam_left": lambda: _make_y_beam(-1),
        "y_beam_right": lambda: _make_y_beam(1),
        "x_beam": _make_x_beam,
        "z_block": _make_z_block,
        "riser_left": lambda: _make_riser(-1),
        "riser_right": lambda: _make_riser(1),
        "z_interface": _make_z_interface,
    }

    if include_base:
        def _make_base():
            base_half = base_size / 2
            base_z = base_ext / 2
            base_parts = []
            for sy in (-1, 1):
                base_parts.append(
                    gmsh.model.occ.addBox(
                        -base_half + base_ext,
                        sy * (base_half - base_ext / 2) - base_ext / 2,
                        base_z - base_ext / 2,
                        base_size - 2 * base_ext,
                        base_ext,
                        base_ext,
                    )
                )
            for sx in (-1, 1):
                base_parts.append(
                    gmsh.model.occ.addBox(
                        sx * (base_half - base_ext / 2) - base_ext / 2,
                        -base_half + base_ext,
                        base_z - base_ext / 2,
                        base_ext,
                        base_size - 2 * base_ext,
                        base_ext,
                    )
                )
            for sx in (-1, 1):
                for sy in (-1, 1):
                    base_parts.append(
                        gmsh.model.occ.addBox(
                            sx * (base_half - base_ext / 2) - base_ext / 2,
                            sy * (base_half - base_ext / 2) - base_ext / 2,
                            base_z - base_ext / 2,
                            base_ext,
                            base_ext,
                            base_ext,
                        )
                    )
            gmsh.model.occ.synchronize()
            base_fuse, _ = gmsh.model.occ.fuse(
                [(3, base_parts[0])],
                [(3, tag) for tag in base_parts[1:]],
                removeObject=True,
                removeTool=True,
            )
            for dim, tag in base_fuse:
                if dim == 3:
                    return tag
            raise RuntimeError("Base fuse failed")

        parts["base_frame"] = _make_base

    def _tag_faces_on_plane(vol_tag: int, axis: int, value: float, name: str):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, vol_tag)
        if axis == 0:
            ents = gmsh.model.getEntitiesInBoundingBox(
                value - surface_tol, ymin - surface_tol, zmin - surface_tol,
                value + surface_tol, ymax + surface_tol, zmax + surface_tol, 2
            )
        elif axis == 1:
            ents = gmsh.model.getEntitiesInBoundingBox(
                xmin - surface_tol, value - surface_tol, zmin - surface_tol,
                xmax + surface_tol, value + surface_tol, zmax + surface_tol, 2
            )
        else:
            ents = gmsh.model.getEntitiesInBoundingBox(
                xmin - surface_tol, ymin - surface_tol, value - surface_tol,
                xmax + surface_tol, ymax + surface_tol, value + surface_tol, 2
            )
        faces = [tag for dim, tag in ents if dim == 2]
        if faces:
            phys = gmsh.model.addPhysicalGroup(2, faces)
            gmsh.model.setPhysicalName(2, phys, name)

    outputs = {}
    for name, make_fn in parts.items():
        path = output_dir / f"{name}.msh"
        def _tagger(vol_tag: int, part_name=name):
            # Tag volume
            phys = gmsh.model.addPhysicalGroup(3, [vol_tag])
            gmsh.model.setPhysicalName(3, phys, part_name)

            # Tag coupling faces
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, vol_tag)
            if part_name == "x_beam" and couple_x_beam_ends:
                _tag_faces_on_plane(vol_tag, 0, xmin, "x_beam_end_left")
                _tag_faces_on_plane(vol_tag, 0, xmax, "x_beam_end_right")
            if part_name == "riser_left" and couple_riser_inner_faces:
                _tag_faces_on_plane(vol_tag, 0, xmax, "riser_left_inner")
            if part_name == "riser_right" and couple_riser_inner_faces:
                _tag_faces_on_plane(vol_tag, 0, xmin, "riser_right_inner")

        try:
            outputs[name] = _mesh_single_volume(make_fn, path, mesh_size_min, mesh_size_max, _tagger)
        except Exception as exc:
            raise RuntimeError(f"Meshing failed for part '{name}': {exc}") from exc

    return outputs


def generate_ttc450_combined_mesh(
    output_msh: Path,
    mesh_size_min: float = 8.0,
    mesh_size_max: float = 25.0,
    beam_length: float = 600.0,
    base_size: float = 600.0,
    base_ext: float = 20.0,
    riser_thk: float = 6.0,
    z_beam_length: float = 150.0,
    profile: str = "cshape",
    interface_depth: float = 20.0,
    interface_gap: float = 0.0,
    include_base: bool = False,
    surface_tol: float = 1e-3,
) -> Path:
    """Generate a combined mesh with multiple disconnected volumes and facet tags.

    This is intended for MPC-based coupling in FEniCSx (all parts in one mesh).
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        width = CBEAM_40X80["width"]
        height = CBEAM_40X80["height"]
        outer_wall = CBEAM_40X80["outer_wall"]
        inner_wall = CBEAM_40X80["inner_wall"]
        slot_depth = CBEAM_40X80["slot_depth"]
        slot_width = CBEAM_40X80["slot_width"]

        y_beam_x = base_size / 2 - width / 2
        y_beam_z = base_ext + height / 2
        x_beam_z = y_beam_z + 100

        # Build volumes
        vol_tags = {}

        def _add_y_beam(name: str, x_sign: int):
            profile_tag = _add_cbeam_profile_surface(
                width=width,
                height=height,
                outer_wall=outer_wall,
                inner_wall=inner_wall,
                slot_depth=slot_depth,
                slot_width=slot_width,
                profile=profile,
                slot_dir=1 if x_sign > 0 else -1,
            )
            gmsh.model.occ.translate([(2, profile_tag)], x_sign * y_beam_x, 0, y_beam_z)
            vol = _extrude_profile(profile_tag, 0, beam_length, 0)
            vol_tags[name] = vol

        def _add_x_beam(name: str):
            profile_tag = _add_cbeam_profile_surface(
                width=width,
                height=height,
                outer_wall=outer_wall,
                inner_wall=inner_wall,
                slot_depth=slot_depth,
                slot_width=slot_width,
                profile=profile,
                slot_dir=-1,
            )
            gmsh.model.occ.rotate([(2, profile_tag)], 0, 0, 0, 0, 0, 1, math.pi / 2)
            gmsh.model.occ.translate([(2, profile_tag)], 0, 0, x_beam_z)
            vol = _extrude_profile(profile_tag, beam_length, 0, 0)
            vol_tags[name] = vol

        def _add_riser(name: str, x_sign: int):
            vol = gmsh.model.occ.addBox(
                x_sign * (beam_length / 2 + interface_gap + riser_thk / 2) - riser_thk / 2,
                -width / 2,
                x_beam_z - height / 2,
                riser_thk,
                width,
                height,
            )
            vol_tags[name] = vol

        def _add_z_block(name: str):
            vol = gmsh.model.occ.addBox(
                -width / 2,
                -(width + interface_depth + 2 * interface_gap) - width / 2,
                x_beam_z - height / 2,
                width,
                width,
                z_beam_length,
            )
            vol_tags[name] = vol

        def _add_z_interface(name: str):
            vol = gmsh.model.occ.addBox(
                -width / 2,
                -width / 2 - interface_gap - interface_depth,
                x_beam_z - height / 2,
                width,
                interface_depth,
                height,
            )
            vol_tags[name] = vol

        _add_y_beam("y_beam_left", -1)
        _add_y_beam("y_beam_right", 1)
        _add_x_beam("x_beam")
        _add_riser("riser_left", -1)
        _add_riser("riser_right", 1)
        _add_z_block("z_block")
        _add_z_interface("z_interface")

        if include_base:
            base_half = base_size / 2
            base_z = base_ext / 2
            base_parts = []
            for sy in (-1, 1):
                base_parts.append(
                    gmsh.model.occ.addBox(
                        -base_half + base_ext,
                        sy * (base_half - base_ext / 2) - base_ext / 2,
                        base_z - base_ext / 2,
                        base_size - 2 * base_ext,
                        base_ext,
                        base_ext,
                    )
                )
            for sx in (-1, 1):
                base_parts.append(
                    gmsh.model.occ.addBox(
                        sx * (base_half - base_ext / 2) - base_ext / 2,
                        -base_half + base_ext,
                        base_z - base_ext / 2,
                        base_ext,
                        base_size - 2 * base_ext,
                        base_ext,
                    )
                )
            for sx in (-1, 1):
                for sy in (-1, 1):
                    base_parts.append(
                        gmsh.model.occ.addBox(
                            sx * (base_half - base_ext / 2) - base_ext / 2,
                            sy * (base_half - base_ext / 2) - base_ext / 2,
                            base_z - base_ext / 2,
                            base_ext,
                            base_ext,
                            base_ext,
                        )
                    )
            gmsh.model.occ.synchronize()
            base_fuse, _ = gmsh.model.occ.fuse(
                [(3, base_parts[0])],
                [(3, tag) for tag in base_parts[1:]],
                removeObject=True,
                removeTool=True,
            )
            for dim, tag in base_fuse:
                if dim == 3:
                    vol_tags["base_frame"] = tag
                    break

        gmsh.model.occ.synchronize()

        # Tag volumes
        for name, tag in vol_tags.items():
            phys = gmsh.model.addPhysicalGroup(3, [tag])
            gmsh.model.setPhysicalName(3, phys, name)

        def _tag_faces_on_plane(vol_tag: int, axis: int, value: float, name: str):
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, vol_tag)
            if axis == 0:
                ents = gmsh.model.getEntitiesInBoundingBox(
                    value - surface_tol, ymin - surface_tol, zmin - surface_tol,
                    value + surface_tol, ymax + surface_tol, zmax + surface_tol, 2
                )
            elif axis == 1:
                ents = gmsh.model.getEntitiesInBoundingBox(
                    xmin - surface_tol, value - surface_tol, zmin - surface_tol,
                    xmax + surface_tol, value + surface_tol, zmax + surface_tol, 2
                )
            else:
                ents = gmsh.model.getEntitiesInBoundingBox(
                    xmin - surface_tol, ymin - surface_tol, value - surface_tol,
                    xmax + surface_tol, ymax + surface_tol, value + surface_tol, 2
                )
            faces = [tag for dim, tag in ents if dim == 2]
            if faces:
                phys = gmsh.model.addPhysicalGroup(2, faces)
                gmsh.model.setPhysicalName(2, phys, name)

        # Interface facet tags
        xb = vol_tags["x_beam"]
        xmin, _, _, xmax, _, _ = gmsh.model.getBoundingBox(3, xb)
        _tag_faces_on_plane(xb, 0, xmin, "x_beam_end_left")
        _tag_faces_on_plane(xb, 0, xmax, "x_beam_end_right")

        rl = vol_tags["riser_left"]
        _, _, _, rlxmax, _, _ = gmsh.model.getBoundingBox(3, rl)
        _tag_faces_on_plane(rl, 0, rlxmax, "riser_left_inner")

        rr = vol_tags["riser_right"]
        rrxmin, _, _, _, _, _ = gmsh.model.getBoundingBox(3, rr)
        _tag_faces_on_plane(rr, 0, rrxmin, "riser_right_inner")

        gmsh.option.setNumber("Mesh.SaveAll", 1)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size_min)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size_max)
        gmsh.model.mesh.generate(3)
        gmsh.write(str(output_msh))
        return output_msh
    finally:
        gmsh.finalize()


def merge_part_meshes(
    input_msh: List[Path],
    output_msh: Path,
) -> Path:
    """Merge multiple .msh files into a single .msh without re-meshing."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        for path in input_msh:
            gmsh.merge(str(path))
        gmsh.write(str(output_msh))
        return output_msh
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
