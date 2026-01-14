#!/usr/bin/env python3
"""Generate hexapod Stewart platform geometry using 4080 C-beam struts.

Uses same 600mm 4080 C-beam profile as TTC450 for direct comparison.
Stewart platform (6-6) configuration with paired joints.
"""

import numpy as np
import gmsh
from pathlib import Path

# Output paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "fem" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Hexapod geometry parameters (mm)
BEAM_LENGTH = 600.0           # Physical beam length (same as TTC450)
STROKE_FRACTION = 0.5         # Carriage position along beam (0.5 = mid-stroke)
STRUT_LENGTH = BEAM_LENGTH * STROKE_FRACTION  # Effective base-to-platform distance

BASE_RADIUS = 300.0           # Base joint circle radius
PLATFORM_RADIUS = 120.0       # Moving platform joint circle radius

# Joint pair angles - symmetric design with minimal twist
# Using same angle for both creates pairs of nearly-parallel struts
PAIR_HALF_ANGLE = 20.0        # Half-angle between joints in each pair (degrees)
PLATFORM_ROTATION = 30.0      # Rotation of platform pairs relative to base (degrees)
                              # 30° gives moderate twist; 0° would be no twist but singular

# Base and platform plate dimensions
BASE_THICKNESS = 20.0         # mm
PLATFORM_THICKNESS = 15.0     # mm

# 4080 C-beam profile (same as TTC450)
PROFILE_WIDTH = 80.0          # mm (4080 = 40x80)
PROFILE_HEIGHT = 40.0         # mm
WALL_THICKNESS = 1.5          # mm

# Mesh parameters
MESH_SIZE_MIN = 8.0
MESH_SIZE_MAX = 20.0


def calculate_joint_positions():
    """Calculate base and platform joint positions for 6-6 Stewart platform.

    Symmetric design: same pair angles on base and platform, with platform
    rotated by PLATFORM_ROTATION degrees. Direct pairing (base[i] to platform[i])
    creates 3 pairs of nearly-parallel struts.
    """
    base_joints = []
    platform_joints = []

    # Base joints: 3 pairs, 120° apart, each pair separated by 2*PAIR_HALF_ANGLE
    for i in range(3):
        base_angle = i * 120.0  # degrees
        for sign in [-1, 1]:
            angle = np.radians(base_angle + sign * PAIR_HALF_ANGLE)
            x = BASE_RADIUS * np.cos(angle)
            y = BASE_RADIUS * np.sin(angle)
            base_joints.append((x, y, 0.0))

    # Platform joints: same pattern but rotated and smaller radius
    for i in range(3):
        platform_angle = i * 120.0 + PLATFORM_ROTATION
        for sign in [-1, 1]:
            angle = np.radians(platform_angle + sign * PAIR_HALF_ANGLE)
            x = PLATFORM_RADIUS * np.cos(angle)
            y = PLATFORM_RADIUS * np.sin(angle)
            platform_joints.append((x, y, 0.0))  # Z will be set by platform height

    return np.array(base_joints), np.array(platform_joints)


def calculate_platform_height(base_joints, platform_joints, strut_length):
    """Calculate platform Z height such that average strut length matches target."""
    # Direct pairing: strut i connects base[i] to platform[i]

    # Start with estimate based on strut angle
    z_estimate = strut_length * np.cos(np.radians(30))

    # Iterate to find exact height
    for _ in range(10):
        total_length = 0
        for i in range(6):
            bi = base_joints[i]
            pi = platform_joints[i].copy()  # Direct pairing
            pi[2] = z_estimate

            length = np.linalg.norm(pi - bi)
            total_length += length

        avg_length = total_length / 6
        # Adjust height to match target strut length
        z_estimate *= strut_length / avg_length

    return z_estimate


def create_cbeam_along_axis(start, end, tag_base, name_prefix):
    """Create a 4080 C-beam profile along an arbitrary axis.

    The C-beam is oriented with the 80mm dimension along the strut axis,
    40mm cross-section perpendicular to it.
    """
    # Vector along strut
    axis = np.array(end) - np.array(start)
    length = np.linalg.norm(axis)
    axis_unit = axis / length

    # Create local coordinate system
    # Z-local = strut axis
    # Need to find perpendicular X and Y axes

    # Find a vector not parallel to axis
    if abs(axis_unit[2]) < 0.9:
        up = np.array([0, 0, 1])
    else:
        up = np.array([1, 0, 0])

    # X-local = perpendicular to axis
    x_local = np.cross(up, axis_unit)
    x_local = x_local / np.linalg.norm(x_local)

    # Y-local = perpendicular to both
    y_local = np.cross(axis_unit, x_local)

    # Profile dimensions (40mm width in x-local, centered)
    hw = PROFILE_HEIGHT / 2  # Half width (20mm)

    # Create hollow rectangular profile at start point
    # Outer rectangle
    outer_points = []
    w, h = PROFILE_HEIGHT, PROFILE_WIDTH  # 40x80, but 80 is along strut

    # For a C-beam along arbitrary axis, we create the profile perpendicular
    # The "length" of the beam is along the strut axis

    # Actually, for simplicity, let's create the C-beam as a simplified
    # hollow rectangular tube (the C-channel detail is less important for
    # axial loading in a hexapod)

    # Outer box corners in local coords (profile is 40x40 for simplicity in cross-section)
    # Using 40x40 cross-section to keep struts from intersecting
    profile_size = PROFILE_HEIGHT  # 40mm square cross-section
    hs = profile_size / 2

    corners_local = [
        (-hs, -hs),
        (hs, -hs),
        (hs, hs),
        (-hs, hs),
    ]

    # Transform to global coords at start point
    def local_to_global(lx, ly, z_along_axis):
        """Convert local coords to global."""
        p = np.array(start) + z_along_axis * axis_unit + lx * x_local + ly * y_local
        return tuple(p)

    # Create outer box points at start and end
    gmsh.model.occ.addPoint(*local_to_global(-hs, -hs, 0), tag=tag_base)
    gmsh.model.occ.addPoint(*local_to_global(hs, -hs, 0), tag=tag_base+1)
    gmsh.model.occ.addPoint(*local_to_global(hs, hs, 0), tag=tag_base+2)
    gmsh.model.occ.addPoint(*local_to_global(-hs, hs, 0), tag=tag_base+3)

    gmsh.model.occ.addPoint(*local_to_global(-hs, -hs, length), tag=tag_base+4)
    gmsh.model.occ.addPoint(*local_to_global(hs, -hs, length), tag=tag_base+5)
    gmsh.model.occ.addPoint(*local_to_global(hs, hs, length), tag=tag_base+6)
    gmsh.model.occ.addPoint(*local_to_global(-hs, hs, length), tag=tag_base+7)

    # Inner box (hollow)
    t = WALL_THICKNESS
    inner_hs = hs - t

    gmsh.model.occ.addPoint(*local_to_global(-inner_hs, -inner_hs, 0), tag=tag_base+8)
    gmsh.model.occ.addPoint(*local_to_global(inner_hs, -inner_hs, 0), tag=tag_base+9)
    gmsh.model.occ.addPoint(*local_to_global(inner_hs, inner_hs, 0), tag=tag_base+10)
    gmsh.model.occ.addPoint(*local_to_global(-inner_hs, inner_hs, 0), tag=tag_base+11)

    gmsh.model.occ.addPoint(*local_to_global(-inner_hs, -inner_hs, length), tag=tag_base+12)
    gmsh.model.occ.addPoint(*local_to_global(inner_hs, -inner_hs, length), tag=tag_base+13)
    gmsh.model.occ.addPoint(*local_to_global(inner_hs, inner_hs, length), tag=tag_base+14)
    gmsh.model.occ.addPoint(*local_to_global(-inner_hs, inner_hs, length), tag=tag_base+15)

    # This is getting complex - let's use OCC primitives instead
    return None  # Will use simpler approach below


def create_hexapod_geometry():
    """Create complete hexapod geometry using OCC primitives."""

    gmsh.initialize()
    gmsh.model.add("hexapod")

    # Calculate joint positions
    base_joints, platform_joints = calculate_joint_positions()

    # Calculate platform height for target strut length
    platform_z = calculate_platform_height(base_joints, platform_joints, STRUT_LENGTH)

    print(f"Hexapod geometry:")
    print(f"  Beam length: {BEAM_LENGTH} mm (4080 C-beam)")
    print(f"  Stroke position: {STROKE_FRACTION*100:.0f}% ({STRUT_LENGTH:.0f}mm effective)")
    print(f"  Base radius: {BASE_RADIUS} mm")
    print(f"  Platform radius: {PLATFORM_RADIUS} mm")
    print(f"  Platform height: {platform_z:.1f} mm")

    # Offset platform joints to correct Z height
    platform_joints_3d = platform_joints.copy()
    platform_joints_3d[:, 2] = platform_z

    # Create base plate (hexagonal, thick)
    base_plate = gmsh.model.occ.addCylinder(0, 0, -BASE_THICKNESS, 0, 0, BASE_THICKNESS,
                                             BASE_RADIUS + 50)

    # Create platform plate (smaller hexagon)
    platform_plate = gmsh.model.occ.addCylinder(0, 0, platform_z, 0, 0, PLATFORM_THICKNESS,
                                                 PLATFORM_RADIUS + 30)

    # Create 6 struts as hollow rectangular tubes (4080 C-beam profile)
    strut_volumes = []
    strut_width = PROFILE_HEIGHT   # 40mm (short dimension)
    strut_depth = PROFILE_WIDTH    # 80mm (long dimension, along radial direction)
    wall_t = WALL_THICKNESS

    print(f"\nStrut geometry (4080 C-beam: {strut_width}x{strut_depth}mm):")
    for i in range(6):
        # Direct pairing: base[i] connects to platform[i]
        bi = base_joints[i]
        pi = platform_joints_3d[i]

        # Strut vector
        strut_vec = pi - bi
        strut_len = np.linalg.norm(strut_vec)
        strut_dir = strut_vec / strut_len

        # Strut angle from vertical
        angle_from_vert = np.degrees(np.arccos(abs(strut_dir[2])))

        print(f"  Strut {i+1}: length={strut_len:.1f}mm, angle={angle_from_vert:.1f}° from vertical")

        # Create strut as a box rotated to align with strut axis
        # Use OCC extrusion along the strut direction

        # Create cross-section at base joint
        # Find local coordinate system for the strut
        # x_loc points radially outward, y_loc is tangential
        if abs(strut_dir[2]) < 0.99:
            up = np.array([0, 0, 1])
        else:
            up = np.array([1, 0, 0])

        x_loc = np.cross(up, strut_dir)
        x_loc = x_loc / np.linalg.norm(x_loc)
        y_loc = np.cross(strut_dir, x_loc)

        # Outer profile points - 40mm (width) x 80mm (depth)
        # depth is along x_loc (radial), width is along y_loc (tangential)
        hw = strut_width / 2   # 20mm half-width
        hd = strut_depth / 2   # 40mm half-depth

        def to_global(lx, ly, along=0):
            return bi + along * strut_dir + lx * x_loc + ly * y_loc

        # Create outer wire at base (80mm in x_loc direction, 40mm in y_loc)
        p1 = gmsh.model.occ.addPoint(*to_global(-hd, -hw))
        p2 = gmsh.model.occ.addPoint(*to_global(hd, -hw))
        p3 = gmsh.model.occ.addPoint(*to_global(hd, hw))
        p4 = gmsh.model.occ.addPoint(*to_global(-hd, hw))

        l1 = gmsh.model.occ.addLine(p1, p2)
        l2 = gmsh.model.occ.addLine(p2, p3)
        l3 = gmsh.model.occ.addLine(p3, p4)
        l4 = gmsh.model.occ.addLine(p4, p1)

        outer_wire = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])

        # Inner profile points (hollow) - same wall thickness all around
        ihw = hw - wall_t  # inner half-width
        ihd = hd - wall_t  # inner half-depth

        p5 = gmsh.model.occ.addPoint(*to_global(-ihd, -ihw))
        p6 = gmsh.model.occ.addPoint(*to_global(ihd, -ihw))
        p7 = gmsh.model.occ.addPoint(*to_global(ihd, ihw))
        p8 = gmsh.model.occ.addPoint(*to_global(-ihd, ihw))

        l5 = gmsh.model.occ.addLine(p5, p6)
        l6 = gmsh.model.occ.addLine(p6, p7)
        l7 = gmsh.model.occ.addLine(p7, p8)
        l8 = gmsh.model.occ.addLine(p8, p5)

        inner_wire = gmsh.model.occ.addCurveLoop([l5, l6, l7, l8])

        # Create hollow cross-section
        cross_section = gmsh.model.occ.addPlaneSurface([outer_wire, inner_wire])

        # Extrude along strut direction
        extrusion = gmsh.model.occ.extrude([(2, cross_section)],
                                            strut_len * strut_dir[0],
                                            strut_len * strut_dir[1],
                                            strut_len * strut_dir[2])

        # Find the volume from extrusion
        for item in extrusion:
            if item[0] == 3:  # Volume
                strut_volumes.append(item[1])
                break

    # Fuse all geometry
    all_volumes = [(3, base_plate), (3, platform_plate)] + [(3, v) for v in strut_volumes]

    print(f"\nFusing {len(all_volumes)} volumes...")

    if len(all_volumes) > 1:
        fused, _ = gmsh.model.occ.fuse([all_volumes[0]], all_volumes[1:])
    else:
        fused = all_volumes

    gmsh.model.occ.synchronize()

    # Add physical group for the volume
    volumes = gmsh.model.getEntities(3)
    if volumes:
        vol_tags = [v[1] for v in volumes]
        gmsh.model.addPhysicalGroup(3, vol_tags, tag=1)
        gmsh.model.setPhysicalName(3, 1, "hexapod")
        print(f"Added physical group with {len(vol_tags)} volume(s)")

    # Mesh
    print(f"\nMeshing with element size {MESH_SIZE_MIN}-{MESH_SIZE_MAX} mm...")
    gmsh.option.setNumber("Mesh.MeshSizeMin", MESH_SIZE_MIN)
    gmsh.option.setNumber("Mesh.MeshSizeMax", MESH_SIZE_MAX)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay

    gmsh.model.mesh.generate(3)

    # Get mesh statistics
    nodes = gmsh.model.mesh.getNodes()
    n_nodes = len(nodes[0])

    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(3)
    n_elements = sum(len(tags) for tags in elem_tags)

    print(f"  Nodes: {n_nodes}")
    print(f"  Elements: {n_elements}")

    # Save mesh
    output_path = RESULTS_DIR / "hexapod.msh"
    gmsh.write(str(output_path))
    print(f"\nSaved: {output_path}")

    # Also save geometry info for analysis
    info = {
        'base_joints': base_joints.tolist(),
        'platform_joints': platform_joints_3d.tolist(),
        'platform_z': platform_z,
        'strut_length': STRUT_LENGTH,
        'base_radius': BASE_RADIUS,
        'platform_radius': PLATFORM_RADIUS,
    }

    import json
    info_path = RESULTS_DIR / "hexapod_geometry.json"
    with open(info_path, 'w') as f:
        json.dump(info, f, indent=2)
    print(f"Saved: {info_path}")

    gmsh.finalize()

    return output_path, info


if __name__ == "__main__":
    create_hexapod_geometry()
