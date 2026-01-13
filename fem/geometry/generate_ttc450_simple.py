"""Generate TTC450 mesh with hollow 4080 C-beam profiles and correct orientations.

Key orientations for CNC gantry:
- Y-beams (side rails): 80mm vertical (Z), 40mm wide (X), extruded along Y
  C-slot faces inward (toward machine center)
- X-gantry beam: 80mm vertical (Z), 40mm deep (Y), extruded along X
  C-slot faces front (-Y direction)
- Riser plates: connect Y-beam ends to X-beam ends

This version uses simple hollow rectangular tubes to ensure robust meshing.
"""

import gmsh
import math
from pathlib import Path

# Dimensions from config (4080 C-beam profile)
PROFILE_HEIGHT = 80.0  # mm - vertical dimension (Z)
PROFILE_WIDTH = 40.0   # mm - width/depth
WALL_THICKNESS = 3.0   # mm - wall thickness for hollow profile


def create_hollow_rect_profile(width: float, height: float, wall: float) -> int:
    """Create hollow rectangular profile in XZ plane (Y=0).

    Args:
        width: profile width (X direction)
        height: profile height (Z direction)
        wall: wall thickness

    Returns:
        Surface tag
    """
    # Outer rectangle
    outer_pts = [
        (-width/2, 0, -height/2),
        ( width/2, 0, -height/2),
        ( width/2, 0,  height/2),
        (-width/2, 0,  height/2),
    ]
    outer_p = [gmsh.model.occ.addPoint(*p) for p in outer_pts]
    outer_lines = [gmsh.model.occ.addLine(outer_p[i], outer_p[(i+1)%4]) for i in range(4)]
    outer_wire = gmsh.model.occ.addWire(outer_lines)

    # Inner rectangle (cavity)
    inner_w = width - 2*wall
    inner_h = height - 2*wall
    inner_pts = [
        (-inner_w/2, 0, -inner_h/2),
        ( inner_w/2, 0, -inner_h/2),
        ( inner_w/2, 0,  inner_h/2),
        (-inner_w/2, 0,  inner_h/2),
    ]
    inner_p = [gmsh.model.occ.addPoint(*p) for p in inner_pts]
    inner_lines = [gmsh.model.occ.addLine(inner_p[i], inner_p[(i+1)%4]) for i in range(4)]
    inner_wire = gmsh.model.occ.addWire(inner_lines)

    # Create surface with hole
    surface = gmsh.model.occ.addPlaneSurface([outer_wire, inner_wire])
    gmsh.model.occ.synchronize()
    return surface


def create_4080_cbeam_profile(wall: float = 2.0) -> int:
    """Create 4080 C-beam profile in XZ plane (Y=0).

    The C-beam is a 2x4 arrangement of 2020 cores with 2 missing:
    {{1,1,1,1},   <- Top row: 4 cores (Z+)
     {1,0,0,1}}   <- Bottom row: 2 cores at ends, C-slot in middle (Z-)

    Profile is 80mm wide (X) x 40mm tall (Z), centered at origin.
    Each 2020 core is 20mm x 20mm with hollow interior.

    Args:
        wall: wall thickness of each 2020 core

    Returns:
        Surface tag
    """
    core_size = 20.0  # mm - 2020 profile size
    inner_size = core_size - 2*wall  # Inner cavity size

    # Core positions (center of each 2020 in the grid)
    # X positions: -30, -10, +10, +30 (4 across, 20mm spacing)
    # Z positions: +10 (top row), -10 (bottom row)
    core_matrix = [
        [1, 1, 1, 1],  # Top row (Z = +10)
        [1, 0, 0, 1],  # Bottom row (Z = -10), C-slot in middle
    ]

    x_positions = [-30.0, -10.0, 10.0, 30.0]
    z_positions = [10.0, -10.0]  # Top row first, then bottom

    def make_rect_in_xz(x_center: float, z_center: float, width: float, height: float) -> int:
        """Create a rectangle in the XZ plane (Y=0)."""
        pts = [
            gmsh.model.occ.addPoint(x_center - width/2, 0, z_center - height/2),
            gmsh.model.occ.addPoint(x_center + width/2, 0, z_center - height/2),
            gmsh.model.occ.addPoint(x_center + width/2, 0, z_center + height/2),
            gmsh.model.occ.addPoint(x_center - width/2, 0, z_center + height/2),
        ]
        lines = [gmsh.model.occ.addLine(pts[i], pts[(i+1) % 4]) for i in range(4)]
        wire = gmsh.model.occ.addWire(lines)
        return gmsh.model.occ.addPlaneSurface([wire])

    # Create outer boundary of entire profile
    # Build as union of core outer boundaries, then cut inner holes
    outer_rects = []
    inner_holes = []

    for row_idx, row in enumerate(core_matrix):
        z_center = z_positions[row_idx]
        for col_idx, present in enumerate(row):
            if present:
                x_center = x_positions[col_idx]

                # Outer boundary of this core
                outer_rect = make_rect_in_xz(x_center, z_center, core_size, core_size)
                outer_rects.append((2, outer_rect))

                # Inner cavity of this core
                inner_rect = make_rect_in_xz(x_center, z_center, inner_size, inner_size)
                inner_holes.append((2, inner_rect))

    gmsh.model.occ.synchronize()

    # Fuse all outer rectangles to get the C-beam outer boundary
    if len(outer_rects) > 1:
        fused_outer, _ = gmsh.model.occ.fuse([outer_rects[0]], outer_rects[1:])
        gmsh.model.occ.synchronize()
    else:
        fused_outer = outer_rects

    # Cut inner holes from the fused outer shape
    if inner_holes:
        result, _ = gmsh.model.occ.cut(fused_outer, inner_holes)
        gmsh.model.occ.synchronize()

        # Get the resulting surface
        for dim, tag in result:
            if dim == 2:
                return tag

    return -1


def generate_ttc450_hollow_cbeam(
    output_msh: Path,
    beam_length: float = 600.0,
    mesh_size_min: float = 3.0,
    mesh_size_max: float = 12.0,
) -> Path:
    """Generate TTC450 gantry mesh with hollow C-beam profiles.

    Based on TTC450 Pro specifications:
    - Overall dimensions: 742mm x 689mm x 413mm (height)
    - Working area: 460mm x 460mm x 80mm
    - X-gantry mounted high on tall gantry risers

    All beams have:
    - 80mm vertical height (Z direction)
    - 40mm width/depth
    - Hollow interior (rectangular tube)
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    try:
        # Layout dimensions based on TTC450 Pro specs
        # Total height: 413mm, X-gantry near top
        y_beam_x_offset = 280.0  # mm from center (560mm between Y-beam centers)

        # Y-beams sit on base frame (2020 extrusions = 20mm tall)
        base_height = 20.0  # mm - base frame extrusion (2020 profile)
        y_beam_z_center = base_height + PROFILE_HEIGHT/2  # 20 + 40 = 60mm center

        # X-gantry aligned with top of riser plates
        # Riser plate top = base_height + 3*80mm = 20 + 240 = 260mm
        # X-beam top should match, so center = 260 - 40 = 220mm
        x_beam_z_center = 220.0  # mm - aligned with riser plate tops

        riser_thickness = 8.0   # mm (upgraded to 8mm on Pro model)

        all_volumes = []

        # Base frame dimensions (2020 extrusion)
        base_extrusion = 20.0  # mm - 2020 profile
        base_frame_x = 600.0   # mm - base frame X dimension
        base_frame_y = 600.0   # mm - base frame Y dimension

        print(f"Creating 2020 base frame: {base_frame_x}mm x {base_frame_y}mm")

        # Front rail (along X, at -Y)
        front_rail = gmsh.model.occ.addBox(
            -base_frame_x/2, -base_frame_y/2, 0,
            base_frame_x, base_extrusion, base_extrusion
        )
        all_volumes.append((3, front_rail))

        # Back rail (along X, at +Y)
        back_rail = gmsh.model.occ.addBox(
            -base_frame_x/2, base_frame_y/2 - base_extrusion, 0,
            base_frame_x, base_extrusion, base_extrusion
        )
        all_volumes.append((3, back_rail))

        # Left rail (along Y, at -X)
        left_rail = gmsh.model.occ.addBox(
            -base_frame_x/2, -base_frame_y/2 + base_extrusion, 0,
            base_extrusion, base_frame_y - 2*base_extrusion, base_extrusion
        )
        all_volumes.append((3, left_rail))

        # Right rail (along Y, at +X)
        right_rail = gmsh.model.occ.addBox(
            base_frame_x/2 - base_extrusion, -base_frame_y/2 + base_extrusion, 0,
            base_extrusion, base_frame_y - 2*base_extrusion, base_extrusion
        )
        all_volumes.append((3, right_rail))

        # Mid-beam (2040 profile along X-axis at center, rotated to lie flat)
        mid_beam_height = 20.0  # mm - 2040 profile rotated, 20mm tall (matches 2020 frame)
        mid_beam_depth = 40.0   # mm - 2040 profile rotated, 40mm deep
        mid_beam = gmsh.model.occ.addBox(
            -base_frame_x/2 + base_extrusion, -mid_beam_depth/2, 0,
            base_frame_x - 2*base_extrusion, mid_beam_depth, mid_beam_height
        )
        all_volumes.append((3, mid_beam))

        print(f"  Base frame rails added (4 x 2020 extrusions + 1 x 2040 mid-beam)")

        print(f"Creating Y-beams: 4080 C-beam profile, extruded {beam_length}mm along Y")

        # LEFT Y-BEAM
        # C-beam profile: 80mm wide (X in profile) x 40mm tall (Z in profile)
        # For Y-beam: need 40mm wide (X) x 80mm tall (Z), C-slot facing -X (outward)
        surf_left = create_4080_cbeam_profile(wall=1.5)
        # Rotate +90° around Y axis to get 40mm(X) x 80mm(Z), C-slot facing -X
        gmsh.model.occ.rotate([(2, surf_left)], 0, 0, 0, 0, 1, 0, math.pi/2)
        # Translate profile to left Y-beam position before extrusion
        gmsh.model.occ.translate([(2, surf_left)], -y_beam_x_offset, -beam_length/2, y_beam_z_center)
        gmsh.model.occ.synchronize()
        # Extrude along +Y
        entities = gmsh.model.occ.extrude([(2, surf_left)], 0, beam_length, 0)
        for dim, tag in entities:
            if dim == 3:
                all_volumes.append((3, tag))
                print(f"  Left Y-beam: volume tag {tag}")
                break

        # RIGHT Y-BEAM (C-slot facing +X, outward)
        surf_right = create_4080_cbeam_profile(wall=1.5)
        # Rotate -90° around Y axis to get 40mm(X) x 80mm(Z), C-slot facing +X
        gmsh.model.occ.rotate([(2, surf_right)], 0, 0, 0, 0, 1, 0, -math.pi/2)
        gmsh.model.occ.translate([(2, surf_right)], y_beam_x_offset, -beam_length/2, y_beam_z_center)
        gmsh.model.occ.synchronize()
        entities = gmsh.model.occ.extrude([(2, surf_right)], 0, beam_length, 0)
        for dim, tag in entities:
            if dim == 3:
                all_volumes.append((3, tag))
                print(f"  Right Y-beam: volume tag {tag}")
                break

        print(f"Creating X-gantry beam: 4080 C-beam profile, extruded {beam_length}mm along X")

        # X-GANTRY BEAM
        # C-beam profile: 80mm wide x 40mm tall
        # For X-beam: need 40mm deep (Y) x 80mm tall (Z), C-slot facing +Y (back)
        surf_x = create_4080_cbeam_profile(wall=1.5)
        # First rotate +90° around Y to get proper height orientation with inverted C-slot
        gmsh.model.occ.rotate([(2, surf_x)], 0, 0, 0, 0, 1, 0, math.pi/2)
        # Then rotate 90° around Z to orient profile for X extrusion
        gmsh.model.occ.rotate([(2, surf_x)], 0, 0, 0, 0, 0, 1, math.pi/2)
        # Translate to X-beam position (shifted +40mm in Y to align with riser plate tops)
        x_beam_y_offset = 40.0  # mm - align with sheared riser plate tops
        gmsh.model.occ.translate([(2, surf_x)], -beam_length/2, x_beam_y_offset, x_beam_z_center)
        gmsh.model.occ.synchronize()
        # Extrude along +X
        entities = gmsh.model.occ.extrude([(2, surf_x)], beam_length, 0, 0)
        for dim, tag in entities:
            if dim == 3:
                all_volumes.append((3, tag))
                print(f"  X-beam: volume tag {tag}")
                break

        print("Creating gantry riser plates (angled design with 23° lean)")

        # Riser plate parameters
        plate_thickness = 8.0   # mm (X direction)
        plate_depth = 80.0      # mm (Y direction)
        square_height = 80.0    # mm - height of top and bottom squares
        para_height = 80.0      # mm - height of parallelogram section
        lean_angle = 23.0       # degrees

        # Calculate parallelogram offset
        lean_offset = para_height * math.tan(math.radians(lean_angle))  # ~34mm

        # Z positions
        plate_z_bottom = y_beam_z_center - PROFILE_HEIGHT/2  # Align with bottom of Y-beams

        def create_riser_plate(x_pos: float, side: str) -> list:
            """Create angled riser plate at given X position with Y-shear.

            Args:
                x_pos: X position of plate (outer edge)
                side: 'left' or 'right' to determine plate orientation

            Returns:
                List of volume tags
            """
            z0 = plate_z_bottom
            z1 = z0 + square_height
            z2 = z1 + para_height
            z3 = z2 + square_height

            y_center = 0.0
            y_offset = lean_offset  # Shear in +Y direction

            if side == 'left':
                x0 = x_pos - plate_thickness
            else:  # right
                x0 = x_pos

            volumes = []

            # Bottom square (simple box)
            bottom_box = gmsh.model.occ.addBox(x0, y_center - plate_depth/2, z0,
                                                plate_thickness, plate_depth, square_height)
            volumes.append(bottom_box)

            # Parallelogram section: create profile at z1, extrude along angled path
            # Profile is a rectangle in XY plane at z = z1
            profile_pts = [
                gmsh.model.occ.addPoint(x0, y_center - plate_depth/2, z1),
                gmsh.model.occ.addPoint(x0 + plate_thickness, y_center - plate_depth/2, z1),
                gmsh.model.occ.addPoint(x0 + plate_thickness, y_center + plate_depth/2, z1),
                gmsh.model.occ.addPoint(x0, y_center + plate_depth/2, z1),
            ]
            profile_lines = [
                gmsh.model.occ.addLine(profile_pts[i], profile_pts[(i+1) % 4])
                for i in range(4)
            ]
            profile_wire = gmsh.model.occ.addWire(profile_lines)
            profile_surf = gmsh.model.occ.addPlaneSurface([profile_wire])
            gmsh.model.occ.synchronize()

            # Extrude along angled direction (0, y_offset, para_height) for Y-shear
            para_entities = gmsh.model.occ.extrude([(2, profile_surf)], 0, y_offset, para_height)
            for dim, tag in para_entities:
                if dim == 3:
                    volumes.append(tag)
                    break

            # Top square (simple box, shifted in +Y by y_offset)
            top_box = gmsh.model.occ.addBox(x0, y_center - plate_depth/2 + y_offset, z2,
                                            plate_thickness, plate_depth, square_height)
            volumes.append(top_box)

            return volumes

        # LEFT RISER PLATE (at X = -beam_length/2)
        left_riser_vols = create_riser_plate(-beam_length/2, 'left')
        for v in left_riser_vols:
            all_volumes.append((3, v))
        total_height = 2*square_height + para_height
        print(f"  Left riser plate: {plate_thickness}mm x {plate_depth}mm x {total_height}mm (23° Y-shear)")

        # RIGHT RISER PLATE (at X = +beam_length/2)
        right_riser_vols = create_riser_plate(beam_length/2, 'right')
        for v in right_riser_vols:
            all_volumes.append((3, v))
        print(f"  Right riser plate: {plate_thickness}mm x {plate_depth}mm x {total_height}mm (23° Y-shear)")

        gmsh.model.occ.synchronize()

        print(f"Fusing {len(all_volumes)} volumes into connected mesh...")

        # Fuse all volumes for connected mesh
        if len(all_volumes) > 1:
            fused, _ = gmsh.model.occ.fuse([all_volumes[0]], all_volumes[1:],
                                          removeObject=True, removeTool=True)
            gmsh.model.occ.synchronize()
            print(f"  Fused result: {len(fused)} volumes")

            # Tag physical group
            vol_tags = [tag for dim, tag in fused if dim == 3]
            if vol_tags:
                phys = gmsh.model.addPhysicalGroup(3, vol_tags)
                gmsh.model.setPhysicalName(3, phys, "gantry_assembly")

        # Mesh settings
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size_min)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size_max)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

        print("Generating 3D mesh...")
        gmsh.model.mesh.generate(3)

        # Get mesh statistics
        node_tags, _, _ = gmsh.model.mesh.getNodes()
        n_nodes = len(node_tags)

        elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim=3)
        n_elements = sum(len(tags) for tags in elem_tags)

        print(f"Mesh complete: {n_nodes} nodes, {n_elements} elements")

        # Write outputs
        gmsh.write(str(output_msh))
        print(f"MSH: {output_msh}")

        vtk_path = output_msh.with_suffix('.vtk')
        gmsh.write(str(vtk_path))
        print(f"VTK: {vtk_path}")

        return output_msh

    finally:
        gmsh.finalize()


if __name__ == "__main__":
    output_dir = Path(__file__).parent.parent / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_msh = output_dir / "ttc450_hollow.msh"
    generate_ttc450_hollow_cbeam(output_msh)
