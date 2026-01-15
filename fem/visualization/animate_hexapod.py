#!/usr/bin/env python3
"""Animate hexapod Stewart platform performing a helical bore operation.

Generates frames showing the platform following a spiral toolpath,
demonstrating the parallel kinematic motion.
"""

import numpy as np
import pyvista as pv
from pathlib import Path
import json

# Paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "fem" / "results"
IMAGES_DIR = PROJECT_ROOT / "docs" / "images"
FRAMES_DIR = IMAGES_DIR / "hexapod_frames"

# Ensure output directories exist
FRAMES_DIR.mkdir(parents=True, exist_ok=True)

# Use off-screen rendering
pv.OFF_SCREEN = True

# Hexapod geometry (must match generate_hexapod.py)
BEAM_LENGTH = 600.0
BASE_RADIUS = 300.0
PLATFORM_RADIUS = 120.0
PAIR_HALF_ANGLE = 20.0
PLATFORM_ROTATION = 30.0
BASE_PILLAR_HEIGHT = 600.0
BASE_PILLAR_RADIUS = 350.0
PLATFORM_THICKNESS = 15.0

# Strut cross-section for visualization
STRUT_WIDTH = 40.0
STRUT_DEPTH = 80.0

# Animation parameters
HELIX_RADIUS = 80.0       # mm - bore radius
HELIX_DEPTH = 180.0       # mm - total depth of cut (shows strut extension)
HELIX_REVOLUTIONS = 4     # number of spiral revolutions
N_FRAMES = 150            # frames for full animation
FPS = 30                  # frames per second for GIF


def calculate_base_joints():
    """Calculate base joint positions (same as generate_hexapod.py)."""
    joints = []
    for i in range(3):
        base_angle = i * 120.0
        for sign in [-1, 1]:
            angle = np.radians(base_angle + sign * PAIR_HALF_ANGLE)
            x = BASE_RADIUS * np.cos(angle)
            y = BASE_RADIUS * np.sin(angle)
            joints.append(np.array([x, y, 0.0]))
    return joints


def calculate_platform_joints_neutral():
    """Calculate platform joint positions at neutral (home) position."""
    joints = []
    for i in range(3):
        platform_angle = i * 120.0 + PLATFORM_ROTATION
        for sign in [-1, 1]:
            angle = np.radians(platform_angle + sign * PAIR_HALF_ANGLE)
            x = PLATFORM_RADIUS * np.cos(angle)
            y = PLATFORM_RADIUS * np.sin(angle)
            joints.append(np.array([x, y, 0.0]))  # Z=0 in platform frame
    return joints


def calculate_neutral_height(base_joints, platform_joints_local, strut_length=300.0):
    """Calculate platform Z height for given strut length."""
    # Iteratively find Z such that average strut length matches target
    z_estimate = strut_length * np.cos(np.radians(30))

    for _ in range(10):
        total_length = 0
        for i in range(6):
            bi = base_joints[i]
            pi = platform_joints_local[i].copy()
            pi[2] = z_estimate
            length = np.linalg.norm(pi - bi)
            total_length += length
        avg_length = total_length / 6
        z_estimate *= strut_length / avg_length

    return z_estimate


def inverse_kinematics(base_joints, platform_joints_local, platform_pose):
    """
    Calculate strut lengths for a given platform pose.

    Args:
        base_joints: List of 6 base joint positions (fixed)
        platform_joints_local: List of 6 platform joint positions in platform frame
        platform_pose: (x, y, z, roll, pitch, yaw) - platform position and orientation

    Returns:
        List of 6 strut lengths, List of 6 platform joint world positions
    """
    x, y, z, roll, pitch, yaw = platform_pose

    # Rotation matrix from Euler angles (XYZ convention)
    cr, sr = np.cos(roll), np.sin(roll)
    cp, sp = np.cos(pitch), np.sin(pitch)
    cy, sy = np.cos(yaw), np.sin(yaw)

    R = np.array([
        [cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr],
        [sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr],
        [-sp, cp*sr, cp*cr]
    ])

    platform_center = np.array([x, y, z])

    strut_lengths = []
    platform_joints_world = []

    for i in range(6):
        # Transform platform joint to world frame
        pj_world = platform_center + R @ platform_joints_local[i]
        platform_joints_world.append(pj_world)

        # Strut length is distance from base to platform joint
        strut_vec = pj_world - base_joints[i]
        strut_lengths.append(np.linalg.norm(strut_vec))

    return strut_lengths, platform_joints_world


def generate_helix_toolpath(center, radius, depth, revolutions, n_points, direction='up'):
    """
    Generate a helical toolpath for boring operation.

    Args:
        center: (x, y, z) center of helix at start
        radius: radius of the helix
        depth: total depth of travel (always positive)
        revolutions: number of complete revolutions
        n_points: number of points on path
        direction: 'up' for lifting out of bore, 'down' for cutting into bore

    Returns:
        List of (x, y, z) points
    """
    points = []
    cx, cy, cz = center

    for i in range(n_points):
        t = i / (n_points - 1)  # 0 to 1
        angle = t * revolutions * 2 * np.pi
        z_offset = t * depth

        x = cx + radius * np.cos(angle)
        y = cy + radius * np.sin(angle)

        if direction == 'up':
            z = cz + z_offset  # Positive - lifting upward out of bore
        else:
            z = cz - z_offset  # Negative - cutting downward into bore

        points.append((x, y, z))

    return points


def create_hexagon_points(radius, z, rotation_deg=0):
    """Create vertices of a hexagon."""
    points = []
    rot = np.radians(rotation_deg)
    for i in range(6):
        angle = rot + i * np.pi / 3
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        points.append([x, y, z])
    return np.array(points)


def create_hexagonal_prism_mesh(radius, height, z_bottom, wall_thickness):
    """Create a hollow hexagonal prism mesh for visualization."""
    # Outer hexagon vertices at bottom and top
    outer_bottom = create_hexagon_points(radius, z_bottom)
    outer_top = create_hexagon_points(radius, z_bottom + height)

    # Inner hexagon (for hollow)
    inner_radius = radius - wall_thickness
    inner_bottom = create_hexagon_points(inner_radius, z_bottom)
    inner_top = create_hexagon_points(inner_radius, z_bottom + height - wall_thickness)

    # Create mesh using PyVista
    # For simplicity, create as a solid and show with transparency
    all_points = np.vstack([outer_bottom, outer_top])

    # Create faces for the prism
    faces = []
    # Bottom face
    faces.append([6, 0, 1, 2, 3, 4, 5])
    # Top face
    faces.append([6, 6, 7, 8, 9, 10, 11])
    # Side faces
    for i in range(6):
        j = (i + 1) % 6
        faces.append([4, i, j, j+6, i+6])

    faces_flat = []
    for f in faces:
        faces_flat.extend(f)

    mesh = pv.PolyData(all_points, faces=faces_flat)
    return mesh


def create_strut_mesh(start, end, width=40, depth=80, base_joint_pos=None):
    """Create a box mesh representing a strut with fixed orientation.

    The wide face (80mm) remains vertical and faces toward the Z-axis.
    y_local is kept horizontal to prevent any roll/twist appearance.
    """
    # Direction vector
    direction = end - start
    length = np.linalg.norm(direction)
    if length < 1e-6:
        return None
    direction = direction / length

    # y_local must be:
    # 1. Horizontal (no Z component) - keeps wide face vertical
    # 2. Perpendicular to strut direction
    # 3. Pointing toward center (Z-axis)
    #
    # The only horizontal vectors perpendicular to direction are ±(-dy, dx, 0)
    y_local = np.array([-direction[1], direction[0], 0])
    y_len = np.linalg.norm(y_local)

    if y_len > 0.01:
        y_local = y_local / y_len

        # Choose sign so y_local points inward (toward Z-axis)
        if base_joint_pos is not None:
            inward = -np.array([base_joint_pos[0], base_joint_pos[1], 0])
            inward = inward / np.linalg.norm(inward)
            if np.dot(y_local, inward) < 0:
                y_local = -y_local
    else:
        # Strut is nearly vertical - use inward direction directly
        if base_joint_pos is not None:
            y_local = -np.array([base_joint_pos[0], base_joint_pos[1], 0])
            y_local = y_local / np.linalg.norm(y_local)
        else:
            y_local = np.array([1, 0, 0])

    # x_local completes the right-hand coordinate system
    x_local = np.cross(y_local, direction)
    x_len = np.linalg.norm(x_local)
    if x_len > 1e-6:
        x_local = x_local / x_len
    else:
        x_local = np.array([0, 0, 1])

    # Create box vertices
    hw, hd = width/2, depth/2
    corners_local = [
        [-hd, -hw], [hd, -hw], [hd, hw], [-hd, hw]
    ]

    vertices = []
    for z_frac in [0, 1]:
        pos = start + z_frac * length * direction
        for lx, ly in corners_local:
            v = pos + lx * x_local + ly * y_local
            vertices.append(v)

    vertices = np.array(vertices)

    # Create faces
    faces = [
        [4, 0, 1, 2, 3],  # bottom
        [4, 4, 5, 6, 7],  # top
        [4, 0, 1, 5, 4],  # side 1
        [4, 1, 2, 6, 5],  # side 2
        [4, 2, 3, 7, 6],  # side 3
        [4, 3, 0, 4, 7],  # side 4
    ]
    faces_flat = [item for sublist in faces for item in sublist]

    return pv.PolyData(vertices, faces=faces_flat)


def create_ball_joint(center, radius=15):
    """Create a ball joint (sphere) mesh."""
    return pv.Sphere(radius=radius, center=center)


def create_universal_joint(center, axis, radius=20, height=25):
    """Create a universal/bearing joint (cylinder) mesh."""
    # Cylinder aligned with strut axis
    cyl = pv.Cylinder(center=center, direction=axis, radius=radius, height=height)
    return cyl


def create_platform_mesh(joints, thickness=15, radius=150):
    """Create a platform mesh as a disk."""
    center = np.mean(joints, axis=0)
    center[2] = joints[0][2]  # Use Z from joints

    disk = pv.Disc(center=center, inner=0, outer=radius, normal=(0, 0, 1))
    return disk


def render_frame(frame_idx, base_joints, platform_joints_local, neutral_z,
                 toolpath, output_dir):
    """Render a single frame of the animation."""

    # Get tool position for this frame
    tool_pos = toolpath[frame_idx]

    # Platform pose: tool position with zero rotation
    # The platform center is offset from tool by platform thickness
    platform_x = tool_pos[0]
    platform_y = tool_pos[1]
    platform_z = tool_pos[2] + PLATFORM_THICKNESS / 2

    platform_pose = (platform_x, platform_y, platform_z, 0, 0, 0)

    # Calculate strut positions
    strut_lengths, platform_joints_world = inverse_kinematics(
        base_joints, platform_joints_local, platform_pose
    )

    # Create plotter
    plotter = pv.Plotter(off_screen=True, window_size=[1280, 720])
    plotter.background_color = 'white'

    # Add base pillar (hexagonal)
    pillar = create_hexagonal_prism_mesh(
        BASE_PILLAR_RADIUS, BASE_PILLAR_HEIGHT, -BASE_PILLAR_HEIGHT, 6.0
    )
    plotter.add_mesh(pillar, color='lightgray', opacity=0.5, show_edges=True)

    # Add struts with joints
    for i in range(6):
        bi = base_joints[i]
        pi = np.array(platform_joints_world[i])

        # Calculate full beam extent
        strut_vec = pi - bi
        strut_len = np.linalg.norm(strut_vec)
        strut_dir = strut_vec / strut_len

        # Beam extends 300mm below base joint
        beam_bottom = bi - 300.0 * strut_dir
        beam_top = pi  # Platform at top

        # Create strut with fixed orientation (no axial rotation)
        strut_mesh = create_strut_mesh(beam_bottom, beam_top, STRUT_WIDTH, STRUT_DEPTH,
                                        base_joint_pos=bi)
        if strut_mesh:
            color = 'steelblue' if i % 2 == 0 else 'royalblue'
            plotter.add_mesh(strut_mesh, color=color, opacity=1.0)

        # Add ball joint at platform end (top)
        ball_joint = create_ball_joint(pi, radius=18)
        plotter.add_mesh(ball_joint, color='silver', opacity=1.0)

        # Add universal joint at base end
        # Position slightly above the base attachment point
        universal_pos = bi + 20 * strut_dir
        universal_joint = create_universal_joint(universal_pos, strut_dir, radius=22, height=30)
        plotter.add_mesh(universal_joint, color='darkgray', opacity=1.0)

    # Add platform
    platform_mesh = create_platform_mesh(platform_joints_world, PLATFORM_THICKNESS,
                                          PLATFORM_RADIUS + 30)
    plotter.add_mesh(platform_mesh, color='orange', opacity=0.9)

    # Add tool point
    tool_sphere = pv.Sphere(radius=10, center=tool_pos)
    plotter.add_mesh(tool_sphere, color='red')

    # Add toolpath trail (points up to current frame)
    if frame_idx > 0:
        trail_points = np.array(toolpath[:frame_idx+1])
        trail = pv.PolyData(trail_points)
        plotter.add_mesh(trail, color='red', point_size=3, render_points_as_spheres=True)

    # Add title
    progress = (frame_idx + 1) / len(toolpath) * 100
    plotter.add_text(f"Hexapod Helical Bore\n"
                     f"Tool: ({tool_pos[0]:.1f}, {tool_pos[1]:.1f}, {tool_pos[2]:.1f}) mm\n"
                     f"Progress: {progress:.0f}%",
                     position='upper_left', font_size=12, color='black')

    # Set camera - fixed position, no rotation
    cx, cy, cz = 0, 0, -100
    dist = 1500
    elevation = 20
    azimuth = 45

    angle_z = np.radians(elevation)
    angle_xy = np.radians(azimuth)

    cam_x = cx + dist * np.cos(angle_xy) * np.cos(angle_z)
    cam_y = cy + dist * np.sin(angle_xy) * np.cos(angle_z)
    cam_z = cz + dist * np.sin(angle_z)

    plotter.camera_position = [
        (cam_x, cam_y, cam_z),
        (cx, cy, cz),
        (0, 0, 1)
    ]

    # Save frame
    output_path = output_dir / f"frame_{frame_idx:04d}.png"
    plotter.screenshot(str(output_path))
    plotter.close()

    return output_path


def create_gif(frames_dir, output_path, fps=30):
    """Create animated GIF from frames."""
    import imageio

    # Get all frame files
    frame_files = sorted(frames_dir.glob("frame_*.png"))

    if not frame_files:
        print("No frames found!")
        return

    # Read frames
    images = []
    for f in frame_files:
        images.append(imageio.imread(f))

    # Write GIF
    duration = 1.0 / fps
    imageio.mimsave(output_path, images, duration=duration, loop=0)
    print(f"Saved animation: {output_path}")


def main():
    print("=" * 60)
    print("Generating Hexapod Helical Bore Animation")
    print("=" * 60)

    # Calculate geometry
    base_joints = calculate_base_joints()
    platform_joints_local = calculate_platform_joints_neutral()
    neutral_z = calculate_neutral_height(base_joints, platform_joints_local, 300.0)

    print(f"Neutral platform height: {neutral_z:.1f} mm")

    # Offset platform joints to include neutral Z
    for pj in platform_joints_local:
        pj[2] = 0  # Keep at 0 in local frame

    # Generate helical toolpath
    # Start with platform LOW (struts at ~2/3 extension) and lift UP (to ~1/3 extension)
    # At neutral (50% stroke), effective strut length is ~300mm
    # We want to show struts going from 2/3 extension (400mm) to 1/3 extension (200mm)
    # This means platform starts LOW and ends HIGH
    #
    # Start position: below neutral by ~60% of helix depth
    # End position: above neutral by ~40% of helix depth
    helix_start_z = neutral_z - HELIX_DEPTH * 0.6 + PLATFORM_THICKNESS

    toolpath = generate_helix_toolpath(
        center=(0, 0, helix_start_z),
        radius=HELIX_RADIUS,
        depth=HELIX_DEPTH,
        revolutions=HELIX_REVOLUTIONS,
        n_points=N_FRAMES,
        direction='up'  # Spiral upward - lifting out of bore
    )

    print(f"Toolpath: {len(toolpath)} points")
    print(f"  Start: {toolpath[0]}")
    print(f"  End: {toolpath[-1]}")
    print(f"  Helix radius: {HELIX_RADIUS} mm")
    print(f"  Helix depth: {HELIX_DEPTH} mm")

    # Render frames
    print(f"\nRendering {N_FRAMES} frames...")

    for i in range(N_FRAMES):
        if i % 10 == 0:
            print(f"  Frame {i+1}/{N_FRAMES}")
        render_frame(i, base_joints, platform_joints_local, neutral_z,
                     toolpath, FRAMES_DIR)

    print(f"\nFrames saved to: {FRAMES_DIR}")

    # Create GIF
    print("\nCreating animated GIF...")
    gif_path = IMAGES_DIR / "hexapod_helix.gif"
    create_gif(FRAMES_DIR, gif_path, FPS)

    print("\n" + "=" * 60)
    print("Animation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
