"""Configuration for FEM analysis of CNC gantry system."""

from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
SCAD_DIR = PROJECT_ROOT / "design"
OUTPUT_DIR = PROJECT_ROOT / "fem" / "results"

# Physical constants
GRAVITY = 9.81  # m/s^2

# =============================================================================
# Material Properties
# =============================================================================

MATERIALS = {
    'aluminum_6061_t6': {
        'name': 'Aluminum 6061-T6',
        'E': 69.0e9,      # Pa (Young's modulus)
        'nu': 0.33,       # Poisson's ratio
        'rho': 2700.0,    # kg/m^3 (density)
        'yield': 276e6,   # Pa (yield strength)
        'notes': 'C-Beam extrusion material',
    },
    'aluminum_6063_t5': {
        'name': 'Aluminum 6063-T5',
        'E': 68.3e9,      # Pa
        'nu': 0.33,
        'rho': 2700.0,
        'yield': 145e6,   # Pa (lower yield than 6061)
        'notes': 'Some V-slot profiles use this softer alloy',
    },
    'steel_a36': {
        'name': 'Steel A36 (Mild)',
        'E': 200e9,       # Pa
        'nu': 0.26,
        'rho': 7850.0,    # kg/m^3
        'yield': 250e6,   # Pa
        'notes': 'Y-axis extension plates (laser cut)',
    },
    'steel_1018': {
        'name': 'Steel 1018 (Mild)',
        'E': 205e9,       # Pa
        'nu': 0.29,
        'rho': 7870.0,
        'yield': 370e6,   # Pa
        'notes': 'Alternative mild steel',
    },
    'steel_52100': {
        'name': 'Steel 52100 (Bearing)',
        'E': 210e9,       # Pa
        'nu': 0.30,
        'rho': 7810.0,
        'yield': 2034e6,  # Pa (hardened)
        'notes': 'Linear rail material (hardened)',
    },
}

# =============================================================================
# C-Beam 40x80 Extrusion Geometry
# =============================================================================

CBEAM_40X80 = {
    'width': 40.0,          # mm (X dimension)
    'height': 80.0,         # mm (Z dimension)
    'Iy': 53.16e4,          # mm^4 (strong axis - about Y, bending in XZ plane)
    'Ix': 11.22e4,          # mm^4 (weak axis - about X)
    'area': 742.0,          # mm^2 (cross-sectional area)
    'mass_per_m': 2.00,     # kg/m (linear mass)
    # Simplified geometry for FEM modeling
    'outer_wall': 3.0,      # mm (outer wall thickness)
    'inner_wall': 2.0,      # mm (internal web thickness)
    'slot_depth': 6.0,      # mm (V-slot depth)
    'slot_width': 11.0,     # mm (slot opening width)
}

# =============================================================================
# HGR20 Linear Rail Geometry
# =============================================================================

HGR20_RAIL = {
    'width': 20.0,          # mm (rail width)
    'height': 27.0,         # mm (rail height from mounting surface)
    'rail_area': 350.0,     # mm^2 (estimated cross-sectional area)
    'mass_per_m': 2.75,     # kg/m (approximate)
    # Carriage dimensions (HGH20CA)
    'carriage_length': 63.0,    # mm
    'carriage_width': 44.0,     # mm
    'carriage_height': 30.0,    # mm
    # Mounting
    'bolt_spacing': 60.0,   # mm (mounting hole pitch)
    'bolt_size': 5.0,       # mm (M5 mounting bolts)
}

# =============================================================================
# RM1605 Ball Screw
# =============================================================================

RM1605_SCREW = {
    'shaft_diameter': 16.0,     # mm
    'lead': 5.0,                # mm/rev
    'precision': 0.05,          # mm (C7 grade)
    'nut_length': 40.0,         # mm (approximate)
    'nut_diameter': 28.0,       # mm (approximate)
}

# =============================================================================
# NEMA-23 Stepper Motor
# =============================================================================

NEMA23_MOTOR = {
    'frame_size': 57.0,     # mm (57mm = NEMA-23)
    'length': 56.0,         # mm (typical for 1.2Nm)
    'shaft_diameter': 6.35, # mm (1/4 inch)
    'torque': 1.2,          # N.m (holding torque)
    'mass': 0.55,           # kg (approximate)
}

# =============================================================================
# CNC Gantry Assembly Parameters
# =============================================================================

@dataclass
class GantryConfig:
    """Configuration for CNC gantry assembly."""
    # Axis lengths
    x_travel: float = 600.0     # mm (X-axis gantry length)
    y_travel: float = 600.0     # mm (Y-axis travel)

    # X-gantry hybrid beam
    x_beam_length: float = 600.0    # mm
    x_rail_spacing: float = 50.0    # mm (distance between HGR20 rail centers)
    x_rail_count: int = 2           # Number of rails on X-gantry

    # Y-axis rails
    y_beam_length: float = 600.0    # mm
    y_rail_count: int = 1           # Single rail per Y-axis

    # Y-axis extension plate (steel I-beam)
    extension_web_height: float = 150.0     # mm (vertical extent)
    extension_flange_width: float = 120.0   # mm (top/bottom flanges)
    extension_depth: float = 100.0          # mm (Y direction)
    extension_thickness: float = 8.0        # mm (plate thickness)

    # Spindle/router mount position (relative to gantry center)
    spindle_x_offset: float = 0.0   # mm from gantry center
    spindle_z_offset: float = -50.0 # mm below rail surface


# Default configuration
DEFAULT_GANTRY = GantryConfig()

# =============================================================================
# Cutting Load Cases
# =============================================================================

# Force components in N at tool tip
CUTTING_LOADS = {
    'heavy_x': {
        'name': 'Heavy X-direction roughing',
        'Fx': 200.0,    # N (feed direction)
        'Fy': 100.0,    # N (cross-feed)
        'Fz': 150.0,    # N (axial/thrust)
        'description': 'Aggressive aluminum slotting in X direction',
    },
    'heavy_y': {
        'name': 'Heavy Y-direction roughing',
        'Fx': 100.0,
        'Fy': 200.0,
        'Fz': 150.0,
        'description': 'Aggressive aluminum slotting in Y direction',
    },
    'plunge': {
        'name': 'Plunge cutting',
        'Fx': 50.0,
        'Fy': 50.0,
        'Fz': 300.0,
        'description': 'Vertical plunge into material',
    },
    'worst_case': {
        'name': 'Worst-case combined',
        'Fx': 200.0,
        'Fy': 200.0,
        'Fz': 300.0,
        'description': 'Maximum expected forces in all directions',
    },
    'light': {
        'name': 'Light finishing',
        'Fx': 30.0,
        'Fy': 30.0,
        'Fz': 50.0,
        'description': 'Finishing pass on aluminum',
    },
}

# =============================================================================
# Mesh Refinement Zones
# =============================================================================

def get_x_gantry_refinement(config: GantryConfig = DEFAULT_GANTRY) -> List[Tuple[float, float, float, float, float]]:
    """Get mesh refinement zones for X-gantry analysis.

    Returns:
        List of (x, y, z, radius, element_size) tuples
    """
    L = config.x_beam_length
    rail_z = config.x_rail_spacing / 2

    return [
        # End mounting zones (high stress at supports)
        (L/2 - 25, 0, 0, 30.0, 1.0),      # Right end
        (-L/2 + 25, 0, 0, 30.0, 1.0),     # Left end
        # Rail mounting surfaces
        (0, -CBEAM_40X80['width']/2, rail_z, 40.0, 1.5),   # Upper rail
        (0, -CBEAM_40X80['width']/2, -rail_z, 40.0, 1.5),  # Lower rail
        # Center (max deflection zone)
        (0, 0, 0, 50.0, 2.0),
    ]


def get_extension_refinement(config: GantryConfig = DEFAULT_GANTRY) -> List[Tuple[float, float, float, float, float]]:
    """Get mesh refinement zones for Y-extension plate analysis.

    Returns:
        List of (x, y, z, radius, element_size) tuples
    """
    h = config.extension_web_height
    t = config.extension_thickness

    return [
        # Top/bottom flange junctions (weld stress concentration)
        (0, 0, h/2 - t, 15.0, 0.5),   # Top junction
        (0, 0, -h/2 + t, 15.0, 0.5),  # Bottom junction
        # End cap junctions
        (0, config.extension_depth/2 - t, 0, 15.0, 0.5),
        (0, -config.extension_depth/2 + t, 0, 15.0, 0.5),
    ]

# =============================================================================
# Analysis Settings
# =============================================================================

DEFAULT_MESH_SIZE_MIN = 1.0     # mm
DEFAULT_MESH_SIZE_MAX = 3.0     # mm
N_MODES_DEFAULT = 10            # Number of modal frequencies to compute

# Chatter analysis parameters
SPINDLE_RPM_RANGE = (8000, 24000)   # Typical router spindle range
FLUTE_COUNTS = [1, 2, 3, 4]         # Common endmill flute counts

# =============================================================================
# Numerical Tolerances
# =============================================================================

BOUNDARY_TOLERANCE_MM = 5.0     # Distance from edge for BC region selection
COORD_TOLERANCE_MM = 0.1        # Coordinate comparison epsilon
DEFAULT_LOAD_RADIUS_MM = 10.0   # Point load distribution radius for FEM

# Spindle/tool parameters for point mass
# TTC450 Pro: 80W spindle ~0.5kg (est), 500W spindle = 1.4kg (from CNX review)
SPINDLE_MASS_KG = 1.4           # 500W spindle mass (kg)

# =============================================================================
# Analytical Validation
# =============================================================================

def analytical_beam_deflection(
    P: float,           # Applied load (N)
    L: float,           # Beam length (mm)
    E: float,           # Young's modulus (Pa)
    I: float,           # Moment of inertia (mm^4)
    support: str = 'simply_supported'
) -> float:
    """Calculate analytical beam deflection for validation.

    Args:
        P: Point load at center (N)
        L: Beam span (mm)
        E: Young's modulus (Pa)
        I: Moment of inertia (mm^4)
        support: 'simply_supported' or 'fixed_fixed'

    Returns:
        Maximum deflection in mm
    """
    # Convert E from Pa to N/mm^2 (MPa)
    E_mm = E / 1e6

    if support == 'simply_supported':
        # delta = P * L^3 / (48 * E * I)
        return P * L**3 / (48 * E_mm * I)
    elif support == 'fixed_fixed':
        # delta = P * L^3 / (192 * E * I)
        return P * L**3 / (192 * E_mm * I)
    else:
        raise ValueError(f"Unknown support type: {support}")


def analytical_beam_frequency(
    L: float,           # Beam length (mm)
    E: float,           # Young's modulus (Pa)
    I: float,           # Moment of inertia (mm^4)
    rho: float,         # Density (kg/m^3)
    A: float,           # Cross-sectional area (mm^2)
    mode: int = 1,
    support: str = 'simply_supported'
) -> float:
    """Calculate analytical beam natural frequency.

    Args:
        L: Beam span (mm)
        E: Young's modulus (Pa)
        I: Moment of inertia (mm^4)
        rho: Density (kg/m^3)
        A: Cross-sectional area (mm^2)
        mode: Mode number (1, 2, 3, ...)
        support: 'simply_supported' or 'fixed_fixed'

    Returns:
        Natural frequency in Hz
    """
    import math

    # Convert units to SI
    L_m = L / 1000      # mm to m
    I_m4 = I / 1e12     # mm^4 to m^4
    A_m2 = A / 1e6      # mm^2 to m^2

    # Beta values for different modes and support conditions
    if support == 'simply_supported':
        beta_n_L = mode * math.pi
    elif support == 'fixed_fixed':
        # Approximate values for fixed-fixed beam
        if mode == 1:
            beta_n_L = 4.730
        elif mode == 2:
            beta_n_L = 7.853
        elif mode == 3:
            beta_n_L = 10.996
        else:
            beta_n_L = (mode + 0.5) * math.pi
    else:
        raise ValueError(f"Unknown support type: {support}")

    # Natural frequency: f_n = (beta_n)^2 / (2*pi*L^2) * sqrt(E*I / (rho*A))
    omega_n = (beta_n_L / L_m)**2 * math.sqrt(E * I_m4 / (rho * A_m2))
    f_n = omega_n / (2 * math.pi)

    return f_n


# =============================================================================
# Component Registry
# =============================================================================

COMPONENTS = {
    'c_beam': {
        'scad_module': 'c_beam',
        'scad_file': 'components/c_beam_40x80.scad',
        'material': 'aluminum_6061_t6',
    },
    'hgr20_rail': {
        'scad_module': 'hgr20_rail',
        'scad_file': 'components/hgr20_rail.scad',
        'material': 'steel_52100',
    },
    'x_gantry': {
        'scad_module': 'x_gantry_hybrid',
        'scad_file': 'components/x_gantry.scad',
        'material': 'multi',  # Requires special handling
    },
    'y_axis': {
        'scad_module': 'y_axis_assembly',
        'scad_file': 'components/y_axis.scad',
        'material': 'multi',
    },
    'y_extension': {
        'scad_module': 'y_extension_plate',
        'scad_file': 'components/y_extension.scad',
        'material': 'steel_a36',
    },
    'ttc450_pro': {
        'scad_module': 'export_ttc450_pro',
        'scad_file': 'ttc450_pro.scad',
        'material': 'multi',  # Full machine with mixed materials
    },
}


# =============================================================================
# Machine Configurations
# =============================================================================

@dataclass
class MachineConfig:
    """Configuration for a complete CNC machine."""
    name: str
    # Working envelope
    x_travel: float         # mm
    y_travel: float         # mm
    z_travel: float         # mm
    # Frame dimensions
    frame_outer_x: float    # mm - overall X dimension
    frame_outer_y: float    # mm - overall Y dimension
    frame_height: float     # mm - overall height
    # Structural parameters
    extrusion_size: float   # mm - frame extrusion size (e.g., 40 for 4040)
    gantry_plate_thickness: float  # mm
    # Drive system
    leadscrew_type: str     # 'T8', 'RM1605', etc.
    leadscrew_diameter: float  # mm
    # Spindle
    spindle_mass: float     # kg
    spindle_diameter: float # mm


# TwoTrees TTC450 Pro - Baseline reference machine
# Specs from CNX-Software review (May 2023)
TTC450_PRO = MachineConfig(
    name='TwoTrees TTC450 Pro',
    x_travel=450.0,         # mm (corrected from 460)
    y_travel=450.0,         # mm (corrected from 460)
    z_travel=80.0,
    frame_outer_x=742.0,
    frame_outer_y=689.0,
    frame_height=413.0,
    extrusion_size=40.0,
    gantry_plate_thickness=8.0,
    leadscrew_type='T8',    # T8 5-start, 2mm pitch, 4mm lead
    leadscrew_diameter=8.0,
    spindle_mass=1.4,       # 500W spindle (from CNX review)
    spindle_diameter=52.0,
)

# 4040 Aluminum Extrusion (used in TTC450 frame)
EXTRUSION_4040 = {
    'size': 40.0,           # mm (square profile)
    'area': 540.0,          # mm² (typical V-slot)
    'Iy': 11.22e4,          # mm⁴ (approximate)
    'mass_per_m': 1.46,     # kg/m
}

# T8 Lead Screw (ACME-style trapezoidal)
T8_LEADSCREW = {
    'major_diameter': 8.0,  # mm
    'root_diameter': 6.5,   # mm
    'pitch': 4.0,           # mm (TTC450 Pro)
    'area': 33.2,           # mm² (at root)
    'material': 'steel_1018',
}
