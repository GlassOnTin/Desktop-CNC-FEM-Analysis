"""Export OpenSCAD components to STL for FEM analysis."""

import subprocess
from pathlib import Path
from typing import Optional, Dict, Any

from ..config import (
    SCAD_DIR, OUTPUT_DIR, COMPONENTS,
    CBEAM_40X80, HGR20_RAIL, DEFAULT_GANTRY
)


def export_component(
    component_name: str,
    output_dir: Optional[Path] = None,
    **params
) -> Path:
    """Export a single component from OpenSCAD to STL.

    Args:
        component_name: One of 'c_beam', 'hgr20_rail', 'x_gantry', 'y_extension'
        output_dir: Directory for STL output (defaults to fem/results)
        **params: Additional parameters passed to OpenSCAD module

    Returns:
        Path to exported STL file
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get component info
    if component_name not in COMPONENTS:
        available = ', '.join(COMPONENTS.keys())
        raise ValueError(f"Unknown component '{component_name}'. Available: {available}")

    comp_info = COMPONENTS[component_name]
    scad_file = SCAD_DIR / comp_info['scad_file']
    module_name = comp_info['scad_module']

    # Build parameter string
    param_str = ', '.join(f'{k}={v}' for k, v in params.items())

    # Create wrapper SCAD that renders just the component
    wrapper_content = f'''
// Wrapper to export single component for FEM analysis
use <{scad_file}>;
{module_name}({param_str});
'''

    wrapper_path = output_dir / f"{component_name}_wrapper.scad"
    stl_path = output_dir / f"{component_name}.stl"

    with open(wrapper_path, 'w') as f:
        f.write(wrapper_content)

    # Run OpenSCAD to export STL
    result = subprocess.run(
        [
            "openscad",
            "-o", str(stl_path),
            "--export-format", "binstl",
            str(wrapper_path)
        ],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(f"OpenSCAD export failed: {result.stderr}")

    # Clean up wrapper file
    wrapper_path.unlink()

    print(f"Exported {component_name} to {stl_path}")
    return stl_path


def export_x_gantry(
    output_dir: Optional[Path] = None,
    length: float = 600,
    rail_spacing: float = 50,
    include_rails: bool = True,
    hollow: bool = True
) -> Path:
    """Export X-gantry hybrid beam assembly.

    Args:
        output_dir: Directory for STL output
        length: Gantry length in mm
        rail_spacing: Distance between rail centers in mm
        include_rails: Whether to include HGR20 rails
        hollow: Use hollow C-beam profile for accurate mass/stiffness

    Returns:
        Path to exported STL file
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    scad_file = SCAD_DIR / "components" / "x_gantry.scad"

    suffix = "_hybrid" if include_rails else "_cbeam_only"
    name = f"x_gantry{suffix}"

    wrapper_content = f'''
// Export X-gantry for FEM analysis
use <{scad_file}>;
x_gantry_hybrid(length={length}, rail_spacing={rail_spacing}, include_rails={'true' if include_rails else 'false'}, hollow={'true' if hollow else 'false'});
'''

    wrapper_path = output_dir / f"{name}_wrapper.scad"
    stl_path = output_dir / f"{name}.stl"

    with open(wrapper_path, 'w') as f:
        f.write(wrapper_content)

    result = subprocess.run(
        [
            "openscad",
            "-o", str(stl_path),
            "--export-format", "binstl",
            str(wrapper_path)
        ],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(f"OpenSCAD export failed: {result.stderr}")

    wrapper_path.unlink()
    print(f"Exported {name} to {stl_path}")
    return stl_path


def export_y_extension(
    output_dir: Optional[Path] = None,
    web_height: float = 150,
    flange_width: float = 120,
    depth: float = 100,
    thickness: float = 8
) -> Path:
    """Export Y-axis extension plate.

    Args:
        output_dir: Directory for STL output
        web_height: I-beam web height in mm
        flange_width: I-beam flange width in mm
        depth: I-beam depth (Y direction) in mm
        thickness: Plate thickness in mm

    Returns:
        Path to exported STL file
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    scad_file = SCAD_DIR / "components" / "y_extension.scad"

    wrapper_content = f'''
// Export Y-extension plate for FEM analysis
use <{scad_file}>;
export_y_extension(web_height={web_height}, flange_width={flange_width}, depth={depth}, thickness={thickness});
'''

    wrapper_path = output_dir / "y_extension_wrapper.scad"
    stl_path = output_dir / "y_extension.stl"

    with open(wrapper_path, 'w') as f:
        f.write(wrapper_content)

    result = subprocess.run(
        [
            "openscad",
            "-o", str(stl_path),
            "--export-format", "binstl",
            str(wrapper_path)
        ],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(f"OpenSCAD export failed: {result.stderr}")

    wrapper_path.unlink()
    print(f"Exported y_extension to {stl_path}")
    return stl_path


def export_x_gantry_with_extensions(
    output_dir: Optional[Path] = None,
    gantry_length: float = 600,
    rail_spacing: float = 50,
    extension_web_height: float = 150,
    extension_flange_width: float = 120,
    extension_depth: float = 100,
    extension_thickness: float = 8,
    hollow: bool = True
) -> Path:
    """Export X-gantry with attached Y-extension plates.

    This creates the complete structural subsystem for analyzing
    gantry deflection under cutting loads.

    Args:
        output_dir: Directory for STL output
        gantry_length: X-gantry length in mm
        rail_spacing: Distance between HGR20 rails in mm
        extension_*: Y-extension plate parameters
        hollow: Use hollow C-beam profile for accurate mass/stiffness

    Returns:
        Path to exported STL file
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Import paths
    x_gantry_scad = SCAD_DIR / "components" / "x_gantry.scad"
    y_extension_scad = SCAD_DIR / "components" / "y_extension.scad"

    # Calculate positions
    extension_z = -(extension_web_height / 2 + 40)  # Below gantry
    extension_y_left = -gantry_length / 2 + extension_depth / 2
    extension_y_right = gantry_length / 2 - extension_depth / 2

    wrapper_content = f'''
// Export X-gantry with Y-extension plates for FEM analysis
use <{x_gantry_scad}>;
use <{y_extension_scad}>;

render() {{
    union() {{
        // X-gantry hybrid beam
        x_gantry_hybrid(length={gantry_length}, rail_spacing={rail_spacing}, hollow={'true' if hollow else 'false'});

        // Left extension plate
        translate([0, {extension_y_left}, {extension_z}])
            export_y_extension(
                web_height={extension_web_height},
                flange_width={extension_flange_width},
                depth={extension_depth},
                thickness={extension_thickness}
            );

        // Right extension plate
        translate([0, {extension_y_right}, {extension_z}])
            export_y_extension(
                web_height={extension_web_height},
                flange_width={extension_flange_width},
                depth={extension_depth},
                thickness={extension_thickness}
            );
    }}
}}
'''

    wrapper_path = output_dir / "x_gantry_system_wrapper.scad"
    stl_path = output_dir / "x_gantry_system.stl"

    with open(wrapper_path, 'w') as f:
        f.write(wrapper_content)

    result = subprocess.run(
        [
            "openscad",
            "-o", str(stl_path),
            "--export-format", "binstl",
            str(wrapper_path)
        ],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(f"OpenSCAD export failed: {result.stderr}")

    # Keep wrapper for debugging
    # wrapper_path.unlink()
    print(f"Exported x_gantry_system to {stl_path}")
    return stl_path


def verify_stl_watertight(stl_path: Path) -> bool:
    """Verify that STL mesh is watertight (closed manifold).

    Args:
        stl_path: Path to STL file

    Returns:
        True if mesh is watertight
    """
    try:
        import trimesh
        mesh = trimesh.load(stl_path)
        return mesh.is_watertight
    except Exception as e:
        print(f"Warning: Could not verify mesh: {e}")
        return False


def export_ttc450_pro(
    output_dir: Optional[Path] = None,
    x_position: float = 0.0,
    y_position: float = 0.0,
    hollow: bool = True,
    simplified: bool = True
) -> Path:
    """Export complete TwoTrees TTC450 Pro CNC for FEM analysis.

    This exports the full structural loop including:
    - Base frame (4040 extrusions)
    - Y-axis rails and lead screws
    - Gantry side plates
    - X-gantry beam
    - Z-axis assembly with spindle

    Args:
        output_dir: Directory for STL output
        x_position: Gantry position along Y-axis (mm)
        y_position: X-carriage position along gantry (mm)
        hollow: Use hollow profiles for accurate mass (ignored if simplified=True)
        simplified: Use simplified model guaranteed to be single body for FEM

    Returns:
        Path to exported STL file
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    scad_file = SCAD_DIR / "ttc450_pro.scad"

    if simplified:
        # Use simplified export that guarantees single connected body
        wrapper_content = f'''
// Export TwoTrees TTC450 Pro (simplified) for FEM analysis
use <{scad_file}>;
export_ttc450_pro_simple(x_position={x_position}, y_position={y_position});
'''
    else:
        wrapper_content = f'''
// Export TwoTrees TTC450 Pro for FEM analysis
use <{scad_file}>;
export_ttc450_pro(x_position={x_position}, y_position={y_position}, hollow={'true' if hollow else 'false'});
'''

    wrapper_path = output_dir / "ttc450_pro_wrapper.scad"
    stl_path = output_dir / "ttc450_pro.stl"

    with open(wrapper_path, 'w') as f:
        f.write(wrapper_content)

    print(f"Exporting TTC450 Pro assembly (this may take a while)...")
    result = subprocess.run(
        [
            "openscad",
            "-o", str(stl_path),
            "--export-format", "binstl",
            str(wrapper_path)
        ],
        capture_output=True,
        text=True,
        timeout=300  # 5 minute timeout for complex model
    )

    if result.returncode != 0:
        raise RuntimeError(f"OpenSCAD export failed: {result.stderr}")

    # Keep wrapper for debugging
    # wrapper_path.unlink()
    print(f"Exported ttc450_pro to {stl_path}")
    return stl_path


def export_all_components(output_dir: Optional[Path] = None) -> Dict[str, Any]:
    """Export all structural components to STL.

    Returns:
        Dict mapping component names to export info
    """
    if output_dir is None:
        output_dir = OUTPUT_DIR

    results = {}

    # Export individual components
    exports = [
        ('x_gantry_cbeam', lambda: export_x_gantry(output_dir, include_rails=False)),
        ('x_gantry_hybrid', lambda: export_x_gantry(output_dir, include_rails=True)),
        ('y_extension', lambda: export_y_extension(output_dir)),
        ('x_gantry_system', lambda: export_x_gantry_with_extensions(output_dir)),
        ('ttc450_pro', lambda: export_ttc450_pro(output_dir)),
    ]

    for name, export_fn in exports:
        try:
            stl_path = export_fn()
            is_watertight = verify_stl_watertight(stl_path)
            results[name] = {
                'path': stl_path,
                'watertight': is_watertight
            }
            print(f"  {name}: watertight={is_watertight}")
        except Exception as e:
            print(f"  {name}: FAILED - {e}")
            results[name] = {'path': None, 'error': str(e)}

    return results


if __name__ == "__main__":
    print("Exporting CNC gantry components...")
    export_all_components()
