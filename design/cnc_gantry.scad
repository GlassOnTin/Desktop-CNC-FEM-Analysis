// CNC Gantry Main Assembly
// Complete parametric assembly of CNC router gantry system
//
// Components:
// - X-axis: C-Beam 40x80 + dual HGR20 rails (hybrid I-beam)
// - Y-axes: Dual C-Beam 40x80 + single HGR20 rails
// - Y-extension plates: Laser-welded steel I-beams
//
// Coordinate system:
//   X: Left-right (gantry travel)
//   Y: Front-back (table direction)
//   Z: Up-down (spindle direction)

use <components/c_beam_40x80.scad>
use <components/hgr20_rail.scad>
use <components/x_gantry.scad>
use <components/y_extension.scad>

// Machine parameters
x_travel = 600;         // mm
y_travel = 600;         // mm
y_axis_spacing = 500;   // mm (distance between Y-axis rails)

// Extension plate parameters
extension_web_height = 150;
extension_flange_width = 120;
extension_depth = 100;
extension_thickness = 8;

// Position parameters
x_gantry_z = 200;       // mm (height of X-gantry above Y-rails)

// Full CNC gantry assembly
module cnc_gantry_assembly() {
    // Left Y-axis rail (C-Beam + HGR20)
    translate([0, -y_axis_spacing/2, 0])
        rotate([0, 0, 90])
        color("silver") {
            c_beam(y_travel);
            translate([0, 40/2, 80/2 + 27/2])
                hgr20_rail(y_travel);
        }

    // Right Y-axis rail
    translate([0, y_axis_spacing/2, 0])
        rotate([0, 0, 90])
        color("silver") {
            c_beam(y_travel);
            translate([0, 40/2, 80/2 + 27/2])
                hgr20_rail(y_travel);
        }

    // X-gantry (at center position)
    translate([0, 0, x_gantry_z])
        color("silver")
        x_gantry_hybrid(x_travel);

    // Left Y-extension plate
    translate([0, -y_axis_spacing/2 + 40/2 + extension_depth/2, x_gantry_z - extension_web_height/2 - 40])
        color("gray")
        y_extension_plate(extension_web_height, extension_flange_width, extension_depth, extension_thickness, false);

    // Right Y-extension plate
    translate([0, y_axis_spacing/2 - 40/2 - extension_depth/2, x_gantry_z - extension_web_height/2 - 40])
        color("gray")
        y_extension_plate(extension_web_height, extension_flange_width, extension_depth, extension_thickness, false);
}

// X-gantry subsystem for isolated FEM analysis
module x_gantry_subsystem(x_pos = 0) {
    translate([x_pos, 0, 0]) {
        // X-gantry beam
        x_gantry_hybrid(x_travel);

        // Left extension plate (attached)
        translate([0, -x_travel/2 + extension_depth/2, -extension_web_height/2 - 40])
            y_extension_plate(extension_web_height, extension_flange_width, extension_depth, extension_thickness, false);

        // Right extension plate (attached)
        translate([0, x_travel/2 - extension_depth/2, -extension_web_height/2 - 40])
            y_extension_plate(extension_web_height, extension_flange_width, extension_depth, extension_thickness, false);
    }
}

// Export modules for FEM
module export_full_assembly() {
    cnc_gantry_assembly();
}

module export_x_gantry_with_extensions() {
    x_gantry_subsystem();
}

// Test render
if ($preview) {
    cnc_gantry_assembly();
}
