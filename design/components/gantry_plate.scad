// Gantry Side Plates for TwoTrees TTC450 Pro style CNC
// 8mm thick aluminum plates connecting Y-carriages to X-gantry beam
//
// These plates are critical structural elements that transfer load
// from the X-gantry to the Y-axis carriages.
//
// TTC450 Pro: Upgraded to 8mm thick for added rigidity

// Default dimensions (estimated from TTC450 Pro photos)
default_height = 150;       // mm - determines Z clearance
default_width = 100;        // mm - Y direction
default_thickness = 8;      // mm - TTC450 Pro uses 8mm (upgraded from 6mm)

// Material: 6061-T6 Aluminum (same as frame)

// Basic gantry side plate
// Simplified solid plate for FEM analysis
// Real plates have mounting holes and weight reduction cutouts
module gantry_side_plate(
    height = default_height,
    width = default_width,
    thickness = default_thickness
) {
    // Plate oriented in YZ plane
    // X is thickness direction
    // Y is width (along Y-axis travel)
    // Z is height (vertical clearance)
    cube([thickness, width, height], center=true);
}

// Plate with stiffening flange
// More realistic representation of actual plate structure
module gantry_side_plate_flanged(
    height = default_height,
    width = default_width,
    thickness = default_thickness,
    flange_depth = 20,      // mm - how far flange extends
    flange_thickness = 4    // mm - flange material thickness
) {
    union() {
        // Main vertical plate
        cube([thickness, width, height], center=true);

        // Top flange (for X-beam mounting)
        translate([flange_depth/2, 0, height/2 - flange_thickness/2])
            cube([flange_depth, width, flange_thickness], center=true);

        // Bottom flange (optional, for Y-carriage mounting)
        translate([flange_depth/2, 0, -height/2 + flange_thickness/2])
            cube([flange_depth, width, flange_thickness], center=true);
    }
}

// Plate with weight reduction (lightening holes)
// For visualization - FEM uses solid version
module gantry_side_plate_lightened(
    height = default_height,
    width = default_width,
    thickness = default_thickness,
    hole_diameter = 30,
    n_holes = 2
) {
    difference() {
        gantry_side_plate(height, width, thickness);

        // Lightening holes (vertical array)
        hole_spacing = (height - 2*hole_diameter) / (n_holes);
        for (i = [0:n_holes-1]) {
            z_pos = -height/2 + hole_diameter + i * hole_spacing + hole_spacing/2;
            translate([0, 0, z_pos])
                rotate([0, 90, 0])
                cylinder(d=hole_diameter, h=thickness+2, center=true, $fn=32);
        }
    }
}

// Paired gantry plates (left and right)
// Positioned for TTC450-style machine
module gantry_plates_pair(
    height = default_height,
    width = default_width,
    thickness = default_thickness,
    y_spacing = 600,        // Distance between Y-rails (outer_y - extrusion_size)
    x_offset = 0            // Position along Y-axis travel
) {
    // Left plate
    translate([x_offset, -y_spacing/2, 0])
        gantry_side_plate(height, width, thickness);

    // Right plate
    translate([x_offset, y_spacing/2, 0])
        gantry_side_plate(height, width, thickness);
}

// Export module for FEM
module export_gantry_plate(height = 150, width = 100, thickness = 8) {
    render() {
        gantry_side_plate(height, width, thickness);
    }
}

module export_gantry_plates_pair(height = 150, width = 100, thickness = 8, y_spacing = 600) {
    render() {
        gantry_plates_pair(height, width, thickness, y_spacing);
    }
}

// Test render
if ($preview) {
    gantry_plates_pair(y_spacing = 600);
}
