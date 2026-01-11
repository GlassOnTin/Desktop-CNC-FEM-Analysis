// Z-Axis Assembly for TwoTrees TTC450 Pro style CNC
// Vertical motion system with T8 lead screw and spindle mount
//
// TTC450 Pro specifications:
//   Z travel: 80 mm
//   Spindle: 775 motor (52mm diameter) or optional 500W spindle
//   Drive: T8 lead screw with optical axis guide rails

// Import dependencies
use <t8_leadscrew.scad>

// Z-axis dimensions (from TTC450 Pro)
default_z_travel = 80;              // mm
default_plate_width = 100;          // mm - X direction
default_plate_height = 150;         // mm - Z direction
default_plate_thickness = 8;        // mm

// Spindle dimensions
spindle_775_diameter = 52;          // mm - 775 motor
spindle_775_length = 70;            // mm - body length (approximate)
spindle_mount_offset = -30;         // mm - below top of Z-plate

// Linear guide dimensions (optical axis style)
guide_rod_diameter = 8;             // mm
guide_rod_spacing = 60;             // mm - horizontal spacing

// Z-axis carriage plate
// This is the moving plate that carries the spindle
module z_carriage_plate(
    width = default_plate_width,
    height = default_plate_height,
    thickness = default_plate_thickness
) {
    // Plate in YZ plane, X is thickness
    cube([thickness, width, height], center=true);
}

// Spindle mount (simplified for FEM)
// Represents the spindle as a cylinder attached to Z-plate
module spindle_mount(
    diameter = spindle_775_diameter,
    length = spindle_775_length,
    mount_offset = spindle_mount_offset
) {
    // Spindle body (oriented along Z, pointing down)
    translate([0, 0, mount_offset - length/2])
        cylinder(d=diameter, h=length, center=true, $fn=32);
}

// Spindle mount clamp (aluminum bracket)
module spindle_clamp(
    spindle_diameter = spindle_775_diameter,
    clamp_thickness = 15,
    clamp_width = 60,
    clamp_height = 50
) {
    difference() {
        // Clamp body
        cube([clamp_thickness, clamp_width, clamp_height], center=true);

        // Spindle bore
        rotate([0, 90, 0])
            cylinder(d=spindle_diameter + 1, h=clamp_thickness + 2, center=true, $fn=32);
    }
}

// Complete Z-axis carriage assembly
// Includes plate, linear bearings positions, and spindle mount
module z_carriage_assembly(
    plate_width = default_plate_width,
    plate_height = default_plate_height,
    plate_thickness = default_plate_thickness,
    include_spindle = true,
    spindle_diameter = spindle_775_diameter
) {
    union() {
        // Z carriage plate
        z_carriage_plate(plate_width, plate_height, plate_thickness);

        // Spindle clamp (extends forward from plate)
        translate([plate_thickness/2 + 7.5, 0, spindle_mount_offset])
            spindle_clamp(spindle_diameter);

        // Spindle (for visualization and mass)
        if (include_spindle) {
            translate([plate_thickness/2 + 15, 0, 0])
                spindle_mount(spindle_diameter);
        }

        // Linear bearing blocks (simplified)
        // These connect to guide rods
        for (y = [guide_rod_spacing/2, -guide_rod_spacing/2]) {
            // Upper bearing
            translate([-plate_thickness/2 - 10, y, plate_height/4])
                cube([20, 20, 30], center=true);
            // Lower bearing
            translate([-plate_thickness/2 - 10, y, -plate_height/4])
                cube([20, 20, 30], center=true);
        }
    }
}

// Z-axis frame (fixed part attached to X-gantry)
// Includes guide rods and lead screw mounting
module z_axis_frame(
    travel = default_z_travel,
    plate_height = default_plate_height,
    frame_thickness = 8
) {
    total_height = travel + plate_height + 20;  // Extra for bearing mounts

    union() {
        // Main mounting plate (attaches to X-gantry carriage)
        translate([0, 0, total_height/2])
            cube([frame_thickness, 80, total_height], center=true);

        // Guide rods (optical axis style)
        for (y = [guide_rod_spacing/2, -guide_rod_spacing/2]) {
            translate([-20, y, total_height/2])
                rotate([0, 0, 0])
                cylinder(d=guide_rod_diameter, h=total_height - 20, center=true, $fn=16);
        }

        // Z-axis lead screw (T8)
        translate([-20, 0, total_height/2])
            rotate([90, 0, 0])
            rotate([0, 90, 0])
            t8_leadscrew(total_height - 40);
    }
}

// Complete Z-axis assembly (frame + carriage)
// For FEM: models the complete Z-axis structural loop
module z_axis_complete(
    travel = default_z_travel,
    z_position = 0,         // 0 = fully up, -travel = fully down
    include_spindle = true
) {
    // Z-axis frame (fixed to X-carriage)
    z_axis_frame(travel);

    // Z carriage (moves along frame)
    translate([-30, 0, z_position + default_plate_height/2 + 10])
        z_carriage_assembly(include_spindle=include_spindle);
}

// Export module for FEM analysis
// Simplified single-body export for meshing
module export_z_axis(travel = 80, include_spindle = true) {
    render() {
        z_axis_complete(travel, z_position=0, include_spindle=include_spindle);
    }
}

// Export just the carriage for separate analysis
module export_z_carriage(include_spindle = true) {
    render() {
        z_carriage_assembly(include_spindle=include_spindle);
    }
}

// Test render
if ($preview) {
    z_axis_complete(include_spindle=true);
}
