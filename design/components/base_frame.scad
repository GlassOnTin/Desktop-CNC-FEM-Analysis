// Base Frame for TwoTrees TTC450 Pro style CNC
// Uses 4040 aluminum extrusion (40x40mm V-slot)
//
// Reference: TTC450 Pro overall dimensions 742 x 689 x 413 mm
// Working area: 460 x 460 x 80 mm

// 4040 V-slot extrusion properties
extrusion_size = 40;        // mm (square profile)
extrusion_area = 540;       // mm² (typical 4040 V-slot)
extrusion_Iy = 11.22e4;     // mm⁴ (approximate)

// Default frame dimensions (TTC450 Pro style)
default_outer_x = 742;      // mm - front to back
default_outer_y = 689;      // mm - left to right
default_frame_height = 40;  // mm - single extrusion height

// Simplified 4040 V-slot profile (solid for FEM)
// Real profile has internal channels, but solid approximation works for stiffness
module extrusion_4040_profile() {
    square([extrusion_size, extrusion_size], center=true);
}

// Hollow 4040 profile for better mass accuracy
// Wall thickness calculated to match real extrusion area (540 mm²)
// Area = 40² - (40-2t)² = 540 → t ≈ 3.75 mm
module extrusion_4040_profile_hollow() {
    t = 3.75;  // Effective wall thickness
    difference() {
        square([extrusion_size, extrusion_size], center=true);
        square([extrusion_size - 2*t, extrusion_size - 2*t], center=true);
    }
}

// Single 4040 extrusion beam
// Extruded along X-axis by default
module extrusion_4040(length, hollow=true) {
    rotate([90, 0, 90])
        linear_extrude(height=length, center=true)
            if (hollow) {
                extrusion_4040_profile_hollow();
            } else {
                extrusion_4040_profile();
            }
}

// Rectangular base frame from 4040 extrusions
// Frame lies in XY plane with Z=0 at bottom
//
// Coordinate system:
//   X: Front-back (Y-axis travel direction)
//   Y: Left-right (X-gantry spans this)
//   Z: Vertical
module base_frame_4040(
    outer_x = default_outer_x,
    outer_y = default_outer_y,
    height = default_frame_height,
    hollow = true
) {
    // Calculate inner lengths (extrusions butt together at corners)
    inner_x = outer_x - 2 * extrusion_size;
    inner_y = outer_y - 2 * extrusion_size;

    // Frame is centered at origin in X-Y, bottom at Z=0
    translate([0, 0, height/2]) {
        // Front rail (along Y)
        translate([outer_x/2 - extrusion_size/2, 0, 0])
            rotate([0, 0, 90])
            extrusion_4040(inner_y, hollow);

        // Back rail (along Y)
        translate([-(outer_x/2 - extrusion_size/2), 0, 0])
            rotate([0, 0, 90])
            extrusion_4040(inner_y, hollow);

        // Left rail (along X)
        translate([0, -(outer_y/2 - extrusion_size/2), 0])
            extrusion_4040(inner_x, hollow);

        // Right rail (along X)
        translate([0, outer_y/2 - extrusion_size/2, 0])
            extrusion_4040(inner_x, hollow);

        // Corner blocks (ensure watertight mesh)
        for (x = [outer_x/2 - extrusion_size/2, -(outer_x/2 - extrusion_size/2)])
            for (y = [outer_y/2 - extrusion_size/2, -(outer_y/2 - extrusion_size/2)])
                translate([x, y, 0])
                    cube([extrusion_size, extrusion_size, extrusion_size], center=true);
    }
}

// Y-axis rail mounts on the frame sides
// These are the extrusions that the gantry rides on
module y_axis_rails(
    outer_x = default_outer_x,
    outer_y = default_outer_y,
    rail_length = 500,
    base_height = default_frame_height,
    hollow = true
) {
    // Y-axis rails sit on top of the base frame
    rail_z = base_height + extrusion_size/2;

    // Left Y-rail (along X direction)
    translate([0, -(outer_y/2 - extrusion_size/2), rail_z])
        extrusion_4040(rail_length, hollow);

    // Right Y-rail (along X direction)
    translate([0, outer_y/2 - extrusion_size/2, rail_z])
        extrusion_4040(rail_length, hollow);
}

// Complete base assembly with Y-rails
module base_with_y_rails(
    outer_x = default_outer_x,
    outer_y = default_outer_y,
    y_rail_length = 500,
    hollow = true
) {
    // Base frame
    base_frame_4040(outer_x, outer_y, hollow=hollow);

    // Y-axis rails on top
    y_axis_rails(outer_x, outer_y, y_rail_length, hollow=hollow);
}

// Export module for FEM analysis
module export_base_frame(outer_x = 742, outer_y = 689, hollow = true) {
    render() {
        base_frame_4040(outer_x, outer_y, hollow=hollow);
    }
}

module export_base_with_y_rails(outer_x = 742, outer_y = 689, y_rail_length = 500, hollow = true) {
    render() {
        base_with_y_rails(outer_x, outer_y, y_rail_length, hollow=hollow);
    }
}

// Test render
if ($preview) {
    base_with_y_rails(hollow=true);
}
