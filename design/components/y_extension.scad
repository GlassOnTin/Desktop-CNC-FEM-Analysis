// Y-Axis Extension Plate
// Laser-welded mild steel I-beam structure
// Critical component for gantry stiffness
//
// Construction:
// - 8mm mild steel plates (A36)
// - Laser-welded assembly
// - I-beam cross-section with end caps
// - Top surface bolts to X-gantry
// - Bottom surface attaches to Y-axis linear rail carriage

// Default parameters
default_web_height = 150;       // mm (vertical extent of I-beam)
default_flange_width = 120;     // mm (top/bottom flange width)
default_depth = 100;            // mm (Y-direction, front-to-back)
default_thickness = 8;          // mm (plate thickness)

// Mounting parameters
carriage_bolt_pattern_x = 32;   // mm (HGH20CA mounting pattern)
carriage_bolt_pattern_y = 60;   // mm
gantry_bolt_pattern_x = 80;     // mm (to C-beam mounting)
gantry_bolt_pattern_y = 60;     // mm

// Y-Extension plate I-beam structure
// Coordinate system:
//   X: Across I-beam web (extension width direction)
//   Y: Along I-beam depth (front-to-back)
//   Z: Vertical (I-beam height)
module y_extension_plate(
    web_height = default_web_height,
    flange_width = default_flange_width,
    depth = default_depth,
    thickness = default_thickness,
    include_mounting_holes = true
) {
    difference() {
        union() {
            // Top flange (bolts to X-gantry)
            translate([0, 0, web_height/2 - thickness/2])
                cube([flange_width, depth, thickness], center=true);

            // Bottom flange (attaches to Y-axis carriage)
            translate([0, 0, -web_height/2 + thickness/2])
                cube([flange_width, depth, thickness], center=true);

            // Central web
            cube([thickness, depth, web_height], center=true);

            // Front end cap
            translate([0, depth/2 - thickness/2, 0])
                cube([flange_width, thickness, web_height], center=true);

            // Rear end cap
            translate([0, -depth/2 + thickness/2, 0])
                cube([flange_width, thickness, web_height], center=true);
        }

        if (include_mounting_holes) {
            // Top flange mounting holes (for X-gantry attachment)
            for (x = [-gantry_bolt_pattern_x/2, gantry_bolt_pattern_x/2]) {
                for (y = [-gantry_bolt_pattern_y/2, gantry_bolt_pattern_y/2]) {
                    translate([x, y, web_height/2 - thickness - 1])
                        cylinder(h=thickness+2, d=6.5, $fn=24);
                }
            }

            // Bottom flange mounting holes (for carriage attachment)
            for (x = [-carriage_bolt_pattern_x/2, carriage_bolt_pattern_x/2]) {
                for (y = [-carriage_bolt_pattern_y/2, carriage_bolt_pattern_y/2]) {
                    translate([x, y, -web_height/2 - 1])
                        cylinder(h=thickness+2, d=5.5, $fn=24);
                }
            }
        }
    }
}

// Lightweight version with material removed from web
module y_extension_plate_lightened(
    web_height = default_web_height,
    flange_width = default_flange_width,
    depth = default_depth,
    thickness = default_thickness
) {
    difference() {
        y_extension_plate(web_height, flange_width, depth, thickness, false);

        // Lightening holes in end caps
        for (y_sign = [-1, 1]) {
            y_pos = y_sign * (depth/2 - thickness/2);
            for (z = [-web_height/4, web_height/4]) {
                translate([0, y_pos, z])
                    rotate([90, 0, 0])
                    cylinder(h=thickness+2, d=30, $fn=36, center=true);
            }
        }
    }
}

// Version optimized for stiffness-to-weight ratio
module y_extension_plate_optimized(
    web_height = default_web_height,
    flange_width = default_flange_width,
    depth = default_depth,
    thickness = default_thickness
) {
    // Thicker flanges, thinner web
    flange_t = thickness * 1.5;
    web_t = thickness * 0.75;

    union() {
        // Top flange (thicker)
        translate([0, 0, web_height/2 - flange_t/2])
            cube([flange_width, depth, flange_t], center=true);

        // Bottom flange (thicker)
        translate([0, 0, -web_height/2 + flange_t/2])
            cube([flange_width, depth, flange_t], center=true);

        // Central web (thinner)
        cube([web_t, depth, web_height - 2*flange_t], center=true);

        // End caps (match web thickness)
        translate([0, depth/2 - web_t/2, 0])
            cube([flange_width, web_t, web_height], center=true);
        translate([0, -depth/2 + web_t/2, 0])
            cube([flange_width, web_t, web_height], center=true);
    }
}

// Export module for FEM analysis
module export_y_extension(
    web_height = 150,
    flange_width = 120,
    depth = 100,
    thickness = 8
) {
    y_extension_plate(web_height, flange_width, depth, thickness, false);
}

// Test render
if ($preview) {
    color("gray")
    y_extension_plate();
}
