// T8 Lead Screw Component for TTC450 Pro style CNC
// Standard ACME-style trapezoidal lead screw
//
// T8 specifications:
//   Major diameter: 8 mm
//   Pitch: 2mm or 8mm (single or 4-start)
//   TTC450 Pro uses 4mm pitch (2-start likely)
//
// For FEM: modeled as simple cylinder
// Lead screws contribute stiffness but are much more compliant than ball screws

// T8 lead screw dimensions
t8_diameter = 8.0;          // mm - major thread diameter
t8_root_diameter = 6.5;     // mm - root diameter (approximate)
t8_pitch = 4.0;             // mm - TTC450 Pro spec
t8_area = 3.14159 * (t8_root_diameter/2)^2;  // ~33 mm²

// Material: Typically stainless steel or carbon steel
// E = 200 GPa for steel

// Simple cylindrical lead screw (for FEM)
// Uses root diameter for conservative stiffness estimate
module t8_leadscrew(length = 500, use_root = true) {
    d = use_root ? t8_root_diameter : t8_diameter;
    // Oriented along X-axis (for Y-axis travel)
    rotate([0, 90, 0])
        cylinder(d=d, h=length, center=true, $fn=24);
}

// Lead screw with end bearings
// More complete representation including bearing blocks
module t8_leadscrew_assembly(
    length = 500,
    bearing_block_size = 30,
    bearing_block_thickness = 15
) {
    union() {
        // Lead screw shaft
        t8_leadscrew(length);

        // Fixed end bearing block
        translate([length/2 - bearing_block_thickness/2, 0, 0])
            cube([bearing_block_thickness, bearing_block_size, bearing_block_size], center=true);

        // Free end bearing block (motor side)
        translate([-(length/2 - bearing_block_thickness/2), 0, 0])
            cube([bearing_block_thickness, bearing_block_size, bearing_block_size], center=true);
    }
}

// Nut (brass or POM anti-backlash)
// For FEM, this connects the carriage to the lead screw
module t8_nut(
    length = 20,        // mm - nut length
    flange_od = 22,     // mm - flange outer diameter
    flange_thickness = 4
) {
    union() {
        // Nut body
        rotate([0, 90, 0])
            cylinder(d=16, h=length, center=true, $fn=24);

        // Mounting flange
        translate([-(length/2 - flange_thickness/2), 0, 0])
            rotate([0, 90, 0])
            cylinder(d=flange_od, h=flange_thickness, center=true, $fn=24);
    }
}

// Dual Y-axis lead screws (TTC450 Pro configuration)
// Both Y-axes are driven by separate lead screws
module dual_y_leadscrews(
    length = 500,
    y_spacing = 649,    // Distance between Y-rails (outer_y - extrusion)
    z_offset = 60       // Height above base frame
) {
    // Left Y-axis lead screw
    translate([0, -y_spacing/2, z_offset])
        t8_leadscrew(length);

    // Right Y-axis lead screw
    translate([0, y_spacing/2, z_offset])
        t8_leadscrew(length);
}

// Calculate lead screw axial stiffness
// k = A*E/L where A = cross-sectional area, E = elastic modulus, L = length
// This is useful for understanding the compliance contribution
function leadscrew_stiffness(length, E_steel = 200e3) =
    t8_area * E_steel / length;  // N/mm

// Export module for FEM
module export_t8_leadscrew(length = 500) {
    render() {
        t8_leadscrew(length);
    }
}

module export_t8_leadscrew_assembly(length = 500) {
    render() {
        t8_leadscrew_assembly(length);
    }
}

// Test render
if ($preview) {
    dual_y_leadscrews();
}
