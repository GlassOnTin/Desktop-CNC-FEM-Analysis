// HGR20 Linear Rail Profile
// HIWIN HGR20 series linear guideway rail
//
// Reference dimensions:
// Width: 20mm
// Height: 27mm (from mounting surface to top of rail)
// Material: Hardened steel (52100 or equivalent)

// Rail dimensions
rail_width = 20;       // mm
rail_height = 27;      // mm
rail_base_height = 8;  // mm (base thickness below raceway)
raceway_height = 19;   // mm (height of profiled section)
raceway_width = 16;    // mm (width at raceway level)

// Mounting holes
mount_hole_diameter = 5.5;  // mm (for M5 bolts)
mount_hole_spacing = 60;    // mm (typical pitch)
mount_counterbore_dia = 9;  // mm
mount_counterbore_depth = 5; // mm

// HGR20 profile (2D cross-section in XZ plane)
// Rail sits on XY plane, extends up in Z
module hgr20_profile() {
    union() {
        // Base section
        translate([0, rail_base_height/2])
            square([rail_width, rail_base_height], center=true);

        // Raceway section (tapered profile)
        translate([0, rail_base_height + raceway_height/2])
            polygon([
                [-rail_width/2, -raceway_height/2],
                [rail_width/2, -raceway_height/2],
                [raceway_width/2, raceway_height/2],
                [-raceway_width/2, raceway_height/2]
            ]);
    }
}

// Simplified rectangular profile for FEM
module hgr20_profile_simple() {
    translate([0, rail_height/2])
        square([rail_width, rail_height], center=true);
}

// HGR20 rail extrusion along X axis
// Rail mounted on XY plane, extends in +Z direction
module hgr20_rail(length=600, simplified=true) {
    if (simplified) {
        // Simple rectangular profile for FEM meshing
        rotate([90, 0, 90])
            linear_extrude(height=length, center=true)
                hgr20_profile_simple();
    } else {
        // Detailed profile with mounting holes
        difference() {
            rotate([90, 0, 90])
                linear_extrude(height=length, center=true)
                    hgr20_profile();

            // Mounting holes along length
            n_holes = floor(length / mount_hole_spacing);
            for (i = [0:n_holes]) {
                x_pos = -length/2 + mount_hole_spacing/2 + i * mount_hole_spacing;
                if (x_pos < length/2 - 10) {
                    translate([x_pos, 0, -1])
                        cylinder(h=rail_base_height+2, d=mount_hole_diameter, $fn=24);
                    translate([x_pos, 0, rail_base_height - mount_counterbore_depth])
                        cylinder(h=mount_counterbore_depth+1, d=mount_counterbore_dia, $fn=24);
                }
            }
        }
    }
}

// Test render
if ($preview) {
    hgr20_rail(length=100, simplified=false);
}
