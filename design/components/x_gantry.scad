// X-Gantry Hybrid Beam Assembly
// Combines C-Beam 40x80 with dual HGR20 linear rails
// Creates hybrid aluminum/steel I-beam for maximum stiffness
//
// Configuration:
// - C-Beam: 40x80mm aluminum extrusion, C-slot facing forward (-Y)
// - Two HGR20 rails: Mounted on back face (+Y side) of C-Beam
// - Rail spacing: Adjustable, default 50mm center-to-center

use <c_beam_40x80.scad>
use <hgr20_rail.scad>

// Default parameters
default_length = 600;           // mm (X-axis travel)
default_rail_spacing = 50;      // mm (center-to-center)

// C-Beam dimensions (from c_beam_40x80.scad)
cbeam_width = 40;   // Y dimension
cbeam_height = 80;  // Z dimension

// Rail dimensions (from hgr20_rail.scad)
rail_height = 27;   // mm
rail_width = 20;    // mm

// Rail embedding distance for proper boolean union
rail_embed = 3.0;  // mm to embed rail into C-Beam for solid fusion

// X-Gantry hybrid beam assembly
// Coordinate system:
//   X: Along beam length (gantry travel direction)
//   Y: Front-back (C-slot opens toward -Y, rails on +Y)
//   Z: Up-down (vertical)
module x_gantry_hybrid(
    length = default_length,
    rail_spacing = default_rail_spacing,
    include_rails = true,
    hollow = true  // Use hollow C-beam profile for accurate mass/stiffness
) {
    // Use render() to force proper CSG evaluation and create solid mesh
    render() {
        union() {
            // C-Beam extrusion (centered at origin)
            // C-slot opens toward -Y
            c_beam(length, hollow=hollow);

            if (include_rails) {
                // Connecting blocks to fuse rails to C-Beam
                // These ensure proper boolean union by bridging the gap
                for (z_offset = [rail_spacing/2, -rail_spacing/2]) {
                    translate([0, cbeam_width/2 - 5, z_offset])
                        cube([length, 10, rail_height], center=true);
                }

                // Upper HGR20 rail (on back face of C-Beam)
                translate([0, cbeam_width/2 + rail_width/2 - rail_embed, rail_spacing/2])
                    rotate([90, 0, 0])
                    hgr20_rail(length);

                // Lower HGR20 rail
                translate([0, cbeam_width/2 + rail_width/2 - rail_embed, -rail_spacing/2])
                    rotate([90, 0, 0])
                    hgr20_rail(length);
            }
        }
    }
}

// Aluminum-only version for comparison analysis
module x_gantry_aluminum_only(length = default_length, hollow = true) {
    x_gantry_hybrid(length, include_rails=false, hollow=hollow);
}

// Full hybrid version with rails
module x_gantry_full(length = default_length, rail_spacing = default_rail_spacing, hollow = true) {
    x_gantry_hybrid(length, rail_spacing, include_rails=true, hollow=hollow);
}

// Version with spindle mount point marker (for visualization)
module x_gantry_with_load_point(
    length = default_length,
    rail_spacing = default_rail_spacing,
    load_x = 0,
    load_z = -50
) {
    union() {
        x_gantry_hybrid(length, rail_spacing);

        // Load application point marker (for visualization only)
        translate([load_x, -cbeam_width/2, load_z])
            color("red")
            sphere(d=10, $fn=24);
    }
}

// Export modules for FEM analysis
// These generate watertight meshes suitable for meshing

module export_x_gantry_hybrid(length = 600, rail_spacing = 50, hollow = true) {
    // Unified solid for single-material FEM approximation
    x_gantry_hybrid(length, rail_spacing, hollow=hollow);
}

module export_x_gantry_cbeam_only(length = 600, hollow = true) {
    // C-Beam alone for baseline comparison
    c_beam(length, hollow=hollow);
}

// Test render
if ($preview) {
    x_gantry_with_load_point(length=200);
}
