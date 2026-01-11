// TwoTrees TTC450 Pro CNC Router - Complete Assembly
// Full machine model for FEM structural analysis
//
// Specifications (from manufacturer):
//   Overall dimensions: 742 x 689 x 413 mm
//   Working area: 460 x 460 x 80 mm
//   Weight: 15.45 kg
//   Frame: 4040 aluminum extrusion
//   Lead screws: T8 (8mm dia, 4mm pitch)
//   Gantry plates: 8mm aluminum (Pro upgrade)
//   Spindle: 775 motor (76W) or optional 500W
//
// Coordinate System:
//   X: Front-back (Y-axis travel direction) - machine "Y"
//   Y: Left-right (X-gantry travel direction) - machine "X"
//   Z: Vertical (up/down)
//
// Note: This follows FEM convention where X is the primary analysis axis.
// Machine coordinates may differ (X/Y swapped).

// Import all components
use <components/base_frame.scad>
use <components/gantry_plate.scad>
use <components/t8_leadscrew.scad>
use <components/z_axis.scad>
use <components/c_beam_40x80.scad>

// TTC450 Pro dimensions (from CNX-Software review)
// https://www.cnx-software.com/2023/05/01/review-of-twotrees-ttc-450-cnc-router-machine-with-80w-and-500w-spindles/
ttc450_outer_x = 742;           // mm - frame front-back (estimated from work area + margins)
ttc450_outer_y = 689;           // mm - frame left-right (estimated)
ttc450_overall_height = 413;    // mm - total height

ttc450_working_x = 450;         // mm - Y travel (corrected from 460)
ttc450_working_y = 450;         // mm - X travel (corrected from 460)
ttc450_working_z = 80;          // mm - Z travel

// Lead screw specs: T8 5-start, 2mm pitch, 4mm lead
t8_leadscrew_dia = 8;           // mm
t8_leadscrew_lead = 4;          // mm per revolution

// Derived dimensions
frame_extrusion = 40;           // mm - 4040 extrusion size
y_rail_spacing = ttc450_outer_y - frame_extrusion;  // 649 mm
base_frame_height = frame_extrusion;  // Single layer base

// Gantry dimensions
gantry_beam_length = ttc450_outer_y;  // Spans full width
gantry_plate_height = 150;      // mm - Z clearance above base
gantry_plate_width = 100;       // mm - Y direction
gantry_plate_thickness = 8;     // mm - TTC450 Pro upgrade

// Position calculations
// NOTE: For FEM meshing, all components must physically overlap
// Positions adjusted to ensure connected geometry
base_top_z = base_frame_height;
y_rail_z = base_top_z + frame_extrusion/2;
gantry_plate_bottom_z = y_rail_z + frame_extrusion/2;  // No gap - directly on rails
gantry_beam_z = gantry_plate_bottom_z + gantry_plate_height/2;
z_axis_z = gantry_beam_z + 40;  // Above X-beam

// Connector dimensions (for ensuring mesh connectivity)
connector_overlap = 5;  // mm overlap for unions

// Complete TTC450 Pro assembly
// All components are connected with overlapping geometry for FEM meshing
module ttc450_pro_assembly(
    x_position = 0,         // Gantry position along Y-axis (mm)
    y_position = 0,         // X-carriage position along gantry (mm)
    z_position = 0,         // Z-axis position (0=up, -80=down)
    include_spindle = true,
    hollow = true           // Use hollow profiles for accurate mass
) {
    union() {
        // === BASE FRAME ===
        // 4040 extrusion rectangular frame
        base_frame_4040(ttc450_outer_x, ttc450_outer_y, hollow=hollow);

        // === Y-AXIS RAILS ===
        // Two parallel rails on top of base frame
        // Extended down slightly to overlap with base frame corners
        translate([0, 0, y_rail_z])
        for (side = [-1, 1]) {
            translate([0, side * y_rail_spacing/2, 0]) {
                extrusion_4040(ttc450_working_x + 40, hollow);
                // Vertical connector to base frame
                translate([0, 0, -frame_extrusion/2 - connector_overlap/2])
                    cube([60, frame_extrusion, connector_overlap + frame_extrusion], center=true);
            }
        }

        // === Y-AXIS LEAD SCREWS ===
        // T8 lead screws parallel to Y-rails
        // Connected to Y-rails via bearing blocks
        translate([0, 0, y_rail_z - 15])
        for (side = [-1, 1]) {
            translate([0, side * (y_rail_spacing/2 - 20), 0]) {
                t8_leadscrew(ttc450_working_x + 20);
                // Bearing block connecting to Y-rail
                translate([0, 20 * side, 15])
                    cube([30, 25, 30], center=true);
            }
        }

        // === GANTRY ASSEMBLY ===
        // Moves along Y-axis rails
        translate([x_position, 0, 0]) {

            // Gantry side plates (left and right)
            // Extended down to overlap with Y-rails
            translate([0, 0, gantry_plate_bottom_z + gantry_plate_height/2])
            for (side = [-1, 1]) {
                translate([0, side * y_rail_spacing/2, 0]) {
                    gantry_side_plate(gantry_plate_height, gantry_plate_width, gantry_plate_thickness);
                    // Carriage block connecting to Y-rail below
                    translate([0, 0, -gantry_plate_height/2 - connector_overlap/2])
                        cube([gantry_plate_width, frame_extrusion + 10, connector_overlap + frame_extrusion/2], center=true);
                }
            }

            // X-Gantry beam (spans between plates)
            // Using 4080 profile
            translate([0, 0, gantry_beam_z])
                rotate([0, 0, 90])
                c_beam(gantry_beam_length - 2*gantry_plate_thickness, hollow=hollow);

            // Connectors between X-beam ends and gantry plates
            for (side = [-1, 1]) {
                translate([0, side * (y_rail_spacing/2 - gantry_plate_thickness), gantry_beam_z])
                    cube([60, gantry_plate_thickness + 10, 80], center=true);
            }

            // === Z-AXIS ASSEMBLY ===
            // Moves along X-gantry beam
            // NOTE: Must overlap with X-beam for single mesh
            translate([0, y_position, 0]) {
                // Carriage mounting block - sits directly on X-beam
                // This creates the connection to the gantry
                translate([0, 0, gantry_beam_z])
                    cube([60, 70, 90], center=true);

                // Z-carriage plate (vertical motion plate)
                translate([40, 0, gantry_beam_z + 60])
                    cube([gantry_plate_thickness, 80, 120], center=true);

                // Spindle mount bracket
                translate([55, 0, gantry_beam_z + 50])
                    cube([20, 60, 100], center=true);

                // Spindle body (775 motor or 500W)
                if (include_spindle) {
                    translate([70, 0, gantry_beam_z + 40 + z_position])
                        cylinder(d=52, h=70, center=true, $fn=32);
                }
            }
        }
    }
}

// Simplified structural model for FEM
// Removes small details, ensures watertight mesh
module ttc450_pro_fem(
    x_position = 0,
    y_position = 0,
    z_position = 0,
    include_spindle = true,
    hollow = true
) {
    render() {
        ttc450_pro_assembly(x_position, y_position, z_position, include_spindle, hollow);
    }
}

// Export module for FEM analysis
// Single watertight mesh of complete machine
module export_ttc450_pro(
    x_position = 0,
    y_position = 0,
    hollow = true
) {
    render() {
        ttc450_pro_assembly(
            x_position = x_position,
            y_position = y_position,
            z_position = 0,
            include_spindle = true,
            hollow = hollow
        );
    }
}

// OpenBuilds C-Beam 4080 profile (from downloaded STL)
// Realistic profile with internal V-slots and webs
// Cross-sectional area: 810.5 mm² (vs ~535 mm² for simplified U-channel)
//
// STL orientation (centered at origin):
//   X: -40 to +40 mm (80mm - profile height)
//   Y: ~-17 to +23 mm (40mm - profile width)
//   Z: -300 to +300 mm (600mm length)
//
// C-beam opens toward +Y direction

// Import the 600mm C-beam for X-gantry and Y-axis beams
module c_beam_4080_600mm() {
    import("components/c_beam_4080_600mm.stl");
}

// Simplified FEM export - uses real OpenBuilds C-beam profile
// C-beams for both X-gantry and Y-axis (600mm each)
// Bed width: 600 - 2×80 = 440mm (C-beam takes 80mm width on each side)
module export_ttc450_pro_simple(
    x_position = 0,
    y_position = 0
) {
    // C-beam dimensions
    c_beam_length = 600;        // mm - length of each C-beam
    c_beam_height = 80;         // mm - profile height (vertical when installed)
    c_beam_width = 40;          // mm - profile width

    // Bed dimensions (derived from C-beam layout)
    bed_width = c_beam_length - 2 * c_beam_height;  // 600 - 160 = 440mm
    bed_depth = c_beam_length - 2 * c_beam_height;  // 600 - 160 = 440mm

    // Frame dimensions
    ext = frame_extrusion;  // 40mm base extrusion size

    // Y-axis beam positions (left and right sides)
    y_beam_y = c_beam_length/2 - c_beam_height/2;  // 260mm from center

    // Heights
    base_z = ext/2;                    // Base frame center height
    y_beam_z = ext + c_beam_height/2;  // Y-beam center height (on top of base)
    x_beam_z = y_beam_z + 100;         // X-beam center height (above Y-beams)

    render() {
        union() {
            // === BASE FRAME (4040 extrusions) ===
            // Simple rectangular frame under the Y-axis beams
            base_x = c_beam_length + 100;  // Slightly larger than beam span
            base_y = c_beam_length + 100;

            // Front/back rails (along Y)
            for (sx = [-1, 1])
                translate([sx * (base_x/2 - ext/2), 0, base_z])
                    cube([ext, base_y - 2*ext, ext], center=true);
            // Left/right rails (along X)
            for (sy = [-1, 1])
                translate([0, sy * (base_y/2 - ext/2), base_z])
                    cube([base_x - 2*ext, ext, ext], center=true);
            // Corner blocks
            for (sx = [-1, 1])
                for (sy = [-1, 1])
                    translate([sx*(base_x/2 - ext/2), sy*(base_y/2 - ext/2), base_z])
                        cube([ext, ext, ext], center=true);

            // === Y-AXIS C-BEAMS (600mm, run along X direction) ===
            // Left Y-beam
            translate([0, -y_beam_y, y_beam_z])
                rotate([0, 90, 0])  // Rotate so length is along X
                rotate([0, 0, 90])  // C-opening faces inward (+Y)
                c_beam_4080_600mm();

            // Right Y-beam
            translate([0, y_beam_y, y_beam_z])
                rotate([0, 90, 0])  // Rotate so length is along X
                rotate([0, 0, -90]) // C-opening faces inward (-Y)
                c_beam_4080_600mm();

            // Risers connecting base to Y-beams
            for (sx = [-1, 1])
                for (sy = [-1, 1])
                    translate([sx * (c_beam_length/2 - ext), sy * y_beam_y, (ext + y_beam_z)/2])
                        cube([ext + 10, c_beam_width + 10, y_beam_z - ext + 10], center=true);

            // === GANTRY ASSEMBLY (moves along Y-axis beams) ===
            translate([x_position, 0, 0]) {

                // Gantry side plates (connect to Y-beam carriages)
                for (sy = [-1, 1])
                    translate([0, sy * y_beam_y, y_beam_z + c_beam_height/2 + 60])
                        cube([80, 10, 120], center=true);

                // Carriage blocks on Y-beams
                for (sy = [-1, 1])
                    translate([0, sy * y_beam_y, y_beam_z])
                        cube([80, c_beam_width + 20, c_beam_height + 10], center=true);

                // === X-GANTRY C-BEAM (600mm, spans between Y-beams) ===
                translate([0, 0, x_beam_z])
                    rotate([90, 0, 0])   // Length along Y
                    rotate([0, 0, 180])  // C-opening faces front (+X)
                    c_beam_4080_600mm();

                // End blocks connecting X-beam to gantry plates
                for (sy = [-1, 1])
                    translate([0, sy * (y_beam_y - 20), x_beam_z])
                        cube([c_beam_height + 20, 50, c_beam_height + 10], center=true);

                // === Z-AXIS CARRIAGE (rides on X-beam) ===
                translate([0, y_position, 0]) {
                    // Carriage block wrapping X-beam
                    translate([0, 0, x_beam_z])
                        cube([c_beam_height + 40, 90, c_beam_height + 20], center=true);

                    // Z-axis plate
                    translate([c_beam_height/2 + 20, 0, x_beam_z + 50])
                        cube([12, 100, 140], center=true);

                    // Spindle mount
                    translate([c_beam_height/2 + 40, 0, x_beam_z + 40])
                        cube([30, 70, 120], center=true);

                    // 500W Spindle (52mm dia, 1.4kg)
                    translate([c_beam_height/2 + 65, 0, x_beam_z + 10])
                        cylinder(d=52, h=100, center=true, $fn=32);
                }
            }
        }
    }
}

// Export just the moving gantry (for isolated analysis)
module export_ttc450_gantry(hollow = true) {
    render() {
        union() {
            // Gantry plates
            for (side = [-1, 1]) {
                translate([0, side * y_rail_spacing/2, gantry_plate_height/2])
                    gantry_side_plate(gantry_plate_height, gantry_plate_width, gantry_plate_thickness);
            }

            // X-beam
            translate([0, 0, gantry_plate_height/2])
                rotate([0, 0, 90])
                c_beam(gantry_beam_length - 2*gantry_plate_thickness, hollow=hollow);
        }
    }
}

// Calculate approximate mass
// For verification against 15.45 kg spec
function ttc450_estimated_mass() =
    // Base frame: ~3 kg
    // Y-rails: ~1.5 kg
    // Lead screws: ~0.5 kg
    // Gantry plates: ~1 kg
    // X-beam: ~1.5 kg
    // Z-axis: ~2 kg
    // Spindle: ~1 kg
    // Hardware: ~3 kg
    // Total: ~13.5 kg (reasonable vs 15.45 kg spec)
    13.5;

// Test render
if ($preview) {
    // Show machine at center position
    ttc450_pro_assembly(
        x_position = 0,
        y_position = 0,
        z_position = 0,
        include_spindle = true,
        hollow = true
    );

    // Show working envelope (for visualization)
    %translate([0, 0, gantry_beam_z + 60])
        cube([ttc450_working_x, ttc450_working_y, ttc450_working_z], center=true);
}
