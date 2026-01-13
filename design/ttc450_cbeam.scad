// TTC450 Pro - Null Baseline with Real C-beam Profiles
// Uses actual OpenBuilds C-beam 4080 STL for accurate cross-section
//
// Coordinate System (global):
//   X: left-right (gantry slot normal direction)
//   Y: front-back (gantry length direction)
//   Z: vertical (up)
//
// C-beam STL orientation (as imported):
//   Length: Z axis (600mm, from -300 to +300)
//   Height: X axis (80mm, from -40 to +40)
//   Width: Y axis (40mm, from -17.4 to +22.6, center at ~2.6mm)
//
// The C-beam slot opens toward +Y

// Key dimensions (mm)
C_BEAM_LENGTH = 600;
C_BEAM_HEIGHT = 80;
C_BEAM_WIDTH = 40;
BASE_EXT = 20;     // 2020 base frame
BASE_SIZE = 600;   // 600x600 square frame
RISER_THK = 6;     // 6mm riser plates

// Positions
Y_BEAM_X = BASE_SIZE / 2 - C_BEAM_WIDTH / 2;  // ±280mm from center (rail spacing in X)
BASE_Z = BASE_EXT / 2;                        // Base frame center at Z=10
Y_BEAM_Z = BASE_EXT + C_BEAM_HEIGHT / 2;      // Y-beam center at Z=60
X_BEAM_Z = Y_BEAM_Z + 100;                    // Gantry center height
Z_BEAM_LENGTH = 150;                          // Z-axis C-beam length

// Mesh-friendly overlap (for watertight union)
OVERLAP = 2;

// Import C-beam STL
module c_beam_stl() {
    import("components/c_beam_4080_600mm.stl");
}

// C-beam oriented for Y-axis rails (length along Y, slot facing outward ±X)
// After rotation: length along Y, height along Z, width along X (slot normal ±X)
module y_axis_cbeam(slot_direction=1) {
    // Original: length=Z, height=X, width=Y
    // Need: length=Y, height=Z, width=X
    rotate([90, 0, 0])      // Z -> Y
        rotate([0, -90, 0]) // X -> Z, Y -> X
            // Flip if slot should face -X instead of +X
            rotate([0, 0, slot_direction > 0 ? 0 : 180])
            c_beam_stl();
}

// C-beam oriented for X-gantry (length along X, slot facing -Y)
// After rotation: length along X, height along Z, width along Y (slot normal -Y)
module x_axis_cbeam() {
    // Original: length=Z, height=X, width=Y
    // Need: length=X, height=Z, width=Y, slot normal -Y
    rotate([0, 90, 0])    // Rotate length from Z to X
        rotate([0, 0, 180])  // Flip slot to face -Y
            c_beam_stl();
}

// C-beam oriented for Z-axis (length along Z, slot facing -Y)
module z_axis_cbeam(length=Z_BEAM_LENGTH) {
    // Original: length=Z, height=X, width=Y, slot normal +Y
    // Rotate so slot faces -Y
    scale([1, 1, length / C_BEAM_LENGTH])
        rotate([0, 0, 180])
            c_beam_stl();
}

// 2020 base frame (600x600 square)
module base_frame_2020(size=BASE_SIZE, ext=BASE_EXT, z=BASE_Z) {
    half = size / 2;
    // Rails along X (front/back)
    for (sy = [-1, 1])
        translate([0, sy * (half - ext/2), z])
            cube([size - 2*ext, ext, ext], center=true);
    // Rails along Y (left/right)
    for (sx = [-1, 1])
        translate([sx * (half - ext/2), 0, z])
            cube([ext, size - 2*ext, ext], center=true);
    // Corner blocks
    for (sx = [-1, 1])
        for (sy = [-1, 1])
            translate([sx * (half - ext/2), sy * (half - ext/2), z])
                cube([ext, ext, ext], center=true);
}

// Interface cuboid (flat-faced connector for mesh robustness)
module interface_block(size=[40, 40, 40], center=[0,0,0]) {
    translate(center)
        cube(size, center=true);
}

// Main assembly module
module ttc450_cbeam_assembly(x_pos=0, y_pos=0) {
    render() {
        union() {
            // === BASE FRAME ===
            base_frame_2020();

            // Vertical posts (connect base to Y-beams)
            post_bottom = BASE_Z + BASE_EXT/2 - OVERLAP;
            post_top = Y_BEAM_Z - C_BEAM_HEIGHT/2 + OVERLAP;
            post_height = post_top - post_bottom;
            for (sx = [-1, 1])
                for (sy = [-1, 1])
                    translate([sx * (C_BEAM_LENGTH/2 - BASE_EXT), sy * (C_BEAM_LENGTH/2 - BASE_EXT), (post_bottom + post_top)/2])
                        cube([BASE_EXT, BASE_EXT, post_height], center=true);

            // === Y-AXIS C-BEAMS ===
            // Left Y-beam (slot faces -Y, outward)
            translate([-Y_BEAM_X, 0, Y_BEAM_Z])
                y_axis_cbeam(slot_direction=-1);

            // Right Y-beam (slot faces +Y, outward)
            translate([Y_BEAM_X, 0, Y_BEAM_Z])
                y_axis_cbeam(slot_direction=1);

            // === GANTRY ASSEMBLY ===
            translate([x_pos, 0, 0]) {
                // Riser plates outside Y-beams (6mm)
                // Plates sit outside the Y-rails and attach to gantry ends
                plate_height = X_BEAM_Z - Y_BEAM_Z + C_BEAM_HEIGHT;
                for (sx = [-1, 1]) {
                    translate([sx * (BASE_SIZE/2 + RISER_THK/2 - OVERLAP), 0, (Y_BEAM_Z + X_BEAM_Z)/2])
                        cube([RISER_THK, C_BEAM_WIDTH, plate_height], center=true);
                }

                // === X-GANTRY C-BEAM ===
                translate([0, 0, X_BEAM_Z])
                    x_axis_cbeam();

                // X-beam end connectors to riser plates (flat-faced cuboids)
                for (sx = [-1, 1])
                    interface_block(
                        size=[RISER_THK + 10, C_BEAM_WIDTH + 10, C_BEAM_HEIGHT + 10],
                        center=[sx * (BASE_SIZE/2), 0, X_BEAM_Z]
                    );

                // === Z-AXIS CARRIAGE ===
                translate([0, y_pos, 0]) {
                    // Z-axis C-beam (150mm) - back face crosses gantry
                    z_center_z = (X_BEAM_Z - C_BEAM_HEIGHT/2 + OVERLAP) + Z_BEAM_LENGTH/2;
                    translate([0, -C_BEAM_WIDTH/2 + OVERLAP, z_center_z])
                        z_axis_cbeam(length=Z_BEAM_LENGTH);

                    // Interface block between gantry and Z-axis beam
                    interface_block(
                        size=[C_BEAM_WIDTH + 10, 60, C_BEAM_HEIGHT + 10],
                        center=[0, -C_BEAM_WIDTH/2 + OVERLAP, X_BEAM_Z]
                    );
                }
            }
        }
    }
}

// Preview
if ($preview) {
    ttc450_cbeam_assembly();
}
