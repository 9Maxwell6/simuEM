// Cube-in-Cube Geometry for Electromagnetic Simulation
// Interior: Conductor region
// Exterior: Surrounding space (air/vacuum)

// ============================================================
// Parameters
// ============================================================
// Inner cube (conductor)
inner_size = 1.0;
inner_x = 0.0;
inner_y = 0.0;
inner_z = 0.0;

// Outer cube (space)
outer_size = 3.0;
outer_x = 0.0;
outer_y = 0.0;
outer_z = 0.0;

// Mesh size
lc_inner = 0.3;   // Finer mesh on conductor
lc_outer = 0.9;   // Coarser mesh in space

// ============================================================
// Inner Cube (Conductor)
// ============================================================
// Half-sizes for convenience
hi = inner_size / 2;

// Inner cube vertices
Point(1) = {inner_x - hi, inner_y - hi, inner_z - hi, lc_inner};
Point(2) = {inner_x + hi, inner_y - hi, inner_z - hi, lc_inner};
Point(3) = {inner_x + hi, inner_y + hi, inner_z - hi, lc_inner};
Point(4) = {inner_x - hi, inner_y + hi, inner_z - hi, lc_inner};
Point(5) = {inner_x - hi, inner_y - hi, inner_z + hi, lc_inner};
Point(6) = {inner_x + hi, inner_y - hi, inner_z + hi, lc_inner};
Point(7) = {inner_x + hi, inner_y + hi, inner_z + hi, lc_inner};
Point(8) = {inner_x - hi, inner_y + hi, inner_z + hi, lc_inner};

// Inner cube edges - bottom face
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Inner cube edges - top face
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Inner cube edges - vertical
Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Inner cube faces
Curve Loop(1) = {1, 2, 3, 4};       // Bottom (-z)
Plane Surface(1) = {1};

Curve Loop(2) = {5, 6, 7, 8};       // Top (+z)
Plane Surface(2) = {2};

Curve Loop(3) = {1, 10, -5, -9};    // Front (-y)
Plane Surface(3) = {3};

Curve Loop(4) = {3, 12, -7, -11};   // Back (+y)
Plane Surface(4) = {4};

Curve Loop(5) = {2, 11, -6, -10};   // Right (+x)
Plane Surface(5) = {5};

Curve Loop(6) = {4, 9, -8, -12};    // Left (-x)
Plane Surface(6) = {6};

// Inner cube surface loop and volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// ============================================================
// Outer Cube (Space boundary)
// ============================================================
ho = outer_size / 2;

// Outer cube vertices
Point(11) = {outer_x - ho, outer_y - ho, outer_z - ho, lc_outer};
Point(12) = {outer_x + ho, outer_y - ho, outer_z - ho, lc_outer};
Point(13) = {outer_x + ho, outer_y + ho, outer_z - ho, lc_outer};
Point(14) = {outer_x - ho, outer_y + ho, outer_z - ho, lc_outer};
Point(15) = {outer_x - ho, outer_y - ho, outer_z + ho, lc_outer};
Point(16) = {outer_x + ho, outer_y - ho, outer_z + ho, lc_outer};
Point(17) = {outer_x + ho, outer_y + ho, outer_z + ho, lc_outer};
Point(18) = {outer_x - ho, outer_y + ho, outer_z + ho, lc_outer};

// Outer cube edges - bottom face
Line(21) = {11, 12};
Line(22) = {12, 13};
Line(23) = {13, 14};
Line(24) = {14, 11};

// Outer cube edges - top face
Line(25) = {15, 16};
Line(26) = {16, 17};
Line(27) = {17, 18};
Line(28) = {18, 15};

// Outer cube edges - vertical
Line(29) = {11, 15};
Line(30) = {12, 16};
Line(31) = {13, 17};
Line(32) = {14, 18};

// Outer cube faces
Curve Loop(11) = {21, 22, 23, 24};    // Bottom (-z)
Plane Surface(11) = {11};

Curve Loop(12) = {25, 26, 27, 28};    // Top (+z)
Plane Surface(12) = {12};

Curve Loop(13) = {21, 30, -25, -29};  // Front (-y)
Plane Surface(13) = {13};

Curve Loop(14) = {23, 32, -27, -31};  // Back (+y)
Plane Surface(14) = {14};

Curve Loop(15) = {22, 31, -26, -30};  // Right (+x)
Plane Surface(15) = {15};

Curve Loop(16) = {24, 29, -28, -32};  // Left (-x)
Plane Surface(16) = {16};

// ============================================================
// Exterior Volume (Space between inner and outer cubes)
// ============================================================
// Surface loop for outer boundary
Surface Loop(2) = {11, 12, 13, 14, 15, 16};

// Exterior volume: outer shell minus inner cube (hole)
Volume(2) = {2, 1};  // {outer surface loop, inner surface loop as hole}

// ============================================================
// Physical Groups (for FEM solver)
// ============================================================
// Interior region - conductor
Physical Volume("Conductor", 1) = {1};

// Exterior region - surrounding space
Physical Volume("Space", 2) = {2};

// Interface surface between conductor and space
Physical Surface("Interface", 10) = {1, 2, 3, 4, 5, 6};

// Outer boundary (for boundary conditions)
Physical Surface("OuterBoundary", 20) = {11, 12, 13, 14, 15, 16};


Mesh.MshFileVersion = 4.5;
//+
Show "*";
