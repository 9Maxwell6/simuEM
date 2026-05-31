// Square-in-Square Geometry for Electromagnetic Simulation (2D)
// Interior: Conductor region
// Exterior: Surrounding space (air/vacuum)
// 2D analog of the original cube-in-cube geometry.

// ============================================================
// Parameters
// ============================================================
// Inner square (conductor)
inner_size = 1.0;
inner_x = 0.0;
inner_y = 0.0;

// Outer square (space)
outer_size = 3.0;
outer_x = 0.0;
outer_y = 0.0;

// Mesh size
lc_inner = 0.3;   // Finer mesh on conductor
lc_outer = 0.3;   // Coarser mesh in space

// ============================================================
// Inner Square (Conductor)
// ============================================================
// Half-size for convenience
hi = inner_size / 2;

// Inner square vertices
Point(1) = {inner_x - hi, inner_y - hi, 0, lc_inner};
Point(2) = {inner_x + hi, inner_y - hi, 0, lc_inner};
Point(3) = {inner_x + hi, inner_y + hi, 0, lc_inner};
Point(4) = {inner_x - hi, inner_y + hi, 0, lc_inner};

// Inner square edges
Line(1) = {1, 2};   // Bottom (-y)
Line(2) = {2, 3};   // Right  (+x)
Line(3) = {3, 4};   // Top    (+y)
Line(4) = {4, 1};   // Left   (-x)

// Inner square surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// ============================================================
// Outer Square (Space boundary)
// ============================================================
ho = outer_size / 2;

// Outer square vertices
Point(5) = {outer_x - ho, outer_y - ho, 0, lc_outer};
Point(6) = {outer_x + ho, outer_y - ho, 0, lc_outer};
Point(7) = {outer_x + ho, outer_y + ho, 0, lc_outer};
Point(8) = {outer_x - ho, outer_y + ho, 0, lc_outer};

// Outer square edges
Line(5) = {5, 6};   // Bottom (-y)
Line(6) = {6, 7};   // Right  (+x)
Line(7) = {7, 8};   // Top    (+y)
Line(8) = {8, 5};   // Left   (-x)

// Outer square boundary loop
Curve Loop(2) = {5, 6, 7, 8};

// ============================================================
// Exterior Surface (Space between inner and outer squares)
// ============================================================
// Exterior surface: outer square minus inner square (hole)
Plane Surface(2) = {2, 1};  // {outer loop, inner loop as hole}

// ============================================================
// Physical Groups (for FEM solver)
// ============================================================
// Interior region - conductor
Physical Surface("Conductor 1", 1) = {1};

// Exterior region - surrounding space
Physical Surface("Insulating region", 2) = {2};

// Interface curves between conductor and space
Physical Curve("Conductor Boundary 1", 10) = {1, 2, 3, 4};

// Outer boundary (for boundary conditions)
Physical Curve("True Boundary", 20) = {5, 6, 7, 8};


Mesh.MshFileVersion = 4.5;
//+
Show "*";
