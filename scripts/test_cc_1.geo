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
lc_outer = 0.4;   // Coarser mesh in space

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
Point(9 ) = {outer_x - ho, outer_y - ho, outer_z - ho, lc_outer};
Point(10) = {outer_x + ho, outer_y - ho, outer_z - ho, lc_outer};
Point(11) = {outer_x + ho, outer_y + ho, outer_z - ho, lc_outer};
Point(12) = {outer_x - ho, outer_y + ho, outer_z - ho, lc_outer};
Point(13) = {outer_x - ho, outer_y - ho, outer_z + ho, lc_outer};
Point(14) = {outer_x + ho, outer_y - ho, outer_z + ho, lc_outer};
Point(15) = {outer_x + ho, outer_y + ho, outer_z + ho, lc_outer};
Point(16) = {outer_x - ho, outer_y + ho, outer_z + ho, lc_outer};

// Outer cube edges - bottom face
Line(13) = {9 , 10};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12,  9};

// Outer cube edges - top face
Line(17) = {13, 14};
Line(18) = {14, 15};
Line(19) = {15, 16};
Line(20) = {16, 13};

// Outer cube edges - vertical
Line(21) = {9 , 13};
Line(22) = {10, 14};
Line(23) = {11, 15};
Line(24) = {12, 16};

// Outer cube faces
Curve Loop(7) = {13, 14, 15, 16};    // Bottom (-z)
Plane Surface(7) = {7};

Curve Loop(8) = {17, 18, 19, 20};    // Top (+z)
Plane Surface(8) = {8};

Curve Loop(9) = {13, 22, -17, -21};  // Front (-y)
Plane Surface(9) = {9};

Curve Loop(10) = {15, 24, -19, -23};  // Back (+y)
Plane Surface(10) = {10};

Curve Loop(11) = {14, 23, -18, -22};  // Right (+x)
Plane Surface(11) = {11};

Curve Loop(12) = {16, 21, -20, -24};  // Left (-x)
Plane Surface(12) = {12};

// ============================================================
// Exterior Volume (Space between inner and outer cubes)
// ============================================================
// Surface loop for outer boundary
Surface Loop(2) = {7, 8, 9, 10, 11, 12};

// Exterior volume: outer shell minus inner cube (hole)
Volume(2) = {2, 1};  // {outer surface loop, inner surface loop as hole}

// ============================================================
// Physical Groups (for FEM solver)
// ============================================================
// Interior region - conductor
Physical Volume("Conductor 1", 1) = {1};

// Exterior region - surrounding space
Physical Volume("Insulating region", 2) = {2};

// Interface surface between conductor and space
Physical Surface("Conductor Boundary 1", 10) = {1, 2, 3, 4, 5, 6};

// Outer boundary (for boundary conditions)
Physical Surface("True Boundary", 20) = {7, 8, 9, 10, 11, 12};


// ============================================================
// Mesh refinement near all 24 edges
// ============================================================

// Distance from Dirichlet surfaces (outer faces here; add {1:6} if inner is also Dirichlet)
Field[1] = Distance;
Field[1].SurfacesList = {7, 8, 9, 10, 11, 12};
Field[1].Sampling = 20;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.15;   // near boundary
Field[2].SizeMax = 0.4;    // bulk (= lc_outer)
Field[2].DistMin = 0.0;
Field[2].DistMax = 0.5;    // half the gap

Background Field = 2;

Mesh.Algorithm3D   = 10;   // HXT: more isotropic tets, far less prone to boundary-flat tets
Mesh.Optimize      = 1;
Mesh.OptimizeNetgen = 1;   // this step specifically kills sliver tets sitting on surfaces

Mesh.MshFileVersion = 2.2;
//+
Show "*";
