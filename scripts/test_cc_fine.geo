// Cube-in-Cube Geometry for Electromagnetic Simulation
// Interior: Conductor region
// Exterior: Surrounding space (air/vacuum)
//
// Mesh refinement strategy:
//   - VERY FINE mesh on the surface of the 1x1x1 inner cube
//   - Coarser mesh inside the inner cube (away from its surface)
//   - Coarser mesh outside the inner cube (in the surrounding space)
// Implemented with a Distance + Threshold size field.

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

// Mesh size targets
lc_fine   = 0.03;   // VERY fine: at/near the inner cube surface
lc_coarse = 0.1;   // Coarser: inside the inner cube and far outside

// Transition distances from the inner cube surface
dist_min = 0.0;     // Up to this distance from surface -> lc_fine
dist_max = 0.3;     // Beyond this distance from surface -> lc_coarse
                    // (smooth linear interpolation in between)

// The Point() lc values below are mostly fallbacks; the size field
// will override them in the region where the field applies.

// ============================================================
// Inner Cube (Conductor)
// ============================================================
hi = inner_size / 2;

// Inner cube vertices - tag with the fine size as a fallback
Point(1) = {inner_x - hi, inner_y - hi, inner_z - hi, lc_fine};
Point(2) = {inner_x + hi, inner_y - hi, inner_z - hi, lc_fine};
Point(3) = {inner_x + hi, inner_y + hi, inner_z - hi, lc_fine};
Point(4) = {inner_x - hi, inner_y + hi, inner_z - hi, lc_fine};
Point(5) = {inner_x - hi, inner_y - hi, inner_z + hi, lc_fine};
Point(6) = {inner_x + hi, inner_y - hi, inner_z + hi, lc_fine};
Point(7) = {inner_x + hi, inner_y + hi, inner_z + hi, lc_fine};
Point(8) = {inner_x - hi, inner_y + hi, inner_z + hi, lc_fine};

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

// Outer cube vertices - tag with the coarse size as a fallback
Point(9 ) = {outer_x - ho, outer_y - ho, outer_z - ho, lc_coarse};
Point(10) = {outer_x + ho, outer_y - ho, outer_z - ho, lc_coarse};
Point(11) = {outer_x + ho, outer_y + ho, outer_z - ho, lc_coarse};
Point(12) = {outer_x - ho, outer_y + ho, outer_z - ho, lc_coarse};
Point(13) = {outer_x - ho, outer_y - ho, outer_z + ho, lc_coarse};
Point(14) = {outer_x + ho, outer_y - ho, outer_z + ho, lc_coarse};
Point(15) = {outer_x + ho, outer_y + ho, outer_z + ho, lc_coarse};
Point(16) = {outer_x - ho, outer_y + ho, outer_z + ho, lc_coarse};

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
Surface Loop(2) = {7, 8, 9, 10, 11, 12};
Volume(2) = {2, 1};  // outer shell minus inner cube (hole)

// ============================================================
// Mesh size field: refine near the inner cube surface
// ============================================================
// Field 1: distance from any point in the domain to the 6 inner cube faces
Field[1] = Distance;
Field[1].SurfacesList = {1, 2, 3, 4, 5, 6};
Field[1].Sampling = 100;  // sampling points per surface for distance evaluation

// Field 2: map distance -> target mesh size
//   distance <= DistMin            -> SizeMin (very fine, at the surface)
//   distance >= DistMax            -> SizeMax (coarse, far away)
//   DistMin < distance < DistMax   -> smooth linear interpolation
Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_fine;
Field[2].SizeMax  = lc_coarse;
Field[2].DistMin  = dist_min;
Field[2].DistMax  = dist_max;

// Use field 2 as the global background mesh size
Background Field = 2;

// Make the field authoritative (don't let point sizes or boundary
// extension override the gradient we just defined).
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;

// ============================================================
// Physical Groups (for FEM solver)
// ============================================================
Physical Volume("Conductor 1", 1) = {1};
Physical Volume("Insulating region", 2) = {2};
Physical Surface("Conductor Boundary 1", 10) = {1, 2, 3, 4, 5, 6};
Physical Surface("True Boundary", 20) = {7, 8, 9, 10, 11, 12};


Mesh.MshFileVersion = 4.5;
//+
Show "*";
