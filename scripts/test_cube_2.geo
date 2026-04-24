// Cube Geometry for test
// Interior: space region

// ============================================================
// Parameters
// ============================================================
// cube
inner_size = 3.0;
inner_x = 0.0;
inner_y = 0.0;
inner_z = 0.0;


// Mesh size
lc_inner = 0.2500;   

// ============================================================
// Cube 
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


Physical Volume("space", 1) = {1};


Physical Surface("True Boundary 1", 10) = {1, 2, 3, 4, 5, 6};



Mesh.MshFileVersion = 2.2;
//+
Show "*";
