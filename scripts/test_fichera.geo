// =====================================================================
// Case 3 - Reentrant-corner (Fichera) conductor in cube domain
//   conductor : [-0.5, 0.5]^3  minus the corner octant [0, 0.5]^3
//   domain    : [-1.5, 1.5]^3
//
// The conductor boundary has 9 faces: the six (notched) outer faces in
// the planes x,y,z = +-0.5 plus the three reentrant faces in the planes
// x,y,z = 0. All belong to physical surface 10.
//
// Physical tags:
//   Volume   1   "conductor"
//   Volume   100 "air"              (includes the notch octant)
//   Surface  10  "conductor_boundary"   -> essential BC  T x n = 0
//   Surface  20  "outer_boundary"       -> essential BC  Omega = 0
//
// Usage:
//   gmsh -3 case3_fichera_in_cube.geo -setnumber h 0.10 -o case3_h0p10.msh
// =====================================================================
SetFactory("OpenCASCADE");

If(!Exists(h))  h = 0.20;  EndIf   // override on CLI:  -setnumber h <value>
If(!Exists(rc)) rc = 1.0;  EndIf   // conductor refinement: size h/rc inside, override: -setnumber rc <value>

// ----- geometry ------------------------------------------------------
Box(1) = {-1.5, -1.5, -1.5,  3, 3, 3};       // domain
Box(2) = {-0.5, -0.5, -0.5,  1, 1, 1};       // full unit cube
Box(3) = { 0.0,  0.0,  0.0,  0.5, 0.5, 0.5}; // corner octant to remove

c() = BooleanDifference{ Volume{2}; Delete; }{ Volume{3}; Delete; };  // Fichera solid
out() = BooleanFragments{ Volume{1}; Delete; }{ Volume{c()}; Delete; };

// ----- entity identification ----------------------------------------
// conductor bbox is exactly [-0.5,0.5]^3 ; the air volume (which owns the
// notch octant) has the full domain bbox, so it is excluded automatically
eps = 1e-3;
cond() = Volume In BoundingBox{-0.5-eps, -0.5-eps, -0.5-eps,  0.5+eps, 0.5+eps, 0.5+eps};
vAll() = Volume In BoundingBox{-1.5-eps, -1.5-eps, -1.5-eps,  1.5+eps, 1.5+eps, 1.5+eps};
air()  = vAll();
air() -= cond();

// all 9 conductor faces (outer notched + 3 reentrant) lie within +-0.5
scond() = Surface In BoundingBox{-0.5-eps, -0.5-eps, -0.5-eps,  0.5+eps, 0.5+eps, 0.5+eps};

L = 1.5;
souter()  = Surface In BoundingBox{-L-eps, -L-eps, -L-eps,  -L+eps,  L+eps,  L+eps};
souter() += Surface In BoundingBox{ L-eps, -L-eps, -L-eps,   L+eps,  L+eps,  L+eps};
souter() += Surface In BoundingBox{-L-eps, -L-eps, -L-eps,   L+eps, -L+eps,  L+eps};
souter() += Surface In BoundingBox{-L-eps,  L-eps, -L-eps,   L+eps,  L+eps,  L+eps};
souter() += Surface In BoundingBox{-L-eps, -L-eps, -L-eps,   L+eps,  L+eps, -L+eps};
souter() += Surface In BoundingBox{-L-eps, -L-eps,  L-eps,   L+eps,  L+eps,  L+eps};

Physical Volume("conductor", 1)            = {cond()};
Physical Volume("air", 100)                = {air()};
Physical Surface("conductor_boundary", 10) = {scond()};
Physical Surface("outer_boundary", 20)     = {souter()};

// ----- mesh size: uniform h, optionally h/rc near the conductor ------
// (the field box covers the notch too - harmless slight refinement of air)
Field[1] = Box;
Field[1].VIn  = h / rc;
Field[1].VOut = h;
Field[1].XMin = -0.5; Field[1].XMax = 0.5;
Field[1].YMin = -0.5; Field[1].YMax = 0.5;
Field[1].ZMin = -0.5; Field[1].ZMax = 0.5;
Field[1].Thickness = 0.5;
Background Field = 1;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeMax = h;
Mesh.ElementOrder = 1;
Mesh.Optimize = 1;
