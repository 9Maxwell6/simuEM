// =====================================================================
//  sphere_with_torus.geo
//  Sphere (diameter 7) with a torus centred at the origin.
//      torus ring (major) diameter = 3  -> major radius 1.5
//      torus tube (minor) diameter = 1  -> minor radius 0.5
//  Torus is a separate, conformal region inside the sphere.
// =====================================================================
SetFactory("OpenCASCADE");

// -------------------- parameters (edit here) -------------------------
sphere_dia     = 7;     // outer sphere diameter
torus_ring_dia = 3;     // diameter of the torus centre-ring  (major)
torus_tube_dia = 1;     // diameter of the tube cross-section (minor)

Rs   = sphere_dia/2;        // 3.5
Rmaj = torus_ring_dia/2;    // 1.5
Rmin = torus_tube_dia/2;    // 0.5

lc_torus  = 0.088388348;   // element size near/inside the torus
lc_sphere = 0.088388348;   // element size in the bulk of the sphere
eps = 1e-3;

// -------------------- geometry ---------------------------------------
Torus(1)  = {0, 0, 0, Rmaj, Rmin};   // axis = z, ring in the xy-plane
Sphere(2) = {0, 0, 0, Rs};

// Fragment -> shared interface, conformal mesh, two volumes
vol() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; };

// -------------------- tag the regions robustly -----------------------
// The torus fits inside this small box; the surrounding region does not.
qxy = (Rmaj + Rmin + Rs) / 2;   // = 2.75, between torus extent (2.0) and Rs (3.5)
qz  = (Rmin        + Rs) / 2;   // = 2.0

torus_vol() = Volume  In BoundingBox{ -qxy, -qxy, -qz,  qxy, qxy, qz };
tsurf()     = Surface In BoundingBox{ -qxy, -qxy, -qz,  qxy, qxy, qz };

If (vol(0) == torus_vol(0))
  matrix_vol = vol(1);
Else
  matrix_vol = vol(0);
EndIf

// -------------------- physical groups --------------------------------
Physical Volume ("Conductor 1",  1) = { torus_vol(0) };
Physical Volume ("Insulating region", 2) = { matrix_vol  };

// outer sphere boundary = matrix boundary that is not the torus interface
mb() = Abs( Boundary{ Volume{matrix_vol}; } );
outer() = {};
For i In { 0 : #mb()-1 }
  hit = 0;
  For j In { 0 : #tsurf()-1 }
    If (mb(i) == tsurf(j))
      hit = 1;
    EndIf
  EndFor
  If (hit == 0)
    outer() += mb(i);
  EndIf
EndFor
Physical Surface ("True Boundary",    3) = outer();
Physical Surface ("Conductor Boundary 1", 4) = tsurf();

// -------------------- mesh size control ------------------------------
Field[1] = Distance;
Field[1].SurfacesList = { tsurf() };
Field[1].Sampling = 200;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = lc_torus;
Field[2].SizeMax  = lc_sphere;
Field[2].DistMin  = 0.2;
Field[2].DistMax  = 1.5;

Background Field = 2;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;
Mesh.Algorithm   = 6;   // 2D: Frontal-Delaunay
Mesh.Algorithm3D = 1;   // 3D: Delaunay