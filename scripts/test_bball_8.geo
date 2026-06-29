
SetFactory("OpenCASCADE");

If(!Exists(h))  h = 0.03125;  EndIf   
If(!Exists(rc)) rc = 1;  EndIf  

// ----- geometry ------------------------------------------------------
Sphere(1) = {0, 0, 0, 1.5};   // domain
Sphere(2) = {0, 0, 0, 0.5};   // conductor

out() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; };

// ----- entity identification ----------------------------------------
eps = 1e-3;
cond() = Volume In BoundingBox{-0.5-eps, -0.5-eps, -0.5-eps,  0.5+eps, 0.5+eps, 0.5+eps};
vAll() = Volume In BoundingBox{-1.5-eps, -1.5-eps, -1.5-eps,  1.5+eps, 1.5+eps, 1.5+eps};
air()  = vAll();
air() -= cond();

scond() = Surface In BoundingBox{-0.5-eps, -0.5-eps, -0.5-eps,  0.5+eps, 0.5+eps, 0.5+eps};
sAll()  = Surface In BoundingBox{-1.5-eps, -1.5-eps, -1.5-eps,  1.5+eps, 1.5+eps, 1.5+eps};
souter()  = sAll();
souter() -= scond();

Physical Volume("Conductor 1", 1)            = {cond()};
Physical Volume("Insulating region", 100)                = {air()};
Physical Surface("Conductor Boundary 1", 10) = {scond()};
Physical Surface("True Boundary", 20)     = {souter()};

// ----- mesh size: uniform h, optionally h/rc in & near the conductor -
Field[1] = Ball;
Field[1].VIn  = h / rc;
Field[1].VOut = h;
Field[1].Radius  = 0.5;
Field[1].Thickness = 0.5;
Field[1].XCenter = 0; Field[1].YCenter = 0; Field[1].ZCenter = 0;
Background Field = 1;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeMax = h;
Mesh.ElementOrder = 1;
Mesh.Optimize = 1;
