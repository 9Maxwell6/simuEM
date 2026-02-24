rcloque = 1.5;
hcloque = 0.;//laisser à 0 (la cloque est faite dans le mehser Z7)
htgo = 0.005;
hsubstrat = 1;
hceramique = 0.150;
hbc = 0.05;
rayon_TMF = 5.5;
l_substrat_x = Pi/2 * 5.5;
l_substrat_z = 15.;
l_rafinement = 0.3;
l_rafinement_cloque = l_rafinement;
l_rafinement_fissure = l_rafinement;
lc_rafinement = htgo;
nelts_tgo = 2;
n_layers_revolution = Round(1.5*Pi / 2 / (lc_rafinement)); // number of layers for the revolution mesh
//for the heat interface meshing
n_rough = 10;
step_rough = (rcloque - l_rafinement) / n_rough;

// parameters for the crack mesh refinement
lcmin = htgo*1.5;
lcmax = hsubstrat/3;
dmin = 0;
dmax = rcloque;
SetFactory("OpenCASCADE");

Point(1) = {0.,-hbc-htgo,0.};
Point(2) = {0.,-hsubstrat-hbc-htgo,0.};
Point(3)={0.,-htgo,0.};
Point(4)={0.,0.,0.};
Point(5)={0.,hceramique+hcloque,0.};

Line(1)= {2,1};
Line(2)={1,3};
Line(3)={3,4};
Line(4)={4,5};
line_id = 3;
For i In {1:n_rough}
  //Printf("'%g'",line_id);
  vector[] = Extrude{step_rough,0.,0.}{Line{line_id};};
  line_id = vector[0];
EndFor

Extrude {rcloque-l_rafinement_cloque,0,0}{Line{1,2,4};}
Coherence;

Extrude {l_rafinement_cloque,0,0}{Line{37,42,40,45};}
Extrude {l_rafinement_fissure,0,0}{Line{48,50,52,54};}

Translate {l_substrat_x, 0, 0} { Duplicata{ Line{1,41,3,44}; } }
Extrude {0,0,-l_substrat_z}{Line{64,65,66,67};}


Extrude { {0,1,0} , {0,0,0} , Pi/2 } { 
   Surface{1:15,17:19,21}; }
Extrude { {0,1,0} , {0,0,0} , Pi/2 } {
   Surface{16,20};Layers{n_layers_revolution};Recombine;}
Coherence;
Extrude {-l_substrat_x,0,0}{Surface{22:25};}
Coherence;


// BooleanFragments{Volume{13:16};Delete;}{Volume{1:12};}
Transfinite Curve{160,197,199} = nelts_tgo+1;
Transfinite Curve{172,177} = Round(l_rafinement_cloque/lc_rafinement);
Transfinite Curve{185,190} = Round(l_rafinement_fissure/lc_rafinement);
Transfinite Surface{121,124};
Recombine Surface{121,124};
//MESH REFINEMENT ARROUND THE CRACK
// Field arround the crack surface
Field[1] = Distance;
Field[1].SurfacesList = {106,115};
Field[1].Sampling = 1000;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lcmin;
Field[2].SizeMax = lcmax;
Field[2].DistMin = dmin;
Field[2].DistMax = dmax;

//céramique
Field[3] = Constant;
Field[3].VIn = hceramique/2;
Field[3].VolumesList = {13,16,19,25};

//céramique cloque
Field[6] = Constant;
Field[6].VIn = hceramique/4;
Field[6].VolumesList = {13};

//bond coat
Field[4] = Constant;
Field[4].VIn = hbc/2;
Field[4].VolumesList = {12,15,18,23};

//TGO
Field[5] = Constant;
Field[5].VIn = hbc/2;
Field[5].VolumesList = {1:10,20,21,24};


Field[99] = Min;
Field[99].FieldsList = {2,3,4,5,6};
Background Field = 99;


Mesh.Algorithm = 6;
Mesh.Algorithm3D = 1;
Mesh.OptimizePyramids = 1;
Mesh 3;
Save "gmsh_mesh3D_5_optpy.msh";
