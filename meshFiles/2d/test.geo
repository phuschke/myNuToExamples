// Commentaires

Mesh.Partitioner = 1;
lc = 0.3;   // longueur caractéristique

// Définition des Points

Point(1) = {-1,1,0,lc};
Point(2) = {1,1,0,lc};
Point(3) = {1,-1,0,lc};
Point(4) = {-1,-1,0,lc};

// Définition des Lignes

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// Définition de la Surface

Line Loop(5) = {1,2,3,4};
Physical Line(1) = {1,2,3,4};
//Physical Line(2) = {1,3};
//Physical Line(1) = {2,4};

Plane Surface(7) = {5};
Physical Surface(8) = {7};

Transfinite Line{1,2,3,4}=31;
Transfinite Surface(7)={1,2,3,4};

Recombine Surface(7)=0;
