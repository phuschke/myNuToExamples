//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
Transfinite Line{1} = 4 Using Progression 2;
Transfinite Line{3} = 4 Using Progression 0.5;
//+
Line Loop(5) = {3, 4, 1, 2};
//+
Plane Surface(6) = {5};
Transfinite Surface{6};
Recombine Surface{6};
