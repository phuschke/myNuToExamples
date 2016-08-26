
a = 9;
h = 2;
b = 5;
c = 10;
delta = 0.1;
meshSize = 9;

Point(1) = {0, 0, 0, meshSize};
Point(2) = {a, 0, 0, meshSize};
Point(3) = {a, h, 0, meshSize};
Point(4) = {0, h, 0, meshSize};

Line(1) = {3, 4};
Line(2) = {1, 4};
Line(3) = {2, 1};
Line(4) = {3, 2};
	
Transfinite Line {1,3} = 4;
Transfinite Line {2,4} = 3;

Point(5) = {3, 1, 0, meshSize};
Point(6) = {6, 1, 0, meshSize};

Line(5) = {5, 6};

Line Loop(18) = {1, -2, -3, -4};
Plane Surface(18) = {18};
Line {5} In Surface {18};
Recombine Surface {18};
Transfinite Surface {18};


Physical Surface(777) = {18};
Physical Line(999) = {5};