height 		= 10;
width 		= 10;
meshSize 	= 2;

Point(1) = {0, 		0, 			0, meshSize};
Point(2) = {width, 	0, 			0, meshSize};
Point(3) = {width, 	height, 	0, meshSize};
Point(4) = {0, 		height, 	0, meshSize};

Line(1) = {3, 4};
Line(2) = {1, 4};
Line(3) = {2, 1};
Line(4) = {3, 2};
	
Point(5) = {4, 4, 0, meshSize};
Point(6) = {6, 6, 0, meshSize};
Line(5) = {5, 6};
Transfinite Line {5} = 2;

Line Loop(18) = {1, -2, -3, -4};
Plane Surface(18) = {18};
//Recombine Surface {18};

Point {5,6} In Surface {18};
Line {5} In Surface {18};

Physical Surface(222) = {18};
Physical Line(999) = {5,6};
