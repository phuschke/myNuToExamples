height 		= 10;
width 		= 10;
meshSize 	= 1;

Point(1) = {0, 		0, 			0, meshSize};
Point(2) = {width, 	0, 			0, meshSize};
Point(3) = {width, 	height, 	0, meshSize};
Point(4) = {0, 		height, 	0, meshSize};

Line(1) = {3, 4};
Line(2) = {1, 4};
Line(3) = {2, 1};
Line(4) = {3, 2};
	
Point(5) = {2, 4, 0, meshSize};
Point(6) = {4, 4, 0, meshSize};
Line(5) = {5, 6};

Point(7) = {2, 6, 0, meshSize};
Point(8) = {4, 6, 0, meshSize};
Line(6) = {7, 8};

Point(9) = {2, 8, 0, meshSize};
Point(10) = {4, 8, 0, meshSize};
Line(7) = {9, 10};

Point(11) = {2, 2, 0, meshSize};
Point(12) = {4, 2, 0, meshSize};
Line(8) = {11, 12};



Point(13) = {6, 4, 0, meshSize};
Point(14) = {8, 4, 0, meshSize};
Line(9) = {13, 14};

Point(15) = {6, 6, 0, meshSize};
Point(16) = {8, 6, 0, meshSize};
Line(10) = {15, 16};

Point(17) = {6, 8, 0, meshSize};
Point(18) = {8, 8, 0, meshSize};
Line(11) = {17, 18};

Point(19) = {6, 2, 0, meshSize};
Point(20) = {8, 2, 0, meshSize};
Line(12) = {19, 20};










Transfinite Line {5:12} = 3;

Line Loop(18) = {1, -2, -3, -4};
Plane Surface(18) = {18};
//Recombine Surface {18};
//Transfinite Surface {18};
Line {5:12} In Surface {18};

Physical Surface(222) = {18};
Physical Line(999) = {5:12};
