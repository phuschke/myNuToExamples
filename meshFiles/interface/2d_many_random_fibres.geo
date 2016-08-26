

height 		= 10;
width 		= 10;
meshSize 	= 0.5;

Point(1) = {0, 		0, 			0, meshSize};
Point(2) = {width, 	0, 			0, meshSize};
Point(3) = {width, 	height, 	0, meshSize};
Point(4) = {0, 		height, 	0, meshSize};

Line(1) = {3, 4};
Line(2) = {1, 4};
Line(3) = {2, 1};
Line(4) = {3, 2};
	
//	 x1		  y1        x2       y2
//   2.03627, 4.67884   3.01564, 4.47672
//   5.23051, 1.19749   6.23004, 1.22808
//   4.90781, 8.06813   5.88419, 8.28420
//   6.15767, 0.67973   7.15266, 0.77976
//   7.49337, 3.30111   8.46214, 3.05318

//   4.51094, 5.92580   5.49133, 5.72877
//   6.69514, 6.65277   7.69336, 6.71237
//   0.31250, 6.27203   1.30920, 6.35318
//   6.39456, 4.60771   7.34705, 4.30314
//   6.56055, 5.24565   7.55190, 5.11441

// fibre 1
Point(5) = {2.03627, 4.67884, 0, meshSize};
Point(6) = {3.01564, 4.47672, 0, meshSize};
Line(5) = {5, 6};

// fibre 2
Point(7) = {5.23051, 1.19749, 0, meshSize};
Point(8) = {6.23004, 1.22808, 0, meshSize};
Line(6) = {7, 8};

// fibre 3
Point(9) = {4.90781, 8.06813, 0, meshSize};
Point(10) = {5.88419, 8.28420, 0, meshSize};
Line(7) = {9, 10};

// fibre 4
Point(11) = {6.15767, 0.67973, 0, meshSize};
Point(12) = {7.15266, 0.77976, 0, meshSize};
Line(8) = {11, 12};


// fibre 5
Point(13) = {7.49337, 3.30111, 0, meshSize};
Point(14) = {8.46214, 3.05318, 0, meshSize};
Line(9) = {13, 14};

// fibre 6
Point(15) = {4.51094, 5.92580, 0, meshSize};
Point(16) = {5.49133, 5.72877, 0, meshSize};
Line(10) = {15, 16};

// fibre 7
Point(17) = {6.69514, 6.65277, 0, meshSize};
Point(18) = {7.69336, 6.71237, 0, meshSize};
Line(11) = {17, 18};

// fibre 8
Point(19) = {0.31250, 6.27203, 0, meshSize};
Point(20) = {1.30920, 6.35318, 0, meshSize};
Line(12) = {19, 20};




Line Loop(18) = {1, -2, -3, -4};
Plane Surface(18) = {18};
//Recombine Surface {18};
//Transfinite Surface {18};
Line {5:12} In Surface {18};

Physical Surface(222) = {18};
Physical Line(999) = {5:12};
