meshSize       =0.05;
offset = 0.01;
// create points
Point(1)  = {0      ,     0    	,   0,  meshSize};
Point(2)  = {1      ,     0   	,   0,  meshSize};
Point(3)  = {1      ,     1   	,   0,  meshSize};
Point(4)  = {0   	  ,	    1   	,   0,  meshSize};

Point(5)  = {0.5   	  ,	    0.0    	,   0,  meshSize};
Point(6)  = {1.0   	  ,	    0.5   	,   0,  meshSize*0.05};
Point(7)  = {0.5   	  ,	    1.0   	,   0,  meshSize};
Point(8)  = {0.0   	  ,	    0.5   	,   0,  meshSize*0.05};

Point(9)  = {0.5   	  ,	    0.5   	,   0,  meshSize*0.05};

Point(10)  = {0.0   	  ,	    0.5-offset   	,   0,  meshSize*0.05};
Point(11)  = {0.5   	  ,	    0.5-offset   	,   0,  meshSize*0.05};
Point(12)  = {1.0   	  ,	    0.5-offset   	,   0,  meshSize*0.05};

Point(13)  = {0.0   	  ,	    0.5+offset   	,   0,  meshSize*0.05};
Point(14)  = {0.5   	  ,	    0.5+offset   	,   0,  meshSize*0.05};
Point(15)  = {1.0   	  ,	    0.5+offset   	,   0,  meshSize*0.05};

Line(2) = {1, 5};
Line(3) = {5, 2};
Line(4) = {2, 12};
Line(5) = {12, 6};
Line(6) = {6, 15};
Line(7) = {15, 3};
Line(8) = {3, 7};
Line(9) = {7, 4};
Line(10) = {4, 13};
Line(11) = {13, 8};
Line(12) = {8, 10};
Line(13) = {10, 1};
Line(14) = {5, 11};
Line(15) = {11, 10};
Line(16) = {8, 9};
Line(17) = {9, 6};
Line(18) = {11, 12};
Line(19) = {11, 9};
Line(20) = {9, 14};
Line(21) = {14, 15};
Line(22) = {14, 7};
Line(23) = {13, 14};
Line Loop(24) = {13, 2, 14, 15};
Plane Surface(25) = {24};
Line Loop(26) = {14, 18, -4, -3};
Plane Surface(27) = {26};
Line Loop(28) = {18, 5, -17, -19};
Plane Surface(29) = {28};
Line Loop(30) = {16, -19, 15, -12};
Plane Surface(31) = {30};
Line Loop(32) = {17, 6, -21, -20};
Plane Surface(33) = {32};
Line Loop(34) = {16, 20, -23, 11};
Plane Surface(35) = {34};
Line Loop(36) = {21, 7, 8, -22};
Plane Surface(37) = {36};
Line Loop(38) = {23, 22, 9, 10};
Plane Surface(39) = {38};


Transfinite Surface {31,29,33,35};
Physical Surface(999) = {31,29,33,35,25,27,37,39};
