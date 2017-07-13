meshSize       =0.01;
x_start = 0;
y_start = 0;
eps = 0.01;
factor = 1;

// create points
Point(1)  = {0      ,     0    	,   0,  meshSize};
Point(2)  = {1      ,     0   	,   0,  meshSize};
Point(3)  = {1      ,     1   	,   0,  meshSize};
Point(4)  = {0   	  ,	    1   	,   0,  meshSize};

Point(5)  = {0   	  ,	    0.5+eps   	,   0,  meshSize};
Point(6)  = {0.5   	  ,	    0.5+eps   	,   0,  factor*meshSize};
Point(7)  = {0.5   	  ,	    0.5-eps   	,   0,  factor*meshSize};
Point(8)  = {0   	  ,	    0.5-eps   	,   0,  meshSize};

Point(9)  = {0.5   	  ,	    0   	,   0,  meshSize};

Point(11)  = {0.5   	  ,	    1   	,   0,  meshSize};

Point(12)  = {1   	  ,	    0.5+eps   	,   0,  factor*meshSize};
Point(13)  = {1   	  ,	    0.5-eps   	,   0,  factor*meshSize};
//+
Line(1) = {1, 9};
//+
Line(2) = {9, 7};
//+
Line(3) = {7, 8};
//+
Line(4) = {8, 1};
//+
Line(5) = {9, 2};
//+
Line(6) = {2, 13};
//+
Line(7) = {13, 7};
//+
Line(8) = {7, 6};
//+
Line(9) = {6, 5};
//+
Line(10) = {13, 12};
//+
Line(11) = {12, 3};
//+
Line(12) = {3, 11};
//+
Line(13) = {11, 4};
//+
Line(14) = {4, 5};
//+
Line(15) = {6, 12};
//+
Line(16) = {6, 11};
//+
Line Loop(17) = {1, 2, 3, 4};
//+
Plane Surface(18) = {17};
//+
Line Loop(19) = {2, -7, -6, -5};
//+
Plane Surface(20) = {-19};
//+
Line Loop(21) = {8, 15, -10, 7};
//+
Plane Surface(22) = {-21};
//+
Line Loop(23) = {16, -12, -11, -15};
//+
Plane Surface(24) = {-23};
//+
Line Loop(25) = {13, 14, -9, 16};
//+
Plane Surface(26) = {25};

Recombine Surface {18,20,22,24,26};
Transfinite Surface {18,20,22,24,26};
Physical Surface(999) = {18,20,22,24,26};
