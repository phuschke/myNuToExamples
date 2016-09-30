meshSize       =0.2;
x_start = 0;
y_start = 0;
eps = 0.00001;
factor = 0.1;

// create points
Point(1)  = {0      ,     0    	,   0,  meshSize};
Point(2)  = {1      ,     0   	,   0,  meshSize};
Point(3)  = {1      ,     1   	,   0,  meshSize};
Point(4)  = {0   	  ,	    1   	,   0,  meshSize};

Point(5)  = {0   	  ,	    0.5+eps   	,   0,  meshSize};
Point(6)  = {0.5   	  ,	    0.5   	,   0,  factor*meshSize};
Point(7)  = {0   	  ,	    0.5-eps   	,   0,  meshSize};

Point(8)  = {0.5   	  ,	    0   	,   0,  meshSize};
Point(9)  = {1   	  ,	    0.5   	,   0,  factor*meshSize};
Point(10)  = {0.5   	  ,	    1   	,   0,  meshSize};

Line(1) = {1, 8};
Line(2) = {8, 6};

Line(3) = {6, 7};
Line(4) = {7, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
Line(7) = {8, 2};
Line(8) = {2, 9};
Line(9) = {9, 6};
Line Loop(10) = {9, -2, 7, 8};
Plane Surface(11) = {10};
Line(12) = {9, 3};
Line(13) = {3, 10};
Line(14) = {10, 6};
Line Loop(15) = {14, -9, 12, 13};
Plane Surface(16) = {15};
Line(17) = {10, 4};
Line(18) = {4, 5};
Line(19) = {5, 6};
Line Loop(20) = {18, 19, -14, 17};
Plane Surface(21) = {20};

Physical Surface(999) = {6,11,16,21};
//Transfinite Line{2} = 10 Using Progression 1.1;
//Transfinite Line{8,12} = 10 Using Progression 1.1;
//Transfinite Line{9} = 100;
