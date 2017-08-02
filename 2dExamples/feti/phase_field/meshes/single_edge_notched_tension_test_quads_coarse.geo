meshSize       =0.05;
offset = 0.01;
// create points
Point(1)  = {0      ,     0    	,   0,  meshSize};
Point(2)  = {1      ,     0   	,   0,  meshSize};
Point(3)  = {1      ,     1   	,   0,  meshSize};
Point(4)  = {0   	  ,	    1   	,   0,  meshSize};

Point(5)  = {0   	  ,	    0.500001   	,   0,  meshSize};
Point(6)  = {0.5   	  ,	    0.5   	,   0,  meshSize};
Point(7)  = {0   	  ,	    0.499999   	,   0,  meshSize};

Point(8)  = {0   	  ,	    0.52   	,   0,  meshSize};
Point(9)  = {0.5   	  ,	    0.52   	,   0,  meshSize};
Point(10)  = {1   	  ,	    0.52   	,   0,  meshSize};

Point(11)  = {0   	  ,	    0.48   	,   0,  meshSize};
Point(12)  = {0.5   	  ,	    0.48   	,   0,  meshSize};
Point(13)  = {1   	  ,	    0.48   	,   0,  meshSize};

Point(14)  = {0.5   	  ,	    1   	,   0,  meshSize};
Point(15)  = {0.5   	  ,	    0   	,   0,  meshSize};
//+
Line(1) = {1, 15};
//+
Line(2) = {15, 2};
//+
Line(3) = {4, 14};
//+
Line(4) = {14, 3};
//+
Line(5) = {8, 9};
//+
Line(6) = {9, 10};
//+
Line(7) = {11, 12};
//+
Line(8) = {12, 13};
//+
Line(9) = {5, 6};
//+
Line(10) = {7, 6};
//+
Line(11) = {1, 11};
//+
Line(12) = {15, 12};
//+
Line(13) = {2, 13};
//+
Line(15) = {12, 6};
//+
Line(16) = {6, 9};
//+
Line(17) = {11, 7};
//+
Line(18) = {5, 8};
//+
Line(19) = {8, 4};
//+
Line(20) = {9, 14};
//+
Line(21) = {10, 3};
//+
Point(16) = {1, 0.5, 0, 1.0};
//+
Line(22) = {13, 16};
//+
Line(23) = {16, 10};
//+
Line(24) = {6, 16};
//+
Line Loop(25) = {19, 3, -20, -5};
//+
Plane Surface(26) = {25};
//+
Line Loop(27) = {20, 4, -21, -6};
//+
Plane Surface(28) = {27};
//+
Line Loop(29) = {6, -23, -24, 16};
//+
Plane Surface(30) = {29};
//+
Line Loop(31) = {15, 24, -22, -8};
//+
Plane Surface(32) = {31};
//+
Line Loop(33) = {18, 5, -16, -9};
//+
Plane Surface(34) = {33};
//+
Line Loop(35) = {10, -15, -7, 17};
//+
Plane Surface(36) = {35};
//+
Line Loop(37) = {11, 7, -12, -1};
//+
Plane Surface(38) = {37};
//+
Line Loop(39) = {12, 8, -13, -2};
//+
Plane Surface(40) = {39};

Transfinite Line {19,20,21} = 20 Using Progression 1.1;
Transfinite Line {11,12,13} = 20 Using Progression 1/1.1;

Transfinite Line {17,15,22} = 4 Using Progression 1/1.01;
Transfinite Line {18,16,23} = 4 Using Progression 1.01;

Transfinite Line {2,4,6,8,24} = 50 Using Progression 1.01;
Transfinite Line {1,3,5,7,9,10} = 50 Using Progression 1/1.01;

Recombine Surface {26,28,30,32,34,36,38,40};
Transfinite Surface {26,28,30,32,34,36,38,40};

Physical Surface(9999) = {26,28,30,32,34,36,38,40};
