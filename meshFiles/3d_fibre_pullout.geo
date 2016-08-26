meshSize       = 5;

x_max = 10;
y_max = 10;
z_max = 10;

// create points
p1 = newp; Point(p1)  = { 0			,   	0   	,   0		,  meshSize}; 
p2 = newp; Point(p2)  = { x_max		,   	0   	,   0		,  meshSize}; 
p3 = newp; Point(p3)  = { x_max		,   	y_max   ,   0		,  meshSize}; 
p4 = newp; Point(p4)  = { 0			,		y_max   ,   0		,  meshSize}; 

p5 = newp; Point(p5)  = { 0			,   	0   	,   z_max	,  meshSize}; 
p6 = newp; Point(p6)  = { x_max		,   	0   	,   z_max	,  meshSize}; 
p7 = newp; Point(p7)  = { x_max		,   	y_max   ,   z_max	,  meshSize}; 
p8 = newp; Point(p8)  = { 0			,		y_max   ,   z_max	,  meshSize}; 


// create lines
l1  = newl; Line(l1)  = {p1,p2};
l2  = newl; Line(l2)  = {p2,p3};
l3  = newl; Line(l3)  = {p3,p4};
l4  = newl; Line(l4)  = {p4,p1};
l5  = newl; Line(l5)  = {p5,p6};
l6  = newl; Line(l6)  = {p6,p7};
l7  = newl; Line(l7)  = {p7,p8};
l8  = newl; Line(l8)  = {p8,p5};
lv1  = newl; Line(lv1)  = {p1,p5};
lv2  = newl; Line(lv2)  = {p2,p6};
lv3  = newl; Line(lv3)  = {p3,p7};
lv4  = newl; Line(lv4)  = {p4,p8};

Line Loop(14) = {1, 10, -5, -9};
Plane Surface(14) = {14};
Line Loop(16) = {4, 9, -8, -12};
Plane Surface(16) = {16};
Line Loop(18) = {8, 5, 6, 7};
Plane Surface(18) = {18};
Line Loop(20) = {6, -11, -2, 10};
Plane Surface(20) = {20};
Line Loop(22) = {1, 2, 3, 4};
Plane Surface(22) = {22};
Line Loop(24) = {11, 7, -12, -3};
Plane Surface(24) = {24};
Surface Loop(26) = {20, 18, 16, 22, 14, 24};
Volume(26) = {26};


Physical Volume (999) = {26};