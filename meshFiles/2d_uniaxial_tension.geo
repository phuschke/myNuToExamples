meshSize       = 10;
notchHeight    = 2;
notchWidth     = 4;
//*********************************************************************************************
Function CreateNotchedRectangle
// create points
Point(1)  = {x_origin       	        ,   y_origin    ,   0,  meshSize}; 
Point(2)  = {x_origin +100              ,   y_origin    ,   0,  meshSize}; 
Point(3)  = {x_origin +150              ,   y_origin +notchHeight    ,   0,  meshSize};
Point(4)  = {x_origin +200              ,   y_origin    ,   0,  meshSize}; 

Point(5)  = {x_end              		,   y_origin	,   0,  meshSize}; 
Point(6)  = {x_end              		,   y_end	,   0,  meshSize}; 

Point(7)  = {x_end -100		            ,   y_end    ,   0,  meshSize}; 
Point(8)  = {x_end -150		            ,   y_end - notchHeight    ,   0,  meshSize}; 
Point(9)  = {x_end -200		            ,   y_end    ,   0,  meshSize}; 
Point(10) = {x_origin	                ,   y_end    ,   0,  meshSize}; 

// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,5};
l5  = newreg; Line(l5)  = {5,6};
l6  = newreg; Line(l6)  = {6,7};
l7  = newreg; Line(l7)  = {7,8};
l8  = newreg; Line(l8)  = {8,9};
l9  = newreg; Line(l9)  = {9,10};
l10  = newreg; Line(l10)  = {10,1};

loop1 = newll; Line Loop(loop1) = {l1:l10};

plane1 = news; Plane Surface(plane1) = {loop1};



Return

// first rectangle
x_origin  =  0;
y_origin  =  0;
x_end     =  300;
y_end     =  20;
index     =    0;
Call CreateNotchedRectangle;


Physical Surface(999) = {12};



