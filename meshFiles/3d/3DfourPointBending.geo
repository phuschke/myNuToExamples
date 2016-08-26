meshSize       = 20;
notchHeight    = 20;
notchWidth     = 20;
//*********************************************************************************************
Function CreateNotchedRectangle
// create points
Point(1)  = {x_origin       	        ,   y_origin    ,   0,  meshSize}; 
Point(2)  = {x_origin + 25              ,   y_origin    ,   0,  meshSize}; 
Point(3)  = {x_origin + 175             ,   y_origin    ,   0,  meshSize}; 
Point(4)  = {x_end/2. -notchWidth/2    	,	y_origin    ,   0,  meshSize}; 

Point(5)  = {x_origin       	        ,   notchHeight    ,   0,  meshSize}; 
Point(6)  = {x_origin + 25              ,   notchHeight    ,   0,  meshSize}; 
Point(7)  = {x_origin + 175             ,   notchHeight    ,   0,  meshSize}; 
Point(8)  = {x_end/2.    	,	notchHeight    ,   0,  meshSize}; 

Point(10)  = {x_origin       	        ,   y_end    ,   0,  meshSize}; 
Point(11)  = {x_origin + 25              ,   y_end    ,   0,  meshSize}; 
Point(12)  = {x_origin + 175             ,   y_end    ,   0,  meshSize}; 
Point(13)  = {x_end/2    	,	y_end    ,   0,  meshSize}; 


// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};

l4  = newreg; Line(l4)  = {5,6};
l5  = newreg; Line(l5)  = {6,7};
l6  = newreg; Line(l6)  = {7,8};


l8  = newreg; Line(l8)  = {10,11};
l9 	= newreg; Line(l9) 	= {11,12};
l10 = newreg; Line(l10) = {12,13};


// create inner vertical lines
lv1 = newreg; Line(lv1) = {1,5};
lv2 = newreg; Line(lv2) = {2,6};
lv3 = newreg; Line(lv3) = {3,7};
lv4 = newreg; Line(lv4) = {4,8};
lv5 = newreg; Line(lv5) = {5,10};
lv6 = newreg; Line(lv6) = {6,11};
lv7 = newreg; Line(lv7) = {7,12};
lv8 = newreg; Line(lv8) = {8,13};


loop1 = newll; Line Loop(loop1) = {l1,lv2,-l4,-lv1};
loop2 = newll; Line Loop(loop2) = {l2,lv3,-l5,-lv2};
loop3 = newll; Line Loop(loop3) = {l3,lv4,-l6,-lv3};

loop4 = newll; Line Loop(loop4) = {l4,lv6,-l8,-lv5};
loop5 = newll; Line Loop(loop5) = {l5,lv7,-l9,-lv6};
loop6 = newll; Line Loop(loop6) = {l6,lv8,-l10,-lv7};



plane1 = news; Plane Surface(plane1) = {loop1};
plane2 = news; Plane Surface(plane2) = {loop2};
plane3 = news; Plane Surface(plane3) = {loop3};
plane4 = news; Plane Surface(plane4) = {loop4};
plane5 = news; Plane Surface(plane5) = {loop5};
plane6 = news; Plane Surface(plane6) = {loop6};





Return
// first rectangle
x_origin  =  0;
y_origin  =  0;
x_end     =  500;
y_end     =  100;
index     =    0;
Call CreateNotchedRectangle;

Rotate {{0, 1, 0}, {x_end/2., 0, 0}, Pi} {
  Duplicata { Surface{plane1:plane6}; }
}

Extrude {0,0,50} 
{
Surface{24:29,30,35,40,45,50,55,62}; 
}

Physical Volume(999) = {1:13};



