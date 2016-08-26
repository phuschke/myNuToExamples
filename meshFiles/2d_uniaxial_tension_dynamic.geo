meshSize       = 1;
x_origin  =  0;
y_origin  =  0;
x_end     =  10;
y_end     =  1;

//*********************************************************************************************

// create points
Point(1)  = {x_origin       	        ,   y_origin    ,   0,  meshSize}; 
Point(2)  = {x_end		                ,   y_origin    ,   0,  meshSize}; 
Point(3)  = {x_end		                ,   y_end       ,   0,  meshSize};
Point(4)  = {x_origin                   ,   y_end       ,   0,  meshSize}; 

// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,1};

loop1 = newll; Line Loop(loop1) = {l1:l4};

plane1 = news; Plane Surface(plane1) = {loop1};
Transfinite Surface {plane1};
Recombine Surface {plane1};


Physical Surface(999) = {plane1};



