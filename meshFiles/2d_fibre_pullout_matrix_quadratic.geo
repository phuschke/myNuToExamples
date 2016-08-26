meshSize       = 0.5;
x_start = 0;
y_start = 0;
x_end = 10;
y_end = 10;

// create points
Point(1)  = {x_start    ,   y_start    	,   0,  meshSize}; 
Point(2)  = {x_end    	,   y_start    	,   0,  meshSize}; 
Point(3)  = {x_end    	,   y_end/2. 	,   0,  meshSize}; 
Point(4)  = {x_end    	,   y_end    	,   0,  meshSize}; 
Point(5)  = {x_start   	,	y_end    	,   0,  meshSize}; 


// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,5};
l5  = newreg; Line(l5)  = {5,1};

Point(6) = {x_end/2., y_end/2., 0, meshSize};
l6  = newreg; Line(l6)  = {6,3};





loop1 = newll; Line Loop(loop1) = {l1,l2,l3,l4, l5};

plane1 = news; Plane Surface(plane1) = {loop1};
Line {l6} In Surface {plane1};
Physical Surface(999) = {plane1};
Physical Line(777) = {l6};


