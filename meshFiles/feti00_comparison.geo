meshSize       = 10;
x_start = 0;
y_start = 0;
x_end = 60;
y_end = 10;

// create points
Point(1)  = {x_start    ,   y_start    	,   0,  meshSize}; 
Point(2)  = {x_end    	,   y_start    	,   0,  meshSize}; 
Point(3)  = {x_end    	,   y_end    	,   0,  meshSize}; 
Point(4)  = {x_start   	,	y_end    	,   0,  meshSize}; 



// create lines
l1  = newreg; Line(l1)  = {1,2};

l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,1};

Transfinite Line {l1,l3} = 301;
Transfinite Line {l2,l4} = 101;


loop1 = newll; Line Loop(loop1) = {l1,l2,l3,l4};

plane1 = news; Plane Surface(plane1) = {loop1};

Transfinite Surface {plane1};
Recombine Surface {plane1};
Physical Surface(111) = {plane1};

