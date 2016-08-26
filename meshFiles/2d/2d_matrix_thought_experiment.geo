meshSize       = 2;
x_start = 0;
y_start = 0;
x_end = 120/2.;
y_end = 20/2.;
delta = 5;

// create points
Point(1)  = {x_start    ,   y_start    	,   0,  meshSize};
Point(2)  = {x_end    	,   y_start    	,   0,  meshSize};
Point(3)  = {x_end      ,   y_end - delta   	,   0,  meshSize};
Point(4)  = {x_start   	,	y_end    	,   0,  meshSize};


// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,1};


loop1 = newll; Line Loop(loop1) = {l1,l2,l3,l4};

plane1 = news; Plane Surface(plane1) = {loop1};

Rotate {{0, 1, 0}, {x_end, 0, 0}, Pi} {
  Duplicata { Surface{plane1}; }
}

Rotate {{1, 0, 0}, {x_end, 0, 0}, Pi} {
  Duplicata { Surface{6,7}; }
}


Physical Surface(999) = {6,-7,-12,17};
