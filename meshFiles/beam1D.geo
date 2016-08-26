meshSize       = 10;
x_start = 0;
x_end = 200;

// create points
Point(1)  = {x_start    ,   0    	,   0,  meshSize}; 
Point(2)  = {x_end    	,   0    	,   0,  meshSize}; 



// create lines
l1  = newreg; Line(l1)  = {1,2};

Physical Line(999) = {l1};
