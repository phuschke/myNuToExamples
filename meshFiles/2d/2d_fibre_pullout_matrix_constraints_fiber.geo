meshSize       = 1;

// create points
Point(1)  = {5.222      ,   5.222    	,   0,  meshSize}; 
Point(2)  = {9.222    	,   5.222    	,   0,  meshSize}; 

// create lines
l1  = newreg; Line(l1)  = {1,2};

Physical Line(777) = {l1};



