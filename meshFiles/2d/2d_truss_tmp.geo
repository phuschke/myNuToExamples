meshSize       = 100;
x_start = 88.001;
y_start = 10.01;
x_end 	= 98.001;
y_end 	= 10.01;

// create points
Point(1)  = {x_start    ,   y_start    	,   0,  meshSize}; 
Point(2)  = {x_end    	,   y_end    	,   0,  meshSize}; 



// create lines
l1  = newreg; Line(l1)  = {1,2};
Transfinite Line {l1} = 10;
Physical Line(777) = {l1};
