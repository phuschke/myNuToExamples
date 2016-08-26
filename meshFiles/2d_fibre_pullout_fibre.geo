meshSize       = 100;
x_start = 0;
y_start = 0;
x_end = 15;
y_end = 10;

Point(1)  = {x_start    ,   y_end/2.0    	,   0,  meshSize}; 
Point(2)  = {x_end    	,   y_end/2.0    	,   0,  meshSize}; 

lFibre  = newreg; Line(lFibre)  = {1,2};
Transfinite Line {lFibre} = 4;
Physical Line(777) = {lFibre};


