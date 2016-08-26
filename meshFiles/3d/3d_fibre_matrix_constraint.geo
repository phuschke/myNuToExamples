meshSize       = 40;
x_start = 0;
y_start = 0;
z_start = 0;
x_end = 100;
y_end = 40;
z_end = 40;

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
Transfinite Line {l1} = 3;
Transfinite Line {l2} = 3;
Transfinite Line {l3} = 3;
Transfinite Line {l4} = 3;



loop1 = newll; Line Loop(loop1) = {l1,l2,l3,l4};

plane1 = news; Plane Surface(plane1) = {loop1};
Transfinite Surface {plane1};




Extrude {0,0,z_end} 
{
Surface{plane1}; 
Layers{1};
}

Transfinite Volume {1};

Physical Volume(999) = {1};

Point(100)  = {20    ,   5    	,   5,  meshSize}; 
Point(200)  = {40    	,   5    	,   5,  meshSize}; 