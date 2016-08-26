meshSize       = 10;
x_start = 0;
y_start = 0;
z_start = 0;
x_end = 200;
y_end = 10;
z_end = 10;

// create points
Point(1)  = {x_start    ,   y_start    	,   z_start,  meshSize}; 
Point(2)  = {x_end    	,   y_start    	,   z_start,  meshSize}; 
Point(3)  = {x_end    	,   y_end    	,   z_start,  meshSize}; 
Point(4)  = {x_start   	,	y_end    	,   z_start,  meshSize}; 

Point(5)  = {x_start    ,   y_start    	,   z_end,  meshSize}; 
Point(6)  = {x_end    	,   y_start    	,   z_end,  meshSize}; 
Point(7)  = {x_end    	,   y_end    	,   z_end,  meshSize}; 
Point(8)  = {x_start    ,	y_end    	,   z_end,  meshSize}; 



// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,1};


l5  = newreg; Line(l5)  = {5,6};
l6  = newreg; Line(l6)  = {6,7};
l7  = newreg; Line(l7)  = {7,8};
l8  = newreg; Line(l8)  = {8,5};

// create vertical lines
lv1 = newreg; Line(lv1) = {1,5};
lv2 = newreg; Line(lv2) = {2,6};
lv3 = newreg; Line(lv3) = {3,7};
lv4 = newreg; Line(lv4) = {4,8};


loop1 = newll; Line Loop(loop1) = {l1,l2,l3,l4};
loop2 = newll; Line Loop(loop2) = {l5,l6,l7,l8};
loop3 = newll; Line Loop(loop3) = {l1,lv2,-l5,-lv1};
loop4 = newll; Line Loop(loop4) = {-l3,lv3,l7,-lv4};
loop5 = newll; Line Loop(loop5) = {l2,lv3,-l6,-lv2};
loop6 = newll; Line Loop(loop6) = {l4,lv1,-l8,-lv4};

plane1 = news; Plane Surface(plane1) = {loop1};
plane2 = news; Plane Surface(plane2) = {loop2};
plane3 = news; Plane Surface(plane3) = {loop3};
plane4 = news; Plane Surface(plane4) = {loop4};
plane5 = news; Plane Surface(plane5) = {loop5};
plane6 = news; Plane Surface(plane6) = {loop6};
Transfinite Surface {plane1:plane6};
Recombine Surface {plane1:plane6};
//physS = newreg; Physical Surface(physS) = {plane1:plane6};
Sloop = newreg; Surface Loop(Sloop) = {plane1:plane6};
vol = newreg;Volume(vol) = {Sloop};
Transfinite Volume {vol};
Recombine Volume {vol};
Physical Volume(999) = {vol};
