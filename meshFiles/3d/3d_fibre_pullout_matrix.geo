meshSize       = 100;
x_start = 0;
y_start = 0;
z_start = 0;
x_end = 10;
y_end = 10;
z_end = 10;

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




Extrude {0,0,z_end} 
{
Surface{plane1}; 
}

Point(100)  = {5    ,   5    	,   5,  meshSize}; 
Point(200)  = {6   	,   6    	,   6,  meshSize}; 
Point(300)  = {6+1   	,   6+1    	,   6,  meshSize}; 
Point(400)  = {5+1    ,   5+1    	,   5,  meshSize}; 

l100  = newreg; Line(l100)  = {100,200};
l200  = newreg; Line(l200)  = {200,300};
l300  = newreg; Line(l300)  = {300,400};
l400  = newreg; Line(l400)  = {400,100};

loop100 = newll; Line Loop(loop100) = {l100,l200,l300,l400};

plane100 = news; Plane Surface(plane100) = {loop100};

Surface {plane100} In Volume {1};


l500  = newreg; Line(l500)  = {100,300};
Line {l500} In Surface {plane100};

Physical Volume(999) = {1};

