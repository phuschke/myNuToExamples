// Standard specimen dimensions: 100*40*40
meshSize       	= 2;


//*********************************************************************************************
Function CreateNotchedRectangle
// create points
Point(1)  = {0       	        ,   0    ,   0,  meshSize}; 
Point(2)  = {75	                ,   0    ,   0,  meshSize}; 
Point(3)  = {80	                ,   5    ,   0,  meshSize};
Point(4)  = {85   	            ,   0    ,   0,  meshSize};
Point(5)  = {160   	            ,   0    ,   0,  meshSize};

Point(6)  = {160       		    ,   40    ,   0,  	meshSize}; 
Point(7)  = {80                 ,   40    ,   0,  	meshSize}; 
Point(8)  = {0                  ,   40    ,   0,  	meshSize}; 

// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,5};
l5  = newreg; Line(l5)  = {5,6};
l6  = newreg; Line(l6)  = {6,7};
l7  = newreg; Line(l7)  = {7,8};
l8  = newreg; Line(l8)  = {8,1};


loop1 = newll; Line Loop(loop1) = {l1,l2,l3,l4,l5,l6,l7,l8};



plane1 = news; Plane Surface(plane1) = {loop1};


Return
// first rectangle

Call CreateNotchedRectangle;



Point(100)  = {70       		    ,   20    ,   0,  	meshSize}; 
Point(200)  = {90                   ,   20    ,   0,  	meshSize}; 
lfibre  = newreg; Line(lfibre)  = {100,200};
Line {lfibre} In Surface {plane1};

Physical Line(999) = {lfibre};
Physical Surface(222) = {plane1};