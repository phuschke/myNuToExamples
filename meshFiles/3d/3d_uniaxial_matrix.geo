meshSize       = 1;
x_origin  =  0;
y_origin  =  0;
x_end     =  50;
y_end     =  10;
notch_width = 1;
notch_height = 1;
//*********************************************************************************************

// create points
Point(1)  = {x_origin   	        					,   y_origin    				,   0,  meshSize}; 
Point(2)  = {x_origin + x_end/2 			        		,   y_origin			 	    ,   0,  meshSize};
Point(3)  = {x_origin + x_end/2				        ,   y_end				        ,   0,  meshSize};
Point(4)  = {x_origin  								 	,   y_end         				,   0,  meshSize};



// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,1};

loop1 = newll; Line Loop(loop1) = {l1:l4};


plane1 = news; Plane Surface(plane1) = {loop1};
Transfinite Surface{plane1};

Rotate {{0, 1, 0}, {x_end/2., 0, 0}, Pi} {
  Duplicata { Surface{plane1}; }
}

Transfinite Surface{7};

Extrude {0,0,20} 
{

Surface{plane1,7}; 
Layers{20};
}

Transfinite Volume{1,2};
Physical Volume(999) = {1,2};



