meshSize       = 1;
x_origin  =  0;
y_origin  =  0;
factor = 1;
//*********************************************************************************************

// create points
Point(1)  = {factor * 0   	        							,   factor * 0		    				,   0,  meshSize}; 
Point(2)  = {factor * 10		 			    			    ,   factor * 0					 	    ,   0,  meshSize};
Point(3)  = {factor * 10						        		,   factor * 10					        ,   0,  meshSize};
Point(4)  = {factor * 0										 	,   factor * 10	         				,   0,  meshSize};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};
Transfinite Surface {1};

Extrude {0,0,factor * 10} 
{
Surface{1}; 
Layers {5};
}

Transfinite Volume {1};
Physical Volume(999) = {1};


