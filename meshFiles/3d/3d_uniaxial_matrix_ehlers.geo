	meshSize       = 3;
x_origin  =  0;
y_origin  =  0;
factor = 1;
//*********************************************************************************************

// create points
Point(1)  = {factor * 0   	        							,   factor * 0		    				,   0,  meshSize}; 
Point(2)  = {factor * 70		 			    			    ,   factor * 0					 	    ,   0,  meshSize};
Point(3)  = {factor * 120	 			    			    	,   factor * 0					 	    ,   0,  meshSize};
Point(4)  = {factor * 120	 			    			    	,   factor * 10					 	    ,   0,  meshSize};
Point(5)  = {factor * 70						        		,   factor * 10					        ,   0,  meshSize};
Point(6)  = {factor * 40										,   factor * 20	         				,   0,  meshSize};
Point(7)  = {factor * 0										 	,   factor * 20	         				,   0,  meshSize};


Point(8) = {factor * 50, factor * 12, 0.0, meshSize};
Spline(5) = {5, 8, 6};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
//Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};

Line Loop(1) = {1:7};

Plane Surface(1) = {1};



Rotate {{0, 1, 0}, {factor * 120, 0, 0}, Pi} {
  Duplicata { Surface{1}; }
}

Rotate {{1, 0, 0}, {factor * 120, 0, 0}, Pi} {
  Duplicata { Surface{1,8}; }
}

Extrude {0,0,factor * 40} 
{
Surface{1,8,16,24}; 
}

Physical Volume(999) = {1,2,3,4};


