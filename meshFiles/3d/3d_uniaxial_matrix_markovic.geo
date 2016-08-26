meshSize       = 5;
x_origin  =  0;
y_origin  =  0;
factor = 1;
//*********************************************************************************************

// create points
Point(1)  = {factor * 0   	        							,   factor * 0		    				,   0,  meshSize}; 
Point(2)  = {factor * 20		 			    			    ,   factor * 0					 	    ,   0,  meshSize};

Point(3)  = {factor * 90	 			    			    	,   factor * 45					 	    ,   0,  meshSize};
Point(4)  = {factor * 90	 			    			    	,   factor * 80					 	    ,   0,  meshSize};

Point(5)  = {factor * 0						        			,   factor * 80					        ,   0,  meshSize};


Point(8) = {factor * 90, factor * -32, 0.0, meshSize};
Circle(2) = {2, 8, 3};

Line(1) = {1, 2};
//Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};

Line Loop(1) = {1:5};

Plane Surface(1) = {1};

Rotate {{0, 1, 0}, {factor * 90, 0, 0}, Pi} {
  Duplicata { Surface{1}; }
}

//Rotate {{1, 0, 0}, {factor * 90, 0, 0}, Pi} {
//  Duplicata { Surface{1,8}; }
//}
//
Extrude {0,0,factor * 35} 
{
Surface{1,6}; 
}

Physical Volume(999) = {1,2};


