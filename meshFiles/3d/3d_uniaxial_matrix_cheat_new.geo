meshSize       = 1;
x_origin  =  0;
y_origin  =  0;
factor = 0.5;
//*********************************************************************************************

// create points
Point(1)  = {factor * 0   	        							,   factor * 0		    				,   0,  meshSize}; 
Point(2)  = {factor * 10		 			    			    ,   factor * 0					 	    ,   0,  meshSize};
Point(3)  = {factor * 60	 			    			    	,   factor * 0					 	    ,   0,  meshSize};
Point(4)  = {factor * 60	 			    			    	,   factor * 5					 	    ,   0,  meshSize};
Point(5)  = {factor * 10						        		,   factor * 5					        ,   0,  meshSize};
Point(6)  = {factor * 0										 	,   factor * 5	         				,   0,  meshSize};




Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {2, 5};

Line Loop(1) = {1,7,5,6};
Line Loop(2) = {2,3,4,-7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Transfinite Surface {1,2};


Rotate {{0, 1, 0}, {factor * 60, 0, 0}, Pi} {
  Duplicata { Surface{1,2}; }
}

Rotate {{1, 0, 0}, {factor * 60, factor * 5, 0}, Pi} {
  Duplicata { Surface{1,2,8,13}; }
}

Transfinite Surface {1,2,8,13,17,22,27,32};
Extrude {0,0,factor * 10} 
{
Surface{1,2,8,13,17,22,27,32}; 
Layers {10};
}

Transfinite Volume {1:8};
Physical Volume(999) = {1:8};


Point(21001) = {0, 0, 0, 1.0};
Point(21002) = {120, 0, 0, 1.0};
Point(21003) = {120, 20, 0, 1.0};
Point(21004) = {0, 20, 0, 1.0};
Point(21005) = {0, 20, 40, 1.0};
Point(21006) = {0, 0, 40, 1.0};
Point(21007) = {120, 0, 40, 1.0};
Point(21008) = {120, 20, 40, 1.0};


