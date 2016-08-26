meshSize       = 100;
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
Point(6)  = {factor * 0										 	,   factor * 10	         				,   0,  meshSize};




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


Rotate {{0, 1, 0}, {factor * 120, 0, 0}, Pi} {
  Duplicata { Surface{1,2}; }
}

Rotate {{1, 0, 0}, {factor * 120, 0, 0}, Pi} {
  Duplicata { Surface{1,2,8,13}; }
}

Transfinite Surface {1,2,8,13,17,22,27,32};
Extrude {0,0,factor * 40} 
{
Surface{1,2,8,13,17,22,27,32}; 
Layers {2};
}

Transfinite Volume {1:8};
Physical Volume(999) = {1:8};


