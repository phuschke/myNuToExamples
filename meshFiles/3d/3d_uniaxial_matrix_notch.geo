meshSize       = 2;
x_origin  =  0;
y_origin  =  0;
x_end     =  50;
y_end     =  10;
notch_width = 1;
notch_height = 1;
//*********************************************************************************************

// create points
Point(1)  = {x_origin   	        					,   y_origin    				,   0,  meshSize}; 
Point(2)  = {x_origin + x_end/2 - notch_width/2         ,   y_origin    				,   0,  meshSize};
Point(3)  = {x_origin + x_end/2					        ,   y_origin + notch_height     ,   0,  meshSize};
Point(4)  = {x_origin + x_end/2					        ,   y_end - notch_height        ,   0,  meshSize};
Point(5)  = {x_origin + x_end/2 - notch_width/2         ,   y_end         				,   0,  meshSize};
Point(6)  = {x_origin  								 	,   y_end         				,   0,  meshSize};



// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,5};
l5  = newreg; Line(l5)  = {5,6};
l6  = newreg; Line(l6)  = {6,1};

loop1 = newll; Line Loop(loop1) = {l1:l6};


plane1 = news; Plane Surface(plane1) = {loop1};



//
//Point(7)  = {x_origin - 5          ,   y_end + 5         				,   0,  meshSize};
//Point(8)  = {x_origin - 5          ,   y_origin - 5       				,   0,  meshSize};
//
//l8  = newreg; Line(l8)  = {7,8};
//l9  = newreg; Line(l9)  = {8,1};
//l10  = newreg; Line(l10)  = {6,7};
//loop2 = newll; Line Loop(loop2) = {l8,l9,-l6,l10};
//plane2 = news; Plane Surface(plane2) = {loop2};
//
//


Rotate {{0, 1, 0}, {x_end/2., 0, 0}, Pi} {
  Duplicata { Surface{plane1}; }
}

Extrude {0,0,20} 
{
Surface{8,9};
}

Physical Volume(999) = {1,2};



