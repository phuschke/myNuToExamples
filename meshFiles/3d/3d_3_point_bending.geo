meshSize       = 20;
notchHeight    = 20;
notchWidth     = 20;
length = 0;
height = 0;
depth  = 0;
//*********************************************************************************************
Function CreateCuboid
// create points
p1 = newp; Point(p1)  = {x_origin       	       ,   y_origin    ,   0,  meshSize};
p2 = newp; Point(p2)  = {x_origin + l_support    ,   y_origin    ,   0,  meshSize};
p3 = newp; Point(p3)  = {x_end/2.     	         ,	y_origin    ,   0,  meshSize};
p4 = newp; Point(p4)  = {x_origin       	       ,   y_end    ,   0,  meshSize};
p5 = newp; Point(p5)  = {x_origin + l_support    ,   y_end    ,   0,  meshSize};
p6 = newp; Point(p6)  = {x_end/2.     	         ,	 y_end    ,   0,  meshSize};


// create lines
l1  = newreg; Line(l1)  = {p1,p2};
l2  = newreg; Line(l2)  = {p2,p3};

l4  = newreg; Line(l4)  = {p4,p5};
l5  = newreg; Line(l5)  = {p5,p6};

// create inner vertical lines
lv1 = newreg; Line(lv1) = {p1,p4};
lv2 = newreg; Line(lv2) = {p2,p5};
lv3 = newreg; Line(lv3) = {p3,p6};

loop1 = newll; Line Loop(loop1) = {l1,lv2,-l4,-lv1};
loop2 = newll; Line Loop(loop2) = {l2,lv3,-l5,-lv2};


plane1 = news; Plane Surface(plane1) = {loop1};
plane2 = news; Plane Surface(plane2) = {loop2};




Return

// first rectangle
x_origin  =  0;
y_origin  =  0;
x_end     =  120;
l_support =  30;
y_end     =  40;
z_end     =  40;
index     =   0;
Call CreateCuboid;

Rotate {{0, 1, 0}, {x_end/2. +20, 0, 0}, Pi} {
  Duplicata { Surface{plane1:plane2}; }
}

Transfinite Surface {10,11,17,12};

Extrude {0,0,z_end}
{
Surface{10,11,17,12};
Layers{2};
}

Physical Volume(777) = {1:4};


// first rectangle
x_origin  =  70;
y_origin  =  0;
x_end     =  160;
l_support =  5;
y_end     =  40;
z_end     =  40;
index     =   0;
Call CreateCuboid;

Rotate {{0, 1, 0}, {x_end/2., 0, 0}, Pi} {
  Duplicata { Surface{plane1:plane2}; }
}

Transfinite Surface {787,788,794,789};

Extrude {0,0,z_end}
{
Surface{787,788,794,789};
Layers{2};
}

Physical Volume(999) = {5:8};
