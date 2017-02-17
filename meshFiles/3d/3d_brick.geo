meshSize       = 20;
notchHeight    = 20;
notchWidth     = 20;
length = 0;
height = 0;
depth  = 0;
//*********************************************************************************************
Function CreateCuboid
// create points
p1 = newp; Point(p1)  = {x_origin       	   ,   y_origin    ,   0,  meshSize};
p2 = newp; Point(p2)  = {x_end               ,   y_origin    ,   0,  meshSize};
p3 = newp; Point(p3)  = {x_end     	         ,	 y_end       ,   0,  meshSize};
p4 = newp; Point(p4)  = {x_origin       	   ,   y_end       ,   0,  meshSize};



// create lines
l1  = newl; Line(l1)  = {p1,p2};
l2  = newl; Line(l2)  = {p2,p3};
l3  = newl; Line(l3)  = {p3,p4};
l4  = newl; Line(l4)  = {p4,p1};

loop1 = newll; Line Loop(loop1) = {l1:l4};

plane1 = news; Plane Surface(plane1) = {loop1};

Return

// first rectangle
x_origin  =  0;
y_origin  =  0;
x_end     =  10;
y_end     =  10;
z_end     =  10;
index     =   0;
Call CreateCuboid;


Transfinite Surface {plane1};
Recombine Surface {plane1};

volumeIds[] = Extrude {0,0,z_end}
{
Surface{plane1};
Layers{1};
Recombine;
};

Physical Volume(777) = {volumeIds[]};
