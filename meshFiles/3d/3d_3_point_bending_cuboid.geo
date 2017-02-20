meshSize       = 5;
notchHeight    = 20;
notchWidth     = 20;
length = 0;
height = 0;
depth  = 0;
//*********************************************************************************************
Function CreateRectangle
// create points
p1 = newp; Point(p1)  = {x_origin       	    ,   y_origin               ,   0,  meshSize};
p2 = newp; Point(p2)  = {x_origin + x_length  ,   y_origin               ,   0,  meshSize};
p3 = newp; Point(p3)  = {x_origin + x_length  ,	  y_origin + y_length    ,   0,  meshSize};
p4 = newp; Point(p4)  = {x_origin       	    ,   y_origin + y_length    ,   0,  meshSize};


// create lines
l1  = newl; Line(l1)  = {p1,p2};
l2  = newl; Line(l2)  = {p2,p3};
l3  = newl; Line(l3)  = {p3,p4};
l4  = newl; Line(l4)  = {p4,p1};

loop1 = newll; Line Loop(loop1) = {l1,l2,l3,l4};

plane1 = news; Plane Surface(plane1) = {loop1};

Return

// rectangle 0
x_origin  =  0;
x_length  =  30;
y_origin  =  0;
y_length  =  35;
z_end     =  20;
Call CreateRectangle;
surfaceId[0]     =   plane1;

// rectangle 4
y_origin  =  35;
y_length  =  5;
Call CreateRectangle;
surfaceId[4]     =   plane1;


// rectangle 1
x_origin  =  30;
y_origin  =  0;
x_length  =  20;
y_length  =  35;
index     =   0;
Call CreateRectangle;
surfaceId[1]     =   plane1;

// rectangle 5
y_origin  =  35;
y_length  =  5;
Call CreateRectangle;
surfaceId[5]     =   plane1;

// rectangle 3
meshSize  /= 5;
x_origin  =  60;
y_origin  =  0;
x_length  =  20;
y_length  =  35;
index     =   0;
Call CreateRectangle;
surfaceId[3]     =   plane1;

// rectangle 6
y_origin  =  35;
y_length  =  5;
Call CreateRectangle;
surfaceId[6]     =   plane1;


meshSize  *= 5;
// rectangle 2
x_origin  =  50;
y_origin  =  0;
x_length  =  10;
y_length  =  35;
index     =   0;
Call CreateRectangle;
surfaceId[2]     =   plane1;

// rectangle 7
y_origin  =  35;
y_length  =  5;
Call CreateRectangle;
surfaceId[7]     =   plane1;


listFine[] = Rotate {{0, 1, 0}, {80, 0, 0}, Pi} {
  Duplicata { Surface{surfaceId[3],surfaceId[6]}; }
};

list[] = Rotate {{0, 1, 0}, {80, 0, 0}, Pi} {
  Duplicata { Surface{surfaceId[0],surfaceId[0],surfaceId[1],surfaceId[2],surfaceId[4],surfaceId[5],surfaceId[7]}; }
};

Printf("surface = %g", list[0]);
//Transfinite Surface {surfaceId[0],surfaceId[1],surfaceId[3],list[0],list[1],list[3]};

volumeIds[] = Extrude {0,0,z_end}
{
Surface{list[],surfaceId[],listFine[]};
};

//Physical Volume(999) = {volumeIds[]};

Physical Volume(999) = {1:9,11:12,14};
Physical Volume(777) = {10,13,16,15};
