meshSize       = 20;
notchHeight    = 20;
notchWidth     = 20;
length = 0;
height = 0;
depth  = 0;
factor = 5;
//*********************************************************************************************
Function CreateRectangle
// create points
p1 = newp; Point(p1)  = {x_origin       ,   y_origin               ,   0,  meshSize};
p2 = newp; Point(p2)  = {x_origin + 30  ,   y_origin               ,   0,  meshSize};
p3 = newp; Point(p3)  = {x_origin + 50  ,	  y_origin               ,   0,  meshSize};
p4 = newp; Point(p4)  = {x_origin + 60  ,   y_origin               ,   0,  meshSize/factor};
p5 = newp; Point(p5)  = {x_origin + 80  ,   y_origin               ,   0,  meshSize/factor};

pv1 = newp; Point(pv1)  = {x_origin       ,   y_origin + y_length               ,   0,  meshSize};
pv2 = newp; Point(pv2)  = {x_origin + 30  ,   y_origin + y_length               ,   0,  meshSize};
pv3 = newp; Point(pv3)  = {x_origin + 50  ,	  y_origin + y_length               ,   0,  meshSize};
pv4 = newp; Point(pv4)  = {x_origin + 60  ,   y_origin + y_length               ,   0,  meshSize/factor};
pv5 = newp; Point(pv5)  = {x_origin + 80  ,   y_origin + y_length               ,   0,  meshSize/factor};

// create lines
l1  = newl; Line(l1)  = {p1,p2};
l2  = newl; Line(l2)  = {p2,p3};
l3  = newl; Line(l3)  = {p3,p4};
l4  = newl; Line(l4)  = {p4,p5};
l5  = newl; Line(l5)  = {p5,pv5};
l6  = newl; Line(l6)  = {pv5,pv4};
l7  = newl; Line(l7)  = {pv4,pv3};
l8  = newl; Line(l8)  = {pv3,pv2};
l9  = newl; Line(l9)  = {pv2,pv1};
l10  = newl; Line(l10)  = {pv1,p1};

loop1 = newll; Line Loop(loop1) = {l1:l10};

plane1 = news; Plane Surface(plane1) = {loop1};

Return

// rectangle 0
x_origin  =  0;
y_origin  =  0;
x_length  =  30;
y_length  =  40;
z_end     =  40;
Call CreateRectangle;

list[] = Rotate {{0, 1, 0}, {80, 0, 0}, Pi} {
  Duplicata { Surface{plane1}; }
};


volumeIds[] = Extrude {0,0,z_end}
{
Surface{list[],plane1};

};

Physical Volume(999) = {volumeIds[]};
