meshSize       = 5;

// create points
Point(1)  = {0,     0,    0,  meshSize};
Point(2)  = {250,   0,    0,  meshSize};
Point(3)  = {250,   250,   0,  meshSize};
Point(4)  = {470,   250,   0,  meshSize};

Point(5)  = {500,   250,   0,  meshSize};
Point(6)  = {500,   500,   0,  meshSize};
Point(7)  = {470,   500,   0,  meshSize};
Point(8)  = {250,   500,   0,  meshSize};
Point(9)  = {0,     500,   0,  meshSize};
Point(10) = {0,     250,   0,  meshSize};

// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,5};
l5  = newreg; Line(l5)  = {5,6};
l6  = newreg; Line(l6)  = {6,7};
l7  = newreg; Line(l7)  = {7,8};
l8  = newreg; Line(l8)  = {8,9};
l9  = newreg; Line(l9)  = {9,10};
l10 = newreg; Line(l10)  = {10,1};

lBetween1 = newreg; Line(lBetween1)  = {3,10};
lBetween2 = newreg; Line(lBetween2)  = {4,7};
lBetween3 = newreg; Line(lBetween3)  = {3,8};

loop1 = newll; Line Loop(loop1) = {l1,l2,lBetween1,l10};
plane1 = news; Plane Surface(plane1) = {loop1};

loop2 = newll; Line Loop(loop2) = {l3,lBetween2,l7,-lBetween3};
plane2 = news; Plane Surface(plane2) = {loop2};
//
loop3  = newll; Line Loop(loop3) = {l4,l5,l6,-lBetween2};
plane3 = news; Plane Surface(plane3) = {loop3};

loop4  = newll; Line Loop(loop4) = {-lBetween1,lBetween3,l8,l9};
plane4 = news; Plane Surface(plane4) = {loop4};

//

Transfinite Surface {plane1,plane2, plane3,plane4};
Recombine Surface {plane1,plane2, plane3,plane4};
Physical Surface(999) = {plane1,plane2, plane3,plane4};
