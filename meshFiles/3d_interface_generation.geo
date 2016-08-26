meshSize       = 0.5;
Geometry.ExtrudeReturnLateralEntities = 2;
// create points
Point(1)  = { 0	,   0   ,   0,  meshSize}; 
Point(2)  = { 1	,   0   ,   0,  meshSize}; 
Point(3)  = { 1	,   1   ,   0,  meshSize}; 
Point(4)  = { 0	,	1   ,   0,  meshSize}; 



// create lines
l1  = newreg; Line(l1)  = {1,2};
l2  = newreg; Line(l2)  = {2,3};
l3  = newreg; Line(l3)  = {3,4};
l4  = newreg; Line(l4)  = {4,1};



loop1 = newll; Line Loop(loop1) = {l1,l2,l3,l4};


plane1 = news; Plane Surface(plane1) = {loop1};


Extrude {0,0,1} 
{
Surface{plane1}; 
}


delta = 0.05;
Point(101)  = { 0.2	,   0.2   ,   0.2,  meshSize}; 
Point(102)  = { 0.4	,   0.2   ,   0.2,  meshSize}; 
Point(103)  = { 0.4	,   0.2+delta   ,   0.2,  meshSize}; 
Point(104)  = { 0.2	,	0.2+delta   ,   0.2,  meshSize}; 



// create lines
l101  = newreg; Line(l101)  = {101,102};
l102  = newreg; Line(l102)  = {102,103};
l103  = newreg; Line(l103)  = {103,104};
l104  = newreg; Line(l104)  = {104,101};

loop101 = newll; Line Loop(loop101) = {l101,l102,l103,l104};

plane101 = news; Plane Surface(plane101) = {loop101};

out[] = Extrude {0,0,delta} 
{
Surface{plane101}; 
};

Printf("plane101 = %g", plane101);
Printf("0 = %g", out[0]);
Printf("1 = %g", out[1]);
Printf("2 = %g", out[2]);
Printf("3 = %g", out[3]);
Printf("4 = %g", out[4]);
Printf("5 = %g", out[5]);

Physical Volume(999) = {plane1,-out[1]};




Delete {
  Volume{1, 2};
}


Surface Loop(1000) = {23, 6, 15, 19, 28, 27};
Surface Loop(1001) = {51, 34, 43, 47, 56, 55};
Volume(1002) = {1000, 1001};

Translate {0,0,delta}
{
	Surface{plane101};
}

Translate {0,0,delta}
{
	Surface{plane101};
}
