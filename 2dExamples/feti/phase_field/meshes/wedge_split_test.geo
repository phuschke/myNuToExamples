meshSize       =0.005;

x0 = 0.0;
x1 = 1.0;

y0 = 0.0;
y1 = 0.6;

// create points
Point(1)  = {x0      ,     0.1   	,   0,  meshSize};
Point(2)  = {0.2      ,     0.1   	,   0,  meshSize};
Point(3)  = {0.25      ,     y0   	,   0,  meshSize};
Point(5)  = {x1   	  ,	    y1   	,   0,  meshSize};
Point(6)  = {x0   	  ,	    y1   	,   0,  meshSize};
Point(8)  = {x1   	  ,	    0.1   	,   0,  meshSize};
Point(10)  = {0.2   	  ,	    y1   	,   0,  meshSize};
Point(11)  = {x1   	  ,	    y0   	,   0,  meshSize};


//Rotate {{1, 0, 0}, {0, 0, 0}, Pi} {
//  Duplicata { Surface{17,19,21,23,25}; }
//}
//

//
//Transfinite Line{1:4} = 20 Using Progression 1.1;
//Transfinite Line{5:7} = 20 Using Progression 1.05;
//Transfinite Line{10,12,15,37} = 50 Using Progression 1.05;
//Transfinite Line{39} = 50 Using Progression 1/1.05;
//
//Transfinite Line{27,38} = 20 Using Progression 1.1;
//Transfinite Line{29,34} = 20 Using Progression 1/1.1;
//Transfinite Line{42,48} = 20 Using Progression 1.05;
//Transfinite Line{44} = 20 Using Progression 1/1.05;
//
//Physical Surface(9999) = {17,19,21,23,25,26,31,36,41,46};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 8};
//+
Line(3) = {1, 6};
//+
Line(4) = {2, 10};
//+
Line(5) = {8, 5};
//+
Line(6) = {11, 8};
//+
Line(7) = {3, 2};
//+
Line(8) = {3, 11};
//+
Line(9) = {6, 10};
//+
Line(10) = {10, 5};
//+
Line Loop(11) = {2, -6, -8, 7};
//+
Plane Surface(12) = {11};
//+
Line Loop(13) = {4, 10, -5, -2};
//+
Plane Surface(14) = {13};
//+
Line Loop(15) = {3, 9, -4, -1};
//+
Plane Surface(16) = {15};

Rotate {{1, 0, 0}, {0, 0, 0}, Pi} {
  Duplicata { Surface{12,14,16}; }
}

Transfinite Surface{12,14,16,17,22,27};
Recombine Surface{12,14,16,17,22,27};
Physical Surface(9999) = {12,14,16,17,22,27};

Transfinite Line{2,8,10,18,24} = 50 Using Progression 1.02;

Transfinite Line{7,6,21} = 20 Using Progression 1.02;
Transfinite Line{19} = 20 Using Progression 1/1.02;

Transfinite Line{3,4,5} = 10 Using Progression 1.1;
Transfinite Line{23,28} = 10 Using Progression 1.1;
Transfinite Line{25} = 10 Using Progression 1/1.1;

Transfinite Line{1,9,29} = 8 Using Progression 1/1.1;
Transfinite Line{31} = 8 Using Progression 1.1;



Mesh.Partitioner              =   1;            // 1= chaco, 2=metis
Mesh.ColorCarousel            =   3;              // color by partition
Mesh.NbPartitions             =   4; // number of partitions
Mesh.IgnorePartitionBoundary  =   0; // 0 = no, 1 = yes
Mesh.RecombineAll             =   1;
Mesh.SurfaceFaces             =   1;
