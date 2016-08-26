Function CreateRectangle

  // create lines
  p1 = newp; Point(p1) = {x,      y,       0,  meshSize};
  p2 = newp; Point(p2) = {x+lx,   y,       0,  meshSize};
  p3 = newp; Point(p3) = {x+lx,   y+ly,    0,  meshSize};
  p4 = newp; Point(p4) = {x   ,   y+ly,    0,  meshSize};

  // create lines
  l1  = newreg; Line(l1)  = {p1,p2};
  l2  = newreg; Line(l2)  = {p2,p3};
  l3  = newreg; Line(l3)  = {p3,p4};
  l4  = newreg; Line(l4)  = {p4,p1};

  Transfinite Line {l1,l3} = 9+1;
  Transfinite Line {l2,l4} = 3+1;

  loop1 = newll; Line Loop(loop1) = {l1,l2,l3,l4};

  plane1 = news; Plane Surface(plane1) = {loop1};
  Transfinite Surface {plane1};

  Physical Surface(4000+partition)  = plane1;

Return

// 100X boundary
// 200X master interface
// 300X slave interface
// 400X subdomain mesh

meshSize       = 10;
x  = 0;
y  = 0;
lx = 20;
ly = 10;
num_partitions = 3;

// left subdomain
partition = 1;
Call CreateRectangle;
Physical Line(1000+partition)     = {l2};
Physical Line(2000+partition)     = {l4};

// inner subdomains
For partition In {2:num_partitions-1}
x += lx;
Call CreateRectangle;
Physical Line(2000+partition)     = {l4};
Physical Line(3000+partition)     = {l2};
EndFor

// right subdomain
partition = num_partitions;
x += lx;
Call CreateRectangle;
Physical Line(3000+partition)     = {l2};


Mesh.Partitioner              =   1;            // 1= chaco, 2=metis
Mesh.ColorCarousel            =   3;              // color by partition
Mesh.MshFilePartitioned       =   1;              // Create physicals by partiion
Mesh.NbPartitions             =   num_partitions; // number of partitions
Mesh.IgnorePartitionBoundary  =   0; // 0 = no, 1 = yes
Mesh.RecombineAll             =   1;
Mesh.SurfaceFaces             =   1;
