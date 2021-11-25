//------------------------------------------------------------
n1=10; // division number of 1/4 circle 

r1= 0.005; // radius of sphere
x2= 0.000;
y2= 0.000;
z2= 0.011; // center of sphere2
x3= 0.000;
y3= 0.000;
z3=-0.011; // center of sphere3
//------------------------------------------------------------

pp2=Pi/2.0;
w1=2.0/Sqrt(3.0)*r1;
h1=0.5*w1;
s1=0.5*Sqrt(2.0)*w1;

// basic point and line
Point(1)={0,0,0};  
Point(2)={s1,0,-h1};
Point(3)={0,s1,-h1};
Point(4)={0,s1,h1};
Point(5)={s1,0,h1};
Point(6)={x2,y2,z2};
Point(7)={x3,y3,z3};

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};

// sphere1
i0[]={1,2,3};
i1[]=Rotate{{0,0,1},{0,0,0},pp2*1.0} {Duplicata{Line{1,2,3};}};
i2[]=Rotate{{0,0,1},{0,0,0},pp2*2.0} {Duplicata{Line{1,2,3};}};
i3[]=Rotate{{0,0,1},{0,0,0},pp2*3.0} {Duplicata{Line{1,2,3};}};
Line Loop(1)={ i0[0], i0[1], i0[2],-i3[1]};  Ruled Surface(1)={1} In Sphere{1};
Line Loop(2)={ i1[0], i1[1], i1[2],-i0[1]};  Ruled Surface(2)={2} In Sphere{1};
Line Loop(3)={ i2[0], i2[1], i2[2],-i1[1]};  Ruled Surface(3)={3} In Sphere{1};
Line Loop(4)={ i3[0], i3[1], i3[2],-i2[1]};  Ruled Surface(4)={4} In Sphere{1};
Line Loop(5)={-i0[2],-i1[2],-i2[2],-i3[2]};  Ruled Surface(5)={5} In Sphere{1};
Line Loop(6)={-i0[0],-i1[0],-i2[0],-i3[0]};  Ruled Surface(6)={6} In Sphere{1};

// sphere2
j0[]=Translate{x2,y2,z2}{Duplicata{Line{1,2,3};}};
j1[]=Translate{x2,y2,z2}{Duplicata{Line{i1[0],i1[1],i1[2]};}};
j2[]=Translate{x2,y2,z2}{Duplicata{Line{i2[0],i2[1],i2[2]};}};
j3[]=Translate{x2,y2,z2}{Duplicata{Line{i3[0],i3[1],i3[2]};}};
Line Loop( 7)={ j0[0], j0[1], j0[2],-j3[1]};  Ruled Surface( 7)={ 7} In Sphere{6};
Line Loop( 8)={ j1[0], j1[1], j1[2],-j0[1]};  Ruled Surface( 8)={ 8} In Sphere{6};
Line Loop( 9)={ j2[0], j2[1], j2[2],-j1[1]};  Ruled Surface( 9)={ 9} In Sphere{6};
Line Loop(10)={ j3[0], j3[1], j3[2],-j2[1]};  Ruled Surface(10)={10} In Sphere{6};
Line Loop(11)={-j0[2],-j1[2],-j2[2],-j3[2]};  Ruled Surface(11)={11} In Sphere{6};
Line Loop(12)={-j0[0],-j1[0],-j2[0],-j3[0]};  Ruled Surface(12)={12} In Sphere{6};

// sphere3
k0[]=Translate{x3,y3,z3}{Duplicata{Line{1,2,3};}};
k1[]=Translate{x3,y3,z3}{Duplicata{Line{i1[0],i1[1],i1[2]};}};
k2[]=Translate{x3,y3,z3}{Duplicata{Line{i2[0],i2[1],i2[2]};}};
k3[]=Translate{x3,y3,z3}{Duplicata{Line{i3[0],i3[1],i3[2]};}};
Line Loop(13)={ k0[0], k0[1], k0[2],-k3[1]};  Ruled Surface(13)={13} In Sphere{7};
Line Loop(14)={ k1[0], k1[1], k1[2],-k0[1]};  Ruled Surface(14)={14} In Sphere{7};
Line Loop(15)={ k2[0], k2[1], k2[2],-k1[1]};  Ruled Surface(15)={15} In Sphere{7};
Line Loop(16)={ k3[0], k3[1], k3[2],-k2[1]};  Ruled Surface(16)={16} In Sphere{7};
Line Loop(17)={-k0[2],-k1[2],-k2[2],-k3[2]};  Ruled Surface(17)={17} In Sphere{7};
Line Loop(18)={-k0[0],-k1[0],-k2[0],-k3[0]};  Ruled Surface(18)={18} In Sphere{7};


// define domain
Physical Surface(1)={ 1, 2, 3, 4, 5, 6}; // Domain 1
Physical Surface(2)={ 7, 8, 9,10,11,12}; // Domain 2
Physical Surface(3)={13,14,15,16,17,18}; // Domain 3 
Physical Surface(99)={ -1, -2, -3, -4, -5, -6,
                       -7, -8, -9,-10,-11,-12,
                      -13,-14,-15,-16,-17,-18}; // Domain 0 ( Open region )

// transfinite interpolation 
Transfinite Line{i0[0],i0[1],i0[2],
                   i1[0],i1[1],i1[2],
                   i2[0],i2[1],i2[2],
                   i3[0],i3[1],i3[2]}=n1+1;
Transfinite Surface{1,2,3,4,5,6};

Transfinite Line{j0[0],j0[1],j0[2],
                 j1[0],j1[1],j1[2],
                 j2[0],j2[1],j2[2],
                 j3[0],j3[1],j3[2]}=n1+1;
Transfinite Surface{7,8,9,10,11,12};

Transfinite Line{k0[0],k0[1],k0[2],
                 k1[0],k1[1],k1[2],
                 k2[0],k2[1],k2[2],
                 k3[0],k3[1],k3[2]}=n1+1;
Transfinite Surface{13,14,15,16,17,18};

// combine triangular element 
Recombine Surface{1,2,3,4,5,6};
Recombine Surface{7,8,9,10,11,12};
Recombine Surface{13,14,15,16,17,18};

