//------------------------------------------------------------
n1=9; // division number of 1/4 inner circle 
n2=12; // division number of 1/4 outer circle 

r1=0.005; // radius of inner sphere
r2=0.008; // radius of outer sphere  
//------------------------------------------------------------

pp2=Pi/2.0;
w1=2.0/Sqrt(3.0)*r1;
w2=2.0/Sqrt(3.0)*r2;
h1=0.5*w1;
h2=0.5*w2;
s1=0.5*Sqrt(2.0)*w1;
s2=0.5*Sqrt(2.0)*w2;

// basic point and line
Point(1)={0,0,0};  
Point(2)={s1,0,-h1};
Point(3)={0,s1,-h1};
Point(4)={0,s1,h1};
Point(5)={s1,0,h1};
Point(6)={s2,0,-h2};
Point(7)={0,s2,-h2};
Point(8)={0,s2,h2};
Point(9)={s2,0,h2};

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={6,1,7};
Circle(5)={7,1,8};
Circle(6)={8,1,9};

// inner sphere
i0[]={1,2,3};
i1[]=Rotate{{0,0,1},{0,0,0},pp2*1.0} {Duplicata{Line{1,2,3};}};
i2[]=Rotate{{0,0,1},{0,0,0},pp2*2.0} {Duplicata{Line{1,2,3};}};
i3[]=Rotate{{0,0,1},{0,0,0},pp2*3.0} {Duplicata{Line{1,2,3};}};
Line Loop(1)={i0[0],i0[1],i0[2],-i3[1]};     Ruled Surface(1)={1} In Sphere{1};
Line Loop(2)={i1[0],i1[1],i1[2],-i0[1]};     Ruled Surface(2)={2} In Sphere{1};
Line Loop(3)={i2[0],i2[1],i2[2],-i1[1]};     Ruled Surface(3)={3} In Sphere{1};
Line Loop(4)={i3[0],i3[1],i3[2],-i2[1]};     Ruled Surface(4)={4} In Sphere{1};
Line Loop(5)={-i0[2],-i1[2],-i2[2],-i3[2]};  Ruled Surface(5)={5} In Sphere{1};
Line Loop(6)={-i0[0],-i1[0],-i2[0],-i3[0]};  Ruled Surface(6)={6} In Sphere{1};

// outer sphere
o0[]={4,5,6};
o1[]=Rotate{{0,0,1},{0,0,0},pp2*1.0} {Duplicata{Line{4,5,6};}};
o2[]=Rotate{{0,0,1},{0,0,0},pp2*2.0} {Duplicata{Line{4,5,6};}};
o3[]=Rotate{{0,0,1},{0,0,0},pp2*3.0} {Duplicata{Line{4,5,6};}};
Line Loop( 7)={o0[0],o0[1],o0[2],-o3[1]};     Ruled Surface( 7)={ 7} In Sphere{1};
Line Loop( 8)={o1[0],o1[1],o1[2],-o0[1]};     Ruled Surface( 8)={ 8} In Sphere{1};
Line Loop( 9)={o2[0],o2[1],o2[2],-o1[1]};     Ruled Surface( 9)={ 9} In Sphere{1};
Line Loop(10)={o3[0],o3[1],o3[2],-o2[1]};     Ruled Surface(10)={10} In Sphere{1};
Line Loop(11)={-o0[2],-o1[2],-o2[2],-o3[2]};  Ruled Surface(11)={11} In Sphere{1};
Line Loop(12)={-o0[0],-o1[0],-o2[0],-o3[0]};  Ruled Surface(12)={12} In Sphere{1};

// define domain
Physical Surface(1)={1,2,3,4,5,6}; // Domain 1 
Physical Surface(2)={-1,-2,-3,-4,-5,-6,7,8,9,10,11,12}; // Domain 2 
Physical Surface(99)={-7,-8,-9,-10,-11,-12}; // Domain 0



// transfinite interpolation 
Transfinite Line{i0[0],i0[1],i0[2],
		   i1[0],i1[1],i1[2],
		   i2[0],i2[1],i2[2],
		   i3[0],i3[1],i3[2]}=n1+1;
Transfinite Line{o0[0],o0[1],o0[2],
		   o1[0],o1[1],o1[2],
		   o2[0],o2[1],o2[2],
		   o3[0],o3[1],o3[2]}=n2+1;
Transfinite Surface{1,2,3,4,5,6,7,8,9,10,11,12};

// combine triangle element 
Recombine Surface{1,2,3,4,5,6,7,8,9,10,11,12};
