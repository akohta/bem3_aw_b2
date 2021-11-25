//------------------------------------------------------------
n1=10; // division number of 1/4 circle 

r1=0.005; // radius of sphere
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

Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};

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

// define domain
Physical Surface(1)={1,2,3,4,5,6}; // Domain 1 
Physical Surface(99)={-1,-2,-3,-4,-5,-6}; // Domain 0


// transfinite interpolation 
Transfinite Line{i0[0],i0[1],i0[2],
		   i1[0],i1[1],i1[2],
		   i2[0],i2[1],i2[2],
		   i3[0],i3[1],i3[2]}=n1+1;
Transfinite Surface{1,2,3,4,5,6};

// combine triangle element 
Recombine Surface{1,2,3,4,5,6};

