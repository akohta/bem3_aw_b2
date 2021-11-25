//---------------------------
n=10; // basic sampling number 
r=0.005; // radius of sphere
//---------------------------
// const
pp2=Pi/2.0;

// point
Point( 1)={0,0,0};
Point( 2)={r,0,0};
Point( 3)={0,0,r};
Point( 4)={0,0,-r};

// line
Circle( 1)={3,1,2};
Circle( 2)={2,1,4};
Transfinite Line{1,2}=n;

// extrude
o0[]=Extrude{{0,0,1},{0,0,0},pp2}{Line{1}; Layers{n}; Recombine;};
o1[]=Extrude{{0,0,1},{0,0,0},pp2}{Line{o0[0]}; Layers{n}; Recombine;};
o2[]=Extrude{{0,0,1},{0,0,0},pp2}{Line{o1[0]}; Layers{n}; Recombine;};
o3[]=Extrude{{0,0,1},{0,0,0},pp2}{Line{o2[0]}; Layers{n}; Recombine;};
o4[]=Extrude{{0,0,1},{0,0,0},pp2}{Line{2}; Layers{n}; Recombine;};
o5[]=Extrude{{0,0,1},{0,0,0},pp2}{Line{o4[0]}; Layers{n}; Recombine;};
o6[]=Extrude{{0,0,1},{0,0,0},pp2}{Line{o5[0]}; Layers{n}; Recombine;};
o7[]=Extrude{{0,0,1},{0,0,0},pp2}{Line{o6[0]}; Layers{n}; Recombine;};

Physical Surface( 1)={o0[1],o1[1],o2[1],o3[1],o4[1],o5[1],o6[1],o7[1]}; // Domain 1
Physical Surface(99)={-o0[1],-o1[1],-o2[1],-o3[1],-o4[1],-o5[1],-o6[1],-o7[1]}; // Domain 0 (opend)


