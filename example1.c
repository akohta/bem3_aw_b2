#include "bem3_aw_b2.h"

int main(int argc,char *argv[])
{
  MOBJ mo;
  double complex p,v[3];
  double F[3],N[3],r[3],rc[3];
  int type,m;

  mo_dat_read(argv[1],&mo);      // read datafile outputed by aw_d3b2_bv_solver
  mo_print_data(&mo);            // print data 

  r[0]=-0.03;                      // set x-coordinate
  r[1]= 0.01;                      // set y-coordinate
  r[2]= 0.002;                      // set z-coordinate 
  type=0;                        // selsect 4 point Gauss-Legendre 
  mo_pv_t(&p,v,r,type,&mo);      // calclation of total field ( add incident field to scattered field )
  //mo_pv_t_dbieq(&p,v,r,type,&mo);// for far-field 
  printf("Total electromagnetic field at r=( % g,% g,% g )\n",r[0],r[1],r[2]);
  printf("type=%d setting ( 4 point Gauss-Legendre )\n",type);
  printf("p  = % 15.14e %+15.14e I \n",creal(p),cimag(p));
  printf("v_x= % 15.14e %+15.14e I \n",creal(v[0]),cimag(v[0])); 
  printf("v_y= % 15.14e %+15.14e I \n",creal(v[1]),cimag(v[1]));
  printf("v_z= % 15.14e %+15.14e I \n",creal(v[2]),cimag(v[2]));
  type=1;                        // selsect 9 point(quadrangular element) or 7 point(triangular element) Gauss-Legendre 
  mo_pv_t(&p,v,r,type,&mo);      // calclation of total field ( add incident field to scattered field )
  //mo_pv_t_dbieq(&p,v,r,type,&mo);// for far-field 
  printf("Total electromagnetic field at r=( % g,% g,% g )\n",r[0],r[1],r[2]);
  printf("type=%d setting ( 9 point or 7 point Gauss-Legendre )\n",type);
  printf("p  = % 15.14e %+15.14e I \n",creal(p),cimag(p));
  printf("v_x= % 15.14e %+15.14e I \n",creal(v[0]),cimag(v[0])); 
  printf("v_y= % 15.14e %+15.14e I \n",creal(v[1]),cimag(v[1]));
  printf("v_z= % 15.14e %+15.14e I \n",creal(v[2]),cimag(v[2]));  

  rc[0]=0.0;                  // set x coordinate of rotation center 
  rc[1]=0.0;                  // set y coordinate of rotation center 
  rc[2]=0.0;                  // set z coordinate of rotation center 
  printf("\nRadiation force and torque\n");
  printf("type=0 setting ( 4 point Gauss-Legendre ) \n"); 
  type=0;                     // select 4 point Gauss-Legendre
  for(m=0;m<mo.N;m++){
    printf("object id = %d\n",m);
    mo_force_FN(F,N,rc,m,type,&mo);
    printf("F = (% 15.14e,% 15.14e,% 15.14e)\n",F[0],F[1],F[2]);
    printf("N = (% 15.14e,% 15.14e,% 15.14e), center of rotation (% g,% g,% g)\n",N[0],N[1],N[2],rc[0],rc[1],rc[2]); 
  }
 
  printf("type=1 setting ( 9 or 7 point Gauss-Legendre ) \n");
  type=1;                     // select 9 or 7 point Gauss-Legendre
  for(m=0;m<mo.N;m++){
    printf("objcet id = %d\n",m);
    mo_force_FN(F,N,rc,m,type,&mo);
    printf("F = (% 15.14e,% 15.14e,% 15.14e)\n",F[0],F[1],F[2]);
    printf("N = (% 15.14e,% 15.14e,% 15.14e), center of rotation (% g,% g,% g)\n",N[0],N[1],N[2],rc[0],rc[1],rc[2]); 
  }
  
  mo_finalize(&mo);
  return 0;
}
