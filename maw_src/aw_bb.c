#include "aw_bb.h"

void read_data_abb(char *fname, Abb *pw)
{
  FILE *fp;
  char buf[256]="";
  double td1,td2;
  
  if((fp=fopen(fname,"rt"))==NULL){ 
    printf("failed to read %s. no beams defined. Exit...\n",fname);
    exit(1);
  }
  
  fgets(buf,256,fp);  fgets(buf,256,fp);
  fscanf(fp,"%lf",&td1);    pw->rho0=td1;
  fscanf(fp,"%lf",&td1);    pw->c0  =td1;
  fscanf(fp,"%lf\n",&td1);  pw->f   =td1;
  fgets(buf,256,fp);
  fscanf(fp,"%lf",&td1);
  fscanf(fp,"%lf",&td2);    pw->p=td1+I*td2;
  fscanf(fp,"%lf",&td1);    pw->d_angle=td1;
  fscanf(fp,"%lf",&td1);    pw->tv[0]=td1;
  fscanf(fp,"%lf",&td1);    pw->tv[1]=td1;
  fscanf(fp,"%lf",&td1);    pw->tv[2]=td1;
  fscanf(fp,"%lf",&td1);    pw->theta=td1;
  fscanf(fp,"%lf",&td1);    pw->psi  =td1;
  
  fclose(fp); 
  
}

void print_data_abb(Abb *pw)
{
  double omega,alpha,r0,k0,kr,kz,absP,j1a;
  
  printf("-- Bessel beam --\n");
  printf("medium density                  [kg/m^3] : %15.14g\n",pw->rho0);
  printf("speed of sound in the medium       [m/s] : %15.14g\n",pw->c0);
  printf("frequency                           [Hz] : %15.14g\n",pw->f);
  printf("sound pressure amplitude            [Pa] :%7.6g+%7.6gI\n",creal(pw->p),cimag(pw->p));
  printf("deflection angle                   [rad] : %15.14g\n",pw->d_angle);
  printf("x-component of translation vector    [m] : %15.14g\n",pw->tv[0]);
  printf("y-component of translation vector    [m] : %15.14g\n",pw->tv[1]);
  printf("z-component of translation vector    [m] : %15.14g\n",pw->tv[2]);
  printf("rotation parameter theta           [rad] : %15.14g\n",pw->theta);
  printf("rotation parameter psi             [rad] : %15.14g\n",pw->psi);
  
  printf("-- additional information --\n");
  printf("wavelength                           [m] : % 15.14g\n",pw->c0/pw->f);
  omega=2.0*M_PI*pw->f;
  k0=omega/pw->c0;
  kr=k0*sin(pw->d_angle);
  kz=k0*cos(pw->d_angle);
  alpha=gsl_sf_bessel_zero_J0(1);
  r0=alpha/kr;
  printf("radius of the center spot            [m] : %15.14g\n",r0);
  j1a=gsl_sf_bessel_J1(alpha);
  absP=0.5*M_PI*alpha*alpha*j1a*j1a*kz/(omega*pw->rho0*kr*kr)*creal(pw->p*conj(pw->p));
  printf("energy passing through the spot      [W] : %15.14g\n",absP);
  printf("x-component of the energy flow       [W] : %15.14g\n",absP*sin(pw->theta)*cos(pw->psi));
  printf("y-component of the energy flow       [W] : %15.14g\n",absP*sin(pw->theta)*sin(pw->psi));
  printf("z-component of the energy flow       [W] : %15.14g\n",absP*cos(pw->theta));
  printf("\n");
}

void setup_abb(Abb *pw)
{
  double ct,st,cp,sp;
 
  pw->wd.omega=2.0*M_PI*pw->f;
  pw->wd.K0=pw->c0*pw->c0*pw->rho0;
  pw->wd.k1=I*pw->wd.omega/pw->wd.K0;
  pw->wd.k2=I*pw->wd.omega*pw->rho0;
  pw->wd.k0=pw->wd.omega/pw->c0;
  pw->wd.lambda0=pw->c0/pw->f;
  pw->wd.kr=pw->wd.k0*sin(pw->d_angle);
  pw->wd.kz=pw->wd.k0*cos(pw->d_angle);
  
  ct=cos(pw->theta);
  st=sin(pw->theta);
  cp=cos(pw->psi);
  sp=sin(pw->psi);
  
  pw->wd.R[0*3+0]=ct*cp*cp+sp*sp;  pw->wd.R[0*3+1]=sp*cp*(ct-1.0);  pw->wd.R[0*3+2]=st*cp;
  pw->wd.R[1*3+0]=sp*cp*(ct-1.0);  pw->wd.R[1*3+1]=ct*sp*sp+cp*cp;  pw->wd.R[1*3+2]=st*sp;
  pw->wd.R[2*3+0]=-st*cp;          pw->wd.R[2*3+1]=-st*sp;          pw->wd.R[2*3+2]=ct;
}

void calc_abb_pv(double complex *p,double complex *v,double *r,Abb *pw)
{
  void abb_pv_axz(double complex *p,double complex *v,double *r,Abb *pw);
  
  double complex vt[3];
  double *R,rc[3],r0[3];
  int i,j;
  
  R=pw->wd.R;

  for(i=0;i<3;i++) rc[i]=r[i]-pw->tv[i];
  r0[0]= R[0*3+0]*rc[0]+R[0*3+1]*rc[1]-R[0*3+2]*rc[2];
  r0[1]= R[1*3+0]*rc[0]+R[1*3+1]*rc[1]-R[1*3+2]*rc[2];
  r0[2]=-R[2*3+0]*rc[0]-R[2*3+1]*rc[1]+R[2*3+2]*rc[2]; 
  
  abb_pv_axz(p,vt,r0,pw);
  
  for(i=0;i<3;i++){
    v[i]=0.0;
    for(j=0;j<3;j++) v[i]+=R[i*3+j]*vt[j];
  } 
}

void calc_abb_dpdn(double complex *p,double complex *dpdn, double *r ,double *n, Abb *pw)
{
  void abb_dpdn_axz(double complex *p,double complex *dpdn,double *r,double *n,Abb *pw);
  
  double *R,rc[3],r0[3],n0[3];
  int i;

  R=pw->wd.R;

  for(i=0;i<3;i++) rc[i]=r[i]-pw->tv[i];
  r0[0]= R[0*3+0]*rc[0]+R[0*3+1]*rc[1]-R[0*3+2]*rc[2];
  r0[1]= R[1*3+0]*rc[0]+R[1*3+1]*rc[1]-R[1*3+2]*rc[2];
  r0[2]=-R[2*3+0]*rc[0]-R[2*3+1]*rc[1]+R[2*3+2]*rc[2]; 
  
  n0[0]= R[0*3+0]*n[0]+R[0*3+1]*n[1]-R[0*3+2]*n[2];
  n0[1]= R[1*3+0]*n[0]+R[1*3+1]*n[1]-R[1*3+2]*n[2];
  n0[2]=-R[2*3+0]*n[0]-R[2*3+1]*n[1]+R[2*3+2]*n[2]; 
  
  abb_dpdn_axz(p,dpdn,r0,n0,pw); 
}

///////////////////////////////////////////////////////////////////////
void abb_pv_axz(double complex *p,double complex *v,double *r,Abb *pw)
{
  double complex pce,ik2,cv;
  double arg,rr,j0,j1;
  
  arg=pw->wd.kz*r[2];
  pce=pw->p*(cos(arg)+I*sin(arg));
  ik2=1.0/pw->wd.k2;
  
  rr=sqrt(r[0]*r[0]+r[1]*r[1]);
  if(rr!=0.0){
    j0=gsl_sf_bessel_J0(pw->wd.kr*rr);
    j1=gsl_sf_bessel_J1(pw->wd.kr*rr);
    cv=pw->wd.kr*ik2/rr;
  
    *p=j0*pce;
    v[0]=-cv*r[0]*j1*pce;
    v[1]=-cv*r[1]*j1*pce;
    v[2]=I*pw->wd.kz*ik2*j0*pce;
  }
  else {
    cv=0.5*pw->wd.kr*pw->wd.kr*ik2;
    *p=pce;
    v[0]=-cv*r[0]*pce;
    v[1]=-cv*r[1]*pce;
    v[2]=I*pw->wd.kz*ik2*pce;
  }
}

void abb_dpdn_axz(double complex *p,double complex *dpdn,double *r,double *n,Abb *pw)
{
  double complex pce,cv,dpdx,dpdy,dpdz;
  double arg,rr,j0,j1;
  
  arg=pw->wd.kz*r[2];
  pce=pw->p*(cos(arg)+I*sin(arg));
  
  rr=sqrt(r[0]*r[0]+r[1]*r[1]);
  if(rr!=0.0){
    j0=gsl_sf_bessel_J0(pw->wd.kr*rr);
    j1=gsl_sf_bessel_J1(pw->wd.kr*rr);
    cv=pw->wd.kr/rr;
  
    *p=j0*pce;
    dpdx=-cv*r[0]*j1*pce;
    dpdy=-cv*r[1]*j1*pce;
    dpdz=I*pw->wd.kz*j0*pce;
    *dpdn=dpdx*n[0]+dpdy*n[1]+dpdz*n[2];
  }
  else {
    *p=pce;
    dpdz=I*pw->wd.kz*pce;
    *dpdn=dpdz*n[2];
  } 
}
