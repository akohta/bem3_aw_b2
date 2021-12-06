#include "aw_pw.h"

void read_data_apw(char *fname, Apw *pw)
{
  FILE *fp;
  char buf[256]="";
  double td1,td2;
  
  if((fp=fopen(fname,"rt"))==NULL){ 
    printf("failed to read %s. no beams defined. Exit...\n",fname);
    exit(1);
  }  
  
  if(fgets(buf,256,fp)==NULL){
    printf("aw_pw.c, read_data_apw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("aw_pw.c, read_data_apw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_pw.c, read_data_apw(), failed to read the rho0. exit...\n");
    exit(1);
  }
  pw->rho0=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_pw.c, read_data_apw(), failed to read the c0. exit...\n");
    exit(1);
  }
  pw->c0  =td1;
  if(fscanf(fp,"%lf\n",&td1)!=1){
    printf("aw_pw.c, read_data_apw(), failed to read the f. exit...\n");
    exit(1);
  }
  pw->f   =td1;
  if(fgets(buf,256,fp)==NULL){
    printf("aw_pw.c, read_data_apw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_pw.c, read_data_apw(), failed to read the real(p). exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&td2)!=1){
    printf("aw_pw.c, read_data_apw(), failed to read the imag(p). exit...\n");
    exit(1);
  }
  pw->p=td1+I*td2;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_pw.c, read_data_apw(), failed to read the tv[0]. exit...\n");
    exit(1);
  }
  pw->tv[0]=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_pw.c, read_data_apw(), failed to read the tv[1]. exit...\n");
    exit(1);
  }
  pw->tv[1]=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_pw.c, read_data_apw(), failed to read the tv[2]. exit...\n");
    exit(1);
  }
  pw->tv[2]=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_pw.c, read_data_apw(), failed to read the theta. exit...\n");
    exit(1);
  }
  pw->theta=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_pw.c, read_data_apw(), failed to read the psi. exit...\n");
    exit(1);
  }
  pw->psi  =td1;
  
  fclose(fp);
}

void print_data_apw(Apw *pw)
{
  double absI;
  
  printf("-- planewave --\n");
  printf("medium density                   [kg/m^3] : %15.14g\n",pw->rho0);
  printf("speed of sound in the medium        [m/s] : %15.14g\n",pw->c0);
  printf("frequency                            [Hz] : %15.14g\n",pw->f);
  printf("sound pressure amplitude             [Pa] :%7.6g+%7.6gI\n",creal(pw->p),cimag(pw->p));
  printf("x-component of translation vector     [m] : %15.14g\n",pw->tv[0]);
  printf("y-component of translation vector     [m] : %15.14g\n",pw->tv[1]);
  printf("z-component of translation vector     [m] : %15.14g\n",pw->tv[2]);
  printf("rotation parameter theta            [rad] : %15.14g\n",pw->theta);
  printf("rotation parameter psi              [rad] : %15.14g\n",pw->psi);
  printf("-- additional information --\n");
  absI=0.5*creal(pw->p*conj(pw->p))/(pw->rho0*pw->c0);
  printf("wavelength                            [m] : %15.14g\n",pw->c0/pw->f);
  printf("absolute value of sound intensity [W/m^2] : %15.14g\n",absI);
  printf("x-component of sound intensity    [W/m^2] : %15.14g\n",absI*sin(pw->theta)*cos(pw->psi));
  printf("y-component of sound intensity    [W/m^2] : %15.14g\n",absI*sin(pw->theta)*sin(pw->psi));
  printf("z-component of sound intensity    [W/m^2] : %15.14g\n",absI*cos(pw->theta));
  printf("\n");
}

void setup_apw(Apw *pw)
{
  double ct,st,cp,sp;

  pw->wd.omega=2.0*M_PI*pw->f;
  pw->wd.K0=pw->c0*pw->c0*pw->rho0;
  pw->wd.k1=I*pw->wd.omega/pw->wd.K0;
  pw->wd.k2=I*pw->wd.omega*pw->rho0;
  pw->wd.k0=pw->wd.omega/pw->c0;
  pw->wd.lambda0=pw->c0/pw->f;

  ct=cos(pw->theta);
  st=sin(pw->theta);
  cp=cos(pw->psi);
  sp=sin(pw->psi);
  
  pw->wd.R[0*3+0]=ct*cp*cp+sp*sp;  pw->wd.R[0*3+1]=sp*cp*(ct-1.0);  pw->wd.R[0*3+2]=st*cp;
  pw->wd.R[1*3+0]=sp*cp*(ct-1.0);  pw->wd.R[1*3+1]=ct*sp*sp+cp*cp;  pw->wd.R[1*3+2]=st*sp;
  pw->wd.R[2*3+0]=-st*cp;          pw->wd.R[2*3+1]=-st*sp;          pw->wd.R[2*3+2]=ct;
}

void calc_apw_pv(double complex *p,double complex *v,double *r,Apw *pw)
{
  void apw_pv_axz(double complex *p,double complex *v,double *r,Apw *pw);
  
  double complex vt[3];
  double *R,rc[3],r0[3];
  int i,j;
  
  R=pw->wd.R;

  for(i=0;i<3;i++) rc[i]=r[i]-pw->tv[i];
  r0[0]= R[0*3+0]*rc[0]+R[0*3+1]*rc[1]-R[0*3+2]*rc[2];
  r0[1]= R[1*3+0]*rc[0]+R[1*3+1]*rc[1]-R[1*3+2]*rc[2];
  r0[2]=-R[2*3+0]*rc[0]-R[2*3+1]*rc[1]+R[2*3+2]*rc[2]; 
  
  apw_pv_axz(p,vt,r0,pw);
  
  for(i=0;i<3;i++){
    v[i]=0.0;
    for(j=0;j<3;j++) v[i]+=R[i*3+j]*vt[j];
  }
}

void calc_apw_dpdn(double complex *p,double complex *dpdn, double *r ,double *n, Apw *pw)
{
  void apw_dpdn_axz(double complex *p,double complex *dpdn,double *r,double *n,Apw *pw);
  
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
  
  apw_dpdn_axz(p,dpdn,r0,n0,pw);
}


//////////////////////////////////////////////////////////////////////
void apw_pv_axz(double complex *p,double complex *v,double *r,Apw *pw)
{
  double complex ce;
  double arg;
  
  arg=pw->wd.k0*r[2];
  ce=cos(arg)+I*sin(arg);
  
  *p=pw->p*ce;
  v[0]=0.0;
  v[1]=0.0;
  v[2]=(*p)/(pw->rho0*pw->c0);
}

void apw_dpdn_axz(double complex *p,double complex *dpdn,double *r,double *n,Apw *pw)
{
  double complex ce,dpdz;
  double arg;
  
  arg=pw->wd.k0*r[2];
  ce=cos(arg)+I*sin(arg);
  
  *p=pw->p*ce;
  dpdz=(*p)*I*pw->wd.k0;
  *dpdn=dpdz*n[2];
}
