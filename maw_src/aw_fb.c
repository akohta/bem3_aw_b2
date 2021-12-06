#include "aw_fb.h"

void read_data_afb(char *fname, Afb *pw)
{
  FILE *fp;
  char buf[256]="";
  double td1,td2;
  int ti;
  
  if((fp=fopen(fname,"rt"))==NULL){ 
    printf("failed to read %s. no beams defined. Exit...\n",fname);
    exit(1);
  }
  
  if(fgets(buf,256,fp)==NULL){
    printf("aw_fb.c, read_data_apw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("aw_fb.c, read_data_apw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the rho0. exit...\n");
    exit(1);
  }
  pw->rho0=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the c0. exit...\n");
    exit(1);
  }
  pw->c0  =td1;
  if(fscanf(fp,"%lf\n",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the f. exit...\n");
    exit(1);
  }
  pw->f   =td1;
  if(fgets(buf,256,fp)==NULL){
    printf("aw_fb.c, read_data_apw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the real(p). exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&td2)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the imag(p). exit...\n");
    exit(1);  
  }
  pw->p=td1+I*td2;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the B. exit...\n");
    exit(1);
  }
  pw->B=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the F. exit...\n");
    exit(1);
  }
  pw->F=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the tv[0]. exit...\n");
    exit(1);
  }
  pw->tv[0]=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the tv[1]. exit...\n");
    exit(1);
  }
  pw->tv[1]=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the tv[2]. exit...\n");
    exit(1);
  }
  pw->tv[2]=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the theta. exit...\n");
    exit(1);
  }
  pw->theta=td1;
  if(fscanf(fp,"%lf",&td1)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the psi. exit...\n");
    exit(1);
  }
  pw->psi  =td1;
  if(fscanf(fp,"%d",&ti)!=1){
    printf("aw_fb.c, read_data_apw(), failed to read the nn. exit...\n");
    exit(1);
  }
  pw->nn=ti;
  
  fclose(fp); 
}

void print_data_afb(Afb *pw)
{
  double Pz,omega,k0;
  
  printf("-- focused beam --\n");
  printf("medium density                       [kg/m^3] : %15.14g\n",pw->rho0);
  printf("speed of sound in the medium            [m/s] : %15.14g\n",pw->c0);
  printf("frequency                                [Hz] : %15.14g\n",pw->f);
  printf("sound pressure amplitude                 [Pa] :%7.6g+%7.6gI\n",creal(pw->p),cimag(pw->p));
  printf("radius of focused transducer              [m] : %15.14g\n",pw->B);
  printf("curvature radius of the transducer        [m] : %15.14g\n",pw->F);
  printf("x-component of translation vector         [m] : %15.14g\n",pw->tv[0]);
  printf("y-component of translation vector         [m] : %15.14g\n",pw->tv[1]);
  printf("z-component of translation vector         [m] : %15.14g\n",pw->tv[2]);
  printf("rotation parameter theta                [rad] : %15.14g\n",pw->theta);
  printf("rotation parameter psi                  [rad] : %15.14g\n",pw->psi);
  printf("sampling number for Gauss-Legendre quadrature : %15d\n",pw->nn);
  
  printf("-- additional information --\n");
  printf("wavelength                                [m] : %15.14g\n",pw->c0/pw->f);
  omega=2.0*M_PI*pw->f;
  k0=omega/pw->c0;
  Pz=2.0*M_PI*M_PI*creal(pw->p*conj(pw->p))*k0*pw->F/(omega*pw->rho0)*(pw->F-sqrt(pw->F*pw->F-pw->B*pw->B));
  printf("energy passing through a plane orthogonal to the acoustic axis [W]\n");
  printf("                                              : %15.14g\n",Pz);
  printf("\n");
  
}

void setup_afb(Afb *pw)
{
  double ct,st,cp,sp,*xt,h,apd;
  int i,j,s,nnc,ss;

  pw->wd.omega=2.0*M_PI*pw->f;
  pw->wd.K0=pw->c0*pw->c0*pw->rho0;
  pw->wd.k1=I*pw->wd.omega/pw->wd.K0;
  pw->wd.k2=I*pw->wd.omega*pw->rho0;
  pw->wd.k0=pw->wd.omega/pw->c0;
  pw->wd.lambda0=pw->c0/pw->f;
  if(pw->F <= pw->B){
    printf("focused transducer radius must be smaller than curvature radius. Exit...\n");
    exit(1);
  }
  pw->wd.alpha=asin(pw->B/pw->F);

  ct=cos(pw->theta);
  st=sin(pw->theta);
  cp=cos(pw->psi);
  sp=sin(pw->psi);
  pw->wd.R[0*3+0]=ct*cp*cp+sp*sp;  pw->wd.R[0*3+1]=sp*cp*(ct-1.0);  pw->wd.R[0*3+2]=st*cp;
  pw->wd.R[1*3+0]=sp*cp*(ct-1.0);  pw->wd.R[1*3+1]=ct*sp*sp+cp*cp;  pw->wd.R[1*3+2]=st*sp;
  pw->wd.R[2*3+0]=-st*cp;          pw->wd.R[2*3+1]=-st*sp;          pw->wd.R[2*3+2]=ct;
  
  // gauss-legendre
  pw->wd.ct=(double *)m_alloc2(pw->nn,sizeof(double),"aw_fb.c,setup_afb(),pw->wd.ct");
  pw->wd.st=(double *)m_alloc2(pw->nn,sizeof(double),"aw_fb.c,setup_afb(),pw->wd.st");
  pw->wd.wt=(double *)m_alloc2(pw->nn,sizeof(double),"aw_fb.c,setup_afb(),pw->wd.wt");
  xt=(double *)m_alloc2(pw->nn,sizeof(double),"aw_fb.c,setup_afb(),pw->wd.xt");
  gauleg( 0.0,pw->wd.alpha,xt,pw->wd.wt,pw->nn);
  for(i=0;i<pw->nn;i++){
    pw->wd.st[i]=sin(xt[i]);
    pw->wd.ct[i]=cos(xt[i]);
    pw->wd.wt[i]*=pw->wd.st[i];
  } 
  free(xt);
  
  // trapezoid
  nnc=(2<<(TRAP_MNUM))+1;
  pw->wd.cp=(double *)m_alloc2(nnc,sizeof(double),"aw_fb.c,setup_afb(),pw->wd.cp");
  pw->wd.sp=(double *)m_alloc2(nnc,sizeof(double),"aw_fb.c,setup_afb(),pw->wd.sp");
  h=2.0*M_PI;
  pw->wd.cp[0]=-1.0;
  pw->wd.sp[0]= 0.0;
  j=1;
  ss=1;
  for(i=0;i<=TRAP_MNUM;i++){
    h*=0.5;
    for(s=1;s<=ss;s++){
      apd=-M_PI+(double)(2*s-1)*h;
      pw->wd.cp[j]=cos(apd);
      pw->wd.sp[j]=sin(apd);
      j++;
    }
    ss*=2;
  }
}

void free_afb(Afb *pw)
{
  free(pw->wd.st);  free(pw->wd.ct);  free(pw->wd.wt);
  free(pw->wd.cp);  free(pw->wd.sp);
}

int calc_afb_pv(double complex *p,double complex *v,double *r,Afb *pw)
{
  int check_region(double *r0,Afb *pw);
  void afb_pv_axz(double complex *p,double complex *v,double *r,Afb *pw);
  
  double complex vt[3];
  double *R,rc[3],r0[3];
  int i,j,ret;
  
  R=pw->wd.R;

  for(i=0;i<3;i++) rc[i]=r[i]-pw->tv[i];
  r0[0]= R[0*3+0]*rc[0]+R[0*3+1]*rc[1]-R[0*3+2]*rc[2];
  r0[1]= R[1*3+0]*rc[0]+R[1*3+1]*rc[1]-R[1*3+2]*rc[2];
  r0[2]=-R[2*3+0]*rc[0]-R[2*3+1]*rc[1]+R[2*3+2]*rc[2]; 

  ret=check_region(r0,pw);
  if(ret==0){
    afb_pv_axz(p,vt,r0,pw);
    for(i=0;i<3;i++){
      v[i]=0.0;
      for(j=0;j<3;j++) v[i]+=R[i*3+j]*vt[j];
    }
  }
  else {
    *p=0.0;
    v[0]=0.0;
    v[1]=0.0;
    v[2]=0.0;
  } 
  return ret;
}

int calc_afb_dpdn(double complex *p,double complex *dpdn, double *r ,double *n, Afb *pw)
{
  int check_region(double *r0,Afb *pw);
  void afb_dpdn_axz(double complex *p,double complex *dpdn,double *r,double *n,Afb *pw);
  
  double *R,rc[3],r0[3],n0[3];
  int i,ret;

  R=pw->wd.R;

  for(i=0;i<3;i++) rc[i]=r[i]-pw->tv[i];
  r0[0]= R[0*3+0]*rc[0]+R[0*3+1]*rc[1]-R[0*3+2]*rc[2];
  r0[1]= R[1*3+0]*rc[0]+R[1*3+1]*rc[1]-R[1*3+2]*rc[2];
  r0[2]=-R[2*3+0]*rc[0]-R[2*3+1]*rc[1]+R[2*3+2]*rc[2]; 
  
  n0[0]= R[0*3+0]*n[0]+R[0*3+1]*n[1]-R[0*3+2]*n[2];
  n0[1]= R[1*3+0]*n[0]+R[1*3+1]*n[1]-R[1*3+2]*n[2];
  n0[2]=-R[2*3+0]*n[0]-R[2*3+1]*n[1]+R[2*3+2]*n[2]; 
  
  ret=check_region(r0,pw);
  if(ret==0) afb_dpdn_axz(p,dpdn,r0,n0,pw); 
  else {
    *p=0.0;
    *dpdn=0.0;
  }
  return ret;
}

///////////////////////////////////////////////////////////////////////
int check_region(double *r0,Afb *pw)
{
  double dr,tr,th,rxy2,r,i_rxy,cp,sp,rd[3],d2;
  
  rxy2=r0[0]*r0[0]+r0[1]*r0[1];
  dr=0.5*TD_TH;
  tr=acos((2.0*pw->F*pw->F-dr*dr)/(2.0*pw->F*pw->F));
  th=atan2(sqrt(r0[0]*r0[0]+r0[1]*r0[1]),r0[2]);
  
  if(fabs(th)<= M_PI-pw->wd.alpha-tr) return 0; // external of transducer region 
  else{
    r=sqrt(rxy2+r0[2]*r0[2]);
    if(r <= pw->F-dr || r >= pw->F+dr ) return 0;
    else if(fabs(th)>= M_PI-pw->wd.alpha) return 1; // internal of transducer region 
    else {
      i_rxy=1.0/sqrt(rxy2);
      cp=r0[0]*i_rxy;
      sp=r0[1]*i_rxy;
      rd[0]=pw->B*cp;
      rd[1]=pw->B*sp;
      rd[2]=-pw->F*cos(pw->wd.alpha);
      d2=pow(r0[0]-rd[0],2)+pow(r0[1]-rd[1],2)+pow(r0[2]-rd[2],2);
      if(d2 >= dr*dr ) return 0;
      else return 1;
    }
  }
}

void afb_pv_axz(double complex *p,double complex *v,double *r,Afb *pw)
{
  void afb_int_psi_pv(double complex *pv,double ct,double st,double *r,Afb *pw);
  
  double complex pv[4]; // pv[0]=p, pv[1]=v_x, pv[2]=v_y, pv[3]=v_z
  double complex c1,c2;
  int j;
  
  c1=pw->p*pw->F*pw->F;
  c2=c1/pw->wd.k2;
  
  *p=0.0;
  v[0]=0.0;  v[1]=0.0;  v[2]=0.0;
  for(j=0;j<pw->nn;j++){
    afb_int_psi_pv(pv,pw->wd.ct[j],pw->wd.st[j],r,pw);
    
    *p+=pv[0]*pw->wd.wt[j];
    v[0]+=pv[1]*pw->wd.wt[j];    v[1]+=pv[2]*pw->wd.wt[j];   v[2]+=pv[3]*pw->wd.wt[j];
  }
  
  *p*=c1;
  v[0]*=c2;
  v[1]*=c2;
  v[2]*=c2;
}

void afb_int_psi_pv(double complex *pv,double ct,double st,double *r,Afb *pw)
{
  void afb_integrand_pv(double complex *pv,double cos_t,double sin_t,double cos_p,double sin_p,double *r,Afb *pw);
  
  double h,cpa,spa,abst,abss;
  double complex Ia[4],It[4],Is[4];
  int i,j,s,sc,cc;
  double f_order,cr_coef;
  double rho0=pw->rho0;
  double i_K0=1.0/pw->wd.K0;
  static int si=0;
  
  cc=0;
  h=2.0*M_PI;
  cpa=pw->wd.cp[0];  spa=pw->wd.sp[0];
  afb_integrand_pv(Ia,ct,st,cpa,spa,r,pw);
  abst=0.0;
  for(i=0;i<4;i++){
    It[i]=h*Ia[i];
    if(i==0) abst+=i_K0*creal(It[i]*conj(It[i])); // 1/K0*|p|^2
    else     abst+=rho0*creal(It[i]*conj(It[i])); // rho0*|v|^2
  }
  f_order=abst;
  j=1;
  sc=1;
  while(cc<=TRAP_MNUM){
    h*=0.5;
    for(i=0;i<4;i++) Is[i]=0.0;
    for(s=1;s<=sc;s++){
      cpa=pw->wd.cp[j];  spa=pw->wd.sp[j];
      j++;
      afb_integrand_pv(Ia,ct,st,cpa,spa,r,pw);
      for(i=0;i<4;i++) Is[i]+=Ia[i];
    }
    abss=0.0;
    for(i=0;i<4;i++){
      It[i]=0.5*It[i]+Is[i]*h;
      if(i==0) abss+=i_K0*creal(It[i]*conj(It[i]));
      else     abss+=rho0*creal(It[i]*conj(It[i]));
    }
    cr_coef=sqrt((double)sc)*f_order;
    if((cc>TRAP_LNUM && fabs(abst-abss) < TRAP_EPS*cr_coef) || (abss==0.0 && abst==0.0) ){
      break;
    }
    abst=abss;
    
    if(abst>f_order) f_order=abst;
    sc*=2;
    cc++;
  }
  if(cc>=TRAP_MNUM) {
    if(si==0){
      printf("aw_fb.c, afb_int_psi_pv(), reached the trapezoidral rule limit.\n");
      printf("The integral for psi (azimuthal part) is not converged, it may be inaccurate.\n");
      printf("Please increase the defined variable TRAP_MNUM in the aw_fb.h.\n");
      si++;
    }
  }
  for(i=0;i<4;i++) pv[i]=It[i];
}

void afb_integrand_pv(double complex *pv,double cos_t,double sin_t,double cos_p,double sin_p,double *r,Afb *pw)
{
  double complex ce,cv;
  double d,i_d,dr[3],arg;
  
  dr[0]=r[0]-pw->F*sin_t*cos_p;
  dr[1]=r[1]-pw->F*sin_t*sin_p;
  dr[2]=r[2]+pw->F*cos_t;
  d=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
  i_d=1.0/d;
  
  arg=pw->wd.k0*d;
  ce=cos(arg)+I*sin(arg);
  pv[0]=ce/d; // p
  
  cv=(I*pw->wd.k0*d-1.0)*i_d*i_d*pv[0];
  pv[1]=cv*dr[0]; // v_x
  pv[2]=cv*dr[1]; // v_y
  pv[3]=cv*dr[2]; // v_z
}

void afb_dpdn_axz(double complex *p,double complex *dpdn,double *r,double *n,Afb *pw)
{
   void afb_int_psi_pv(double complex *pv,double ct,double st,double *r,Afb *pw);
  
  double complex tpv[4],pv[4]; // pv[0]=p, pv[1]=v_x, pv[2]=v_y, pv[3]=v_z
  double complex c1;
  int j,i;
  
  c1=pw->p*pw->F*pw->F;
  
  for(i=0;i<4;i++) pv[i]=0.0;
  for(j=0;j<pw->nn;j++){
    afb_int_psi_pv(tpv,pw->wd.ct[j],pw->wd.st[j],r,pw);
    
    for(i=0;i<4;i++) pv[i]+=tpv[i]*pw->wd.wt[j];
  }
  
  *p=c1*pv[0];
  *dpdn=c1*(pv[1]*n[0]+pv[2]*n[1]+pv[3]*n[2]);
}
