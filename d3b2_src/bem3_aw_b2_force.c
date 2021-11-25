#include "bem3_aw_b2.h"

int mo_force_FN(double *F,double *N,double *rc,int oid,int type,MOBJ *mo)
{
  void mo_bil_force_4p(double *F,double *N,double *rc,int oid,int s,MOBJ *mo);
  void mo_bil_force_9p(double *F,double *N,double *rc,int oid,int s,MOBJ *mo);
  void mo_lit_force_4p(double *F,double *N,double *rc,int oid,int s,MOBJ *mo);
  void mo_lit_force_7p(double *F,double *N,double *rc,int oid,int s,MOBJ *mo);

  double tf[3],tn[3],Fx,Fy,Fz,Nx,Ny,Nz;
  int t,td,nc;
  
  Fx=0.0;    Fy=0.0;    Fz=0.0;
  Nx=0.0;    Ny=0.0;    Nz=0.0;
  nc=0;

  // object id check
  if(oid<0) return -1;
  if(oid>=mo->N) return -1;

  #pragma omp parallel for schedule(dynamic) reduction(+:Fx,Fy,Fz,Nx,Ny,Nz,nc) private(td,tf,tn)
  for(t=1;t<=mo->md[oid].bd.sb[0].Ne;t++){
    td=mo->md[oid].bd.sb[0].sid[t];
    if( ELT4==check_element_type(td,&(mo->md[oid].bd)) ){
      if(type==0) mo_bil_force_4p(tf,tn,rc,oid,t,mo);
      else        mo_bil_force_9p(tf,tn,rc,oid,t,mo);
    }
    else {
      if(type==0) mo_lit_force_4p(tf,tn,rc,oid,t,mo);
      else        mo_lit_force_7p(tf,tn,rc,oid,t,mo);
    }
    Fx+=tf[0];
    Fy+=tf[1];
    Fz+=tf[2];
    Nx+=tn[0];
    Ny+=tn[1];
    Nz+=tn[2];
    nc++;
  }

  F[0]=Fx;
  F[1]=Fy;
  F[2]=Fz;
  N[0]=Nx;
  N[1]=Ny;
  N[2]=Nz;  

  if(nc!=0) return 0;
  else return -1;
}

/////////////////////////////////////////////////////////////////////////
void mo_bil_force_4p(double *F,double *N,double *rc,int oid,int s,MOBJ *mo)
{
  double complex p,v[3];
  double r[3],w[3],T[9],c1,c2,c3,t1;
  int sd,asd,i,j,l,m;

  c1=0.25/(mo->md[0].aw.c0*mo->md[0].aw.c0*mo->md[0].aw.rho0); // 1/(4*K0)
  c2=0.25*mo->md[0].aw.rho0; // rho0/4
  c3=0.50*mo->md[0].aw.rho0; // rho0/2

  sd=mo->md[oid].bd.sb[0].sid[s];
  asd=abs(sd);
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<4;i++){
    for(j=0;j<3;j++){
      r[j]=mo->md[oid].bd.ren[asd][i][j];
      w[j]=mo->md[oid].bd.wen[asd][i][j];
    }
    mo_pv_t_bd(&p,v,oid,0,s,mo->md[oid].bd.zt_44[i],mo->md[oid].bd.et_44[i],0,mo);
    
    // radiation stress tensor
    t1=c1*creal(p*conj(p))-c2*creal(v[0]*conj(v[0])+v[1]*conj(v[1])+v[2]*conj(v[2]));
    for(l=0;l<3;l++) for(m=0;m<3;m++) T[l*3+m]=c3*creal(v[l]*conj(v[m]));
    T[0*3+0]+=t1;
    T[1*3+1]+=t1;
    T[2*3+2]+=t1; 
    // force
    F[0]+=(T[0*3+0]*w[0]+T[0*3+1]*w[1]+T[0*3+2]*w[2])*mo->md[oid].bd.wt_44[i];
    F[1]+=(T[1*3+0]*w[0]+T[1*3+1]*w[1]+T[1*3+2]*w[2])*mo->md[oid].bd.wt_44[i];
    F[2]+=(T[2*3+0]*w[0]+T[2*3+1]*w[1]+T[2*3+2]*w[2])*mo->md[oid].bd.wt_44[i]; 
    // torque
    N[0]+=( ((r[1]-rc[1])*T[2*3+0]-(r[2]-rc[2])*T[1*3+0])*w[0]
           +((r[1]-rc[1])*T[2*3+1]-(r[2]-rc[2])*T[1*3+1])*w[1]
           +((r[1]-rc[1])*T[2*3+2]-(r[2]-rc[2])*T[1*3+2])*w[2])*mo->md[oid].bd.wt_44[i];
    N[1]+=( ((r[2]-rc[2])*T[0*3+0]-(r[0]-rc[0])*T[2*3+0])*w[0]
           +((r[2]-rc[2])*T[0*3+1]-(r[0]-rc[0])*T[2*3+1])*w[1]
           +((r[2]-rc[2])*T[0*3+2]-(r[0]-rc[0])*T[2*3+2])*w[2])*mo->md[oid].bd.wt_44[i];
    N[2]+=( ((r[0]-rc[0])*T[1*3+0]-(r[1]-rc[1])*T[0*3+0])*w[0]
           +((r[0]-rc[0])*T[1*3+1]-(r[1]-rc[1])*T[0*3+1])*w[1]
           +((r[0]-rc[0])*T[1*3+2]-(r[1]-rc[1])*T[0*3+2])*w[2])*mo->md[oid].bd.wt_44[i];
  }
}

void mo_bil_force_9p(double *F,double *N,double *rc,int oid,int s,MOBJ *mo)
{
  double complex p,v[3];
  double r[3],w[3],T[9],c1,c2,c3,t1,cr[3][4],cw[3][3];
  int sd,asd,i,j,l,m;

  c1=0.25/(mo->md[0].aw.c0*mo->md[0].aw.c0*mo->md[0].aw.rho0); // 1/(4*K0)
  c2=0.25*mo->md[0].aw.rho0; // rho0/4
  c3=0.50*mo->md[0].aw.rho0; // rho0/2

  sd=mo->md[oid].bd.sb[0].sid[s];
  asd=abs(sd);
  bil_copy_elem_const_rw(cr,cw,asd,&(mo->md[oid].bd));
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<9;i++){
    bil_rw_zeta_eta(r,w,mo->md[oid].bd.zt_49[i],mo->md[oid].bd.et_49[i],cr,cw);
    mo_pv_t_bd(&p,v,oid,0,s,mo->md[oid].bd.zt_49[i],mo->md[oid].bd.et_49[i],1,mo);
    // radiation stress tensor
    t1=c1*creal(p*conj(p))-c2*creal(v[0]*conj(v[0])+v[1]*conj(v[1])+v[2]*conj(v[2]));
    for(l=0;l<3;l++) for(m=0;m<3;m++) T[l*3+m]=c3*creal(v[l]*conj(v[m]));
    T[0*3+0]+=t1;
    T[1*3+1]+=t1;
    T[2*3+2]+=t1; 
    // force
    F[0]+=(T[0*3+0]*w[0]+T[0*3+1]*w[1]+T[0*3+2]*w[2])*mo->md[oid].bd.wt_49[i];
    F[1]+=(T[1*3+0]*w[0]+T[1*3+1]*w[1]+T[1*3+2]*w[2])*mo->md[oid].bd.wt_49[i];
    F[2]+=(T[2*3+0]*w[0]+T[2*3+1]*w[1]+T[2*3+2]*w[2])*mo->md[oid].bd.wt_49[i]; 
    // torque
    N[0]+=( ((r[1]-rc[1])*T[2*3+0]-(r[2]-rc[2])*T[1*3+0])*w[0]
           +((r[1]-rc[1])*T[2*3+1]-(r[2]-rc[2])*T[1*3+1])*w[1]
           +((r[1]-rc[1])*T[2*3+2]-(r[2]-rc[2])*T[1*3+2])*w[2])*mo->md[oid].bd.wt_49[i];
    N[1]+=( ((r[2]-rc[2])*T[0*3+0]-(r[0]-rc[0])*T[2*3+0])*w[0]
           +((r[2]-rc[2])*T[0*3+1]-(r[0]-rc[0])*T[2*3+1])*w[1]
           +((r[2]-rc[2])*T[0*3+2]-(r[0]-rc[0])*T[2*3+2])*w[2])*mo->md[oid].bd.wt_49[i];
    N[2]+=( ((r[0]-rc[0])*T[1*3+0]-(r[1]-rc[1])*T[0*3+0])*w[0]
           +((r[0]-rc[0])*T[1*3+1]-(r[1]-rc[1])*T[0*3+1])*w[1]
           +((r[0]-rc[0])*T[1*3+2]-(r[1]-rc[1])*T[0*3+2])*w[2])*mo->md[oid].bd.wt_49[i];
  }
}

void mo_lit_force_4p(double *F,double *N,double *rc,int oid,int s,MOBJ *mo)
{
  double complex p,v[3];
  double r[3],w[3],t1,T[9],c1,c2,c3;
  int i,j,sd,asd,l,m;
  
  c1=0.25/(mo->md[0].aw.c0*mo->md[0].aw.c0*mo->md[0].aw.rho0); // 1/(4*K0)
  c2=0.25*mo->md[0].aw.rho0; // rho0/4
  c3=0.50*mo->md[0].aw.rho0; // rho0/2
  
  sd=mo->md[oid].bd.sb[0].sid[s];
  asd=abs(sd);

  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<4;i++){
    for(j=0;j<3;j++){
      r[j]=mo->md[oid].bd.ren[asd][i][j];
      w[j]=mo->md[oid].bd.wen[asd][i][j];
    }
    mo_pv_t_bd(&p,v,oid,0,s,mo->md[oid].bd.zt_34[i],mo->md[oid].bd.et_34[i],0,mo);
    
    // radiation stress tensor
    t1=c1*creal(p*conj(p))-c2*creal(v[0]*conj(v[0])+v[1]*conj(v[1])+v[2]*conj(v[2]));
    for(l=0;l<3;l++) for(m=0;m<3;m++) T[l*3+m]=c3*creal(v[l]*conj(v[m]));
    T[0*3+0]+=t1;
    T[1*3+1]+=t1;
    T[2*3+2]+=t1; 
    // force
    F[0]+=(T[0*3+0]*w[0]+T[0*3+1]*w[1]+T[0*3+2]*w[2])*mo->md[oid].bd.wt_34[i];
    F[1]+=(T[1*3+0]*w[0]+T[1*3+1]*w[1]+T[1*3+2]*w[2])*mo->md[oid].bd.wt_34[i];
    F[2]+=(T[2*3+0]*w[0]+T[2*3+1]*w[1]+T[2*3+2]*w[2])*mo->md[oid].bd.wt_34[i]; 
    // torque
    N[0]+=( ((r[1]-rc[1])*T[2*3+0]-(r[2]-rc[2])*T[1*3+0])*w[0]
           +((r[1]-rc[1])*T[2*3+1]-(r[2]-rc[2])*T[1*3+1])*w[1]
           +((r[1]-rc[1])*T[2*3+2]-(r[2]-rc[2])*T[1*3+2])*w[2])*mo->md[oid].bd.wt_34[i];
    N[1]+=( ((r[2]-rc[2])*T[0*3+0]-(r[0]-rc[0])*T[2*3+0])*w[0]
           +((r[2]-rc[2])*T[0*3+1]-(r[0]-rc[0])*T[2*3+1])*w[1]
           +((r[2]-rc[2])*T[0*3+2]-(r[0]-rc[0])*T[2*3+2])*w[2])*mo->md[oid].bd.wt_34[i];
    N[2]+=( ((r[0]-rc[0])*T[1*3+0]-(r[1]-rc[1])*T[0*3+0])*w[0]
           +((r[0]-rc[0])*T[1*3+1]-(r[1]-rc[1])*T[0*3+1])*w[1]
           +((r[0]-rc[0])*T[1*3+2]-(r[1]-rc[1])*T[0*3+2])*w[2])*mo->md[oid].bd.wt_34[i];
  }

  for(j=0;j<3;j++){
    F[j]*=0.5;
    N[j]*=0.5;
  }
}

void mo_lit_force_7p(double *F,double *N,double *rc,int oid,int s,MOBJ *mo)
{
  double complex p,v[3];
  double r[3],w[3],cr[3][4],cw[3][3],c1,c2,c3,t1,T[9];
  int i,j,asd,sd,l,m;
  
  c1=0.25/(mo->md[0].aw.c0*mo->md[0].aw.c0*mo->md[0].aw.rho0); // 1/(4*K0)
  c2=0.25*mo->md[0].aw.rho0; // rho0/4
  c3=0.50*mo->md[0].aw.rho0; // rho0/2
  
  sd=mo->md[oid].bd.sb[0].sid[s];
  asd=abs(sd);
  lit_copy_elem_const_rw(cr,cw,asd,&(mo->md[oid].bd));

  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<7;i++){
    lit_rw_zeta_eta(r,w,mo->md[oid].bd.zt_37[i],mo->md[oid].bd.et_37[i],cr,cw);
    mo_pv_t_bd(&p,v,oid,0,s,mo->md[oid].bd.zt_37[i],mo->md[oid].bd.et_37[i],1,mo);
    
    // radiation stress tensor
    t1=c1*creal(p*conj(p))-c2*creal(v[0]*conj(v[0])+v[1]*conj(v[1])+v[2]*conj(v[2]));
    for(l=0;l<3;l++) for(m=0;m<3;m++) T[l*3+m]=c3*creal(v[l]*conj(v[m]));
    T[0*3+0]+=t1;
    T[1*3+1]+=t1;
    T[2*3+2]+=t1; 
    // force
    F[0]+=(T[0*3+0]*w[0]+T[0*3+1]*w[1]+T[0*3+2]*w[2])*mo->md[oid].bd.wt_37[i];
    F[1]+=(T[1*3+0]*w[0]+T[1*3+1]*w[1]+T[1*3+2]*w[2])*mo->md[oid].bd.wt_37[i];
    F[2]+=(T[2*3+0]*w[0]+T[2*3+1]*w[1]+T[2*3+2]*w[2])*mo->md[oid].bd.wt_37[i]; 
    // torque
    N[0]+=( ((r[1]-rc[1])*T[2*3+0]-(r[2]-rc[2])*T[1*3+0])*w[0]
           +((r[1]-rc[1])*T[2*3+1]-(r[2]-rc[2])*T[1*3+1])*w[1]
           +((r[1]-rc[1])*T[2*3+2]-(r[2]-rc[2])*T[1*3+2])*w[2])*mo->md[oid].bd.wt_37[i];
    N[1]+=( ((r[2]-rc[2])*T[0*3+0]-(r[0]-rc[0])*T[2*3+0])*w[0]
           +((r[2]-rc[2])*T[0*3+1]-(r[0]-rc[0])*T[2*3+1])*w[1]
           +((r[2]-rc[2])*T[0*3+2]-(r[0]-rc[0])*T[2*3+2])*w[2])*mo->md[oid].bd.wt_37[i];
    N[2]+=( ((r[0]-rc[0])*T[1*3+0]-(r[1]-rc[1])*T[0*3+0])*w[0]
           +((r[0]-rc[0])*T[1*3+1]-(r[1]-rc[1])*T[0*3+1])*w[1]
           +((r[0]-rc[0])*T[1*3+2]-(r[1]-rc[1])*T[0*3+2])*w[2])*mo->md[oid].bd.wt_37[i];
  }

  for(j=0;j<3;j++){
    F[j]*=0.5;
    N[j]*=0.5;
  }  
}
