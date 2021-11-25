#include "bem3_aw_b2.h"

int domain_id_m_dmda(double *rt,DMDA *ad);


int phi_s_dmda(double complex *phi,double *rt,int type,DMDA *ad)
{
  double complex CC[9],kc;
  double F;
  int did,s,sd,n;

  did=domain_id_m_dmda(rt,ad);

  *phi=0.0;
  F=0.0;
  kc=(double complex)ad->k0[did];

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];
    coef_rt(CC,rt,sd,kc,type,&(ad->bd));

    for(n=0;n<4;n++) *phi+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];

    F+=creal(CC[8]);
  }

  if(did==0) *phi/=1.0+F;
  else *phi/=F;

  return did;
}

int phi_t_dmda(double complex *phi,double *rt,int type,DMDA *ad)
{
  double complex p,v[3];
  int did;

  did=phi_s_dmda(phi,rt,type,ad);
  if(did==0){
    calc_maw_pv(&p,v,rt,&(ad->aw)); // multi_aw.h
    *phi+=-p/ad->aw.k2;
  }

  return did;
}

int phi_i_dmda(double complex *phi,double *rt,int type,DMDA *ad)
{
  double complex p,v[3];

  calc_maw_pv(&p,v,rt,&(ad->aw)); // multi_aw.h
  *phi=-p/ad->aw.k2;
  return domain_id_m_dmda(rt,ad);
}

int p_s_dmda(double complex *p,double *rt,int type,DMDA *ad)
{
  int did;

  did=phi_s_dmda(p,rt,type,ad);
  *p*=-ad->k2[did];
  return did;
}

int p_t_dmda(double complex *p,double *rt,int type,DMDA *ad)
{
  int did;

  did=phi_t_dmda(p,rt,type,ad);
  *p*=-ad->k2[did];
  return did;
}

int p_i_dmda(double complex *p,double *rt,int type,DMDA *ad)
{
  double complex v[3];

  calc_maw_pv(p,v,rt,&(ad->aw)); // multi_aw.h 
  return domain_id_m_dmda(rt,ad);
}

int pv_s_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad)
{
  double complex CC[9];
  double F,cf;
  int did,s,sd,l,n;

  did=domain_id_m_dmda(rt,ad);

  *p=0.0;
  for(l=0;l<3;l++) pv[l]=0.0;

  F=0.0;
  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];
    coef_rt(CC,rt,sd,ad->k0[did],type,&(ad->bd));

    for(n=0;n<4;n++){
      *p+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
      for(l=0;l<3;l++) pv[l]+=CC[n+0]*ad->bd.sb[did].dpv[s][n][l]-CC[n+4]*ad->bd.sb[did].pv[s][n][l];
    }

    F+=creal(CC[8]);
  }

  if(did==0){
    cf=1.0/(1.0+F);
    *p*=cf;
    for(l=0;l<3;l++) pv[l]*=cf;
  }
  else{
    cf=1.0/F;
    *p*=cf;
    for(l=0;l<3;l++) pv[l]*=cf;
  }
  *p*=-ad->k2[did];
  return did;
}

int pv_t_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad)
{
  double complex pi,pvi[3];
  int did,l;

  did=pv_s_dmda(p,pv,rt,type,ad);
  if(did==0){
    calc_maw_pv(&pi,pvi,rt,&(ad->aw)); // multi_aw.h
    *p+=pi;
    for(l=0;l<3;l++) pv[l]+=pvi[l];
  }

  return did;
}

int pv_i_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad)
{
  int did;

  did=domain_id_m_dmda(rt,ad);
  calc_maw_pv(p,pv,rt,&(ad->aw)); // multi_aw.h
  return did;
}


int pv_s_dbieq_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad)
{
  double complex CC[9],dCC[9],kc,P,dP[3];
  double v[3][3],F,dF[3];
  int did,l,sd,s,i,n;

  did=domain_id_m_dmda(rt,ad);
  kc=ad->k0[did];

  F=0.0;
  P=0.0;
  for(i=0;i<3;i++){
    dP[i]=0.0;

    for(l=0;l<3;l++) {
      if(i==l) v[i][l]=1.0;
      else v[i][l]=0.0;
    }
    dF[i]=0.0;
  }

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];

    for(i=0;i<3;i++){
      dcoef_rt(CC,dCC,rt,v[i],sd,kc,type,&(ad->bd));

      for(n=0;n<4;n++) dP[i]+=dCC[n+0]*ad->bd.sb[did].dP[s][n]-dCC[n+4]*ad->bd.sb[did].P[s][n];
      dF[i]+=creal(dCC[8]);
    }

    for(n=0;n<4;n++) P+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
    F+=creal(CC[8]);
  }

  if(did==0){
    P/=1.0+F;
    for(i=0;i<3;i++) dP[i]=(dP[i]-P*dF[i])/(1.0+F);
  }
  else {
    P/=F;
    for(i=0;i<3;i++) dP[i]=(dP[i]-P*dF[i])/F;
  }

  if(fabs(dF[0])>CBD_CDF || fabs(dF[1])>CBD_CDF || fabs(dF[2])>CBD_CDF){
    *p=0.0;
    for(i=0;i<3;i++) pv[i]=0.0;
    return -did;
  }

  *p=-ad->k2[did]*P;
  for(i=0;i<3;i++) pv[i]=-dP[i];

  return did;
}

int pv_t_dbieq_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad)
{
  double complex pi,pvi[3];
  int did,l;

  did=pv_s_dbieq_dmda(p,pv,rt,type,ad);
  if(did<0) return did;

  if(did==0){
    calc_maw_pv(&pi,pvi,rt,&(ad->aw)); // multi_aw.h
    *p+=pi;
    for(l=0;l<3;l++) pv[l]+=pvi[l];
  }

  return did;
}

int pv_i_dbieq_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad)
{
  return pv_i_dmda(p,pv,rt,type,ad);
}

void phi_s_bd(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad)
{
  double complex CC[9];
  double F,rt[3];
  int s,sd,n,td;

  *phi=0.0;
  F=0.0;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];

    coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did]+0.0*I,type,&(ad->bd));

    for(n=0;n<4;n++){
      *phi+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    *phi/=1.0+F;
  }
  else{
    *phi/=F;
  }
}

void phi_t_bd(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad)
{
  double complex CC[9],pi,v[3];
  double F,rt[3];
  int s,sd,n,td;

  *phi=0.0;
  F=0.0;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];

    coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did]+0.0*I,type,&(ad->bd));

    for(n=0;n<4;n++){
      *phi+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    *phi/=1.0+F;
    calc_maw_pv(&pi,v,rt,&(ad->aw)); // multi_aw.h
    *phi+=-pi/ad->aw.k2;
  }
  else{
    *phi/=F;
  }
}

void phi_i_bd(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad)
{
  double complex v[3];
  double rt[3];
  int td;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));
  calc_maw_pv(phi,v,rt,&(ad->aw)); // multi_aw.h
  *phi/=-ad->aw.k2;
}

void pv_s_bd_dmda(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad)
{
  double complex CC[9];
  double F,rt[3],cf;
  int s,sd,l,n,td;

  *p=0.0;
  for(l=0;l<3;l++) pv[l]=0.0;
  F=0.0;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];

    coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did],type,&(ad->bd));

    for(n=0;n<4;n++){
      *p+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
      for(l=0;l<3;l++) pv[l]+=CC[n+0]*ad->bd.sb[did].dpv[s][n][l]-CC[n+4]*ad->bd.sb[did].pv[s][n][l];
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    cf=1.0/(1.0+F);
    *p*=-cf*ad->k2[did];
    for(l=0;l<3;l++) pv[l]*=cf;
  }
  else{
    cf=1.0/F;
    *p*=-cf*ad->k2[did];
    for(l=0;l<3;l++) pv[l]*=cf;
  } 
}

void pv_t_bd_dmda(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad)
{
  double complex CC[9],pi,vi[3];
  double F,rt[3],cf;
  int s,sd,l,n,td;

  *p=0.0;
  for(l=0;l<3;l++) pv[l]=0.0;
  F=0.0;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];

    coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did],type,&(ad->bd));

    for(n=0;n<4;n++){
      *p+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
      for(l=0;l<3;l++) pv[l]+=CC[n+0]*ad->bd.sb[did].dpv[s][n][l]-CC[n+4]*ad->bd.sb[did].pv[s][n][l];
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    cf=1.0/(1.0+F);
    calc_maw_pv(&pi,vi,rt,&(ad->aw)); // multi_aw.h
    *p=-(*p)*cf*ad->k2[did]+pi;
    for(l=0;l<3;l++) pv[l]=cf*pv[l]+vi[l];
  }
  else{
    cf=1.0/F;
    *p*=-cf*ad->k2[did];
    for(l=0;l<3;l++) pv[l]*=cf;
  } 
}

void pv_i_bd_dmda(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad)
{
  int td;
  double rt[3];
  int i;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));

  if(did==0){
    calc_maw_pv(p,pv,rt,&(ad->aw)); // multi_aw.h
  }
  else {
    *p=0.0;
    for(i=0;i<3;i++) pv[i]=0.0;
  }
}
