/*
 * bil_coef.c
 *
 *  Created on: Dec 18, 2018
 *      Author: ohta
 */

#include "d3b1_elem.h"

double complex bil_sG_zeta(double zeta,void *tmp);
double complex bil_sG_eta(double eta,void *tmp);
double complex bil_sH_zeta(double zeta,void *tmp);
double complex bil_sH_eta(double eta,void *tmp);
double bil_sF_zeta(double zeta,void *tmp);
double bil_sF_eta(double eta,void *tmp);

double complex bil_sG_theta(double theta,void *tmp);
double complex bil_sG_r(double r,void *tmp);
double complex bil_sH_theta(double theta,void *tmp);
double complex bil_sH_r(double r,void *tmp);
double bil_sF_theta(double theta,void *tmp);
double bil_sF_r(double r,void *tmp);


void bil_coef_4p(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce;
  double vR[3],aR,i_aR,aW,tD,*vW;
  int sig,i,l;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  CC[8]=0.0;
  for(i=0;i<4;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];
    vW=bd->wen[s][i];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    aW=vabs_d(vW);
    ca=I*k*aR;
    if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
    else ce=cexp(ca);
    tD=(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;

    CC[i  ]=I_FP*ce*i_aR*aW; // G
    CC[i+4]=I_FP*(ca-1.0)*ce*tD; // H
    CC[8]+=bd->wt_44[i]*tD;
  }
  CC[8]*=I_FP;

  if(sig==1) bil_convert_CC(CC);
}

void bil_coef_9p(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np;
  int sig,i,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<4;ep++){
    for(i=0;i<9;i++){
      zeta=bd->zt_49[i];
      eta =bd->et_49[i];
      bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      aR=vabs_d(vR);
      i_aR=1.0/aR;
      aW=vabs_d(vW);
      ca=I*k*aR;
      if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
      else ce=cexp(ca);
      tD=(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;
      Np=bil_Nn(ep,zeta,eta);

      CC[ep  ]+=bd->wt_49[i]*ce*i_aR*aW*Np;
      CC[ep+4]+=bd->wt_49[i]*(ca-1.0)*ce*tD*Np;
      if(ep==0) CC[8]+=bd->wt_49[i]*tD;
    }
    CC[ep  ]*=I_FP;
    CC[ep+4]*=I_FP;
    if(ep==0) CC[8]*=I_FP;
  }

  if(sig==1) bil_convert_CC(CC);
}

void bil_coef_GL(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,tmpg,tmph;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,tmpf;
  int sig,i,j,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<4;ep++){
    for(i=0;i<GLN;i++){
      tmpg=0.0;
      tmph=0.0;
      tmpf=0.0;
      for(j=0;j<GLN;j++){
        zeta=bd->xli[i];
        eta =bd->xli[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        ca=I*k*aR;
        if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
        else ce=cexp(ca);
        tD=(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;
        Np=bil_Nn(ep,zeta,eta);

        tmpg+=bd->wli[j]*ce*i_aR*aW*Np;
        tmph+=bd->wli[j]*(ca-1.0)*ce*tD*Np;
        if(ep==0) tmpf+=bd->wli[j]*tD;
      }
      CC[ep  ]+=bd->wli[i]*tmpg;
      CC[ep+4]+=bd->wli[i]*tmph;
      if(ep==0) CC[8]+=bd->wli[i]*tmpf;
    }
    CC[ep  ]*=I_FP;
    CC[ep+4]*=I_FP;
    if(ep==0) CC[8]*=I_FP;
  }

  if(sig==1) bil_convert_CC(CC);
}

void bil_coef_GH(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,tmpg,tmph;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,tmpf;
  int sig,i,j,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<4;ep++){
    for(i=0;i<GHN;i++){
      tmpg=0.0;
      tmph=0.0;
      tmpf=0.0;
      for(j=0;j<GHN;j++){
        zeta=bd->xhi[i];
        eta =bd->xhi[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        ca=I*k*aR;
        if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
        else ce=cexp(ca);
        tD=(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;
        Np=bil_Nn(ep,zeta,eta);

        tmpg+=bd->whi[j]*ce*i_aR*aW*Np;
        tmph+=bd->whi[j]*(ca-1.0)*ce*tD*Np;
        if(ep==0) tmpf+=bd->whi[j]*tD;
      }
      CC[ep  ]+=bd->whi[i]*tmpg;
      CC[ep+4]+=bd->whi[i]*tmph;
      if(ep==0) CC[8]+=bd->whi[i]*tmpf;
    }
    CC[ep  ]*=I_FP;
    CC[ep+4]*=I_FP;
    if(ep==0) CC[8]*=I_FP;
  }

  if(sig==1) bil_convert_CC(CC);
}

void bil_coef_DE(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double err;
  int i,sig,l;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  bil_copy_elem_const_rw(td.cr,td.cw,s,bd);
  for(l=0;l<3;l++) td.rt[l]=rt[l];

  for(i=0;i<9;i++) CC[i]=0.0; //
  for(i=0;i<4;i++){ // G
    td.nid=i;
    CC[i]=I_FP*deintz(bil_sG_zeta,-1.0,1.0,&td,IEPS,&err);
    if(err<0){ printf("bil_coef.c, DE integration error, bil_coef_DE(), bil_sG_zeta(), nid=%d. Exit...\n",td.nid);  exit(1);  }
  }

  if(bil_check_on_plane(td.rt,td.cr,td.cw)==1){
    for(i=0;i<5;i++)CC[i+4]=0.0; // H,F=0
  }
  else {
    for(i=0;i<4;i++){ // H
      td.nid=i;
      CC[i+4]=I_FP*deintz(bil_sH_zeta,-1.0,1.0,&td,IEPS,&err);
      if(err<0){ printf("bil_coef.c, DE integration error, bil_coef_DE(), bil_sH_zeta(), nid=%d. Exit...\n",td.nid);  exit(1);  }
    }
    // F
    CC[8]=I_FP*deintd(bil_sF_zeta,-1.0,1.0,&td,IEPS,&err);
    if(err<0){ printf("bil_coef.c, DE integration error, bil_coef_DE(), bil_sF_zeta(). Exit...\n");  exit(1);  }
  }

  if(sig==1) bil_convert_CC(CC);
}

void bil_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex tC;
  double r[5],th[5],err,bdC,tD;
  int sig,l,m,i,j;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  td.zeta_t=zeta_t;
  td.eta_t = eta_t;
  if(SW_BDR_BIEQ==0){
    td.gn=GLN;
    td.xi=bd->xli;
    td.wi=bd->wli;
  }
  else {
    td.gn=GHN;
    td.xi=bd->xhi;
    td.wi=bd->whi;
  }
  bil_copy_elem_const_rw(td.cr,td.cw,s,bd);
  bil_calc_node_r_th(r,th,zeta_t,eta_t);

  // G
  for(i=0;i<4;i++){
    td.nid=i;
    CC[i]=0.0;
    for(j=0;j<4;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];

      if(SW_BDT_BIEQ==0){
        tC=0.0;
        for(m=0;m<GHN;m++) tC+=( bil_sG_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +bil_sG_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tC*=0.25*(th[j+1]-th[j]);
        CC[i]+=tC;
      }
      else {
        CC[i]+=deintz(bil_sG_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_coef.c, bil_coef_bd(), bil_sG_theta()! err=%f Exit...\n",err);  exit(1);  }
      }
    }
    CC[i]*=I_FP;
  }

  if(bil_check_plane(td.cr,td.cw)==1){ // element is plane
    for(l=0;l<5;l++)CC[l+4]=0.0;
  }
  else {
    // H,F
    bdC=td.cr[0][1]*td.cw[0][2]+td.cr[1][1]*td.cw[1][2]+td.cr[2][1]*td.cw[2][2];
    for(i=0;i<4;i++){
      td.nid=i;
      CC[i+4]=0.0;
      if(i==0)CC[8]=0.0;
      for(j=0;j<4;j++){
        td.r0=r[j];      td.r1=r[j+1];
        td.th0=th[j];    td.th1=th[j+1];
        // H
        if(SW_BDT_BIEQ==0){
          tC=0.0;
          for(m=0;m<GHN;m++) tC+=( bil_sH_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                  +bil_sH_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
          tC*=0.25*(th[j+1]-th[j]);
          CC[i+4]+=tC;
        }
        else {
          CC[i+4]+=deintz(bil_sH_theta,th[j],th[j+1],&td,IEPS,&err);
          if(err<0.0){ printf("DE integration error bil_coef.c, bil_coef_bd(), bil_sH_theta()! err=%f Exit...\n",err);  exit(1);  }
        }
        // F
        if(i==0){
          if(SW_BDT_BIEQ==0){
            tD=0.0;
            for(m=0;m<GHN;m++) tD+=( bil_sF_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                    +bil_sF_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
            tD*=0.25*(th[j+1]-th[j]);
            CC[8]+=tD;
          }
          else {
            CC[8]+=deintd(bil_sF_theta,th[j],th[j+1],&td,IEPS,&err);
            if(err<0.0){ printf("DE integration error bil_coef.c, bil_coef_bd(), bil_sF_theta()! err=%f Exit...\n",err);  exit(1);  }
          }
        }
      }
      CC[i+4]*=I_FP*bdC;
      if(i==0) CC[8]*=I_FP*bdC;
    }
  }

  if(sig==1) bil_convert_CC(CC);
}

//////////////////////////////////////////////////////////////////////

double complex bil_sG_zeta(double zeta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta=zeta;
  res= deintz(bil_sG_eta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("bil_coef.c, DE integration error, bil_sG_zeta(). Exit...\n");  exit(1);  }
  return res;
}

double complex bil_sG_eta(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ce;
  double r[3],w[3],aR,aW;

  bil_rw_zeta_eta(r,w,td->zeta,eta,td->cr,td->cw);
  aR=vabs_2dm(r,td->rt);
  aW=vabs_d(w);
  if(cimag(td->k)==0.0) ce=cos(creal(td->k)*aR)+I*sin(creal(td->k)*aR);
  else ce=cexp(I*td->k*aR);

  return ce/aR*aW*bil_Nn(td->nid,td->zeta,eta);
}

double complex bil_sH_zeta(double zeta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta=zeta;
  res= deintz(bil_sH_eta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("bil_coef.c, DE integration error, bil_sH_zeta(). err=%15.14e. Exit...\n",err);  exit(1);  }
  return res;
}

double complex bil_sH_eta(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ikr,ce;
  double r[3],w[3],aR,i_aR,vR[3];

  bil_rw_zeta_eta(r,w,td->zeta,eta,td->cr,td->cw);
  vsub_d(vR,r,td->rt);
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  ikr=I*td->k*aR;
  if(cimag(td->k)==0.0) ce=cos(creal(td->k)*aR)+I*sin(creal(td->k)*aR);
  else ce=cexp(ikr);

  return (ikr-1.0)*ce*i_aR*i_aR*i_aR*vdot_d(vR,w)*bil_Nn(td->nid,td->zeta,eta);
}

double bil_sF_zeta(double zeta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double res,err;

  td->zeta=zeta;
  res= deintd(bil_sF_eta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("bil_coef.c, DE integration error, bil_sF_zeta(). Exit...\n");  exit(1);  }
  return res;
}

double bil_sF_eta(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double r[3],w[3],aR,i_aR,vR[3];

  bil_rw_zeta_eta(r,w,td->zeta,eta,td->cr,td->cw);
  vsub_d(vR,r,td->rt);
  aR=vabs_d(vR);
  i_aR=1.0/aR;

  return i_aR*i_aR*i_aR*vdot_d(vR,w);
}

double complex bil_sG_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("integral region error bil_coef.c, bil_sG_theta()! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  ret=0.0;
  for(i=0;i<td->gn;i++) ret+=bil_sG_r(0.5*rm*(td->xi[i]+1.0),td)*td->wi[i];
  ret*=0.5*rm;

  return ret;
}

double complex bil_sG_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex ce;
  double K[3],W[3],aK,aW,zeta,eta;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
          +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
    W[i]=td->cw[i][0]+td->cw[i][1]*zeta+td->cw[i][2]*eta;
  }
  aK=vabs_d(K);
  aW=vabs_d(W);

  if(cimag(td->k)==0.0) ce=cos(creal(td->k)*aK*r)+I*sin(creal(td->k)*aK*r);
  else ce=cexp(I*td->k*aK*r);

  return ce/aK*aW*bil_Nn(td->nid,zeta,eta);
}

double complex bil_sH_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("integral region error bil_coef.c, bil_sH_theta()! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  ret=0.0;
  for(i=0;i<td->gn;i++) ret+=bil_sH_r(0.5*rm*(td->xi[i]+1.0),td)*td->wi[i];
  ret*=0.5*rm;

  return ret;
}

double complex bil_sH_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kR,ce;
  double K[3],aK,zeta,eta;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++) K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
                        +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  aK=vabs_d(K);
  kR=td->k*r*aK;
  if(cimag(td->k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
  else ce=cexp(I*kR);

  return (I*kR-1.0)*ce/(aK*aK*aK)*td->cth*td->sth*bil_Nn(td->nid,zeta,eta);
}

double bil_sF_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double ret;
  double rm;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("integral region error bil_coef.c, bil_sH_theta()! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  ret=0.0;
  for(i=0;i<td->gn;i++) ret+=bil_sF_r(0.5*rm*(td->xi[i]+1.0),td)*td->wi[i];
  ret*=0.5*rm;

  return ret;
}

double bil_sF_r(double r,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double K[3],aK;
  int i;

  for(i=0;i<3;i++) K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
                        +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  aK=vabs_d(K);

  return 1.0/(aK*aK*aK)*td->cth*td->sth;
}



