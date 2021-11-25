/*
 * lit_coef.c
 *
 *  Created on: Dec 18, 2018
 *      Author: ohta
 */

#include "d3b1_elem.h"

void lit_sGL_zeta_eta(double complex *ret,TDATA *td);
double complex lit_sG_zeta_eta(TDATA *td);
double complex lit_sH_zeta_eta(TDATA *td);
double lit_sF_zeta_eta(TDATA *td);

double complex lit_sG_alpha(double alpha,void *tmp);
double complex lit_sG_beta(double beta,void *tmp);
double complex lit_sH_alpha(double alpha,void *tmp);
double complex lit_sH_beta(double beta,void *tmp);
double lit_sF_alpha(double alpha,void *tmp);
double lit_sF_beta(double beta,void *tmp);

double complex lit_sG_theta(double theta,void *tmp);
double complex lit_sG_r(double r,void *tmp);


void lit_coef_4p(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce;
  double vR[3],*vW,aR,i_aR,aW,tD,Nn;
  int sig,i,l;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  for(i=0;i<9;i++) CC[i]=0.0;
  vW=bd->wen[s][0];
  aW=vabs_d(vW);

  for(i=0;i<3;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    ca=I*k*aR;
    if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
    else ce=cexp(ca);
    tD=(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;

    CC[i+0]=bd->wt_34[i]*ce*i_aR;
    CC[i+4]=bd->wt_34[i]*(ca-1.0)*ce*tD;
    CC[8]+=bd->wt_34[i]*tD;
  }

  for(l=0;l<3;l++) vR[l]=bd->ren[s][3][l]-rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  ca=I*k*aR;
  if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
  else ce=cexp(ca);
  tD=(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;
  for(i=0;i<3;i++){
    Nn=lit_Nn(i,bd->zt_34[3],bd->et_34[3]);
    CC[i+0]+=bd->wt_34[3]*ce*i_aR*Nn;
    CC[i+4]+=bd->wt_34[3]*(ca-1.0)*ce*tD*Nn;
    if(i==0) CC[8]+=bd->wt_34[3]*tD;
  }

  for(i=0;i<3;i++){
    CC[i+0]*=I_EP*aW;
    CC[i+4]*=I_EP;
    if(i==0) CC[8]*=I_EP;
  }

  if(sig==1) lit_convert_CC(CC);
}

void lit_coef_7p(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce;
  double cr[3][4],cw[3][3],zeta,eta,vW[3],aW,vR[3],aR,i_aR,tD,Np;
  int sig,i,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  lit_w_zeta_eta(vW,0.0,0.0,cw);
  aW=vabs_d(vW);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<3;ep++){
    for(i=0;i<7;i++){
      zeta=bd->zt_37[i];
      eta =bd->et_37[i];
      lit_r_zeta_eta(vR,zeta,eta,cr);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      aR=vabs_d(vR);
      i_aR=1.0/aR;
      ca=I*k*aR;
      if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
      else ce=cexp(ca);
      tD=(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;
      Np=lit_Nn(ep,zeta,eta);

      CC[ep  ]+=bd->wt_37[i]*ce*i_aR*Np;
      CC[ep+4]+=bd->wt_37[i]*(ca-1.0)*ce*tD*Np;
      if(ep==0) CC[8]+=bd->wt_37[i]*tD;
    }
    CC[ep  ]*=I_EP*aW;
    CC[ep+4]*=I_EP;
    if(ep==0) CC[8]*=I_EP;
  }

  if(sig==1) lit_convert_CC(CC);
}

void lit_coef_GL(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex ret[3],tpj[3],tpi[3];
  double vW[3],a,b,wt;
  int sig,i,j,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vW,0.0,0.0,td.cw);
  td.aW=vabs_d(vW);
  for(l=0;l<3;l++) td.rt[l]=rt[l];

  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GLN;i++){
      for(l=0;l<3;l++) tpi[l]=0.0;
      a=bd->xli[i];
      for(j=0;j<GLN;j++){
        for(l=0;l<3;l++) tpj[l]=0.0;
        b=bd->xli[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        lit_sGL_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        lit_sGL_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        lit_sGL_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];

        for(l=0;l<3;l++) tpi[l]+=tpj[l]*wt*bd->wli[j];
      }
      CC[ep  ]+=tpi[0]*bd->wli[i];
      CC[ep+4]+=tpi[1]*bd->wli[i];
      if(ep==0) CC[8]+=tpi[2]*bd->wli[i];
    }

    CC[ep  ]*=0.0625*SQ3*I_SPSQ3*td.aW;
    CC[ep+4]*=0.0625*SQ3*I_SPSQ3;
    if(ep==0) CC[8]*=0.0625*SQ3*I_SPSQ3;
  }

  if(sig==1) lit_convert_CC(CC);
}

void lit_coef_GH(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex ret[3],tpj[3],tpi[3];
  double vW[3],a,b,wt;
  int sig,i,j,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vW,0.0,0.0,td.cw);
  td.aW=vabs_d(vW);
  for(l=0;l<3;l++) td.rt[l]=rt[l];

  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GHN;i++){
      for(l=0;l<3;l++) tpi[l]=0.0;
      a=bd->xhi[i];
      for(j=0;j<GHN;j++){
        for(l=0;l<3;l++) tpj[l]=0.0;
        b=bd->xhi[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        lit_sGL_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        lit_sGL_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        lit_sGL_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];

        for(l=0;l<3;l++) tpi[l]+=tpj[l]*wt*bd->whi[j];
      }
      CC[ep  ]+=tpi[0]*bd->whi[i];
      CC[ep+4]+=tpi[1]*bd->whi[i];
      if(ep==0) CC[8]+=tpi[2]*bd->whi[i];
    }

    CC[ep  ]*=0.0625*SQ3*I_SPSQ3*td.aW;
    CC[ep+4]*=0.0625*SQ3*I_SPSQ3;
    if(ep==0) CC[8]*=0.0625*SQ3*I_SPSQ3;
  }

  if(sig==1) lit_convert_CC(CC);
}

void lit_coef_DE(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double vW[3],err;
  int sig,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vW,0.0,0.0,td.cw);
  td.aW=vabs_d(vW);
  for(l=0;l<3;l++) td.rt[l]=rt[l];

  // G
  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    CC[ep+0]=0.0625*SQ3*I_SPSQ3*td.aW*deintz(lit_sG_alpha,-1.0,1.0,&td,IEPS,&err);
    if(err<0.0){ printf("lit_coef.c, DE integration error, lit_coef_DE(), lit_sG_alpha(). err=%g Exit...\n",err);  exit(1);  }
  }

  // H,F
  if(lit_check_on_plane(td.rt,td.cr,td.cw)==0){ // rt is not on element plane
    //H
    for(ep=0;ep<3;ep++){
      td.nid=ep;
      CC[ep+4]=0.0625*SQ3*I_SPSQ3*deintz(lit_sH_alpha,-1.0,1.0,&td,IEPS,&err);
      if(err<0.0){ printf("lit_coef.c, DE integration error, lit_coef_DE(), lit_sH_alpha(). err=%g Exit...\n",err);  exit(1);  }
    }
    //F
    CC[8]=0.0625*SQ3*I_SPSQ3*deintd(lit_sF_alpha,-1.0,1.0,&td,IEPS,&err);
    if(err<0.0){ printf("lit_coef.c, DE integration error, lit_coef_DE(), lit_sF_alpha(). err=%g Exit...\n",err);  exit(1);  }
  }
  else {
    for(ep=0;ep<3;ep++)  CC[ep+4]=0.0;
    CC[8]=0.0;
  }

  if(sig==1) lit_convert_CC(CC);
}

void lit_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex tC;
  double r[4],th[4],err,vW[3],aW;
  int sig,i,j,l,m;

  if(s>0){
    sig=0;
  }
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  td.zeta_t=zeta_t;
  td.eta_t =eta_t;
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
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vW,0.0,0.0,td.cw);
  aW=vabs_d(vW);
  lit_calc_node_r_th(r,th,td.zeta_t,td.eta_t);

  for(l=0;l<9;l++) CC[l]=0.0;
  // G
  for(i=0;i<3;i++){
    td.nid=i;
    CC[i]=0.0;
    for(j=0;j<3;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];

      if(SW_BDT_BIEQ==0){
        tC=0.0;
        for(m=0;m<GHN;m++) tC+=( lit_sG_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +lit_sG_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tC*=0.25*(th[j+1]-th[j]);
        CC[i]+=tC;
      }
      else {
        CC[i]+=deintz(lit_sG_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error lit_coef.c, lit_coef_bd(), lit_sG_theta()! err=%f Exit...\n",err);  exit(1);  }
      }
    }
    CC[i]*=I_SPSQ3*aW;
  }

  if(sig==1) lit_convert_CC(CC);
}


//////////////////////////////////////////////////////////////////////////////////////////////
void lit_sGL_zeta_eta(double complex *ret,TDATA *td)
{
  double complex ca,ce;
  double vR[3],aR,i_aR,tD,Np;
  int l;

  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);

  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  ca=I*td->k*aR;
  if(cimag(td->k)==0.0) ce=cos(creal(td->k)*aR)+I*sin(creal(td->k)*aR);
  else ce=cexp(ca);
  tD=(vR[0]*td->cw[0][0]+vR[1]*td->cw[1][0]+vR[2]*td->cw[2][0])*i_aR*i_aR*i_aR;
  Np=lit_Nn(td->nid,td->zeta,td->eta);

  ret[0]=ce*i_aR*Np;
  ret[1]=(ca-1.0)*ce*tD*Np;
  ret[2]=tD;
}

double complex lit_sG_zeta_eta(TDATA *td)
{
  double complex ca,ce;
  double vR[3],aR,i_aR,Np;
  int l;

  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);

  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  ca=I*td->k*aR;
  if(cimag(td->k)==0.0) ce=cos(creal(td->k)*aR)+I*sin(creal(td->k)*aR);
  else ce=cexp(ca);
  Np=lit_Nn(td->nid,td->zeta,td->eta);

  return ce*i_aR*Np;
}

double complex lit_sH_zeta_eta(TDATA *td)
{
  double complex ca,ce;
  double vR[3],aR,i_aR,tD,Np;
  int l;

  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);

  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  ca=I*td->k*aR;
  if(cimag(td->k)==0.0) ce=cos(creal(td->k)*aR)+I*sin(creal(td->k)*aR);
  else ce=cexp(ca);
  tD=(vR[0]*td->cw[0][0]+vR[1]*td->cw[1][0]+vR[2]*td->cw[2][0])*i_aR*i_aR*i_aR;
  Np=lit_Nn(td->nid,td->zeta,td->eta);

  return (ca-1.0)*ce*tD*Np;
}

double lit_sF_zeta_eta(TDATA *td)
{
  double vR[3],aR,i_aR,tD;
  int l;

  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);

  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  tD=(vR[0]*td->cw[0][0]+vR[1]*td->cw[1][0]+vR[2]*td->cw[2][0])*i_aR*i_aR*i_aR;

  return tD;
}



double complex lit_sG_alpha(double alpha,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta_t=alpha;
  res= deintz(lit_sG_beta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("lit_coef.c, DE integration error, lit_sG_alpha(). err=%g. Exit...\n",err);  exit(1);  }
  return res;
}

double complex lit_sG_beta(double b,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double a,wt;

  a=td->zeta_t;
  wt=1.0+0.25*(a+b);

  ret=0.0;
  // point 1
  td->zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
  td->eta =0.125*SQ3*(-a+b);
  ret+=lit_sG_zeta_eta(td);
  // part 2
  td->zeta=0.0625*(-3.0+a-5.0*b-a*b);
  td->eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
  ret+=lit_sG_zeta_eta(td);
  // part 3
  td->zeta=0.0625*(-3.0-5.0*a+b-a*b);
  td->eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
  ret+=lit_sG_zeta_eta(td);

  return ret*wt;
}

double complex lit_sH_alpha(double alpha,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta_t=alpha;
  res= deintz(lit_sH_beta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("lit_coef.c, DE integration error, lit_sH_alpha(). err=%g. Exit...\n",err);  exit(1);  }
  return res;
}

double complex lit_sH_beta(double b,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double a,wt;

  a=td->zeta_t;
  wt=1.0+0.25*(a+b);

  ret=0.0;
  // point 1
  td->zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
  td->eta =0.125*SQ3*(-a+b);
  ret+=lit_sH_zeta_eta(td);
  // part 2
  td->zeta=0.0625*(-3.0+a-5.0*b-a*b);
  td->eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
  ret+=lit_sH_zeta_eta(td);
  // part 3
  td->zeta=0.0625*(-3.0-5.0*a+b-a*b);
  td->eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
  ret+=lit_sH_zeta_eta(td);

  return ret*wt;
}

double lit_sF_alpha(double alpha,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double res,err;

  td->zeta_t=alpha;
  res= deintd(lit_sF_beta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("lit_coef.c, DE integration error, lit_sF_alpha(). err=%g. Exit...\n",err);  exit(1);  }
  return res;
}

double lit_sF_beta(double b,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double ret;
  double a,wt;

  a=td->zeta_t;
  wt=1.0+0.25*(a+b);

  ret=0.0;
  // point 1
  td->zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
  td->eta =0.125*SQ3*(-a+b);
  ret+=lit_sF_zeta_eta(td);
  // part 2
  td->zeta=0.0625*(-3.0+a-5.0*b-a*b);
  td->eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
  ret+=lit_sF_zeta_eta(td);
  // part 3
  td->zeta=0.0625*(-3.0-5.0*a+b-a*b);
  td->eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
  ret+=lit_sF_zeta_eta(td);

  return ret*wt;
}

double complex lit_sG_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("integral region error lit_coef.c, lit_sG_theta()! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  ret=0.0;
  for(i=0;i<td->gn;i++) ret+=lit_sG_r(0.5*rm*(td->xi[i]+1.0),td)*td->wi[i];
  ret*=0.5*rm;

  return ret;
}

double complex lit_sG_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex ce;
  double vK[3],aK,zeta,eta;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++) vK[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth;
  aK=vabs_d(vK);

  if(cimag(td->k)==0.0) ce=cos(creal(td->k)*aK*r)+I*sin(creal(td->k)*aK*r);
  else ce=cexp(I*td->k*aK*r);

  return ce/aK*lit_Nn(td->nid,zeta,eta);
}





