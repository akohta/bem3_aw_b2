/*
 * lit_dcoef_v2.c
 *
 *  Created on: Jan 24, 2019
 *      Author: ohta
 */

#include "d3b1_elem.h"

void lit_sdGL_v2_zeta_eta(double complex *ret,TDATA *td);


void lit_dcoef_v2_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd)
{
  double complex ca,cb,ce,kR;
  double vR[3],*vW,aR,i_aR,aW,rdW,rdv0,Wdv0,Nn,rdv1,Wdv1;
  int sig,i,l;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  for(i=0;i<9;i++){
    CC [i]=0.0;
    dC0[i]=0.0;
    dC1[i]=0.0;
  }
  vW=bd->wen[s][0];
  aW=vabs_d(vW);

  for(i=0;i<3;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    kR=k*aR;
    ca=I*kR;
    cb=ca-1.0;
    if(cimag(k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
    else ce=cexp(ca);
    rdW=vdot_d(vR,vW)*i_aR;
    rdv0=vdot_d(vR,vt0)*i_aR;
    rdv1=vdot_d(vR,vt1)*i_aR;
    Wdv0=vdot_d(vW,vt0);
    Wdv1=vdot_d(vW,vt1);

    CC[i+0]=bd->wt_34[i]*ce*i_aR;
    CC[i+4]=bd->wt_34[i]*(ca-1.0)*ce*rdW*i_aR*i_aR;

    dC0[i+0]=bd->wt_34[i]*cb*ce*i_aR*i_aR*rdv0; // dG/dv
    dC1[i+0]=bd->wt_34[i]*cb*ce*i_aR*i_aR*rdv1; // dG/dv
    dC0[i+4]=bd->wt_34[i]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv0*rdW-cb*Wdv0 ); // dH/dv
    dC1[i+4]=bd->wt_34[i]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv1*rdW-cb*Wdv1 ); // dH/dv
    CC [8]+=bd->wt_34[i]*rdW*i_aR*i_aR; // F
    dC0[8]+=bd->wt_34[i]*i_aR*i_aR*i_aR*( 3.0*rdv0*rdW-Wdv0 ); // dF/dv
    dC1[8]+=bd->wt_34[i]*i_aR*i_aR*i_aR*( 3.0*rdv1*rdW-Wdv1 ); // dF/dv
  }

  for(l=0;l<3;l++) vR[l]=bd->ren[s][3][l]-rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  kR=k*aR;
  ca=I*kR;
  cb=ca-1.0;
  if(cimag(k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
  else ce=cexp(ca);
  rdW=vdot_d(vR,vW)*i_aR;
  rdv0=vdot_d(vR,vt0)*i_aR;
  rdv1=vdot_d(vR,vt1)*i_aR;
  Wdv0=vdot_d(vW,vt0);
  Wdv1=vdot_d(vW,vt1);

  for(i=0;i<3;i++){
    Nn=lit_Nn(i,bd->zt_34[3],bd->et_34[3]);
    CC[i+0]+=bd->wt_34[3]*ce*i_aR*Nn;
    CC[i+4]+=bd->wt_34[3]*(ca-1.0)*ce*rdW*i_aR*i_aR*Nn;
    dC0[i+0]+=bd->wt_34[3]*cb*ce*i_aR*i_aR*rdv0*Nn;
    dC1[i+0]+=bd->wt_34[3]*cb*ce*i_aR*i_aR*rdv1*Nn;
    dC0[i+4]+=bd->wt_34[3]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv0*rdW-cb*Wdv0 )*Nn;
    dC1[i+4]+=bd->wt_34[3]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv1*rdW-cb*Wdv1 )*Nn;
    if(i==0){
      CC [8]+=bd->wt_34[3]*rdW*i_aR*i_aR;
      dC0[8]+=bd->wt_34[3]*i_aR*i_aR*i_aR*( 3.0*rdv0*rdW-Wdv0 );
      dC1[8]+=bd->wt_34[3]*i_aR*i_aR*i_aR*( 3.0*rdv1*rdW-Wdv1 );
    }
  }

  for(i=0;i<3;i++){
    CC[i+0]*=I_EP*aW;
    CC[i+4]*=I_EP;
    dC0[i+0]*=-I_EP*aW;
    dC1[i+0]*=-I_EP*aW;
    dC0[i+4]*= I_EP;
    dC1[i+4]*= I_EP;
    if(i==0){
      CC [8]*=I_EP;
      dC0[8]*=I_EP;
      dC1[8]*=I_EP;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    lit_convert_dCC(dC0);
    lit_convert_dCC(dC1);
  }
}

void lit_dcoef_v2_7p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd)
{
  double complex ca,cb,ce,kR;
  double cr[3][4],cw[3][3],zeta,eta,vW[3],aW,vR[3],aR,i_aR,rdW,rdv0,Wdv0,Np,rdv1,Wdv1;
  int sig,i,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  lit_w_zeta_eta(vW,0.0,0.0,cw);
  aW=vabs_d(vW);

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }
  for(ep=0;ep<3;ep++){
    for(i=0;i<7;i++){
      zeta=bd->zt_37[i];
      eta =bd->et_37[i];
      lit_r_zeta_eta(vR,zeta,eta,cr);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      aR=vabs_d(vR);
      i_aR=1.0/aR;
      kR=k*aR;
      ca=I*kR;
      cb=ca-1.0;
      if(cimag(k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
      else ce=cexp(ca);
      rdW=vdot_d(vR,vW)*i_aR;
      rdv0=vdot_d(vR,vt0)*i_aR;
      rdv1=vdot_d(vR,vt1)*i_aR;
      Wdv0=vdot_d(vW,vt0);
      Wdv1=vdot_d(vW,vt1);
      Np=lit_Nn(ep,zeta,eta);

      CC[ep  ]+=bd->wt_37[i]*ce*i_aR*Np;
      CC[ep+4]+=bd->wt_37[i]*(ca-1.0)*ce*rdW*i_aR*i_aR*Np;
      dC0[ep  ]+=bd->wt_37[i]*cb*ce*i_aR*i_aR*rdv0*Np;
      dC1[ep  ]+=bd->wt_37[i]*cb*ce*i_aR*i_aR*rdv1*Np;
      dC0[ep+4]+=bd->wt_37[i]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv0*rdW-cb*Wdv0 )*Np;
      dC1[ep+4]+=bd->wt_37[i]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv1*rdW-cb*Wdv1 )*Np;
      if(ep==0){
        CC [8]+=bd->wt_37[i]*rdW*i_aR*i_aR;
        dC0[8]+=bd->wt_37[i]*i_aR*i_aR*i_aR*( 3.0*rdv0*rdW-Wdv0 );
        dC1[8]+=bd->wt_37[i]*i_aR*i_aR*i_aR*( 3.0*rdv1*rdW-Wdv1 );
      }
    }
    CC [ep  ]*= I_EP*aW;
    CC [ep+4]*= I_EP;
    dC0[ep  ]*=-I_EP*aW;
    dC1[ep  ]*=-I_EP*aW;
    dC0[ep+4]*= I_EP;
    dC1[ep+4]*= I_EP;
    if(ep==0){
      CC [8]*=I_EP;
      dC0[8]*=I_EP;
      dC1[8]*=I_EP;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    lit_convert_dCC(dC0);
    lit_convert_dCC(dC1);
  }
}

void lit_dcoef_v2_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex ret[9],tpj[9],tpi[9];
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
  for(l=0;l<3;l++){
    td.rt[l]=rt[l];
    td.vt[l]=vt0[l];
    td.vu[l]=vt1[l];
  }

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GLN;i++){
      for(l=0;l<9;l++) tpi[l]=0.0;
      a=bd->xli[i];
      for(j=0;j<GLN;j++){
        for(l=0;l<9;l++) tpj[l]=0.0;
        b=bd->xli[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        lit_sdGL_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        lit_sdGL_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        lit_sdGL_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];

        for(l=0;l<9;l++) tpi[l]+=tpj[l]*wt*bd->wli[j];
      }
      CC [ep  ]+=tpi[0]*bd->wli[i];
      CC [ep+4]+=tpi[1]*bd->wli[i];
      dC0[ep  ]+=tpi[3]*bd->wli[i];
      dC0[ep+4]+=tpi[4]*bd->wli[i];
      dC1[ep  ]+=tpi[6]*bd->wli[i];
      dC1[ep+4]+=tpi[7]*bd->wli[i];
      if(ep==0){
        CC [8]+=tpi[2]*bd->wli[i];
        dC0[8]+=tpi[5]*bd->wli[i];
        dC1[8]+=tpi[8]*bd->wli[i];
      }
    }
    CC[ep  ]*=0.0625*SQ3*I_SPSQ3*td.aW;
    CC[ep+4]*=0.0625*SQ3*I_SPSQ3;
    dC0[ep  ]*=-0.0625*SQ3*I_SPSQ3*td.aW;
    dC0[ep+4]*= 0.0625*SQ3*I_SPSQ3;
    dC1[ep  ]*=-0.0625*SQ3*I_SPSQ3*td.aW;
    dC1[ep+4]*= 0.0625*SQ3*I_SPSQ3;
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      dC0[8]*=0.0625*SQ3*I_SPSQ3;
      dC1[8]*=0.0625*SQ3*I_SPSQ3;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    lit_convert_dCC(dC0);
    lit_convert_dCC(dC1);
  }
}

void lit_dcoef_v2_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex ret[9],tpj[9],tpi[9];
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
  for(l=0;l<3;l++){
    td.rt[l]=rt[l];
    td.vt[l]=vt0[l];
    td.vu[l]=vt1[l];
  }

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GHN;i++){
      for(l=0;l<9;l++) tpi[l]=0.0;
      a=bd->xhi[i];
      for(j=0;j<GHN;j++){
        for(l=0;l<9;l++) tpj[l]=0.0;
        b=bd->xhi[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        lit_sdGL_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        lit_sdGL_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        lit_sdGL_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];

        for(l=0;l<9;l++) tpi[l]+=tpj[l]*wt*bd->whi[j];
      }
      CC [ep  ]+=tpi[0]*bd->whi[i];
      CC [ep+4]+=tpi[1]*bd->whi[i];
      dC0[ep  ]+=tpi[3]*bd->whi[i];
      dC0[ep+4]+=tpi[4]*bd->whi[i];
      dC1[ep  ]+=tpi[6]*bd->whi[i];
      dC1[ep+4]+=tpi[7]*bd->whi[i];
      if(ep==0){
        CC [8]+=tpi[2]*bd->whi[i];
        dC0[8]+=tpi[5]*bd->whi[i];
        dC1[8]+=tpi[8]*bd->whi[i];
      }
    }
    CC[ep  ]*=0.0625*SQ3*I_SPSQ3*td.aW;
    CC[ep+4]*=0.0625*SQ3*I_SPSQ3;
    dC0[ep  ]*=-0.0625*SQ3*I_SPSQ3*td.aW;
    dC0[ep+4]*= 0.0625*SQ3*I_SPSQ3;
    dC1[ep  ]*=-0.0625*SQ3*I_SPSQ3*td.aW;
    dC1[ep+4]*= 0.0625*SQ3*I_SPSQ3;
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      dC0[8]*=0.0625*SQ3*I_SPSQ3;
      dC1[8]*=0.0625*SQ3*I_SPSQ3;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    lit_convert_dCC(dC0);
    lit_convert_dCC(dC1);
  }
}

void lit_dcoef_v2_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd)
{
  lit_coef_DE(CC,rt,s,k,bd);
  lit_dcoef_DE0(dC0,rt,vt0,s,k,bd);
  lit_dcoef_DE0(dC1,rt,vt1,s,k,bd);
}

void lit_dcoef_bd_vt2_PV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  lit_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  lit_dcoef_bd_tz0_PV(dCdtz,zeta_t,eta_t,s,k,bd);
  lit_dcoef_bd_te0_PV(dCdte,zeta_t,eta_t,s,k,bd);
}

void lit_dcoef_bd_vt2_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  lit_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  lit_dcoef_bd_tz0_DV(dCdtz,zeta_t,eta_t,s,k,bd);
  lit_dcoef_bd_te0_DV(dCdte,zeta_t,eta_t,s,k,bd);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
void lit_sdGL_v2_zeta_eta(double complex *ret,TDATA *td)
{
  double complex ca,cb,ce,kR;
  double vR[3],aR,i_aR,rdW,rdv0,Wdv0,rdv1,Wdv1,Np;
  int l;

  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);

  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  kR=td->k*aR;
  ca=I*kR;
  cb=ca-1.0;
  if(cimag(td->k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
  else ce=cexp(ca);
  rdW=(vR[0]*td->cw[0][0]+vR[1]*td->cw[1][0]+vR[2]*td->cw[2][0])*i_aR;
  rdv0=vdot_d(vR,td->vt)*i_aR;
  rdv1=vdot_d(vR,td->vu)*i_aR;
  Wdv0=td->cw[0][0]*td->vt[0]+td->cw[1][0]*td->vt[1]+td->cw[2][0]*td->vt[2];
  Wdv1=td->cw[0][0]*td->vu[0]+td->cw[1][0]*td->vu[1]+td->cw[2][0]*td->vu[2];
  Np=lit_Nn(td->nid,td->zeta,td->eta);

  ret[0]=ce*i_aR*Np; // G
  ret[1]=(ca-1.0)*ce*rdW*i_aR*i_aR*Np; // H
  ret[2]=rdW*i_aR*i_aR; // F
  ret[3]=cb*ce*i_aR*i_aR*rdv0*Np; // dG/dv
  ret[4]=ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv0*rdW-cb*Wdv0 )*Np; // dH/dv
  ret[5]=i_aR*i_aR*i_aR*( 3.0*rdv0*rdW-Wdv0 ); // dF/dv
  ret[6]=cb*ce*i_aR*i_aR*rdv1*Np; // dG/du
  ret[7]=ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv1*rdW-cb*Wdv1 )*Np; // dH/du
  ret[8]=i_aR*i_aR*i_aR*( 3.0*rdv1*rdW-Wdv1 ); // dF/du
}


