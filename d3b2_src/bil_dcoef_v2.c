/*
 * bil_dcoef_v2.c
 *
 *  Created on: Jan 24, 2019
 *      Author: ohta
 */

#include "d3b1_elem.h"


void bil_dcoef_v2_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,cb;
  double vR[3],aR,i_aR,aW,tD,*vW,rdW,rdv0,Wdv0,rdv1,Wdv1;
  int sig,i,l;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  CC[8]=0.0;
  dC0[8]=0.0;
  dC1[8]=0.0;
  for(i=0;i<4;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];
    vW=bd->wen[s][i];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    aW=vabs_d(vW);
    ca=I*k*aR;
    cb=ca-1.0;
    if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
    else ce=cexp(ca);
    rdW=vdot_d(vR,vW)*i_aR;
    tD=rdW*i_aR*i_aR;
    rdv0=vdot_d(vR,vt0)*i_aR;
    Wdv0=vdot_d(vW,vt0);
    rdv1=vdot_d(vR,vt1)*i_aR;
    Wdv1=vdot_d(vW,vt1);

    CC[i  ]=I_FP*ce*i_aR*aW; // G
    CC[i+4]=I_FP*(ca-1.0)*ce*tD; // H

    dC0[i  ]=-I_FP*cb*ce*i_aR*i_aR*rdv0*aW; // dG/dv
    dC1[i  ]=-I_FP*cb*ce*i_aR*i_aR*rdv1*aW; // dG/dv
    dC0[i+4]= I_FP*(ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv0*rdW-cb*Wdv0)); // dH/dv
    dC1[i+4]= I_FP*(ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv1*rdW-cb*Wdv1)); // dH/dv
    CC[8]+=bd->wt_44[i]*tD; // F
    dC0[8]+=bd->wt_44[i]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0); // dF/dv
    dC1[8]+=bd->wt_44[i]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1); // dF/dv
  }

  CC[8]*=I_FP;
  dC0[8]*=I_FP;
  dC1[8]*=I_FP;

  if(sig==1){
    bil_convert_CC(CC);
    bil_convert_dCC(dC0);
    bil_convert_dCC(dC1);
  }

}

void bil_dcoef_v2_9p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd)
{
  double complex ca,cb,ce;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,rdW,rdv0,Wdv0,rdv1,Wdv1;
  int sig,i,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }
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
      cb=ca-1.0;
      if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
      else ce=cexp(ca);
      rdW=vdot_d(vR,vW)*i_aR;
      tD=rdW*i_aR*i_aR;
      rdv0=vdot_d(vR,vt0)*i_aR;
      Wdv0=vdot_d(vW,vt0);
      rdv1=vdot_d(vR,vt1)*i_aR;
      Wdv1=vdot_d(vW,vt1);
      Np=bil_Nn(ep,zeta,eta);

      CC[ep  ]+=bd->wt_49[i]*ce*i_aR*aW*Np;
      CC[ep+4]+=bd->wt_49[i]*(ca-1.0)*ce*tD*Np;

      dC0[ep  ]+=bd->wt_49[i]*cb*ce*i_aR*i_aR*rdv0*aW*Np;
      dC1[ep  ]+=bd->wt_49[i]*cb*ce*i_aR*i_aR*rdv1*aW*Np;
      dC0[ep+4]+=bd->wt_49[i]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv0*rdW-cb*Wdv0)*Np;
      dC1[ep+4]+=bd->wt_49[i]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv1*rdW-cb*Wdv1)*Np;
      if(ep==0){
        CC [8]+=bd->wt_49[i]*tD;
        dC0[8]+=bd->wt_49[i]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0);
        dC1[8]+=bd->wt_49[i]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1);
      }
    }
    CC [ep  ]*= I_FP;
    CC [ep+4]*= I_FP;
    dC0[ep  ]*=-I_FP;
    dC0[ep+4]*= I_FP;
    dC1[ep  ]*=-I_FP;
    dC1[ep+4]*= I_FP;
    if(ep==0){
      CC [8]*=I_FP;
      dC0[8]*=I_FP;
      dC1[8]*=I_FP;
    }
  }

  if(sig==1){
    bil_convert_CC(CC);
    bil_convert_dCC(dC0);
    bil_convert_dCC(dC1);
  }
}

void bil_dcoef_v2_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,cb,tmpdg0,tmpdh0,tmpg,tmph,tmpdg1,tmpdh1;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,tmpf,tmpdf0,rdW,rdv0,Wdv0,tmpdf1,rdv1,Wdv1;
  int sig,i,j,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }
  for(ep=0;ep<4;ep++){
    for(i=0;i<GLN;i++){
      tmpg=0.0;      tmph=0.0;
      tmpdg0=0.0;      tmpdh0=0.0;
      tmpdg1=0.0;      tmpdh1=0.0;
      tmpf =0.0;
      tmpdf0=0.0;
      tmpdf1=0.0;
      for(j=0;j<GLN;j++){
        zeta=bd->xli[i];
        eta =bd->xli[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        ca=I*k*aR;
        cb=ca-1.0;
        if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
        else ce=cexp(ca);
        rdW=vdot_d(vR,vW)*i_aR;
        tD=rdW*i_aR*i_aR;
        rdv0=vdot_d(vR,vt0)*i_aR;
        rdv1=vdot_d(vR,vt1)*i_aR;
        Wdv0=vdot_d(vW,vt0);
        Wdv1=vdot_d(vW,vt1);
        Np=bil_Nn(ep,zeta,eta);

        tmpg+=bd->wli[j]*ce*i_aR*aW*Np;
        tmph+=bd->wli[j]*(ca-1.0)*ce*tD*Np;

        tmpdg0+=bd->wli[j]*cb*ce*i_aR*i_aR*rdv0*aW*Np;
        tmpdg1+=bd->wli[j]*cb*ce*i_aR*i_aR*rdv1*aW*Np;
        tmpdh0+=bd->wli[j]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv0*rdW-cb*Wdv0)*Np;
        tmpdh1+=bd->wli[j]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv1*rdW-cb*Wdv1)*Np;
        if(ep==0){
          tmpf +=bd->wli[j]*tD;
          tmpdf0+=bd->wli[j]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0);
          tmpdf1+=bd->wli[j]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1);
        }
      }
      CC [ep  ]+=bd->wli[i]*tmpg;
      CC [ep+4]+=bd->wli[i]*tmph;
      dC0[ep  ]+=bd->wli[i]*tmpdg0;
      dC1[ep  ]+=bd->wli[i]*tmpdg1;
      dC0[ep+4]+=bd->wli[i]*tmpdh0;
      dC1[ep+4]+=bd->wli[i]*tmpdh1;
      if(ep==0){
        CC [8]+=bd->wli[i]*tmpf;
        dC0[8]+=bd->wli[i]*tmpdf0;
        dC1[8]+=bd->wli[i]*tmpdf1;
      }
    }
    CC [ep  ]*= I_FP;
    CC [ep+4]*= I_FP;
    dC0[ep  ]*=-I_FP;
    dC1[ep  ]*=-I_FP;
    dC0[ep+4]*= I_FP;
    dC1[ep+4]*= I_FP;
    if(ep==0){
      CC [8]*=I_FP;
      dC0[8]*=I_FP;
      dC1[8]*=I_FP;
    }
  }

  if(sig==1){
    bil_convert_CC(CC);
    bil_convert_dCC(dC0);
    bil_convert_dCC(dC1);
  }
}

void bil_dcoef_v2_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,cb,tmpdg0,tmpdh0,tmpg,tmph,tmpdg1,tmpdh1;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,tmpf,tmpdf0,rdW,rdv0,Wdv0,tmpdf1,rdv1,Wdv1;
  int sig,i,j,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }
  for(ep=0;ep<4;ep++){
    for(i=0;i<GHN;i++){
      tmpg=0.0;      tmph=0.0;
      tmpdg0=0.0;      tmpdh0=0.0;
      tmpdg1=0.0;      tmpdh1=0.0;
      tmpf =0.0;
      tmpdf0=0.0;
      tmpdf1=0.0;
      for(j=0;j<GHN;j++){
        zeta=bd->xhi[i];
        eta =bd->xhi[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        ca=I*k*aR;
        cb=ca-1.0;
        if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
        else ce=cexp(ca);
        rdW=vdot_d(vR,vW)*i_aR;
        tD=rdW*i_aR*i_aR;
        rdv0=vdot_d(vR,vt0)*i_aR;
        rdv1=vdot_d(vR,vt1)*i_aR;
        Wdv0=vdot_d(vW,vt0);
        Wdv1=vdot_d(vW,vt1);
        Np=bil_Nn(ep,zeta,eta);

        tmpg+=bd->whi[j]*ce*i_aR*aW*Np;
        tmph+=bd->whi[j]*(ca-1.0)*ce*tD*Np;

        tmpdg0+=bd->whi[j]*cb*ce*i_aR*i_aR*rdv0*aW*Np;
        tmpdg1+=bd->whi[j]*cb*ce*i_aR*i_aR*rdv1*aW*Np;
        tmpdh0+=bd->whi[j]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv0*rdW-cb*Wdv0)*Np;
        tmpdh1+=bd->whi[j]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv1*rdW-cb*Wdv1)*Np;
        if(ep==0){
          tmpf +=bd->whi[j]*tD;
          tmpdf0+=bd->whi[j]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0);
          tmpdf1+=bd->whi[j]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1);
        }
      }
      CC [ep  ]+=bd->whi[i]*tmpg;
      CC [ep+4]+=bd->whi[i]*tmph;
      dC0[ep  ]+=bd->whi[i]*tmpdg0;
      dC1[ep  ]+=bd->whi[i]*tmpdg1;
      dC0[ep+4]+=bd->whi[i]*tmpdh0;
      dC1[ep+4]+=bd->whi[i]*tmpdh1;
      if(ep==0){
        CC [8]+=bd->whi[i]*tmpf;
        dC0[8]+=bd->whi[i]*tmpdf0;
        dC1[8]+=bd->whi[i]*tmpdf1;
      }
    }
    CC [ep  ]*= I_FP;
    CC [ep+4]*= I_FP;
    dC0[ep  ]*=-I_FP;
    dC1[ep  ]*=-I_FP;
    dC0[ep+4]*= I_FP;
    dC1[ep+4]*= I_FP;
    if(ep==0){
      CC [8]*=I_FP;
      dC0[8]*=I_FP;
      dC1[8]*=I_FP;
    }
  }

  if(sig==1){
    bil_convert_CC(CC);
    bil_convert_dCC(dC0);
    bil_convert_dCC(dC1);
  }
}

void bil_dcoef_v2_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd)
{
  bil_coef_DE(CC,rt,s,k,bd);
  bil_dcoef_DE0(dC0,rt,vt0,s,k,bd);
  bil_dcoef_DE0(dC1,rt,vt1,s,k,bd);
}

void bil_dcoef_bd_vt2_PV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  bil_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  bil_dcoef_bd_tz0_PV(dCdtz,zeta_t,eta_t,s,k,bd);
  bil_dcoef_bd_te0_PV(dCdte,zeta_t,eta_t,s,k,bd);
}

void bil_dcoef_bd_vt2_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  bil_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  bil_dcoef_bd_tz0_DV(dCdtz,zeta_t,eta_t,s,k,bd);
  bil_dcoef_bd_te0_DV(dCdte,zeta_t,eta_t,s,k,bd);
}


