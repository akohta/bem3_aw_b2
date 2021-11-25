/*
 * bil_dcoef.c
 *
 *  Created on: Dec 25, 2018
 *      Author: ohta
 */

#include "d3b1_elem.h"

double complex bil_sdG_zeta(double zeta,void *tmp);
double complex bil_sdG_eta(double eta,void *tmp);
double complex bil_sdH_zeta(double zeta,void *tmp);
double complex bil_sdH_eta(double eta,void *tmp);
double bil_sdF_zeta(double zeta,void *tmp);
double bil_sdF_eta(double eta,void *tmp);

double complex bil_sdG_Db_theta(double theta,void *tmp);
double complex bil_sdG_Db_r(double r,void *tmp);
double complex bil_sdG_De_theta(double theta,void *tmp);
double complex bil_sdG_De_r(double r,void *tmp);
double complex bil_sdG_DP_theta(double theta,void *tmp);
double complex bil_sdG_DP_r(double r,void *tmp);

double complex bil_sdHz_Db_theta(double theta,void *tmp);
double complex bil_sdHz_Db_r(double r,void *tmp);
double complex bil_sdHz_De_theta(double theta,void *tmp);
double complex bil_sdHz_De0_r(double r,void *tmp);
double complex bil_sdHz_De1_r(double r,void *tmp);
double complex bil_sdHz_DP_theta(double theta,void *tmp);
double complex bil_sdHz_DP_r(double r,void *tmp);
double bil_sdFz_Db_theta(double theta,void *tmp);
double bil_sdFz_Db_r(double r,void *tmp);
double bil_sdFz_DP_theta(double theta,void *tmp);
double bil_sdFz_DP_r(double r,void *tmp);

double complex bil_sdHe_Db_theta(double theta,void *tmp);
double complex bil_sdHe_Db_r(double r,void *tmp);
double complex bil_sdHe_De_theta(double theta,void *tmp);
double complex bil_sdHe_De0_r(double r,void *tmp);
double complex bil_sdHe_De1_r(double r,void *tmp);
double complex bil_sdHe_DP_theta(double theta,void *tmp);
double complex bil_sdHe_DP_r(double r,void *tmp);
double bil_sdFe_Db_theta(double theta,void *tmp);
double bil_sdFe_Db_r(double r,void *tmp);
double bil_sdFe_DP_theta(double theta,void *tmp);
double bil_sdFe_DP_r(double r,void *tmp);



void bil_dcoef_4p(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,cb;
  double vR[3],aR,i_aR,aW,tD,*vW,rdW,rdv,Wdv;
  int sig,i,l;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  CC[8]=0.0;
  dCC[8]=0.0;
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
    rdv=vdot_d(vR,vt)*i_aR;
    Wdv=vdot_d(vW,vt);

    CC[i  ]=I_FP*ce*i_aR*aW; // G
    CC[i+4]=I_FP*(ca-1.0)*ce*tD; // H

    dCC[i  ]=-I_FP*cb*ce*i_aR*i_aR*rdv*aW; // dG/dv
    dCC[i+4]= I_FP*(ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv*rdW-cb*Wdv)); // dH/dv
    CC[8]+=bd->wt_44[i]*tD; // F
    dCC[8]+=bd->wt_44[i]*i_aR*i_aR*i_aR*(3.0*rdv*rdW-Wdv); // dF/dv
  }

  CC[8]*=I_FP;
  dCC[8]*=I_FP;

  if(sig==1){
    bil_convert_CC(CC);
    bil_convert_dCC(dCC);
  }
}

void bil_dcoef_9p(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
{
  double complex ca,cb,ce;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,rdW,rdv,Wdv;
  int sig,i,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dCC[l]=0.0;
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
      rdv=vdot_d(vR,vt)*i_aR;
      Wdv=vdot_d(vW,vt);
      Np=bil_Nn(ep,zeta,eta);

      CC[ep  ]+=bd->wt_49[i]*ce*i_aR*aW*Np;
      CC[ep+4]+=bd->wt_49[i]*(ca-1.0)*ce*tD*Np;

      dCC[ep  ]+=bd->wt_49[i]*cb*ce*i_aR*i_aR*rdv*aW*Np;
      dCC[ep+4]+=bd->wt_49[i]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv*rdW-cb*Wdv)*Np;
      if(ep==0){
        CC [8]+=bd->wt_49[i]*tD;
        dCC[8]+=bd->wt_49[i]*i_aR*i_aR*i_aR*(3.0*rdv*rdW-Wdv);
      }
    }
    CC [ep  ]*= I_FP;
    CC [ep+4]*= I_FP;
    dCC[ep  ]*=-I_FP;
    dCC[ep+4]*= I_FP;
    if(ep==0){
      CC [8]*=I_FP;
      dCC[8]*=I_FP;
    }
  }

  if(sig==1){
    bil_convert_CC(CC);
    bil_convert_dCC(dCC);
  }
}

void bil_dcoef_GL(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,cb,tmpdg,tmpdh,tmpg,tmph;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,tmpf,tmpdf,rdW,rdv,Wdv;
  int sig,i,j,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dCC[l]=0.0;
  }
  for(ep=0;ep<4;ep++){
    for(i=0;i<GLN;i++){
      tmpg=0.0;
      tmph=0.0;
      tmpdg=0.0;
      tmpdh=0.0;
      tmpf =0.0;
      tmpdf=0.0;
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
        rdv=vdot_d(vR,vt)*i_aR;
        Wdv=vdot_d(vW,vt);
        Np=bil_Nn(ep,zeta,eta);

        tmpg+=bd->wli[j]*ce*i_aR*aW*Np;
        tmph+=bd->wli[j]*(ca-1.0)*ce*tD*Np;

        tmpdg+=bd->wli[j]*cb*ce*i_aR*i_aR*rdv*aW*Np;
        tmpdh+=bd->wli[j]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv*rdW-cb*Wdv)*Np;
        if(ep==0){
          tmpf +=bd->wli[j]*tD;
          tmpdf+=bd->wli[j]*i_aR*i_aR*i_aR*(3.0*rdv*rdW-Wdv);
        }
      }
      CC [ep  ]+=bd->wli[i]*tmpg;
      CC [ep+4]+=bd->wli[i]*tmph;
      dCC[ep  ]+=bd->wli[i]*tmpdg;
      dCC[ep+4]+=bd->wli[i]*tmpdh;
      if(ep==0){
        CC [8]+=bd->wli[i]*tmpf;
        dCC[8]+=bd->wli[i]*tmpdf;
      }
    }
    CC [ep  ]*= I_FP;
    CC [ep+4]*= I_FP;
    dCC[ep  ]*=-I_FP;
    dCC[ep+4]*= I_FP;
    if(ep==0){
      CC [8]*=I_FP;
      dCC[8]*=I_FP;
    }
  }

  if(sig==1){
    bil_convert_CC(CC);
    bil_convert_dCC(dCC);
  }
}

void bil_dcoef_GH(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,cb,tmpdg,tmpdh,tmpg,tmph;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,tmpf,tmpdf,rdW,rdv,Wdv;
  int sig,i,j,l,ep;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dCC[l]=0.0;
  }
  for(ep=0;ep<4;ep++){
    for(i=0;i<GHN;i++){
      tmpg=0.0;
      tmph=0.0;
      tmpdg=0.0;
      tmpdh=0.0;
      tmpf =0.0;
      tmpdf=0.0;
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
        rdv=vdot_d(vR,vt)*i_aR;
        Wdv=vdot_d(vW,vt);
        Np=bil_Nn(ep,zeta,eta);

        tmpg+=bd->whi[j]*ce*i_aR*aW*Np;
        tmph+=bd->whi[j]*(ca-1.0)*ce*tD*Np;

        tmpdg+=bd->whi[j]*cb*ce*i_aR*i_aR*rdv*aW*Np;
        tmpdh+=bd->whi[j]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*rdv*rdW-cb*Wdv)*Np;
        if(ep==0){
          tmpf +=bd->whi[j]*tD;
          tmpdf+=bd->whi[j]*i_aR*i_aR*i_aR*(3.0*rdv*rdW-Wdv);
        }
      }
      CC [ep  ]+=bd->whi[i]*tmpg;
      CC [ep+4]+=bd->whi[i]*tmph;
      dCC[ep  ]+=bd->whi[i]*tmpdg;
      dCC[ep+4]+=bd->whi[i]*tmpdh;
      if(ep==0){
        CC [8]+=bd->whi[i]*tmpf;
        dCC[8]+=bd->whi[i]*tmpdf;
      }
    }
    CC [ep  ]*= I_FP;
    CC [ep+4]*= I_FP;
    dCC[ep  ]*=-I_FP;
    dCC[ep+4]*= I_FP;
    if(ep==0){
      CC [8]*=I_FP;
      dCC[8]*=I_FP;
    }
  }

  if(sig==1){
    bil_convert_CC(CC);
    bil_convert_dCC(dCC);
  }
}

void bil_dcoef_DE(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
{
  bil_coef_DE(CC,rt,s,k,bd);
  bil_dcoef_DE0(dCC,rt,vt,s,k,bd);
}

void bil_dcoef_DE0(double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
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
  for(l=0;l<3;l++){
    td.rt[l]=rt[l];
    td.vt[l]=vt[l];
  }

  for(i=0;i<9;i++) dCC[i]=0.0; //
  for(i=0;i<4;i++){ // dG/dv,dH/dv
    td.nid=i;
    dCC[i]=-I_FP*deintz(bil_sdG_zeta,-1.0,1.0,&td,IEPS,&err);
    if(err<0){ printf("bil_dcoef.c, DE integration error, bil_dcoef_DE(), bil_sdG_zeta(), nid=%d. Exit...\n",td.nid);  exit(1);  }
    dCC[i+4]=I_FP*deintz(bil_sdH_zeta,-1.0,1.0,&td,IEPS,&err);
    if(err<0){ printf("bil_dcoef.c, DE integration error, bil_dcoef_DE(), bil_sdH_zeta(), nid=%d. Exit...\n",td.nid);  exit(1);  }
  }

  // dF/dv
  dCC[8]=I_FP*deintd(bil_sdF_zeta,-1.0,1.0,&td,IEPS,&err);
  if(err<0){ printf("bil_dcoef.c, DE integration error, bil_dcoef_DE(), bil_sdF_zeta(). Exit...\n");  exit(1);  }

  if(sig==1) bil_convert_dCC(dCC);
}

void bil_dcoef_bd_tz_PV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  bil_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  bil_dcoef_bd_tz0_PV(dCC,zeta_t,eta_t,s,k,bd);
}

void bil_dcoef_bd_tz0_PV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex tC;
  double r[5],th[5],err,bdC,i_aT,tD;
  int sig,i,j,l,m;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  td.zeta_t=zeta_t;
  td.eta_t = eta_t;

  td.gn=GLN;
  td.xi=bd->xli;
  td.wi=bd->wli;

  td.ghn=GHN;
  td.xhi=bd->xhi;
  td.whi=bd->whi;

  bil_copy_elem_const_rw(td.cr,td.cw,s,bd);
  bil_calc_node_r_th(r,th,zeta_t,eta_t);

  bil_t_zeta(td.vt,td.zeta_t,td.eta_t,td.cr);
  i_aT=1.0/vabs_d(td.vt);

  for(i=0;i<9;i++) dCC[i]=0.0;

  // dG/dt_zeta
  for(i=0;i<4;i++){
    td.nid=i;
    dCC[i]=0.0;
    // D_bar(epsilon)
    for(j=0;j<4;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];

      if(SW_DBDT_BDE==0){
        tC=0.0;
        for(m=0;m<GHN;m++) tC+=( bil_sdG_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +bil_sdG_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tC*=0.25*(th[j+1]-th[j]);
        dCC[i]+=tC;
      }
      else {
        dCC[i]+=deintz(bil_sdG_Db_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdG_Db_theta()! err=%f Exit...\n",err);  exit(1);  }
      }
    }
    // D_epsilon
    if(SW_DBDT_DE==0){
      tC=0.0;
      for(m=0;m<GHN;m++) tC+=( bil_sdG_De_theta(0.25*M_PI*bd->xhi[m]-0.75*M_PI,&td)
                              +bil_sdG_De_theta(0.25*M_PI*bd->xhi[m]-0.25*M_PI,&td)
                              +bil_sdG_De_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                              +bil_sdG_De_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
      tC*=0.25*M_PI;
      dCC[i]+=I*tC;
    }
    else {
      dCC[i]+=I*deintz(bil_sdG_De_theta,0.0,2.0*M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdG_De_theta()! err=%f Exit...\n",err);  exit(1);  }
    }

    // D_epsilon PV
    if(SW_DBDT_DE==0){
      tC=0.0;
      for(m=0;m<GHN;m++) tC+=( bil_sdG_DP_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                              +bil_sdG_DP_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
      tC*=0.25*M_PI;
      dCC[i]-=tC;
    }
    else {
      dCC[i]-=deintz(bil_sdG_DP_theta,0.0,M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdG_DP_theta()! err=%f Exit...\n",err);  exit(1);  }
    }
    dCC[i]*=-I_FP*i_aT;
  }

  if(bil_check_plane(td.cr,td.cw)==1){ // element is plane
    for(l=0;l<5;l++) dCC[l+4]=0.0;
  }
  else {
    bdC=td.cr[0][1]*td.cw[0][2]+td.cr[1][1]*td.cw[1][2]+td.cr[2][1]*td.cw[2][2];
    // dH/dt_zeta
    for(i=0;i<4;i++){
      td.nid=i;
      dCC[i+4]=0.0;
      // D\bar(epsilon)
      for(j=0;j<4;j++){
        td.r0=r[j];      td.r1=r[j+1];
        td.th0=th[j];    td.th1=th[j+1];

        if(SW_DBDT_BDE==0){
          tC=0.0;
          for(m=0;m<GHN;m++) tC+=( bil_sdHz_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                  +bil_sdHz_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
          tC*=0.25*(th[j+1]-th[j]);
          dCC[i+4]+=tC;
        }
        else {
          dCC[i+4]+=deintz(bil_sdHz_Db_theta,th[j],th[j+1],&td,IEPS,&err);
          if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdHz_Db_theta()! err=%f Exit...\n",err);  exit(1);  }
        }
      }

      // D_epsilon
      if(SW_DBDT_DE==0){
        tC=0.0;
        for(m=0;m<GHN;m++) tC+=( bil_sdHz_De_theta(0.25*M_PI*bd->xhi[m]-0.75*M_PI,&td)
                                +bil_sdHz_De_theta(0.25*M_PI*bd->xhi[m]-0.25*M_PI,&td)
                                +bil_sdHz_De_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                                +bil_sdHz_De_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
        tC*=0.25*M_PI;
        dCC[i+4]+=tC;
      }
      else {
        dCC[i+4]+=deintz(bil_sdHz_De_theta,0.0,M_PI,&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdHz_De_theta()0! err=%f Exit...\n",err);  exit(1);  }
        dCC[i+4]+=deintz(bil_sdHz_De_theta,M_PI,2.0*M_PI,&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdHz_De_theta()1! err=%f Exit...\n",err);  exit(1);  }
      }

      // D_epsilon PV
      if(SW_DBDT_DE==0){
        tC=0.0;
        for(m=0;m<GHN;m++) tC+=( bil_sdHz_DP_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                                +bil_sdHz_DP_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
        tC*=0.25*M_PI;
        dCC[i+4]-=tC;
      }
      else {
        dCC[i+4]-=deintz(bil_sdHz_DP_theta,0.0,M_PI,&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdHz_DP_theta()! err=%f Exit...\n",err);  exit(1);  }
      }

      dCC[i+4]*=bdC*I_FP*i_aT;
    }
    // dF/dt_zeta
    dCC[8]=0.0;
    // D\bar(epsilon)
    for(j=0;j<4;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];

      if(SW_DBDT_BDE==0){
        tD=0.0;
        for(m=0;m<GHN;m++) tD+=( bil_sdFz_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +bil_sdFz_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tD*=0.25*(th[j+1]-th[j]);
        dCC[8]+=tD;
      }
      else {
        dCC[8]+=deintd(bil_sdFz_Db_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdFz_Db_theta()! err=%f Exit...\n",err);  exit(1);  }
      }
    }

    // D_epsilon PV
    if(SW_DBDT_DE==0){
      tD=0.0;
      for(m=0;m<GHN;m++) tD+=( bil_sdFz_DP_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                              +bil_sdFz_DP_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
      tD*=0.25*M_PI;
      dCC[8]+=tD;
    }
    else {
      dCC[8]+=deintd(bil_sdFz_DP_theta,0.0,M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdFz_DP_theta()! err=%f Exit...\n",err);  exit(1);  }
    }

    dCC[8]*=bdC*I_FP*i_aT;
  }

  if(sig==1) bil_convert_dCC(dCC);
}

void bil_dcoef_bd_te_PV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  bil_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  bil_dcoef_bd_te0_PV(dCC,zeta_t,eta_t,s,k,bd);
}

void bil_dcoef_bd_te0_PV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex tC;
  double r[5],th[5],err,bdC,i_aT,tD;
  int sig,i,j,l,m;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  td.zeta_t=zeta_t;
  td.eta_t = eta_t;

  td.gn=GLN;
  td.xi=bd->xli;
  td.wi=bd->wli;

  td.ghn=GHN;
  td.xhi=bd->xhi;
  td.whi=bd->whi;

  bil_copy_elem_const_rw(td.cr,td.cw,s,bd);
  bil_calc_node_r_th(r,th,zeta_t,eta_t);

  bil_t_eta(td.vt,td.zeta_t,td.eta_t,td.cr);
  i_aT=1.0/vabs_d(td.vt);

  for(i=0;i<9;i++) dCC[i]=0.0;

  // dG/dt_eta
  for(i=0;i<4;i++){
    td.nid=i;
    dCC[i]=0.0;
    // D_bar(epsilon)
    for(j=0;j<4;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];

      if(SW_DBDT_BDE==0){
        tC=0.0;
        for(m=0;m<GHN;m++) tC+=( bil_sdG_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +bil_sdG_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tC*=0.25*(th[j+1]-th[j]);
        dCC[i]+=tC;
      }
      else {
        dCC[i]+=deintz(bil_sdG_Db_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdG_Db_theta()! err=%f Exit...\n",err);  exit(1);  }
      }
    }

    // D_epsilon
    if(SW_DBDT_DE==0){
      tC=0.0;
      for(m=0;m<GHN;m++) tC+=( bil_sdG_De_theta(0.25*M_PI*bd->xhi[m]-0.75*M_PI,&td)
                              +bil_sdG_De_theta(0.25*M_PI*bd->xhi[m]-0.25*M_PI,&td)
                              +bil_sdG_De_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                              +bil_sdG_De_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
      tC*=0.25*M_PI;
      dCC[i]+=I*tC;
    }
    else {
      dCC[i]+=I*deintz(bil_sdG_De_theta,0.0,2.0*M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdG_De_theta()! err=%f Exit...\n",err);  exit(1);  }
    }

    // D_epsilon PV
    if(SW_DBDT_DE==0){
      tC=0.0;
      for(m=0;m<GHN;m++) tC+=( bil_sdG_DP_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                              +bil_sdG_DP_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
      tC*=0.25*M_PI;
      dCC[i]-=tC;
    }
    else {
      dCC[i]-=deintz(bil_sdG_DP_theta,0.0,M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_tz_PV(), bil_sdG_DP_theta()! err=%f Exit...\n",err);  exit(1);  }
    }

    dCC[i]*=-I_FP*i_aT;
  }

  if(bil_check_plane(td.cr,td.cw)==1){ // element is plane
    for(l=0;l<5;l++) dCC[l+4]=0.0;
  }
  else {
    bdC=td.cr[0][1]*td.cw[0][2]+td.cr[1][1]*td.cw[1][2]+td.cr[2][1]*td.cw[2][2];
    // dH/dt_eta
    for(i=0;i<4;i++){
      td.nid=i;
      dCC[i+4]=0.0;
      // D\bar(epsilon)
      for(j=0;j<4;j++){
        td.r0=r[j];      td.r1=r[j+1];
        td.th0=th[j];    td.th1=th[j+1];

        if(SW_DBDT_BDE==0){
          tC=0.0;
          for(m=0;m<GHN;m++) tC+=( bil_sdHe_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                  +bil_sdHe_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
          tC*=0.25*(th[j+1]-th[j]);
          dCC[i+4]+=tC;
        }
        else {
          dCC[i+4]+=deintz(bil_sdHe_Db_theta,th[j],th[j+1],&td,IEPS,&err);
          if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_te_PV(), bil_sdHe_Db_theta()! err=%f Exit...\n",err);  exit(1);  }
        }
      }

      // D_epsilon
      if(SW_DBDT_DE==0){
        tC=0.0;
        for(m=0;m<GHN;m++) tC+=( bil_sdHe_De_theta(0.25*M_PI*bd->xhi[m]-0.75*M_PI,&td)
                                +bil_sdHe_De_theta(0.25*M_PI*bd->xhi[m]-0.25*M_PI,&td)
                                +bil_sdHe_De_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                                +bil_sdHe_De_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
        tC*=0.25*M_PI;
        dCC[i+4]+=tC;
      }
      else {
        dCC[i+4]+=deintz(bil_sdHe_De_theta,0.0,M_PI,&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_te_PV(), bil_sdHe_De_theta()0! err=%f Exit...\n",err);  exit(1);  }
        dCC[i+4]+=deintz(bil_sdHe_De_theta,M_PI,2.0*M_PI,&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_te_PV(), bil_sdHe_De_theta()1! err=%f Exit...\n",err);  exit(1);  }
      }

      // D_epsilon PV
      if(SW_DBDT_DE==0){
        tC=0.0;
        for(m=0;m<GHN;m++) tC+=( bil_sdHe_DP_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                                +bil_sdHe_DP_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
        tC*=0.25*M_PI;
        dCC[i+4]-=tC;
      }
      else {
        dCC[i+4]-=deintz(bil_sdHe_DP_theta,0.0,M_PI,&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_te_PV(), bil_sdHe_DP_theta()! err=%f Exit...\n",err);  exit(1);  }
      }

      dCC[i+4]*=bdC*I_FP*i_aT;
    }
    // dF/dt_eta
    dCC[8]=0.0;
    // D\bar(epsilon)
    for(j=0;j<4;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];

      if(SW_DBDT_BDE==0){
        tD=0.0;
        for(m=0;m<GHN;m++) tD+=( bil_sdFe_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +bil_sdFe_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tD*=0.25*(th[j+1]-th[j]);
        dCC[8]+=tD;
      }
      else {
        dCC[8]+=deintd(bil_sdFe_Db_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_te_PV(), bil_sdFe_Db_theta()! err=%f Exit...\n",err);  exit(1);  }
      }
    }

    // D_epsilon PV
    if(SW_DBDT_DE==0){
      tD=0.0;
      for(m=0;m<GHN;m++) tD+=( bil_sdFe_DP_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                              +bil_sdFe_DP_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
      tD*=0.25*M_PI;
      dCC[8]+=tD;
    }
    else {
      dCC[8]+=deintd(bil_sdFe_DP_theta,0.0,M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_dcoef_bd_te_PV(), bil_sdFe_DP_theta()! err=%f Exit...\n",err);  exit(1);  }
    }

    dCC[8]*=bdC*I_FP*i_aT;
  }

  if(sig==1) bil_convert_dCC(dCC);
}

void bil_dcoef_bd_tz_DV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  bil_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  bil_dcoef_bd_tz0_DV(dCC,zeta_t,eta_t,s,k,bd);
}

void bil_dcoef_bd_tz0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex Cp[9],Cm[9];
  double r[5],th[5],vT[3],i_aT;
  int sig,i;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  td.zeta_t=zeta_t;
  td.eta_t =eta_t;
  bil_copy_elem_const_rw(td.cr,td.cw,s,bd);
  bil_calc_node_r_th(r,th,zeta_t,eta_t);

  bil_t_zeta(vT,td.zeta_t,td.eta_t,td.cr);
  i_aT=1.0/vabs_d(vT);

  for(i=0;i<9;i++) dCC[i]=0.0;

  bil_coef_bd(Cp,zeta_t+CDFH,eta_t,s,k,bd);
  bil_coef_bd(Cm,zeta_t-CDFH,eta_t,s,k,bd);
  for(i=0;i<4;i++){
    dCC[i+0]=i_aT*(Cp[i+0]-Cm[i+0])/(2.0*CDFH);
    dCC[i+4]=i_aT*(Cp[i+4]-Cm[i+4])/(2.0*CDFH);
  }
  dCC[8]=i_aT*(Cp[8]-Cm[8])/(2.0*CDFH);

  if(sig==1) bil_convert_dCC(dCC);
}

void bil_dcoef_bd_te_DV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  bil_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  bil_dcoef_bd_te0_DV(dCC,zeta_t,eta_t,s,k,bd);
}

void bil_dcoef_bd_te0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex Cp[9],Cm[9];
  double r[5],th[5],vT[3],i_aT;
  int sig,i;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  td.zeta_t=zeta_t;
  td.eta_t = eta_t;
  bil_copy_elem_const_rw(td.cr,td.cw,s,bd);
  bil_calc_node_r_th(r,th,zeta_t,eta_t);

  bil_t_eta(vT,td.zeta_t,td.eta_t,td.cr);
  i_aT=1.0/vabs_d(vT);

  for(i=0;i<9;i++) dCC[i]=0.0;

  bil_coef_bd(Cp,zeta_t,eta_t+CDFH,s,k,bd);
  bil_coef_bd(Cm,zeta_t,eta_t-CDFH,s,k,bd);
  for(i=0;i<4;i++){
    dCC[i+0]=i_aT*(Cp[i+0]-Cm[i+0])/(2.0*CDFH);
    dCC[i+4]=i_aT*(Cp[i+4]-Cm[i+4])/(2.0*CDFH);
  }
  dCC[8]=i_aT*(Cp[8]-Cm[8])/(2.0*CDFH);

  if(sig==1) bil_convert_dCC(dCC);
}

/////////////////////////////////////////////////////////////////////////////
double complex bil_sdG_zeta(double zeta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta=zeta;
  res= deintz(bil_sdG_eta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("bil_dcoef.c, DE integration error, bil_sdG_zeta(). Exit...\n");  exit(1);  }
  return res;
}

double complex bil_sdG_eta(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ca,cb,ce;
  double r[3],vW[3],vR[3],aR,i_aR,aW,rdv;

  bil_rw_zeta_eta(r,vW,td->zeta,eta,td->cr,td->cw);
  vsub_d(vR,r,td->rt);
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  aW=vabs_d(vW);
  ca=I*td->k*aR;
  cb=ca-1.0;
  if(cimag(td->k)==0.0) ce=cos(creal(td->k)*aR)+I*sin(creal(td->k)*aR);
  else ce=cexp(ca);
  rdv=vdot_d(vR,td->vt)*i_aR;

  return cb*ce*i_aR*i_aR*rdv*aW*bil_Nn(td->nid,td->zeta,eta);
}

double complex bil_sdH_zeta(double zeta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta=zeta;
  res= deintz(bil_sdH_eta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("bil_dcoef.c, DE integration error, bil_sdH_zeta(). err=%15.14e. Exit...\n",err);  exit(1);  }
  return res;
}

double complex bil_sdH_eta(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ca,cb,ce;
  double r[3],vR[3],vW[3],aR,i_aR,rdv,rdW,Wdv;

  bil_rw_zeta_eta(r,vW,td->zeta,eta,td->cr,td->cw);
  vsub_d(vR,r,td->rt);
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  ca=I*td->k*aR;
  cb=ca-1.0;
  if(cimag(td->k)==0.0) ce=cos(creal(td->k)*aR)+I*sin(creal(td->k)*aR);
  else ce=cexp(ca);
  rdW=vdot_d(vR,vW)*i_aR;
  rdv=vdot_d(vR,td->vt)*i_aR;
  Wdv=vdot_d(vW,td->vt);

  return ce*i_aR*i_aR*i_aR*((td->k*td->k*aR*aR+3.0*cb)*rdv*rdW-cb*Wdv)*bil_Nn(td->nid,td->zeta,eta);
}

double bil_sdF_zeta(double zeta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double res,err;

  td->zeta=zeta;
  res= deintd(bil_sdF_eta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("bil_dcoef.c, DE integration error, bil_sdF_zeta(). Exit...\n");  exit(1);  }
  return res;
}

double bil_sdF_eta(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double r[3],w[3],aR,i_aR,vR[3],rdW,rdv,Wdv;

  bil_rw_zeta_eta(r,w,td->zeta,eta,td->cr,td->cw);
  vsub_d(vR,r,td->rt);
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  rdW=vdot_d(vR,w)*i_aR;
  rdv=vdot_d(vR,td->vt)*i_aR;
  Wdv=vdot_d(w,td->vt);

  return i_aR*i_aR*i_aR*(3.0*rdv*rdW-Wdv);
}

double complex bil_sdG_Db_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm,err;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("integral region error bil_dcoef.c, bil_sdG_Db_theta()! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_BDE==0){
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdG_Db_r(0.5*(rm-DEPS)*td->xhi[i]+0.5*(rm+DEPS),td)*td->whi[i];
    ret*=0.5*(rm-DEPS);
  }
  else {
    ret=deintz(bil_sdG_Db_r,DEPS,rm,td,IEPS,&err);
    if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_sdG_Db_theta(), bil_sdG_Db_r()! err=%f Exit...\n",err);  exit(1);  }
  }

  return ret;
}

double complex bil_sdG_Db_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex ca,cb,ce;
  double K[3],W[3],aK,i_aK,aW,zeta,eta,KdT;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
          +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
    W[i]=td->cw[i][0]+td->cw[i][1]*zeta+td->cw[i][2]*eta;
  }
  aK=vabs_d(K);
  i_aK=1.0/aK;
  aW=vabs_d(W);
  ca=I*td->k*aK*r;
  cb=ca-1.0;
  if(cimag(td->k)==0.0) ce=cos(creal(td->k)*aK*r)+I*sin(creal(td->k)*aK*r);
  else ce=cexp(ca);
  KdT=vdot_d(K,td->vt)*i_aK;

  return cb*ce*i_aK*i_aK/r*KdT*bil_Nn(td->nid,zeta,eta)*aW;
}

double complex bil_sdG_De_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  int i;

  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_DE==0){
    ret=0.0;
    for(i=0;i<td->gn;i++) ret+=bil_sdG_De_r(0.5*DEPS*(td->xi[i]+1.0),td)*td->wi[i];
    ret*=0.5*DEPS;
  }
  else {
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdG_De_r(0.5*DEPS*(td->xhi[i]+1.0),td)*td->whi[i];
    ret*=0.5*DEPS;
  }

  return ret;
}

double complex bil_sdG_De_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kR,ca,cb,cc,cs;
  double K[3],W[3],aK,i_aK,aW,zeta,eta,KdT;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
          +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
    W[i]=td->cw[i][0]+td->cw[i][1]*zeta+td->cw[i][2]*eta;
  }
  aK=vabs_d(K);
  i_aK=1.0/aK;
  aW=vabs_d(W);
  kR=td->k*aK*r;
  ca=I*kR;
  cb=ca-1.0;
  if(cimag(td->k)==0.0){
    cc=cos(creal(kR));
    cs=sin(creal(kR));
  }
  else{
    cc=ccos(kR);
    cs=csin(kR);
  }
  KdT=vdot_d(K,td->vt)*i_aK;

  return (cb*cs+kR*cc)*i_aK*i_aK/r*KdT*aW*bil_Nn(td->nid,zeta,eta);
}

double complex bil_sdG_DP_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  int i;

  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_DE==0){
    ret=0.0;
    for(i=0;i<td->gn;i++) ret+=bil_sdG_DP_r(0.5*DEPS*(td->xi[i]+1.0),td)*td->wi[i];
    ret*=0.5*DEPS;
  }
  else {
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdG_DP_r(0.5*DEPS*(td->xhi[i]+1.0),td)*td->whi[i];
    ret*=0.5*DEPS;
  }

  return ret;
}

double complex bil_sdG_DP_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kRp,kRm,ccp,ccm;
  double Kp[3],Wp[3],Km[3],Wm[3],aKp,i_aKp,aWp,aKm,i_aKm,aWm,ztp,etp,ztm,etm,KpdT,KmdT;
  int i;

  ztp= r*td->cth+td->zeta_t;
  etp= r*td->sth+td->eta_t;
  ztm=-r*td->cth+td->zeta_t;
  etm=-r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    Kp[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
          +td->cr[i][3]*( r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
    Wp[i]=td->cw[i][0]+td->cw[i][1]*ztp+td->cw[i][2]*etp;
    Km[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
          +td->cr[i][3]*(-r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
    Wm[i]=td->cw[i][0]+td->cw[i][1]*ztm+td->cw[i][2]*etm;
  }
  aKp=vabs_d(Kp);
  i_aKp=1.0/aKp;
  aWp=vabs_d(Wp);
  kRp=td->k*aKp*r;
  aKm=vabs_d(Km);
  i_aKm=1.0/aKm;
  aWm=vabs_d(Wm);
  kRm=-td->k*aKm*r;
  if(cimag(td->k)==0.0){
    ccp=cos(creal(kRp));
    ccm=cos(creal(kRm));
  }
  else{
    ccp=ccos(kRp);
    ccm=ccos(kRm);
  }
  KpdT=vdot_d(Kp,td->vt)*i_aKp;
  KmdT=vdot_d(Km,td->vt)*i_aKm;

  return 1.0/r*( ccp*i_aKp*i_aKp*KpdT*bil_Nn(td->nid,ztp,etp)*aWp
                -ccm*i_aKm*i_aKm*KmdT*bil_Nn(td->nid,ztm,etm)*aWm);
}

double complex bil_sdHz_Db_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm,err;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("integral region error bil_dcoef.c, bil_sdHz_Db_theta()! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_BDE==0){
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdHz_Db_r(0.5*(rm-DEPS)*td->xhi[i]+0.5*(rm+DEPS),td)*td->whi[i];
    ret*=0.5*(rm-DEPS);
  }
  else {
    ret=deintz(bil_sdHz_Db_r,DEPS,rm,td,IEPS,&err);
    if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_sdHz_Db_theta()! err=%f Exit...\n",err);  exit(1);  }
  }

  return td->sth*ret;
}

double complex bil_sdHz_Db_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kR,ca,cb,ce;
  double K[3],aK,zeta,eta,KdT;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++) K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
                        +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  aK=vabs_d(K);
  kR=td->k*r*aK;
  ca=I*kR;
  cb=ca-1.0;
  if(cimag(td->k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
  else ce=cexp(ca);
  KdT=vdot_d(K,td->vt);

  return ce/(r*aK*aK*aK*aK*aK)*( (kR*kR+3.0*cb)*KdT*td->cth-cb*aK*aK )*bil_Nn(td->nid,zeta,eta);
}

double complex bil_sdHz_De_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret,tr;
  int i;

  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_DE==0){
    tr=0.0;
    for(i=0;i<td->gn;i++) tr+=bil_sdHz_De0_r(0.5*DEPS*(td->xi[i]+1.0),td)*td->wi[i];
    tr*=0.5*DEPS;

    ret=tr;

    tr=0.0;
    for(i=0;i<td->gn;i++) tr+=bil_sdHz_De1_r(0.5*DEPS*(td->xi[i]+1.0),td)*td->wi[i];
    tr*=0.5*DEPS;
  }
  else {
    tr=0.0;
    for(i=0;i<td->ghn;i++) tr+=bil_sdHz_De0_r(0.5*DEPS*(td->xhi[i]+1.0),td)*td->whi[i];
    tr*=0.5*DEPS;

    ret=tr;

    tr=0.0;
    for(i=0;i<td->ghn;i++) tr+=bil_sdHz_De1_r(0.5*DEPS*(td->xhi[i]+1.0),td)*td->whi[i];
    tr*=0.5*DEPS;
  }

  ret-=I*tr;

  return td->sth*ret;
}

double complex bil_sdHz_De0_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kR,ca,cb,ce;
  double K[3],aK,zeta,eta,KdT;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
          +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  }
  aK=vabs_d(K);
  kR=td->k*aK*r;
  ca=I*kR;
  cb=I*td->k;
  if(cimag(td->k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
  else ce=cexp(ca);
  KdT=vdot_d(K,td->vt);

  return ce/(aK*aK*aK*aK)*( (td->k*kR+3.0*cb)*KdT*td->cth-cb*aK*aK )*bil_Nn(td->nid,zeta,eta);
}

double complex bil_sdHz_De1_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kR,cs;
  double K[3],aK,zeta,eta,KdT;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
          +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  }
  aK=vabs_d(K);
  kR=td->k*aK*r;
  if(cimag(td->k)==0.0) cs=sin(creal(kR));
  else cs=csin(kR);
  KdT=vdot_d(K,td->vt);

  return cs/(aK*aK*aK*aK*aK*r)*( 3.0*KdT*td->cth-aK*aK )*bil_Nn(td->nid,zeta,eta);
}

double complex bil_sdHz_DP_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  int i;

  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_DE==0){
    ret=0.0;
    for(i=0;i<td->gn;i++) ret+=bil_sdHz_DP_r(0.5*DEPS*(td->xi[i]+1.0),td)*td->wi[i];
    ret*=0.5*DEPS;
  }
  else {
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdHz_DP_r(0.5*DEPS*(td->xhi[i]+1.0),td)*td->whi[i];
    ret*=0.5*DEPS;
  }

  return td->sth*ret;
}

double complex bil_sdHz_DP_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kRp,kRm,ccp,ccm;
  double Kp[3],Km[3],aKp,aKm,ztp,etp,ztm,etm,KpdT,KmdT;
  int i;

  ztp= r*td->cth+td->zeta_t;
  etp= r*td->sth+td->eta_t;
  ztm=-r*td->cth+td->zeta_t;
  etm=-r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    Kp[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
        +td->cr[i][3]*( r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
    Km[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
        +td->cr[i][3]*(-r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  }
  aKp=vabs_d(Kp);
  kRp= td->k*aKp*r;
  aKm=vabs_d(Km);
  kRm=-td->k*aKm*r;
  if(cimag(td->k)==0.0){
    ccp=cos(creal(kRp));
    ccm=cos(creal(kRm));
  }
  else {
    ccp=ccos(kRp);
    ccm=ccos(kRm);
  }
  KpdT=vdot_d(Kp,td->vt);
  KmdT=vdot_d(Km,td->vt);

  return 1.0/r*( ccp/(aKp*aKp*aKp*aKp*aKp)*(3.0*KpdT*td->cth-aKp*aKp)*bil_Nn(td->nid,ztp,etp)
                -ccm/(aKm*aKm*aKm*aKm*aKm)*(3.0*KmdT*td->cth-aKm*aKm)*bil_Nn(td->nid,ztm,etm) );
}

double bil_sdFz_Db_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double ret,rm,err;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("integral region error bil_dcoef.c, bil_sdFz_Db_theta()! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_BDE==0){
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdFz_Db_r(0.5*(rm-DEPS)*td->xhi[i]+0.5*(rm+DEPS),td)*td->whi[i];
    ret*=0.5*(rm-DEPS);
  }
  else {
    ret=deintd(bil_sdFz_Db_r,DEPS,rm,td,IEPS,&err);
    if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_sdFz_Db_theta(), bil_sdFz_Db_theta(). err=%f Exit...\n",err);  exit(1);  }
  }

  return td->sth*ret;
}

double bil_sdFz_Db_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double K[3],aK,KdT;
  int i;

  for(i=0;i<3;i++) K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
                        +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  aK=vabs_d(K);
  KdT=vdot_d(K,td->vt);

  return 1.0/(r*aK*aK*aK*aK*aK)*( 3.0*KdT*td->cth-aK*aK );
}

double bil_sdFz_DP_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double ret;
  int i;

  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_DE==0){
    ret=0.0;
    for(i=0;i<td->gn;i++) ret+=bil_sdFz_DP_r(0.5*DEPS*(td->xi[i]+1.0),td)*td->wi[i];
    ret*=0.5*DEPS;
  }
  else {
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdFz_DP_r(0.5*DEPS*(td->xhi[i]+1.0),td)*td->whi[i];
    ret*=0.5*DEPS;
  }

   return td->sth*ret;
}

double bil_sdFz_DP_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double Kp[3],Km[3],aKp,aKm,KpdT,KmdT;
  int i;

  for(i=0;i<3;i++){
    Kp[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
        +td->cr[i][3]*( r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
    Km[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
        +td->cr[i][3]*(-r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  }
  aKp=vabs_d(Kp);
  aKm=vabs_d(Km);
  KpdT=vdot_d(Kp,td->vt);
  KmdT=vdot_d(Km,td->vt);

  return 1.0/r*( 1.0/(aKp*aKp*aKp*aKp*aKp)*(3.0*KpdT*td->cth-aKp*aKp)
                -1.0/(aKm*aKm*aKm*aKm*aKm)*(3.0*KmdT*td->cth-aKm*aKm) );
}

double complex bil_sdHe_Db_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm,err;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("integral region error bil_dcoef.c, bil_sdHe_Db_theta()! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_BDE==0){
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdHe_Db_r(0.5*(rm-DEPS)*td->xhi[i]+0.5*(rm+DEPS),td)*td->whi[i];
    ret*=0.5*(rm-DEPS);
  }
  else {
    ret=deintz(bil_sdHe_Db_r,DEPS,rm,td,IEPS,&err);
    if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_sdHe_Db_theta()! err=%f Exit...\n",err);  exit(1);  }
  }

  return td->cth*ret;
}

double complex bil_sdHe_Db_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kR,ca,cb,ce;
  double K[3],aK,zeta,eta,KdT;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++) K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
                        +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  aK=vabs_d(K);
  kR=td->k*r*aK;
  ca=I*kR;
  cb=ca-1.0;
  if(cimag(td->k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
  else ce=cexp(ca);
  KdT=vdot_d(K,td->vt);

  return ce/(r*aK*aK*aK*aK*aK)*( (kR*kR+3.0*cb)*KdT*td->sth-cb*aK*aK )*bil_Nn(td->nid,zeta,eta);
}

double complex bil_sdHe_De_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret,tr;
  int i;

  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_DE==0){
    tr=0.0;
    for(i=0;i<td->gn;i++) tr+=bil_sdHe_De0_r(0.5*DEPS*(td->xi[i]+1.0),td)*td->wi[i];
    tr*=0.5*DEPS;

    ret=tr;

    tr=0.0;
    for(i=0;i<td->gn;i++) tr+=bil_sdHe_De1_r(0.5*DEPS*(td->xi[i]+1.0),td)*td->wi[i];
    tr*=0.5*DEPS;

    ret-=I*tr;
  }
  else {
    tr=0.0;
    for(i=0;i<td->ghn;i++) tr+=bil_sdHe_De0_r(0.5*DEPS*(td->xhi[i]+1.0),td)*td->whi[i];
    tr*=0.5*DEPS;

    ret=tr;

    tr=0.0;
    for(i=0;i<td->ghn;i++) tr+=bil_sdHe_De1_r(0.5*DEPS*(td->xhi[i]+1.0),td)*td->whi[i];
    tr*=0.5*DEPS;

    ret-=I*tr;
  }

  return td->cth*ret;
}

double complex bil_sdHe_De0_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kR,ca,cb,ce;
  double K[3],aK,zeta,eta,KdT;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
          +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  }
  aK=vabs_d(K);
  kR=td->k*aK*r;
  ca=I*kR;
  cb=I*td->k;
  if(cimag(td->k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
  else ce=cexp(ca);
  KdT=vdot_d(K,td->vt);

  return ce/(aK*aK*aK*aK)*( (td->k*kR+3.0*cb)*KdT*td->sth-cb*aK*aK )*bil_Nn(td->nid,zeta,eta);
}

double complex bil_sdHe_De1_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kR,cs;
  double K[3],aK,zeta,eta,KdT;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
          +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  }
  aK=vabs_d(K);
  kR=td->k*aK*r;
  if(cimag(td->k)==0.0) cs=sin(creal(kR));
  else cs=csin(kR);
  KdT=vdot_d(K,td->vt);

  return cs/(aK*aK*aK*aK*aK*r)*( 3.0*KdT*td->sth-aK*aK )*bil_Nn(td->nid,zeta,eta);
}

double complex bil_sdHe_DP_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  int i;

  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_DE==0){
    ret=0.0;
    for(i=0;i<td->gn;i++) ret+=bil_sdHe_DP_r(0.5*DEPS*(td->xi[i]+1.0),td)*td->wi[i];
    ret*=0.5*DEPS;
  }
  else {
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdHe_DP_r(0.5*DEPS*(td->xhi[i]+1.0),td)*td->whi[i];
    ret*=0.5*DEPS;
  }

  return td->cth*ret;
}

double complex bil_sdHe_DP_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex kRp,kRm,ccp,ccm;
  double Kp[3],Km[3],aKp,aKm,ztp,etp,ztm,etm,KpdT,KmdT;
  int i;

  ztp= r*td->cth+td->zeta_t;
  etp= r*td->sth+td->eta_t;
  ztm=-r*td->cth+td->zeta_t;
  etm=-r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    Kp[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
        +td->cr[i][3]*( r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
    Km[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
        +td->cr[i][3]*(-r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  }
  aKp=vabs_d(Kp);
  kRp= td->k*aKp*r;
  aKm=vabs_d(Km);
  kRm=-td->k*aKm*r;
  if(cimag(td->k)==0.0){
    ccp=cos(creal(kRp));
    ccm=cos(creal(kRm));
  }
  else {
    ccp=ccos(kRp);
    ccm=ccos(kRm);
  }
  KpdT=vdot_d(Kp,td->vt);
  KmdT=vdot_d(Km,td->vt);

  return 1.0/r*( ccp/(aKp*aKp*aKp*aKp*aKp)*(3.0*KpdT*td->sth-aKp*aKp)*bil_Nn(td->nid,ztp,etp)
                -ccm/(aKm*aKm*aKm*aKm*aKm)*(3.0*KmdT*td->sth-aKm*aKm)*bil_Nn(td->nid,ztm,etm) );
}

double bil_sdFe_Db_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double ret,rm,err;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("integral region error bil_dcoef.c, bil_sdFe_Db_theta()! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_BDE==0){
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdFe_Db_r(0.5*(rm-DEPS)*td->xhi[i]+0.5*(rm+DEPS),td)*td->whi[i];
    ret*=0.5*(rm-DEPS);
  }
  else {
    ret=deintd(bil_sdFe_Db_r,DEPS,rm,td,IEPS,&err);
    if(err<0.0){ printf("DE integration error bil_dcoef.c, bil_sdFe_Db_theta(), bil_sdFe_Db_theta(). err=%f Exit...\n",err);  exit(1);  }
  }

  return td->cth*ret;
}

double bil_sdFe_Db_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double K[3],aK,KdT;
  int i;

  for(i=0;i<3;i++) K[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
                        +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  aK=vabs_d(K);
  KdT=vdot_d(K,td->vt);

  return 1.0/(r*aK*aK*aK*aK*aK)*( 3.0*KdT*td->sth-aK*aK );
}

double bil_sdFe_DP_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double ret;
  int i;

  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_DE==0){
    ret=0.0;
    for(i=0;i<td->gn;i++) ret+=bil_sdFe_DP_r(0.5*DEPS*(td->xi[i]+1.0),td)*td->wi[i];
    ret*=0.5*DEPS;
  }
  else {
    ret=0.0;
    for(i=0;i<td->ghn;i++) ret+=bil_sdFe_DP_r(0.5*DEPS*(td->xhi[i]+1.0),td)*td->whi[i];
    ret*=0.5*DEPS;
  }

   return td->cth*ret;
}

double bil_sdFe_DP_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double Kp[3],Km[3],aKp,aKm,KpdT,KmdT;
  int i;

  for(i=0;i<3;i++){
    Kp[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
        +td->cr[i][3]*( r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
    Km[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
        +td->cr[i][3]*(-r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth);
  }
  aKp=vabs_d(Kp);
  aKm=vabs_d(Km);
  KpdT=vdot_d(Kp,td->vt);
  KmdT=vdot_d(Km,td->vt);

  return 1.0/r*( 1.0/(aKp*aKp*aKp*aKp*aKp)*(3.0*KpdT*td->sth-aKp*aKp)
                -1.0/(aKm*aKm*aKm*aKm*aKm)*(3.0*KmdT*td->sth-aKm*aKm) );
}

