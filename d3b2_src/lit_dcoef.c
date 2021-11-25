/*
 * lit_dcoef.c
 *
 *  Created on: Jan 3, 2019
 *      Author: ohta
 */
#include "d3b1_elem.h"

void lit_sdGL_zeta_eta(double complex *ret,TDATA *td);
double complex lit_sdG_zeta_eta(TDATA *td);
double complex lit_sdH_zeta_eta(TDATA *td);
double lit_sdF_zeta_eta(TDATA *td);

double complex lit_sdG_alpha(double alpha,void *tmp);
double complex lit_sdG_beta(double beta,void *tmp);
double complex lit_sdH_alpha(double alpha,void *tmp);
double complex lit_sdH_beta(double beta,void *tmp);
double lit_sdF_alpha(double alpha,void *tmp);
double lit_sdF_beta(double beta,void *tmp);

double lit_dNn_pr_cs(int n,double cth,double sth);
double complex lit_sdG_Db_theta(double theta,void *tmp);
double complex lit_sdG_Db_r(double r,void *tmp);
double complex lit_sdG_De_theta(double theta,void *tmp);
double complex lit_sdG_DP_theta(double theta,void *tmp);



void lit_dcoef_4p(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
{
  double complex ca,cb,ce,kR;
  double vR[3],*vW,aR,i_aR,aW,rdW,rdv,Wdv,Nn;
  int sig,i,l;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  for(i=0;i<9;i++){
    CC [i]=0.0;
    dCC[i]=0.0;
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
    rdv=vdot_d(vR,vt)*i_aR;
    Wdv=vdot_d(vW,vt);

    CC[i+0]=bd->wt_34[i]*ce*i_aR;
    CC[i+4]=bd->wt_34[i]*(ca-1.0)*ce*rdW*i_aR*i_aR;

    dCC[i+0]=bd->wt_34[i]*cb*ce*i_aR*i_aR*rdv; // dG/dv
    dCC[i+4]=bd->wt_34[i]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv*rdW-cb*Wdv ); // dH/dv
    CC [8]+=bd->wt_34[i]*rdW*i_aR*i_aR; // F
    dCC[8]+=bd->wt_34[i]*i_aR*i_aR*i_aR*( 3.0*rdv*rdW-Wdv ); // dF/dv
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
  rdv=vdot_d(vR,vt)*i_aR;
  Wdv=vdot_d(vW,vt);

  for(i=0;i<3;i++){
    Nn=lit_Nn(i,bd->zt_34[3],bd->et_34[3]);
    CC[i+0]+=bd->wt_34[3]*ce*i_aR*Nn;
    CC[i+4]+=bd->wt_34[3]*(ca-1.0)*ce*rdW*i_aR*i_aR*Nn;
    dCC[i+0]+=bd->wt_34[3]*cb*ce*i_aR*i_aR*rdv*Nn;
    dCC[i+4]+=bd->wt_34[3]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv*rdW-cb*Wdv )*Nn;
    if(i==0){
      CC [8]+=bd->wt_34[3]*rdW*i_aR*i_aR;
      dCC[8]+=bd->wt_34[3]*i_aR*i_aR*i_aR*( 3.0*rdv*rdW-Wdv );
    }
  }

  for(i=0;i<3;i++){
    CC[i+0]*=I_EP*aW;
    CC[i+4]*=I_EP;
    dCC[i+0]*=-I_EP*aW;
    dCC[i+4]*= I_EP;
    if(i==0){
      CC [8]*=I_EP;
      dCC[8]*=I_EP;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    lit_convert_dCC(dCC);
  }
}

void lit_dcoef_7p(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
{
  double complex ca,cb,ce,kR;
  double cr[3][4],cw[3][3],zeta,eta,vW[3],aW,vR[3],aR,i_aR,rdW,rdv,Wdv,Np;
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
    dCC[l]=0.0;
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
      rdv=vdot_d(vR,vt)*i_aR;
      Wdv=vdot_d(vW,vt);
      Np=lit_Nn(ep,zeta,eta);

      CC[ep  ]+=bd->wt_37[i]*ce*i_aR*Np;
      CC[ep+4]+=bd->wt_37[i]*(ca-1.0)*ce*rdW*i_aR*i_aR*Np;
      dCC[ep  ]+=bd->wt_37[i]*cb*ce*i_aR*i_aR*rdv*Np;
      dCC[ep+4]+=bd->wt_37[i]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv*rdW-cb*Wdv )*Np;
      if(ep==0){
        CC [8]+=bd->wt_37[i]*rdW*i_aR*i_aR;
        dCC[8]+=bd->wt_37[i]*i_aR*i_aR*i_aR*( 3.0*rdv*rdW-Wdv );
      }
    }
    CC [ep  ]*= I_EP*aW;
    CC [ep+4]*= I_EP;
    dCC[ep  ]*=-I_EP*aW;
    dCC[ep+4]*= I_EP;
    if(ep==0){
      CC [8]*=I_EP;
      dCC[8]*=I_EP;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    lit_convert_dCC(dCC);
  }
}

void lit_dcoef_GL(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex ret[6],tpj[6],tpi[6];
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
    td.vt[l]=vt[l];
  }

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dCC[l]=0.0;
  }
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GLN;i++){
      for(l=0;l<6;l++) tpi[l]=0.0;
      a=bd->xli[i];
      for(j=0;j<GLN;j++){
        for(l=0;l<6;l++) tpj[l]=0.0;
        b=bd->xli[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        lit_sdGL_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        lit_sdGL_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        lit_sdGL_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];

        for(l=0;l<6;l++) tpi[l]+=tpj[l]*wt*bd->wli[j];
      }
      CC [ep  ]+=tpi[0]*bd->wli[i];
      CC [ep+4]+=tpi[1]*bd->wli[i];
      dCC[ep  ]+=tpi[3]*bd->wli[i];
      dCC[ep+4]+=tpi[4]*bd->wli[i];
      if(ep==0){
        CC [8]+=tpi[2]*bd->wli[i];
        dCC[8]+=tpi[5]*bd->wli[i];
      }
    }
    CC[ep  ]*=0.0625*SQ3*I_SPSQ3*td.aW;
    CC[ep+4]*=0.0625*SQ3*I_SPSQ3;
    dCC[ep  ]*=-0.0625*SQ3*I_SPSQ3*td.aW;
    dCC[ep+4]*= 0.0625*SQ3*I_SPSQ3;
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      dCC[8]*=0.0625*SQ3*I_SPSQ3;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    lit_convert_dCC(dCC);
  }
}

void lit_dcoef_GH(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex ret[6],tpj[6],tpi[6];
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
    td.vt[l]=vt[l];
  }

  for(l=0;l<9;l++){
    CC [l]=0.0;
    dCC[l]=0.0;
  }
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GHN;i++){
      for(l=0;l<6;l++) tpi[l]=0.0;
      a=bd->xhi[i];
      for(j=0;j<GHN;j++){
        for(l=0;l<6;l++) tpj[l]=0.0;
        b=bd->xhi[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        lit_sdGL_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        lit_sdGL_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        lit_sdGL_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];

        for(l=0;l<6;l++) tpi[l]+=tpj[l]*wt*bd->whi[j];
      }
      CC [ep  ]+=tpi[0]*bd->whi[i];
      CC [ep+4]+=tpi[1]*bd->whi[i];
      dCC[ep  ]+=tpi[3]*bd->whi[i];
      dCC[ep+4]+=tpi[4]*bd->whi[i];
      if(ep==0){
        CC [8]+=tpi[2]*bd->whi[i];
        dCC[8]+=tpi[5]*bd->whi[i];
      }
    }
    CC[ep  ]*=0.0625*SQ3*I_SPSQ3*td.aW;
    CC[ep+4]*=0.0625*SQ3*I_SPSQ3;
    dCC[ep  ]*=-0.0625*SQ3*I_SPSQ3*td.aW;
    dCC[ep+4]*= 0.0625*SQ3*I_SPSQ3;
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      dCC[8]*=0.0625*SQ3*I_SPSQ3;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    lit_convert_dCC(dCC);
  }
}

void lit_dcoef_DE(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
{
  lit_coef_DE(CC,rt,s,k,bd);
  lit_dcoef_DE0(dCC,rt,vt,s,k,bd);
}

void lit_dcoef_DE0(double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd)
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
  for(l=0;l<3;l++){
    td.rt[l]=rt[l];
    td.vt[l]=vt[l];
  }

  for(l=0;l<9;l++) dCC[l]=0.0;
  // dG/dv
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    dCC[ep+0]=-0.0625*SQ3*I_SPSQ3*td.aW*deintz(lit_sdG_alpha,-1.0,1.0,&td,IEPS,&err);
    if(err<0.0){ printf("lit_coef.c, DE integration error, lit_dcoef_DE(), lit_sdG_alpha(). err=%g Exit...\n",err);  exit(1);  }
  }

  if(lit_check_on_plane(td.rt,td.cr,td.cw)==0){ // rt is not on element plane
    // dH/dv
    for(ep=0;ep<3;ep++){
      td.nid=ep;
      dCC[ep+4]=0.0625*SQ3*I_SPSQ3*deintz(lit_sdH_alpha,-1.0,1.0,&td,IEPS,&err);
      if(err<0.0){ printf("lit_dcoef.c, DE integration error, lit_dcoef_DE(), lit_sdH_alpha(). err=%g Exit...\n",err);  exit(1);  }
    }
    // dF/dv
    dCC[8]=0.0625*SQ3*I_SPSQ3*deintd(lit_sdF_alpha,-1.0,1.0,&td,IEPS,&err);
    if(err<0.0){ printf("lit_dcoef.c, DE integration error, lit_dcoef_DE(), lit_sdF_alpha(). err=%g Exit...\n",err);  exit(1);  }
  }
  else {
    for(ep=0;ep<3;ep++)  dCC[ep+4]=0.0;
    dCC[8]=0.0;
  }

  if(sig==1) lit_convert_dCC(dCC);
}

void lit_dcoef_bd_tz_DV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  lit_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  lit_dcoef_bd_tz0_DV(dCC,zeta_t,eta_t,s,k,bd);
}

void lit_dcoef_bd_tz0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double complex Cp[9],Cm[9];
  double vT[3],i_aT,cr[3][4],cw[3][3];
  int sig,i;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  lit_t_zeta(vT,zeta_t,eta_t,cr);
  i_aT=1.0/vabs_d(vT);

  for(i=0;i<9;i++) dCC[i]=0.0;

  lit_coef_bd(Cp,zeta_t+CDFH,eta_t,s,k,bd);
  lit_coef_bd(Cm,zeta_t-CDFH,eta_t,s,k,bd);
  for(i=0;i<3;i++){
    dCC[i+0]=i_aT*(Cp[i+0]-Cm[i+0])/(2.0*CDFH);
    dCC[i+4]=i_aT*(Cp[i+4]-Cm[i+4])/(2.0*CDFH);
  }

  if(sig==1) lit_convert_dCC(dCC);
}

void lit_dcoef_bd_te_DV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  lit_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  lit_dcoef_bd_te0_DV(dCC,zeta_t,eta_t,s,k,bd);
}

void lit_dcoef_bd_te0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double complex Cp[9],Cm[9];
  double vT[3],i_aT,cr[3][4],cw[3][3];
  int sig,i;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  lit_t_eta(vT,zeta_t,eta_t,cr);
  i_aT=1.0/vabs_d(vT);

  for(i=0;i<9;i++) dCC[i]=0.0;

  lit_coef_bd(Cp,zeta_t,eta_t+CDFH,s,k,bd);
  lit_coef_bd(Cm,zeta_t,eta_t-CDFH,s,k,bd);
  for(i=0;i<3;i++){
    dCC[i+0]=i_aT*(Cp[i+0]-Cm[i+0])/(2.0*CDFH);
    dCC[i+4]=i_aT*(Cp[i+4]-Cm[i+4])/(2.0*CDFH);
  }

  if(sig==1) lit_convert_dCC(dCC);
}

void lit_dcoef_bd_tz_PV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  lit_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  lit_dcoef_bd_tz0_PV(dCC,zeta_t,eta_t,s,k,bd);
}

void lit_dcoef_bd_tz0_PV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex tC;
  double r[4],th[4],err,i_aT,vW[3],aW;
  int sig,i,j,m;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  td.zeta_t=zeta_t;
  td.eta_t = eta_t;

  td.gn=GHN;
  td.xi=bd->xhi;
  td.wi=bd->whi;

  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_calc_node_r_th(r,th,zeta_t,eta_t);
  lit_w_zeta_eta(vW,zeta_t,eta_t,td.cw);
  aW=vabs_d(vW);
  lit_t_zeta(td.vt,zeta_t,eta_t,td.cr);
  i_aT=1.0/vabs_d(td.vt);

  for(i=0;i<9;i++) dCC[i]=0.0;

  // dG/dt_zeta
  for(i=0;i<3;i++){
    td.nid=i;
    dCC[i]=0.0;
    // D_bar(epsilon)
    for(j=0;j<3;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];

      if(SW_DBDT_BDE==0){
        tC=0.0;
        for(m=0;m<GHN;m++) tC+=( lit_sdG_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +lit_sdG_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tC*=0.25*(th[j+1]-th[j]);
        dCC[i]+=tC;
      }
      else {
        dCC[i]+=deintz(lit_sdG_Db_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error lit_dcoef.c, lit_dcoef_bd_tz_PV(), lit_sdG_Db_theta()! err=%f Exit...\n",err);  exit(1);  }
      }
    }
    // D_epsilon
    if(SW_DBDT_DE==0){
      tC=0.0;
      for(m=0;m<GHN;m++) tC+=( lit_sdG_De_theta(0.25*M_PI*bd->xhi[m]-0.75*M_PI,&td)
                              +lit_sdG_De_theta(0.25*M_PI*bd->xhi[m]-0.25*M_PI,&td)
                              +lit_sdG_De_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                              +lit_sdG_De_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
      tC*=0.25*M_PI;
      dCC[i]+=I*tC/k;
    }
    else {
      dCC[i]+=I/k*deintz(lit_sdG_De_theta,0.0,0.5*M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error lit_dcoef.c, lit_dcoef_bd_tz_PV(), lit_sdG_De_theta()1! err=%f Exit...\n",err);  exit(1);  }
      dCC[i]+=I/k*deintz(lit_sdG_De_theta,0.5*M_PI,M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error lit_dcoef.c, lit_dcoef_bd_tz_PV(), lit_sdG_De_theta()2! err=%f Exit...\n",err);  exit(1);  }
    }
    // D_epsilon PV
    if(SW_DBDT_DE==0){
      tC=0.0;
      for(m=0;m<GHN;m++) tC+=( lit_sdG_DP_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                              +lit_sdG_DP_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
      tC*=0.25*M_PI;
      dCC[i]-=tC/k;
    }
    else {
      dCC[i]-=1.0/k*deintz(lit_sdG_DP_theta,0.0,M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error lit_dcoef.c, lit_dcoef_bd_tz_PV(), lit_sdG_DP_theta()! err=%f Exit...\n",err);  exit(1);  }
    }

    dCC[i]*=-I_SPSQ3*aW*i_aT;
  }
  if(sig==1) lit_convert_dCC(dCC);
}

void lit_dcoef_bd_te_PV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  lit_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  lit_dcoef_bd_te0_PV(dCC,zeta_t,eta_t,s,k,bd);
}

void lit_dcoef_bd_te0_PV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  TDATA td;
  double complex tC;
  double r[4],th[4],err,i_aT,vW[3],aW;
  int sig,i,j,m;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  td.zeta_t=zeta_t;
  td.eta_t = eta_t;

  td.gn=GHN;
  td.xi=bd->xhi;
  td.wi=bd->whi;

  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_calc_node_r_th(r,th,zeta_t,eta_t);
  lit_w_zeta_eta(vW,zeta_t,eta_t,td.cw);
  aW=vabs_d(vW);
  lit_t_eta(td.vt,zeta_t,eta_t,td.cr);
  i_aT=1.0/vabs_d(td.vt);

  for(i=0;i<9;i++) dCC[i]=0.0;

  // dG/dt_eta
  for(i=0;i<3;i++){
    td.nid=i;
    dCC[i]=0.0;
    // D_bar(epsilon)
    for(j=0;j<3;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];

      if(SW_DBDT_BDE==0){
        tC=0.0;
        for(m=0;m<GHN;m++) tC+=( lit_sdG_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +lit_sdG_Db_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tC*=0.25*(th[j+1]-th[j]);
        dCC[i]+=tC;
      }
      else {
        dCC[i]+=deintz(lit_sdG_Db_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("DE integration error lit_dcoef.c, lit_dcoef_bd_te_PV(), lit_sdG_Db_theta()! err=%f Exit...\n",err);  exit(1);  }
      }
    }
    // D_epsilon
    if(SW_DBDT_DE==0){
      tC=0.0;
      for(m=0;m<GHN;m++) tC+=( lit_sdG_De_theta(0.25*M_PI*bd->xhi[m]-0.75*M_PI,&td)
                              +lit_sdG_De_theta(0.25*M_PI*bd->xhi[m]-0.25*M_PI,&td)
                              +lit_sdG_De_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                              +lit_sdG_De_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
      tC*=0.25*M_PI;
      dCC[i]+=I*tC/k;
    }
    else {
      dCC[i]+=I/k*deintz(lit_sdG_De_theta,0.0,0.5*M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error lit_dcoef.c, lit_dcoef_bd_te_PV(), lit_sdG_De_theta()1! err=%f Exit...\n",err);  exit(1);  }
      dCC[i]+=I/k*deintz(lit_sdG_De_theta,0.5*M_PI,M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error lit_dcoef.c, lit_dcoef_bd_te_PV(), lit_sdG_De_theta()2! err=%f Exit...\n",err);  exit(1);  }
    }

    // D_epsilon PV
    if(SW_DBDT_DE==0){
      tC=0.0;
      for(m=0;m<GHN;m++) tC+=( lit_sdG_DP_theta(0.25*M_PI*bd->xhi[m]+0.25*M_PI,&td)
                              +lit_sdG_DP_theta(0.25*M_PI*bd->xhi[m]+0.75*M_PI,&td) )*bd->whi[m];
      tC*=0.25*M_PI;
      dCC[i]-=tC/k;
    }
    else {
      dCC[i]-=1.0/k*deintz(lit_sdG_DP_theta,0.0,M_PI,&td,IEPS,&err);
      if(err<0.0){ printf("DE integration error lit_dcoef.c, lit_dcoef_bd_te_PV(), lit_sdG_DP_theta()! err=%f Exit...\n",err);  exit(1);  }
    }
    dCC[i]*=-I_SPSQ3*aW*i_aT;
  }

  if(sig==1) lit_convert_dCC(dCC);
}

///////////////////////////////////////////////////////////////////////////
void lit_sdGL_zeta_eta(double complex *ret,TDATA *td)
{
  double complex ca,cb,ce,kR;
  double vR[3],aR,i_aR,rdW,rdv,Wdv,Np;
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
  rdv=vdot_d(vR,td->vt)*i_aR;
  Wdv=td->cw[0][0]*td->vt[0]+td->cw[1][0]*td->vt[1]+td->cw[2][0]*td->vt[2];
  Np=lit_Nn(td->nid,td->zeta,td->eta);

  ret[0]=ce*i_aR*Np; // G
  ret[1]=(ca-1.0)*ce*rdW*i_aR*i_aR*Np; // H
  ret[2]=rdW*i_aR*i_aR; // F
  ret[3]=cb*ce*i_aR*i_aR*rdv*Np; // dG/dv
  ret[4]=ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv*rdW-cb*Wdv )*Np; // dH/dv
  ret[5]=i_aR*i_aR*i_aR*( 3.0*rdv*rdW-Wdv ); // dF/dv
}

double complex lit_sdG_zeta_eta(TDATA *td)
{
  double complex ca,cb,ce,kR;
  double vR[3],aR,i_aR,rdv,Np;
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
  rdv=vdot_d(vR,td->vt)*i_aR;
  Np=lit_Nn(td->nid,td->zeta,td->eta);

  return cb*ce*i_aR*i_aR*rdv*Np;
}

double complex lit_sdH_zeta_eta(TDATA *td)
{
  double complex ca,cb,ce,kR;
  double vR[3],aR,i_aR,rdW,rdv,Wdv,Np;
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
  rdv=vdot_d(vR,td->vt)*i_aR;
  Wdv=td->cw[0][0]*td->vt[0]+td->cw[1][0]*td->vt[1]+td->cw[2][0]*td->vt[2];
  Np=lit_Nn(td->nid,td->zeta,td->eta);

  return ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*rdv*rdW-cb*Wdv )*Np;
}

double lit_sdF_zeta_eta(TDATA *td)
{
  double vR[3],aR,i_aR,rdW,rdv,Wdv;
  int l;

  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);

  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  rdW=(vR[0]*td->cw[0][0]+vR[1]*td->cw[1][0]+vR[2]*td->cw[2][0])*i_aR;
  rdv=vdot_d(vR,td->vt)*i_aR;
  Wdv=td->cw[0][0]*td->vt[0]+td->cw[1][0]*td->vt[1]+td->cw[2][0]*td->vt[2];

  return i_aR*i_aR*i_aR*( 3.0*rdv*rdW-Wdv );
}


double complex lit_sdG_alpha(double alpha,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta_t=alpha;
  res= deintz(lit_sdG_beta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("lit_dcoef.c, DE integration error, lit_sdG_alpha(). err=%g. Exit...\n",err);  exit(1);  }
  return res;
}

double complex lit_sdG_beta(double b,void *tmp)
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
  ret+=lit_sdG_zeta_eta(td);
  // part 2
  td->zeta=0.0625*(-3.0+a-5.0*b-a*b);
  td->eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
  ret+=lit_sdG_zeta_eta(td);
  // part 3
  td->zeta=0.0625*(-3.0-5.0*a+b-a*b);
  td->eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
  ret+=lit_sdG_zeta_eta(td);

  return ret*wt;
}

double complex lit_sdH_alpha(double alpha,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta_t=alpha;
  res= deintz(lit_sdH_beta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("lit_dcoef.c, DE integration error, lit_sdH_alpha(). err=%g. Exit...\n",err);  exit(1);  }
  return res;
}

double complex lit_sdH_beta(double b,void *tmp)
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
  ret+=lit_sdH_zeta_eta(td);
  // part 2
  td->zeta=0.0625*(-3.0+a-5.0*b-a*b);
  td->eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
  ret+=lit_sdH_zeta_eta(td);
  // part 3
  td->zeta=0.0625*(-3.0-5.0*a+b-a*b);
  td->eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
  ret+=lit_sdH_zeta_eta(td);

  return ret*wt;
}

double lit_sdF_alpha(double alpha,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double res,err;

  td->zeta_t=alpha;
  res= deintd(lit_sdF_beta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("lit_dcoef.c, DE integration error, lit_sdF_alpha(). err=%g. Exit...\n",err);  exit(1);  }
  return res;
}

double lit_sdF_beta(double b,void *tmp)
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
  ret+=lit_sdF_zeta_eta(td);
  // part 2
  td->zeta=0.0625*(-3.0+a-5.0*b-a*b);
  td->eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
  ret+=lit_sdF_zeta_eta(td);
  // part 3
  td->zeta=0.0625*(-3.0-5.0*a+b-a*b);
  td->eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
  ret+=lit_sdF_zeta_eta(td);

  return ret*wt;
}

double lit_dNn_pr_cs(int n,double cth,double sth)
{
  double i_a=0.5/P34_N0;
  switch (n){
  case 0: return -2.0/3.0*i_a*(cth+SQ3*sth);
  case 1: return 4.0/3.0*i_a*cth;
  case 2: return -2.0/3.0*i_a*(cth-SQ3*sth);
  default:
    printf("lit_dcoef.c, lit_dNn_pr_cs(), mode number error! n=%d. Exit...\n",n);
    exit(1);
  }
}

double complex lit_sdG_Db_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm,err;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("integral region error lit_dcoef.c, lit_sdG_Db_theta()! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  if(SW_DBDR_BDE==0){
    ret=0.0;
    for(i=0;i<td->gn;i++) ret+=lit_sdG_Db_r(0.5*(rm-DEPS)*td->xi[i]+0.5*(rm+DEPS),td)*td->wi[i];
    ret*=0.5*(rm-DEPS);
  }
  else {
    ret=deintz(lit_sdG_Db_r,DEPS,rm,td,IEPS,&err);
    if(err<0.0){ printf("DE integration error lit_dcoef.c, lit_sdG_Db_theta(), lit_sdG_Db_r()! err=%f Exit...\n",err);  exit(1);  }
  }

  return ret;
}

double complex lit_sdG_Db_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex ca,cb,ce,kR;
  double vK[3],aK,i_aK,zeta,eta,Kdv;
  int i;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++) vK[i]=td->cr[i][1]*td->cth+td->cr[i][2]*td->sth;
  aK=vabs_d(vK);
  i_aK=1.0/aK;

  kR=td->k*aK*r;
  ca=I*kR;
  cb=ca-1.0;
  if(cimag(td->k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
  else ce=cexp(ca);
  Kdv=vdot_d(vK,td->vt)*i_aK;

  return cb*ce/r*i_aK*i_aK*Kdv*lit_Nn(td->nid,zeta,eta);
}

double complex lit_sdG_De_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex keK,cs,cc;
  double cth,sth,vK[3],aK,i_aK,Kdv;
  int i;

  cth=cos(theta);
  sth=sin(theta);

  for(i=0;i<3;i++) vK[i]=td->cr[i][1]*cth+td->cr[i][2]*sth;
  aK=vabs_d(vK);
  i_aK=1.0/aK;
  Kdv=vdot_d(vK,td->vt)*i_aK;

  keK=td->k*DEPS*aK;
  if(cimag(td->k)==0.0) {
    cc=cos(creal(keK));
    cs=sin(creal(keK));
  }
  else {
    cc=ccos(keK);
    cs=csin(keK);
  }

  return lit_dNn_pr_cs(td->nid,cth,sth)*i_aK*i_aK*i_aK*Kdv*((keK+I)*cs+(2.0-I*keK)*cc-2.0);
}

double complex lit_sdG_DP_theta(double theta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex cs;
  double Kdv,i_aK,sth,cth,vK[3],aK;
  int i;

  cth=cos(theta);
  sth=sin(theta);

  for(i=0;i<3;i++) vK[i]=td->cr[i][1]*cth+td->cr[i][2]*sth; // K
  aK=vabs_d(vK); // |K|
  i_aK=1.0/aK;
  Kdv=vdot_d(vK,td->vt)*i_aK;
  if(cimag(td->k)==0.0) cs=sin(creal(td->k)*aK*DEPS);
  else cs=csin(td->k*aK*DEPS);

  return lit_dNn_pr_cs(td->nid,cth,sth)*cs*i_aK*i_aK*i_aK*Kdv;
}




