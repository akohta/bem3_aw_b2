/*
 * d3b1_elem.c
 *
 *  Created on: Dec 14, 2018
 *      Author: ohta
 */

#include "d3b1_elem.h"

// ---- common function -----------------------------------------------------------------------------------
int check_element_type(int seid,BOUD *bd)
{
  int aeid=abs(seid);
  if(bd->ed[aeid][3]==0) return ELT3; // triangular element
  else return ELT4; // quadrangular element
}

void r_bd(double *r,int seid,double zeta,double eta,BOUD *bd)
{
  int i,s;

  s=abs(seid);
  if(ELT4==check_element_type(seid,bd)){
    for(i=0;i<3;i++) r[i]=bd->cr[s][i][0]+bd->cr[s][i][1]*zeta+bd->cr[s][i][2]*eta+bd->cr[s][i][3]*zeta*eta;
  }
  else {
    for(i=0;i<3;i++) r[i]=bd->cr[s][i][0]+bd->cr[s][i][1]*zeta+bd->cr[s][i][2]*eta;
  }
}

void r_bd_node(double *r,int seid,int nid,BOUD *bd)
{
  int l;

  for(l=0;l<3;l++) r[l]=bd->ren[abs(seid)][nid][l];
}

void n_bd_node(double *n,int seid,int nid,BOUD *bd)
{
  int l;

  for(l=0;l<3;l++) n[l]=bd->wen[abs(seid)][nid][l];
  vuni_d(n);
}

void r_tz_bd(double *r,double *vtz,int seid,double zeta,double eta,BOUD *bd)
{
  int i,s;

  s=abs(seid);
  if(ELT4==check_element_type(seid,bd)){
    for(i=0;i<3;i++){
      r[i]=bd->cr[s][i][0]+bd->cr[s][i][1]*zeta+bd->cr[s][i][2]*eta+bd->cr[s][i][3]*zeta*eta;
      vtz[i]=bd->cr[s][i][1]+bd->cr[s][i][3]*eta;
    }
    vuni_d(vtz);
  }
  else {
    for(i=0;i<3;i++){
      r[i]=bd->cr[s][i][0]+bd->cr[s][i][1]*zeta+bd->cr[s][i][2]*eta;
      vtz[i]=bd->cr[s][i][1];
    }
    vuni_d(vtz);
  }
}

void r_te_bd(double *r,double *vte,int seid,double zeta,double eta,BOUD *bd)
{
  int i,s;

  s=abs(seid);
  if(ELT4==check_element_type(seid,bd)){
    for(i=0;i<3;i++){
      r[i]=bd->cr[s][i][0]+bd->cr[s][i][1]*zeta+bd->cr[s][i][2]*eta+bd->cr[s][i][3]*zeta*eta;
      vte[i]=bd->cr[s][i][2]+bd->cr[s][i][3]*zeta;
    }
    vuni_d(vte);
  }
  else {
    for(i=0;i<3;i++){
      r[i]=bd->cr[s][i][0]+bd->cr[s][i][1]*zeta+bd->cr[s][i][2]*eta;
      vte[i]=bd->cr[s][i][2];
    }
    vuni_d(vte);
  }
}

void Tz_bd(double *vTz,int seid,double zeta,double eta,BOUD *bd)
{
  int i,s;

  s=abs(seid);
  if(ELT4==check_element_type(seid,bd)){
    for(i=0;i<3;i++){
      vTz[i]=bd->cr[s][i][1]+bd->cr[s][i][3]*eta;
    }
  }
  else {
    for(i=0;i<3;i++){
      vTz[i]=bd->cr[s][i][1];
    }
  }
}

void Te_bd(double *vTe,int seid,double zeta,double eta,BOUD *bd)
{
  int i,s;

  s=abs(seid);
  if(ELT4==check_element_type(seid,bd)){
    for(i=0;i<3;i++){
      vTe[i]=bd->cr[s][i][2]+bd->cr[s][i][3]*zeta;
    }
  }
  else {
    for(i=0;i<3;i++){
      vTe[i]=bd->cr[s][i][2];
    }
  }
}

void tz_te_bd_node(double *vtz,double *vte,int seid,int nid,BOUD *bd)
{
  int i,s;

  s=abs(seid);
  if(ELT4==check_element_type(seid,bd)){
    for(i=0;i<3;i++){
      vtz[i]=bd->cr[s][i][1]+bd->cr[s][i][3]*bd->et_44[nid];
      vte[i]=bd->cr[s][i][2]+bd->cr[s][i][3]*bd->zt_44[nid];
    }
    vuni_d(vtz);
    vuni_d(vte);
  }
  else {
    for(i=0;i<3;i++){
      vtz[i]=bd->cr[s][i][1];
      vte[i]=bd->cr[s][i][2];
    }
    vuni_d(vtz);
    vuni_d(vte);
  }
}




void coef_rt(double complex *CC,double *rt,int s,double complex k,int type,BOUD *bd)
{

  if(ELT4==check_element_type(s,bd)){
    switch(type){
    case 0: bil_coef_4p(CC,rt,s,k,bd);
    break;
    case 1: bil_coef_9p(CC,rt,s,k,bd);
    break;
    case 2: bil_coef_GL(CC,rt,s,k,bd);
    break;
    case 3: bil_coef_GH(CC,rt,s,k,bd);
    break;
    case 4: bil_coef_DE(CC,rt,s,k,bd);
    break;
    default:
      printf("b_elem.c, coef_rt(), type number error! type=%d is not defined.Exit...\n",type);
      exit(1);
    }
  }
  else { // ELT3
    switch(type){
    case 0: lit_coef_4p(CC,rt,s,k,bd);
    break;
    case 1: lit_coef_7p(CC,rt,s,k,bd);
    break;
    case 2: lit_coef_GL(CC,rt,s,k,bd);
    break;
    case 3: lit_coef_GH(CC,rt,s,k,bd);
    break;
    case 4: lit_coef_DE(CC,rt,s,k,bd);
    break;
    default:
      printf("b_elem.c, coef_rt(), type number error! type=%d is not defined.Exit...\n",type);
      exit(1);
    }
  }
}

void coef_bd(double complex *CC,double *rt,int t,double zeta_t,double eta_t,int s,double complex k,int type,BOUD *bd)
{

  if(ELT4==check_element_type(s,bd)){
    if(t!=s){
      switch(type){
      case 0: bil_coef_4p(CC,rt,s,k,bd);
      break;
      case 1: bil_coef_9p(CC,rt,s,k,bd);
      break;
      case 2: bil_coef_GL(CC,rt,s,k,bd);
      break;
      case 3: bil_coef_GH(CC,rt,s,k,bd);
      break;
      case 4: bil_coef_DE(CC,rt,s,k,bd);
      break;
      default:
        printf("b_elem.c, coef_bd(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }
    }
    else bil_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  }
  else { // ELT3
    if(t!=s){
      switch(type){
      case 0: lit_coef_4p(CC,rt,s,k,bd);
      break;
      case 1: lit_coef_7p(CC,rt,s,k,bd);
      break;
      case 2: lit_coef_GL(CC,rt,s,k,bd);
      break;
      case 3: lit_coef_GH(CC,rt,s,k,bd);
      break;
      case 4: lit_coef_DE(CC,rt,s,k,bd);
      break;
      default:
        printf("b_elem.c, coef_bd(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }
    }
    else lit_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  }
}

void coef_bd_node(double complex *CC,int t,int nid,int s,double complex k,int type,BOUD *bd)
{

  if(ELT4==check_element_type(s,bd)){
    if(t!=s){
      switch(type){
      case 0: bil_coef_4p(CC,bd->ren[abs(t)][nid],s,k,bd);
      break;
      case 1: bil_coef_9p(CC,bd->ren[abs(t)][nid],s,k,bd);
      break;
      case 2: bil_coef_GL(CC,bd->ren[abs(t)][nid],s,k,bd);
      break;
      case 3: bil_coef_GH(CC,bd->ren[abs(t)][nid],s,k,bd);
      break;
      case 4: bil_coef_DE(CC,bd->ren[abs(t)][nid],s,k,bd);
      break;
      default:
        printf("b_elem.c, coef_bd_node(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }
    }
    else bil_coef_bd(CC,bd->zt_44[nid],bd->et_44[nid],s,k,bd);
  }
  else { // ELT3
    if(t!=s){
      switch(type){
      case 0: lit_coef_4p(CC,bd->ren[abs(t)][nid],s,k,bd);
      break;
      case 1: lit_coef_7p(CC,bd->ren[abs(t)][nid],s,k,bd);
      break;
      case 2: lit_coef_GL(CC,bd->ren[abs(t)][nid],s,k,bd);
      break;
      case 3: lit_coef_GH(CC,bd->ren[abs(t)][nid],s,k,bd);
      break;
      case 4: lit_coef_DE(CC,bd->ren[abs(t)][nid],s,k,bd);
      break;
      default:
        printf("b_elem.c, coef_bd_node(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }
    }
    else lit_coef_bd(CC,bd->zt_34[nid],bd->et_34[nid],s,k,bd);
  }
}

void dcoef_rt(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,int type,BOUD *bd)
{

  if(ELT4==check_element_type(s,bd)){
    switch(type){
    case 0: bil_dcoef_4p(CC,dCC,rt,vt,s,k,bd);
    break;
    case 1: bil_dcoef_9p(CC,dCC,rt,vt,s,k,bd);
    break;
    case 2: bil_dcoef_GL(CC,dCC,rt,vt,s,k,bd);
    break;
    case 3: bil_dcoef_GH(CC,dCC,rt,vt,s,k,bd);
    break;
    case 4: bil_dcoef_DE(CC,dCC,rt,vt,s,k,bd);
    break;
    default:
      printf("b_elem.c, dcoef_rt(), type number error! type=%d is not defined.Exit...\n",type);
      exit(1);
    }
    if(type<3 && fabs(creal(dCC[8]))> CBD_SDF) bil_dcoef_GH(CC,dCC,rt,vt,s,k,bd);
  }
  else { // ELT3
    switch(type){
    case 0: lit_dcoef_4p(CC,dCC,rt,vt,s,k,bd);
    break;
    case 1: lit_dcoef_7p(CC,dCC,rt,vt,s,k,bd);
    break;
    case 2: lit_dcoef_GL(CC,dCC,rt,vt,s,k,bd);
    break;
    case 3: lit_dcoef_GH(CC,dCC,rt,vt,s,k,bd);
    break;
    case 4: lit_dcoef_DE(CC,dCC,rt,vt,s,k,bd);
    break;
    default:
      printf("b_elem.c, dcoef_rt(), type number error! type=%d is not defined.Exit...\n",type);
      exit(1);
    }
    if(type<3 && fabs(creal(dCC[8]))> CBD_SDF) lit_dcoef_GH(CC,dCC,rt,vt,s,k,bd);
  }
}

void dcoef_bd_tz(double complex *CC,double complex *dCC,double *rt,double *vtz,int t,double zeta_t,double eta_t,int s,double complex k,int type,BOUD *bd)
{

  if(ELT4==check_element_type(s,bd)){
    if(t!=s){
      switch(type){
      case 0: bil_dcoef_4p(CC,dCC,rt,vtz,s,k,bd);
      break;
      case 1: bil_dcoef_9p(CC,dCC,rt,vtz,s,k,bd);
      break;
      case 2: bil_dcoef_GL(CC,dCC,rt,vtz,s,k,bd);
      break;
      case 3: bil_dcoef_GH(CC,dCC,rt,vtz,s,k,bd);
      break;
      case 4: bil_dcoef_DE(CC,dCC,rt,vtz,s,k,bd);
      break;
      default:
        printf("b_elem.c, dcoef_bd_tz(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }
      if(type<3 && fabs(creal(dCC[8]))> CBD_SDF) bil_dcoef_GH(CC,dCC,rt,vtz,s,k,bd);
    }
    else {
      if(SW_DBIEQ==0) bil_dcoef_bd_tz_PV(CC,dCC,zeta_t,eta_t,s,k,bd);
      else             bil_dcoef_bd_tz_DV(CC,dCC,zeta_t,eta_t,s,k,bd);
    }
  }
  else {
    if(t!=s){
      switch(type){
      case 0: lit_dcoef_4p(CC,dCC,rt,vtz,s,k,bd);
      break;
      case 1: lit_dcoef_7p(CC,dCC,rt,vtz,s,k,bd);
      break;
      case 2: lit_dcoef_GL(CC,dCC,rt,vtz,s,k,bd);
      break;
      case 3: lit_dcoef_GH(CC,dCC,rt,vtz,s,k,bd);
      break;
      case 4: lit_dcoef_DE(CC,dCC,rt,vtz,s,k,bd);
      break;
      default:
        printf("b_elem.c, dcoef_bd_tz(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }
      if(type<3 && fabs(creal(dCC[8]))>CBD_SDF) lit_dcoef_GH(CC,dCC,rt,vtz,s,k,bd);
    }
    else {
      if(SW_DBIEQ==0) lit_dcoef_bd_tz_PV(CC,dCC,zeta_t,eta_t,s,k,bd);
      else             lit_dcoef_bd_tz_DV(CC,dCC,zeta_t,eta_t,s,k,bd);
    }
  }
}

void dcoef_bd_te(double complex *CC,double complex *dCC,double *rt,double *vte,int t,double zeta_t,double eta_t,int s,double complex k,int type,BOUD *bd)
{
  if(ELT4==check_element_type(s,bd)){
    if(t!=s){
      switch(type){
      case 0: bil_dcoef_4p(CC,dCC,rt,vte,s,k,bd);
      break;
      case 1: bil_dcoef_9p(CC,dCC,rt,vte,s,k,bd);
      break;
      case 2: bil_dcoef_GL(CC,dCC,rt,vte,s,k,bd);
      break;
      case 3: bil_dcoef_GH(CC,dCC,rt,vte,s,k,bd);
      break;
      case 4: bil_dcoef_DE(CC,dCC,rt,vte,s,k,bd);
      break;
      default:
        printf("b_elem.c, dcoef_bd_te(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }
      if(type<3 && fabs(creal(dCC[8]))> CBD_SDF) bil_dcoef_GH(CC,dCC,rt,vte,s,k,bd);
    }
    else {
      if(SW_DBIEQ==0) bil_dcoef_bd_te_PV(CC,dCC,zeta_t,eta_t,s,k,bd);
      else             bil_dcoef_bd_te_DV(CC,dCC,zeta_t,eta_t,s,k,bd);
    }
  }
  else {
    if(t!=s){
      switch(type){
      case 0: lit_dcoef_4p(CC,dCC,rt,vte,s,k,bd);
      break;
      case 1: lit_dcoef_7p(CC,dCC,rt,vte,s,k,bd);
      break;
      case 2: lit_dcoef_GL(CC,dCC,rt,vte,s,k,bd);
      break;
      case 3: lit_dcoef_GH(CC,dCC,rt,vte,s,k,bd);
      break;
      case 4: lit_dcoef_DE(CC,dCC,rt,vte,s,k,bd);
      break;
      default:
        printf("b_elem.c, dcoef_bd_te(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }
      if(type<3 && fabs(creal(dCC[8]))>CBD_SDF) lit_dcoef_GH(CC,dCC,rt,vte,s,k,bd);
    }
    else {
      if(SW_DBIEQ==0) lit_dcoef_bd_te_PV(CC,dCC,zeta_t,eta_t,s,k,bd);
      else             lit_dcoef_bd_te_DV(CC,dCC,zeta_t,eta_t,s,k,bd);
    }
  }
}

void dcoef_bd_t2(double complex *CC,double complex *dCdtz,double complex *dCdte,
                  double *rt,double *vtz,double *vte,int t,double zeta_t,double eta_t,int s,double complex k,int type, BOUD *bd)
{
  if(ELT4==check_element_type(s,bd)){
    if(t!=s){
      switch(type){
      case 0: bil_dcoef_v2_4p(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);
      break;
      case 1: bil_dcoef_v2_9p(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);
      break;
      case 2: bil_dcoef_v2_GL(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);
      break;
      case 3: bil_dcoef_v2_GH(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);
      break;
      case 4: bil_dcoef_v2_DE(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);
      break;
      default:
        printf("b_elem.c, dcoef_bd_t2(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }

      if(type<3 && ( fabs(creal(dCdtz[8]))> CBD_SDF || fabs(creal(dCdte[8]))> CBD_SDF ))
        bil_dcoef_v2_GH(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);

    }
    else {
      if(SW_DBIEQ==0) bil_dcoef_bd_vt2_PV(CC,dCdtz,dCdte,zeta_t,eta_t,s,k,bd);
      else             bil_dcoef_bd_vt2_DV(CC,dCdtz,dCdte,zeta_t,eta_t,s,k,bd);
    }
  }
  else {
    if(t!=s){
      switch(type){
      case 0: lit_dcoef_v2_4p(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);
      break;
      case 1: lit_dcoef_v2_7p(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);
      break;
      case 2: lit_dcoef_v2_GL(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);
      break;
      case 3: lit_dcoef_v2_GH(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);
      break;
      case 4: lit_dcoef_v2_DE(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);
      break;
      default:
        printf("b_elem.c, dcoef_bd_t2(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }

      if(type<3 && ( fabs(creal(dCdtz[8]))>CBD_SDF || fabs(creal(dCdte[8]))>CBD_SDF ))
        lit_dcoef_v2_GH(CC,dCdtz,dCdte,rt,vtz,vte,s,k,bd);

    }
    else {
      if(SW_DBIEQ==0) lit_dcoef_bd_vt2_PV(CC,dCdtz,dCdte,zeta_t,eta_t,s,k,bd);
      else             lit_dcoef_bd_vt2_DV(CC,dCdtz,dCdte,zeta_t,eta_t,s,k,bd);
    }
  }
}

void dcoef_bd_node_t2(double complex *CC,double complex *dCdtz,double complex *dCdte,
                  int t,int nid,double *vtz,double *vte,int s,double complex k,int type,BOUD *bd)
{
  if(ELT4==check_element_type(s,bd)){
    if(t!=s){
      switch(type){
      case 0: bil_dcoef_v2_4p(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);
      break;
      case 1: bil_dcoef_v2_9p(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);
      break;
      case 2: bil_dcoef_v2_GL(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);
      break;
      case 3: bil_dcoef_v2_GH(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);
      break;
      case 4: bil_dcoef_v2_DE(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);
      break;
      default:
        printf("b_elem.c, dcoef_bd_t2(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }

      if(type<3 && ( fabs(creal(dCdtz[8]))> CBD_SDF || fabs(creal(dCdte[8]))> CBD_SDF ))
        bil_dcoef_v2_GH(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);

    }
    else {
      if(SW_DBIEQ==0) bil_dcoef_bd_vt2_PV(CC,dCdtz,dCdte,bd->zt_44[nid],bd->et_44[nid],s,k,bd);
      else             bil_dcoef_bd_vt2_DV(CC,dCdtz,dCdte,bd->zt_44[nid],bd->et_44[nid],s,k,bd);
    }
  }
  else {
    if(t!=s){
      switch(type){
      case 0: lit_dcoef_v2_4p(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);
      break;
      case 1: lit_dcoef_v2_7p(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);
      break;
      case 2: lit_dcoef_v2_GL(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);
      break;
      case 3: lit_dcoef_v2_GH(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);
      break;
      case 4: lit_dcoef_v2_DE(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);
      break;
      default:
        printf("b_elem.c, dcoef_bd_t2(), type number error! type=%d is not defined.Exit...\n",type);
        exit(1);
      }

      if(type<3 && ( fabs(creal(dCdtz[8]))>CBD_SDF || fabs(creal(dCdte[8]))>CBD_SDF ))
        lit_dcoef_v2_GH(CC,dCdtz,dCdte,bd->ren[abs(t)][nid],vtz,vte,s,k,bd);

    }
    else {
      if(SW_DBIEQ==0) lit_dcoef_bd_vt2_PV(CC,dCdtz,dCdte,bd->zt_34[nid],bd->et_34[nid],s,k,bd);
      else             lit_dcoef_bd_vt2_DV(CC,dCdtz,dCdte,bd->zt_34[nid],bd->et_34[nid],s,k,bd);
    }
  }
}


// ---- bi-linear element ---------------------------------------------------------------------------------
double bil_Mn(int n,double zeta,double eta)
{
  switch (n){
  case 0:  return 0.25*(1.0-zeta)*(1.0-eta);
  case 1:  return 0.25*(1.0+zeta)*(1.0-eta);
  case 2:  return 0.25*(1.0+zeta)*(1.0+eta);
  case 3:  return 0.25*(1.0-zeta)*(1.0+eta);
  default:
    printf("b_elem.c, bli_Mn(), shape function mode number error! n=%d. Exit...\n",n);
    exit(1);
  }
}

double bil_Nn(int n,double zeta,double eta)
{
  double i_a=1.0/P44_N;
  switch (n){
  case 0: return 0.25*(1.0-zeta*i_a)*(1.0-eta*i_a);
  case 1: return 0.25*(1.0+zeta*i_a)*(1.0-eta*i_a);
  case 2: return 0.25*(1.0+zeta*i_a)*(1.0+eta*i_a);
  case 3: return 0.25*(1.0-zeta*i_a)*(1.0+eta*i_a);
  default:
    printf("b_elem.c, bli_Nn(), interpolate function mode number error! n=%d. Exit...\n",n);
    exit(1);
  }
}

void bil_copy_elem_const_rw(double cr[][4],double cw[][3],int s,BOUD *bd)
{
  int l,m;

  for(l=0;l<3;l++){
    for(m=0;m<4;m++) cr[l][m]=bd->cr[s][l][m];
    for(m=0;m<3;m++) cw[l][m]=bd->cw[s][l][m];
  }
}

void bil_convert_CC(double complex *CC)
{
  // coef H
  CC[4]*=-1.0;
  CC[5]*=-1.0;
  CC[6]*=-1.0;
  CC[7]*=-1.0;
  // coef F
  CC[8]*=-1.0;
}

void bil_convert_dCC(double complex *dCC)
{
  // coef dH/dv
  dCC[4]*=-1.0;
  dCC[5]*=-1.0;
  dCC[6]*=-1.0;
  dCC[7]*=-1.0;
  // coef dF/dv
  dCC[8]*=-1.0;
}

void bil_r_zeta_eta(double *r,double zeta,double eta,double cr[][4])
{
  int i;

  for(i=0;i<3;i++) r[i]=cr[i][0]+cr[i][1]*zeta+cr[i][2]*eta+cr[i][3]*zeta*eta;
}

void bil_w_zeta_eta(double *w,double zeta,double eta,double cw[][3])
{
  int i;

  for(i=0;i<3;i++) w[i]=cw[i][0]+cw[i][1]*zeta+cw[i][2]*eta;
}

void bil_rw_zeta_eta(double *r,double *w,double zeta,double eta,double cr[][4],double cw[][3])
{
  int i;

  for(i=0;i<3;i++){
    r[i]=cr[i][0]+cr[i][1]*zeta+cr[i][2]*eta+cr[i][3]*zeta*eta;
    w[i]=cw[i][0]+cw[i][1]*zeta+cw[i][2]*eta;
  }
}

void bil_t_zeta(double *t,double zeta,double eta,double cr[][4])
{
  int i;

  for(i=0;i<3;i++) t[i]=cr[i][1]+cr[i][3]*eta;
}

void bil_t_eta(double *t,double zeta,double eta,double cr[][4])
{
  int i;

  for(i=0;i<3;i++) t[i]=cr[i][2]+cr[i][3]*zeta;
}

int bil_check_plane(double cr[][4],double cw[][3])
{
  double bdc,ab,aC;

  bdc=cr[0][1]*cw[0][2]+cr[1][1]*cw[1][2]+cr[2][1]*cw[2][2];
  ab=sqrt(cr[0][1]*cr[0][1]+cr[1][1]*cr[1][1]+cr[2][1]*cr[2][1]);
  aC=sqrt(cw[0][2]*cw[0][2]+cw[1][2]*cw[1][2]+cw[2][2]*cw[2][2]);
  if(fabs(bdc/(ab*aC))<GEPS) return 1;
  else return 0;
}

int bil_check_on_plane(double *rt,double cr[][4],double cw[][3])
{
  double aW,r0[3],ar0,dp;
  int l;

  if(bil_check_plane(cr,cw)==0) return 0;
  else { // element is plane
    aW=sqrt(cw[0][0]*cw[0][0]+cw[1][0]*cw[1][0]+cw[2][0]*cw[2][0]);
    for(l=0;l<3;l++) r0[l]=cr[l][0]-rt[l];
    ar0=sqrt(r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2]);
    dp=(r0[0]*cw[0][0]+r0[1]*cw[1][0]+r0[2]*cw[2][0])/(aW*ar0);
    if(fabs(dp)<GEPS) return 1;
    else return 0;
  }
}

void bil_calc_node_r_th(double *r,double *th,double zeta_t,double eta_t)
{
  double ztp,ztm,etp,etm;

  ztp= 1.0-zeta_t;    ztm=-1.0-zeta_t;
  etp= 1.0-eta_t;      etm=-1.0-eta_t;

  r[0]=sqrt(ztm*ztm+etm*etm);    r[1]=sqrt(ztp*ztp+etm*etm);
  r[2]=sqrt(ztp*ztp+etp*etp);    r[3]=sqrt(ztm*ztm+etp*etp);
  r[4]=r[0];

  th[0]=atan2(etm,ztm);    th[1]=atan2(etm,ztp);
  th[2]=atan2(etp,ztp);    th[3]=atan2(etp,ztm);
  th[4]=th[0]+2.0*M_PI;
}





// ---- linear triangular element ------------------------------------------------------------------
double lit_Mn(int n,double zeta,double eta)
{
  switch(n){
  case 0: return 1.0/3.0*(1.0-zeta-sqrt(3.0)*eta);
  case 1: return 1.0/3.0*(1.0+2.0*zeta);
  case 2: return 1.0/3.0*(1.0-zeta+sqrt(3.0)*eta);
  default:
    printf("b_elem.c, lit_Mn(), shape function mode number error! n=%d. Exit...\n",n);
    exit(1);
  }
}

double lit_Nn(int n,double zeta,double eta)
{
  double i_a=0.5/P34_N0;
  switch(n){
  case 0: return 1.0/3.0*(1.0-i_a*(zeta+sqrt(3.0)*eta));
  case 1: return 1.0/3.0*(1.0+2.0*i_a*zeta);
  case 2: return 1.0/3.0*(1.0-i_a*(zeta-sqrt(3.0)*eta));
  default:
    printf("b_elem.c, lit_Mn(), interpolate function mode number error! n=%d. Exit...\n",n);
    exit(1);
  }
}

void lit_copy_elem_const_rw(double cr[][4],double cw[][3],int s,BOUD *bd)
{
  int l;
  for(l=0;l<3;l++){
    cr[l][0]=bd->cr[s][l][0];
    cr[l][1]=bd->cr[s][l][1];
    cr[l][2]=bd->cr[s][l][2];
    cr[l][3]=0.0;

    cw[l][0]=bd->cw[s][l][0];
    cw[l][1]=0.0;
    cw[l][2]=0.0;
  }
}

void lit_convert_CC(double complex *CC)
{
  // H
  CC[0+4]*=-1.0;
  CC[1+4]*=-1.0;
  CC[2+4]*=-1.0;
  // F
  CC[8]*=-1.0;
}

void lit_convert_dCC(double complex *dCC)
{
  // dH/dv
  dCC[0+4]*=-1.0;
  dCC[1+4]*=-1.0;
  dCC[2+4]*=-1.0;
  // dF/dv
  dCC[8]*=-1.0;
}

void lit_r_zeta_eta(double *r,double zeta,double eta,double cr[][4])
{
  int i;

  for(i=0;i<3;i++) r[i]=cr[i][0]+cr[i][1]*zeta+cr[i][2]*eta;
}

void lit_w_zeta_eta(double *w,double zeta,double eta,double cw[][3])
{
  int i;

  for(i=0;i<3;i++) w[i]=cw[i][0];
}

void lit_rw_zeta_eta(double *r,double *w,double zeta,double eta,double cr[][4],double cw[][3])
{
  int i;

  for(i=0;i<3;i++){
    r[i]=cr[i][0]+cr[i][1]*zeta+cr[i][2]*eta;
    w[i]=cw[i][0];
  }
}

void lit_t_zeta(double *t,double zeta,double eta,double cr[][4])
{
  int i;

  for(i=0;i<3;i++) t[i]=cr[i][1];
}

void lit_t_eta (double *t,double zeta,double eta,double cr[][4])
{
  int i;

  for(i=0;i<3;i++) t[i]=cr[i][2];
}

int lit_check_on_plane(double *rt,double cr[][4],double cw[][3])
{
  double vW[3],aW,vR[3],aR,dp;
  int l;

  for(l=0;l<3;l++){
    vW[l]=cw[l][0];
    vR[l]=cr[l][0]-rt[l];
  }
  aW=vabs_d(vW);
  aR=vabs_d(vR);
  dp=(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])/(aW*aR);
  if(fabs(dp)<GEPS) return 1;
  else return 0;
}

void lit_calc_node_r_th(double *r,double *th,double zeta_t,double eta_t)
{
   double ztp,ztm,etp,etm;

   ztp=1.0-zeta_t;    ztm=-0.5-zeta_t;
   etp= 0.5*SQ3-eta_t;    etm=-0.5*SQ3-eta_t;
   r[0]=sqrt(ztm*ztm+etm*etm);    r[1]=sqrt(ztp*ztp+eta_t*eta_t);
   r[2]=sqrt(ztm*ztm+etp*etp);
   r[3]=r[0];

   th[0]=atan2(etm,ztm);    th[1]=atan2(-eta_t,ztp);
   th[2]=atan2(etp,ztm);
   th[3]=th[0]+2.0*M_PI;
}

