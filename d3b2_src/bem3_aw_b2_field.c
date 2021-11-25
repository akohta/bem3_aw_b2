#include "bem3_aw_b2.h"


void mo_object_domain_id(int *oid,int *did,double *rt,MOBJ *mo)
{
  int domain_id_m_dmda(double *rt,DMDA *ad);
  
  int o;

  for(o=0;o<mo->N;o++){
    *did=domain_id_m_dmda(rt,&(mo->md[o]));
    if(*did!=0){
      *oid=o;
      return;
    }
  }

  *oid=-1; // open region
}

int mo_phi_s(double complex *phi,double *rt,int type,MOBJ *mo)
{
  int phi_s_dmda(double complex *phi,double *rt,int type,DMDA *ad); //bem3_aw_b1_field.c
  
  double complex tp;
  int oid,did,m;

  mo_object_domain_id(&oid,&did,rt,mo);
  if(oid<0){ // open region
    *phi=0.0;
    for(m=0;m<mo->N;m++){
      phi_s_dmda(&tp,rt,type,&(mo->md[m]));
      *phi+=tp;
    }
  }
  else{
    phi_s_dmda(phi,rt,type,&(mo->md[oid]));
  }

  return oid;
}

int mo_phi_t(double complex *phi,double *rt,int type,MOBJ *mo)
{
  int phi_s_dmda(double complex *phi,double *rt,int type,DMDA *ad); //bem3_aw_b1_field.c
  int phi_t_dmda(double complex *phi,double *rt,int type,DMDA *ad); //bem3_aw_b1_field.c
  
  double complex tp,p,vp[3];
  int oid,did,m;

  mo_object_domain_id(&oid,&did,rt,mo);
  if(oid<0){ // open region
    calc_maw_pv(&p,vp,rt,&(mo->md[0].aw));
    *phi=-p/mo->md[0].aw.k2;
    for(m=0;m<mo->N;m++){
      phi_s_dmda(&tp,rt,type,&(mo->md[m]));
      *phi+=tp;
    }
  }
  else {
    phi_t_dmda(phi,rt,type,&(mo->md[oid]));
  }

  return oid;
}

int mo_phi_i(double complex *phi,double *rt,int type,MOBJ *mo)
{
  double complex p,vp[3];
  int oid,did;

  mo_object_domain_id(&oid,&did,rt,mo);
  calc_maw_pv(&p,vp,rt,&(mo->md[0].aw));
  *phi=-p/mo->md[0].aw.k2;

  return oid;
}

int mo_p_s(double complex *p,double *rt,int type,MOBJ *mo)
{
  int p_s_dmda(double complex *p,double *rt,int type,DMDA *ad); // bem3_aw_b1_field.c
  
  double complex tp;
  int oid,did,m;

  mo_object_domain_id(&oid,&did,rt,mo);
  if(oid<0){ // open region
    *p=0.0;
    for(m=0;m<mo->N;m++){
      p_s_dmda(&tp,rt,type,&(mo->md[m]));
      *p+=tp;
    }
  }
  else p_s_dmda(p,rt,type,&(mo->md[oid]));

  return oid;
}

int mo_p_t(double complex *p,double *rt,int type,MOBJ *mo)
{
  int p_s_dmda(double complex *p,double *rt,int type,DMDA *ad); // bem3_aw_b1_field.c
  int p_t_dmda(double complex *p,double *rt,int type,DMDA *ad); // bem3_aw_b1_field.c
  
  double complex tp,vp[3];
  int oid,did,m;

  mo_object_domain_id(&oid,&did,rt,mo);
  if(oid<0){ // open region
    calc_maw_pv(p,vp,rt,&(mo->md[0].aw));
    for(m=0;m<mo->N;m++){
      p_s_dmda(&tp,rt,type,&(mo->md[m]));
      *p+=tp;
    }
  }
  else p_t_dmda(p,rt,type,&(mo->md[oid]));

  return oid;  
}

int mo_p_i(double complex *p,double *rt,int type,MOBJ *mo)
{
  double complex vp[3];
  int oid,did;

  mo_object_domain_id(&oid,&did,rt,mo);
  calc_maw_pv(p,vp,rt,&(mo->md[0].aw));

  return oid;
}

int mo_pv_s(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo)
{
  int pv_s_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // bem3_aw_b1_field.c
  
  double complex tp,tv[3];
  int oid,did,m;

  mo_object_domain_id(&oid,&did,rt,mo);
  if(oid<0){ // open region
    *p=0.0;
    pv[0]=0.0;    pv[1]=0.0;    pv[2]=0.0;
    for(m=0;m<mo->N;m++){
      pv_s_dmda(&tp,tv,rt,type,&(mo->md[m]));
      *p+=tp;
      pv[0]+=tv[0];      pv[1]+=tv[1];      pv[2]+=tv[2];
    }
  }
  else pv_s_dmda(p,pv,rt,type,&(mo->md[oid]));

  return oid;
}

int mo_pv_t(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo)
{
  int pv_s_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // bem3_aw_b1_field.c
  int pv_t_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // bem3_aw_b1_field.c
  
  double complex tp,tv[3];
  int oid,did,m;

  mo_object_domain_id(&oid,&did,rt,mo);
  if(oid<0){ // open region
    calc_maw_pv(p,pv,rt,&(mo->md[0].aw));
    for(m=0;m<mo->N;m++){
      pv_s_dmda(&tp,tv,rt,type,&(mo->md[m]));
      *p+=tp;
      pv[0]+=tv[0];      pv[1]+=tv[1];      pv[2]+=tv[2];
    }
  }
  else pv_t_dmda(p,pv,rt,type,&(mo->md[oid]));

  return oid;
}

int mo_pv_i(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo)
{
  int oid,did;

  mo_object_domain_id(&oid,&did,rt,mo);
  calc_maw_pv(p,pv,rt,&(mo->md[0].aw));

  return oid;
}

int mo_pv_s_dbieq(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo)
{
  int pv_s_dbieq_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // bem3_aw_b1_field.c
  
  double complex tp,tv[3];
  int oid,did,m,i;

  mo_object_domain_id(&oid,&did,rt,mo);
  if(oid<0){ // open region
    *p=0.0;
    for(i=0;i<3;i++) pv[i]=0.0;
    for(m=0;m<mo->N;m++){
      pv_s_dbieq_dmda(&tp,tv,rt,type,&(mo->md[m]));
      *p+=tp;
      for(i=0;i<3;i++) pv[i]+=tv[i];
    }
  }
  else pv_s_dbieq_dmda(p,pv,rt,type,&(mo->md[oid]));

  return oid;
}

int mo_pv_t_dbieq(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo)
{
  int pv_s_dbieq_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // bem3_aw_b1_field.c
  int pv_t_dbieq_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // bem3_aw_b1_field.c
  
  double complex tp,tv[3];
  int oid,did,m;

  mo_object_domain_id(&oid,&did,rt,mo);
  if(oid<0){ // open region
    calc_maw_pv(p,pv,rt,&(mo->md[0].aw));
    for(m=0;m<mo->N;m++){
      pv_s_dbieq_dmda(&tp,tv,rt,type,&(mo->md[m]));
      *p+=tp;
      pv[0]+=tv[0];      pv[1]+=tv[1];      pv[2]+=tv[2];
    }
  }
  else pv_t_dbieq_dmda(p,pv,rt,type,&(mo->md[oid]));

  return oid;
}  

int mo_pv_i_dbieq(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo)
{
  int oid,did;

  mo_object_domain_id(&oid,&did,rt,mo);
  calc_maw_pv(p,pv,rt,&(mo->md[0].aw));

  return oid;  
}

void mo_pv_s_bd(double complex *p,double complex *pv,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo)
{
  void pv_s_bd_dmda(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad);
  
  pv_s_bd_dmda(p,pv,did,t,zeta_t,eta_t,type,&(mo->md[oid]));
}

void mo_pv_t_bd(double complex *p,double complex *pv,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo)
{
  void pv_s_bd_dmda(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad);
  int pv_s_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); 
  
  double complex tp,tv[3];
  double rt[3];
  int m;

  pv_s_bd_dmda(p,pv,did,t,zeta_t,eta_t,type,&(mo->md[oid]));
  if(did==0){
    r_bd(rt,mo->md[oid].bd.sb[did].sid[t],zeta_t,eta_t,&(mo->md[oid].bd));
    for(m=0;m<mo->N;m++){
      if(m==oid) continue;
      pv_s_dmda(&tp,tv,rt,type,&(mo->md[m]));
      *p+=tp;
      pv[0]+=tv[0];
      pv[1]+=tv[1];
      pv[2]+=tv[2];
    }
    calc_maw_pv(&tp,tv,rt,&(mo->md[0].aw));
    *p+=tp;
    pv[0]+=tv[0];
    pv[1]+=tv[1];
    pv[2]+=tv[2];
  }
}

void mo_pv_i_bd(double complex *p,double complex *pv,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo)
{
  int pv_s_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); 
  
  double complex tp,tv[3];
  double rt[3];
  int m;
  
  *p=0.0;
  pv[0]=0.0;
  pv[1]=0.0;
  pv[2]=0.0;
  if(did==0){
    r_bd(rt,mo->md[oid].bd.sb[did].sid[t],zeta_t,eta_t,&(mo->md[oid].bd));
    for(m=0;m<mo->N;m++){
      if(m==oid) continue;
      pv_s_dmda(&tp,tv,rt,type,&(mo->md[m]));
      *p+=tp;
      pv[0]+=tv[0];
      pv[1]+=tv[1];
      pv[2]+=tv[2];
    }
    calc_maw_pv(&tp,tv,rt,&(mo->md[0].aw));
    *p+=tp;
    pv[0]+=tv[0];
    pv[1]+=tv[1];
    pv[2]+=tv[2];
  }   
}
