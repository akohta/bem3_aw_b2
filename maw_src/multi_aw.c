#include "multi_aw.h"

void init_maw(Maw *aw)
{
  aw->n_pw=0;  sprintf(aw->fname_pw,"%s","null");
  aw->bd.pw=(Apw *)m_alloc2(0,sizeof(Apw),"multi_aw.c,init_maw(),aw->bd.pw");
  aw->n_bb=0;  sprintf(aw->fname_bb,"%s","null");
  aw->bd.bb=(Abb *)m_alloc2(0,sizeof(Abb),"multi_aw.c,init_maw(),aw->bd.bb");
  aw->n_fb=0;  sprintf(aw->fname_fb,"%s","null");
  aw->bd.fb=(Afb *)m_alloc2(0,sizeof(Afb),"multi_aw.c,init_maw(),aw->bd.fb");
}

void read_data_maw(Maw *aw)
{
  int wn=0;

  wn+=read_data_maw_pw(fn_pw,aw); 
  wn+=read_data_maw_bb(fn_bb,aw); 
  wn+=read_data_maw_fb(fn_fb,aw); 
  
  if(wn==0){   
    printf("read_data_maw() No beam defined, check beam data file : %s, %s, %s. Exit...\n",
                fn_pw,fn_bb,fn_fb);
    exit(1);
  }
}

void print_data_maw(Maw *aw)
{
  int i;

  for(i=0;i<aw->n_pw;i++){
    printf(" \"%s\" No.%02d ",aw->fname_pw,i);    print_data_apw(&(aw->bd.pw[i]));
  }
  for(i=0;i<aw->n_bb;i++){
    printf(" \"%s\" No.%02d ",aw->fname_bb,i);    print_data_abb(&(aw->bd.bb[i]));
  }
  for(i=0;i<aw->n_fb;i++){
    printf(" \"%s\" No.%02d ",aw->fname_fb,i);    print_data_afb(&(aw->bd.fb[i]));
  }
  

}

void setup_maw(Maw *aw)
{
  void check_data_maw(Maw *aw);
  
  int i;
  
  check_data_maw(aw);
  // pw
  for(i=0;i<aw->n_pw;i++)      setup_apw(&(aw->bd.pw[i])); 
  // bb
  for(i=0;i<aw->n_bb;i++)      setup_abb(&(aw->bd.bb[i]));
  // fb
  for(i=0;i<aw->n_fb;i++)      setup_afb(&(aw->bd.fb[i]));
  
}

void free_maw(Maw *aw)
{
  int i;
  
  // pw
  free(aw->bd.pw);  aw->n_pw=0; 
  // bb
  free(aw->bd.bb);  aw->n_bb=0;
  // fb
  for(i=0;i<aw->n_fb;i++)  free_afb(&(aw->bd.fb[i]));
  free(aw->bd.fb);  aw->n_fb=0;

}

int read_data_maw_pw(char *fname,Maw *aw)
{
  FILE *fp;
  char buf[256]="";    int i,nn;  double td1,td2,rho0,c0,f;
   
  if((fp=fopen(fname,"rt"))==NULL) return 0;
  sprintf(aw->fname_pw,"%s",fname);   
  
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_pw(), failed to read the line. exit...\n");
    exit(1);  
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_pw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&rho0)!=1){
    printf("multi_aw.c, read_data_maw_pw(), failed to read the rho0. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&c0)!=1){
    printf("multi_aw.c, read_data_maw_pw(), failed to read the c0. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf\n",&f)!=1){
    printf("multi_aw.c, read_data_maw_pw(), failed to read the f. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_pw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%d\n",&nn)!=1){
    printf("multi_aw.c, read_data_maw_pw(), failed to read the nn. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_pw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(nn==0) return 0;

  free(aw->bd.pw);
  aw->n_pw=nn;     aw->bd.pw=(Apw *)m_alloc2(nn,sizeof(Apw),"multi_aw.c,read_data_maw_pw(),aw->bd.pw");

  for(i=0;i<nn;i++){
    aw->bd.pw[i].rho0=rho0;
    aw->bd.pw[i].c0=c0;
    aw->bd.pw[i].f=f;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_pw(), failed to read the real(p). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&td2)!=1){
      printf("multi_aw.c, read_data_maw_pw(), failed to read the imag(p). exit...\n");
      exit(1);
    }
    aw->bd.pw[i].p=td1+I*td2;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_pw(), failed to read the tv[0]. exit...\n");
      exit(1);
    }
    aw->bd.pw[i].tv[0]=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_pw(), failed to read the tv[1]. exit...\n");
      exit(1);
    }
    aw->bd.pw[i].tv[1]=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_pw(), failed to read the tv[2]. exit...\n");
      exit(1);
    }
    aw->bd.pw[i].tv[2]=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_pw(), failed to read the theta. exit...\n");
      exit(1); 
    }
    aw->bd.pw[i].theta=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_pw(), failed to read the psi. exit...\n");
      exit(1);
    }
    aw->bd.pw[i].psi  =td1;
  }
  fclose(fp);
  return nn;
}

int read_data_maw_bb(char *fname,Maw *aw)
{
  FILE *fp;
  char buf[256]="";    int i,nn;  double td1,td2,rho0,c0,f;
   
  if((fp=fopen(fname,"rt"))==NULL) return 0;
  sprintf(aw->fname_bb,"%s",fname);   
  
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_bb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_bb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&rho0)!=1){
    printf("multi_aw.c, read_data_maw_bb(), failed to read the rho0. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&c0)!=1){
    printf("multi_aw.c, read_data_maw_bb(), failed to read the c0. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf\n",&f)!=1){
    printf("multi_aw.c, read_data_maw_bb(), failed to read the f. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_bb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%d\n",&nn)!=1){
    printf("multi_aw.c, read_data_maw_bb(), failed to read the nn. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_bb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(nn==0) return 0;

  free(aw->bd.bb);
  aw->n_bb=nn;     aw->bd.bb=(Abb *)m_alloc2(nn,sizeof(Abb),"multi_aw.c,read_data_maw_bb(),aw->bd.bb");
  
  for(i=0;i<nn;i++){
    aw->bd.bb[i].rho0=rho0;
    aw->bd.bb[i].c0  =c0;
    aw->bd.bb[i].f   =f;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_bb(), failed to read the real(p). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&td2)!=1){
      printf("multi_aw.c, read_data_maw_bb(), failed to read the imag(p). exit...\n");
      exit(1);
    }
    aw->bd.bb[i].p=td1+I*td2;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_bb(), failed to read the d_angle. exit...\n");
      exit(1);
    }
    aw->bd.bb[i].d_angle=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_bb(), failed to read the tv[0]. exit...\n");
      exit(1);
    }
    aw->bd.bb[i].tv[0]=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_bb(), failed to read the tv[1]. exit...\n");
      exit(1);
    }
    aw->bd.bb[i].tv[1]=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_bb(), failed to read the tv[2]. exit...\n");
      exit(1);
    }
    aw->bd.bb[i].tv[2]=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_bb(), failed to read the theta. exit...\n");
      exit(1);
    }
    aw->bd.bb[i].theta=td1;
    if(fscanf(fp,"%lf\n",&td1)!=1){
      printf("multi_aw.c, read_data_maw_bb(), failed to read the psi. exit...\n");
      exit(1);
    }
    aw->bd.bb[i].psi  =td1;
  }
  fclose(fp);
  return nn;
}

int read_data_maw_fb(char *fname,Maw *aw)
{
  FILE *fp;
  char buf[256]="";    int i,nn,ti;  double td1,td2,rho0,c0,f;
   
  if((fp=fopen(fname,"rt"))==NULL) return 0;
  sprintf(aw->fname_fb,"%s",fname);   
  
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_fb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_fb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&rho0)!=1){
    printf("multi_aw.c, read_data_maw_fb(), failed to read the rho0. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&c0)!=1){
    printf("multi_aw.c, read_data_maw_fb(), failed to read the c0. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf\n",&f)!=1){
    printf("multi_aw.c, read_data_maw_fb(), failed to read the f. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_fb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%d\n",&nn)!=1){
    printf("multi_aw.c, read_data_maw_fb(), failed to read the nn. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_aw.c, read_data_maw_fb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(nn==0) return 0;

  free(aw->bd.fb);
  aw->n_fb=nn;     aw->bd.fb=(Afb *)m_alloc2(nn,sizeof(Afb),"multi_aw.c,read_data_maw_fb(),aw->bd.fb");
  
  for(i=0;i<nn;i++){
    aw->bd.fb[i].rho0=rho0;
    aw->bd.fb[i].c0  =c0;
    aw->bd.fb[i].f   =f;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_fb(), failed to read the real(p). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&td2)!=1){
      printf("multi_aw.c, read_data_maw_fb(), failed to read the imag(p). exit...\n");
      exit(1);
    }
    aw->bd.fb[i].p=td1+I*td2;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_fb(), failed to read the B. exit...\n");
      exit(1);
    }
    aw->bd.fb[i].B=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_fb(), failed to read the F. exit...\n");
      exit(1);
    }
    aw->bd.fb[i].F=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_fb(), failed to read the tv[0]. exit...\n");
      exit(1);
    }
    aw->bd.fb[i].tv[0]=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_fb(), failed to read the tv[1]. exit...\n");
      exit(1);  
    }
    aw->bd.fb[i].tv[1]=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_fb(), failed to read the tv[2]. exit...\n");
      exit(1);
    }
    aw->bd.fb[i].tv[2]=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_fb(), failed to read the theta. exit...\n");
      exit(1);
    }
    aw->bd.fb[i].theta=td1;
    if(fscanf(fp,"%lf",&td1)!=1){
      printf("multi_aw.c, read_data_maw_fb(), failed to read the psi. exit...\n");
      exit(1);
    }
    aw->bd.fb[i].psi  =td1;
    if(fscanf(fp,"%d\n",&ti)!=1){
      printf("multi_aw.c, read_data_maw_fb(), failed to read the nn. exit...\n");
      exit(1);
    }
    aw->bd.fb[i].nn=ti;
  }
  fclose(fp);
  return nn;
}

int calc_maw_pv(double complex *p,double complex *v,double *r,Maw *aw)
{
  double complex tp,tv[3];
  int i,j,ret,ti;
  
  *p=0.0;
  v[0]=0.0;  v[1]=0.0;  v[2]=0.0;
  ret=-1;
  // pw
  for(i=0;i<aw->n_pw;i++){
    calc_apw_pv(&tp,tv,r,&(aw->bd.pw[i]));
    *p+=tp;
    for(j=0;j<3;j++) v[j]+=tv[j];
  }
  // bb
  for(i=0;i<aw->n_bb;i++){
    calc_abb_pv(&tp,tv,r,&(aw->bd.bb[i]));
    *p+=tp;
    for(j=0;j<3;j++) v[j]+=tv[j];
  }
  // fb
  for(i=0;i<aw->n_fb;i++){
    ti=calc_afb_pv(&tp,tv,r,&(aw->bd.fb[i]));
    if(ti==0) ret=0; // outside the transducer
    else{ // inside the transducer 
      ret=i; 
      *p=0.0;
      v[0]=0.0;  v[1]=0.0;  v[2]=0.0;
      return ret;
    }
    *p+=tp;
    for(j=0;j<3;j++) v[j]+=tv[j];
  }
  
  return ret;
}

int calc_maw_dpdn(double complex *p,double complex *dpdn,double *r,double *n,Maw *aw)
{
  double complex tp,tdp;
  int i,ret,ti;
  
  *p=0.0;
  *dpdn=0.0;
  ret=-1;
  // pw
  for(i=0;i<aw->n_pw;i++){
    calc_apw_dpdn(&tp,&tdp,r,n,&(aw->bd.pw[i]));
    *p+=tp;
    *dpdn+=tdp;
  }
  // bb
  for(i=0;i<aw->n_bb;i++){
    calc_abb_dpdn(&tp,&tdp,r,n,&(aw->bd.bb[i]));
    *p+=tp;
    *dpdn+=tdp;
  }
  // fb
  for(i=0;i<aw->n_fb;i++){
    ti=calc_afb_dpdn(&tp,&tdp,r,n,&(aw->bd.fb[i]));
    if(ti==0) ret=0; // outside the transducer
    else{ // inside the transducer 
      ret=i; 
      *p=0.0;
      *dpdn=0.0;
      return ret;
    }
    *p+=tp;
    *dpdn+=tdp;
  }
  
  return ret; 
  
}

////////////////////////////////////////////////////////////////////
void check_data_maw(Maw *aw)
{
  int i;
  int st=0;
  
  aw->rho0=0.0;
  aw->c0=0.0;
  aw->f=0.0;
  
  // pw
  for(i=0;i<aw->n_pw;i++){
    if(st==0){
      aw->rho0=aw->bd.pw[i].rho0;
      aw->c0  =aw->bd.pw[i].c0;
      aw->f   =aw->bd.pw[i].f;
      st++;
    }
    else {
      if(aw->rho0!=aw->bd.pw[i].rho0){ printf("data error pw rho0. Exit\n"); exit(1);}
      if(aw->c0  !=aw->bd.pw[i].c0  ){ printf("data error pw c0. Exit\n"); exit(1);}
      if(aw->f   !=aw->bd.pw[i].f   ){ printf("data error pw f. Exit\n"); exit(1);}
    }
  }
  // bb
  for(i=0;i<aw->n_bb;i++){
    if(st==0){
      aw->rho0=aw->bd.bb[i].rho0;
      aw->c0  =aw->bd.bb[i].c0;
      aw->f   =aw->bd.bb[i].f;
      st++;
    }
    else {
      if(aw->rho0!=aw->bd.bb[i].rho0){ printf("data error bb rho0. Exit\n"); exit(1);}
      if(aw->c0  !=aw->bd.bb[i].c0  ){ printf("data error bb c0. Exit\n"); exit(1);}
      if(aw->f   !=aw->bd.bb[i].f   ){ printf("data error bb f. Exit\n"); exit(1);}
    }
  }
  // fb
  for(i=0;i<aw->n_fb;i++){
    if(st==0){
      aw->rho0=aw->bd.fb[i].rho0;
      aw->c0  =aw->bd.fb[i].c0;
      aw->f   =aw->bd.fb[i].f;
      st++;
    }
    else {
      if(aw->rho0!=aw->bd.fb[i].rho0){ printf("data error fb rho0. Exit\n"); exit(1);}
      if(aw->c0  !=aw->bd.fb[i].c0  ){ printf("data error fb c0. Exit\n"); exit(1);}
      if(aw->f   !=aw->bd.fb[i].f   ){ printf("data error fb f. Exit\n"); exit(1);}
    }
  }
  
  if(st==0){
    printf("no beams defined. check datafile name. Exit...\n");
    exit(0);
  }
  
  aw->lambda0=aw->c0/aw->f;
  aw->k0=2.0*M_PI/aw->lambda0;
  aw->k1=I*2.0*M_PI*aw->f/(aw->c0*aw->c0*aw->rho0);
  aw->k2=I*2.0*M_PI*aw->f*aw->rho0;
}



