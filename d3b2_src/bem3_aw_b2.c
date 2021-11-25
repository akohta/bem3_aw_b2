#include "bem3_aw_b2.h"


void read_dmda2(int argc,char **argv,DMDA *ad)
{
  void read_medium_data(char *med_fn,DMDA *ad); // bem3_aw_b1.c
  void read_mesh_data(char *msh_fn,DMDA *ad);   // bem3_aw_b1.c
  
  // check arguments
  if(argc!=11 && argc!=4){
    printf("This program needs command line arguments as follows.\n");
    printf("%s medium_datafile_name mesh_datafile_name output_datafile_name [rv_x rv_y rv_z theta tr_x tr_y tr_z](optional)\n",argv[0]);
    printf("rv : vector defining rotation axis, theta : rotation angle ( using Rodrigues' rotation formula ), tr : translation vector\n");
    exit(0);
  }
  
  read_data_maw(&(ad->aw)); // multi_aw.h
  read_medium_data(argv[1],ad);
  read_mesh_data(argv[2],ad);
  
  if(argc==11){
    ad->rv[0]=atof(argv[ 4]);
    ad->rv[1]=atof(argv[ 5]);
    ad->rv[2]=atof(argv[ 6]);
    ad->th   =atof(argv[ 7]);
    ad->tv[0]=atof(argv[ 8]);
    ad->tv[1]=atof(argv[ 9]);
    ad->tv[2]=atof(argv[10]);
  }
  else {
    ad->rv[0]=1.0;
    ad->rv[1]=0.0;
    ad->rv[2]=0.0;
    ad->th=0.0;
    ad->tv[0]=0.0;
    ad->tv[1]=0.0;
    ad->tv[2]=0.0;
  }
}

void print_dmda2(DMDA *ad)
{
  void print_medium_data(DMDA *ad); // bem3_aw_b1.c
  void print_mesh_data(DMDA *ad);   // bem3_aw_b1.c

  print_data_maw(&(ad->aw));
  print_medium_data(ad);
  print_mesh_data(ad);

  if(ad->th!=0.0 || vabs_d(ad->tv)!=0.0){
    printf("-- rotation and translation settings --\n");
    if(ad->th!=0.0){
      printf("vector defining rotation axis :(% 8.7g,% 8.7g,% 8.7g)\n",ad->rv[0],ad->rv[1],ad->rv[2]);
      printf("rotation angle           [rad]: %8.7g\n",ad->th);
    }
    if(vabs_d(ad->tv)!=0.0){
      printf("translation vector            :(%8.7g,%8.7g,%8.7g)\n",ad->tv[0],ad->tv[1],ad->tv[2]);
    }
  }
}

void initialize_dmda2(DMDA *ad)
{
  void rotation_translation_obj(double *rv,double th,double *tv,DMDA *ad); // bem3_aw_b1.c
  void init_elem_const(BOUD *bd);   // bem3_aw_b1.c
  void malloc_sub_domain(DMDA *ad); // bem3_aw_b1.c
  void init_sub_domain(DMDA *ad);   // bem3_aw_b1.c
  void init_boundary_data2(DMDA *ad);

  double omega;
  int i;

  // incident field
  setup_maw(&(ad->aw));
  // medium
  ad->rho0[0]=ad->aw.rho0;
  ad->c0[0]=ad->aw.c0;
  ad->k0[0]=ad->aw.k0;
  ad->K0[0]=ad->aw.rho0*ad->aw.c0*ad->aw.c0;
  ad->k1[0]=ad->aw.k1;
  ad->k2[0]=ad->aw.k2;
  for(i=1;i<=ad->MN;i++){
    omega=2.0*M_PI*ad->aw.f;
    ad->k0[i]=omega/ad->c0[i];
    ad->K0[i]=ad->rho0[i]*ad->c0[i]*ad->c0[i];
    ad->k1[i]=I*omega/ad->K0[i];
    ad->k2[i]=I*omega*ad->rho0[i];
  }
  // rotation and translation object
  if(ad->th!=0.0 || vabs_d(ad->tv)!=0.0) rotation_translation_obj(ad->rv,ad->th,ad->tv,ad);  
  // element constant
  init_elem_const(&(ad->bd));
  // sub domain
  malloc_sub_domain(ad);
  init_sub_domain(ad);
  // boundary data
  init_boundary_data2(ad);
}

void finalize_dmda2(DMDA *ad)
{
  void finalize_dmda(DMDA *ad); // bem3_aw_b1.c
  
  finalize_dmda(ad);
}
  
void mo_initialize(int argc,char **argv,MOBJ *mo)
{
  void mo_malloc(MOBJ *mo);
  void dat_read2(char *fname,DMDA *md);
  void rotation_translation_obj(double *rv,double th,double *tv,DMDA *md); // bem3_aw_b1.c
  void init_elem_const(BOUD *bd);   // bem3_aw_b1.c
  void init_boundary_data2(DMDA *md);
  void init_boundary_data_u2(DMDA *md);
  
  FILE *fp;
  int ti,i;
  char tmp[256];
  double rv[3];

  rv[0]=1.0;  rv[1]=0.0;  rv[2]=0.0;

  if(argc!=3){
    printf("This program needs command line argument as follows.\n");
    printf("%s model_setting_filename output_data_filename\n",argv[0]);
    exit(0);
  }

  if((fp=fopen(argv[1],"rt"))==NULL){    printf("bem3_aw_b2.c, mo_initialize(), Can not open the %s file.\n",argv[1]);    exit(1);  }

  fgets(tmp,256,fp);
  fgets(tmp,256,fp);

  fscanf(fp,"%d\n",&ti); // number of objects
  mo->N=ti;
  mo_malloc(mo); // malloc

  fgets(tmp,256,fp);
  for(i=0;i<mo->N;i++){
    fscanf(fp,"%s %lf %lf %lf",mo->mod_fn[i],&(mo->rs[i][0]),&(mo->rs[i][1]),&(mo->rs[i][2]));
    dat_read2(mo->mod_fn[i],&(mo->md[i]));
    rotation_translation_obj(rv,0.0,mo->rs[i],&(mo->md[i]));
    init_elem_const(&(mo->md[i].bd));
    init_boundary_data2(&(mo->md[i]));
    init_boundary_data_u2(&(mo->md[i]));
  }
  fclose(fp);
}

void mo_print_data(MOBJ *mo)
{
  void print_medium_data(DMDA *ad); // bem3_aw_b1.c
  void print_mesh_data(DMDA *ad);   // bem3_aw_b1.c

  int i;

  print_data_maw(&(mo->md[0].aw));

  for(i=0;i<mo->N;i++){
    printf("********** object %d data **********\n",i);
    printf("data file name  : %s\n",mo->mod_fn[i]);
    printf("shift parameter : (%8.7g, %8.7g, %8.7g)\n",mo->rs[i][0],mo->rs[i][1],mo->rs[i][2]);
    print_medium_data(&(mo->md[i]));
    print_mesh_data(&(mo->md[i]));

    if(mo->md[i].th!=0.0 || vabs_d(mo->md[i].tv)!=0.0){
      printf("-- rotation and translation settings --\n");
      if(mo->md[i].th!=0.0){
        printf("vector defining rotation axis :(% 8.7g,% 8.7g,% 8.7g)\n",mo->md[i].rv[0],mo->md[i].rv[1],mo->md[i].rv[2]);
        printf("rotation angle           [rad]: %8.7g\n",mo->md[i].th);
      }
      if(vabs_d(mo->md[i].tv)!=0.0){
        printf("translation vector         [m]:(%8.7g,%8.7g,%8.7g)\n",mo->md[i].tv[0],mo->md[i].tv[1],mo->md[i].tv[2]);
      }
      printf("\n");
    }
  }
}

void mo_finalize(MOBJ *mo)
{
  void finalize_cmd2(CMD *cm);  // bem3_aw_b2_solve_bieq.c
  void finalize_dmda(DMDA *ad); // bem3_aw_b1.c

  int N,i;

  N=mo->N;
  for(i=0;i<N;i++){
    finalize_cmd2(&(mo->md[i].cm));
    finalize_dmda(&(mo->md[i]));
    free(mo->mod_fn[i]);
    free(mo->rs[i]);
  }
  free(mo->md);

  free(mo->mod_fn);
  free(mo->rs);
}

void mo_dat_write(char *fname,MOBJ *mo)
{
  FILE *fp;
  int i,j,d,m;

  if((fp=fopen(fname,"wb"))==NULL){    printf("bem3_aw_b2.c, mo_dat_write(), Failed to create the %s file.\n",fname);    exit(1);  }
  
  fwrite(&(mo->N),sizeof(int),1,fp);
  for(m=0;m<mo->N;m++){
    fwrite(mo->mod_fn[m],sizeof(char),256,fp);
    fwrite(mo->rs[m],sizeof(double),3,fp);
    fwrite(&(mo->md[m]),sizeof(DMDA),1,fp);
    // material def
    fwrite(mo->md[m].rho0,sizeof(double),mo->md[m].MN+1,fp);
    fwrite(mo->md[m].c0,sizeof(double),mo->md[m].MN+1,fp);
    fwrite(mo->md[m].k0,sizeof(double),mo->md[m].MN+1,fp);
    fwrite(mo->md[m].K0,sizeof(double),mo->md[m].MN+1,fp);
    fwrite(mo->md[m].k1,sizeof(double complex),mo->md[m].MN+1,fp);
    fwrite(mo->md[m].k2,sizeof(double complex),mo->md[m].MN+1,fp);
    // beam data
    fwrite(mo->md[m].aw.bd.pw,sizeof(Apw),mo->md[m].aw.n_pw,fp);
    fwrite(mo->md[m].aw.bd.bb,sizeof(Abb),mo->md[m].aw.n_bb,fp);
    fwrite(mo->md[m].aw.bd.fb,sizeof(Afb),mo->md[m].aw.n_fb,fp);
    // BOUD
    for(i=0;i<=mo->md[m].bd.Nn;i++) fwrite(mo->md[m].bd.rn[i],sizeof(double),3,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fwrite(mo->md[m].bd.ed[i],sizeof(int),4,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fwrite(mo->md[m].bd.eni[i],sizeof(int),4,fp);
    fwrite(mo->md[m].bd.md,sizeof(int),mo->md[m].bd.Ne+1,fp);
    fwrite(mo->md[m].bd.sd,sizeof(int),mo->md[m].bd.Ne+1,fp);
    fwrite(mo->md[m].bd.gd,sizeof(int),mo->md[m].bd.Ne+1,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.ren[i][j],sizeof(double),3,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.wen[i][j],sizeof(double),3,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fwrite(mo->md[m].bd. Pi[i],sizeof(double complex),4,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fwrite(mo->md[m].bd.dPi[i],sizeof(double complex),4,fp);
    // sub domain data
    for(d=0;d<=mo->md[m].MN;d++){
      fwrite(&(mo->md[m].bd.sb[d].Ne),sizeof(int),1,fp);
      fwrite(mo->md[m].bd.sb[d].sid,sizeof(int),mo->md[m].bd.sb[d].Ne+1,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) fwrite(mo->md[m].bd.sb[d]. P[i],sizeof(double complex),4,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) fwrite(mo->md[m].bd.sb[d].dP[i],sizeof(double complex),4,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.sb[d].pv [i][j],sizeof(double complex),3,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(mo->md[m].bd.sb[d].dpv[i][j],sizeof(double complex),3,fp);
    }
  }
  fclose(fp);
}

void mo_dat_read(char *fname,MOBJ *mo)
{
  void mo_malloc(MOBJ *mo);
  void malloc_node(BOUD *bd); // bem3_aw_b1.c
  void malloc_elem(BOUD *bd); // bem3_aw_b1.c
  void init_elem_const(BOUD *bd); // bem3_aw_b1.c
  void malloc_sub_domain(DMDA *ad); // bem3_aw_b1.c
  void initalize_cmd2(CMD *cm); // bem3_aw_b2_solve_bieq.c

  FILE *fp;
  int i,j,d,tmp,m;

  if((fp=fopen(fname,"rb"))==NULL){    printf("bem3_aw_b2.c, mo_dat_read(), Failed to open the %s file.\n",fname);    exit(1);  }

  fread(&(mo->N),sizeof(int),1,fp);
  mo_malloc(mo);
  for(m=0;m<mo->N;m++){
    fread(mo->mod_fn[m],sizeof(char),256,fp);
    fread(mo->rs[m],sizeof(double),3,fp);
    fread(&(mo->md[m]),sizeof(DMDA),1,fp);
     // material def
    mo->md[m].rho0=(double *)m_alloc2(mo->md[m].MN+1,sizeof(double),"bem3_aw_b2.c,mo_dat_read(),mo->md[m].rho0");
    mo->md[m].c0  =(double *)m_alloc2(mo->md[m].MN+1,sizeof(double),"bem3_aw_b2.c,mo_dat_read(),mo->md[m].c0");
    mo->md[m].k0  =(double *)m_alloc2(mo->md[m].MN+1,sizeof(double),"bem3_aw_b2.c,mo_dat_read(),mo->md[m].k0");
    mo->md[m].K0  =(double *)m_alloc2(mo->md[m].MN+1,sizeof(double),"bem3_aw_b2.c,mo_dat_read(),mo->md[m].K0");
    mo->md[m].k1=(double complex *)m_alloc2(mo->md[m].MN+1,sizeof(double complex),"bem3_aw_b2.c,mo_dat_read(),mo->md[m].k1");
    mo->md[m].k2=(double complex *)m_alloc2(mo->md[m].MN+1,sizeof(double complex),"bem3_aw_b2.c,mo_dat_read(),mo->md[m].k2");
    fread(mo->md[m].rho0,sizeof(double),mo->md[m].MN+1,fp);
    fread(mo->md[m].c0,sizeof(double),mo->md[m].MN+1,fp);
    fread(mo->md[m].k0,sizeof(double),mo->md[m].MN+1,fp);
    fread(mo->md[m].K0,sizeof(double),mo->md[m].MN+1,fp);
    fread(mo->md[m].k1,sizeof(double complex),mo->md[m].MN+1,fp);
    fread(mo->md[m].k2,sizeof(double complex),mo->md[m].MN+1,fp);
    // beam data
    mo->md[m].aw.bd.pw=(Apw *)m_alloc2(mo->md[m].aw.n_pw,sizeof(Apw),"bem3_aw_b2.c,mo_dat_read(),mo->md[m].aw.bd.pw");
    fread(mo->md[m].aw.bd.pw,sizeof(Apw),mo->md[m].aw.n_pw,fp);
    mo->md[m].aw.bd.bb=(Abb *)m_alloc2(mo->md[m].aw.n_bb,sizeof(Abb),"bem3_aw_b2.c,mo_dat_read(),mo->md[m].aw.bd.bb");
    fread(mo->md[m].aw.bd.bb,sizeof(Abb),mo->md[m].aw.n_bb,fp);
    mo->md[m].aw.bd.fb=(Afb *)m_alloc2(mo->md[m].aw.n_fb,sizeof(Afb),"bem3_aw_b2.c,mo_dat_read(),mo->md[m].aw.bd.fb");
    fread(mo->md[m].aw.bd.fb,sizeof(Afb),mo->md[m].aw.n_fb,fp);
    setup_maw(&(mo->md[m].aw));
    // BOUD
    malloc_node(&(mo->md[m].bd)); // malloc
    malloc_elem(&(mo->md[m].bd)); // malloc
    for(i=0;i<=mo->md[m].bd.Nn;i++) fread(mo->md[m].bd.rn[i],sizeof(double),3,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fread(mo->md[m].bd.ed[i],sizeof(int),4,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fread(mo->md[m].bd.eni[i],sizeof(int),4,fp);
    fread(mo->md[m].bd.md,sizeof(int),mo->md[m].bd.Ne+1,fp);
    fread(mo->md[m].bd.sd,sizeof(int),mo->md[m].bd.Ne+1,fp);
    fread(mo->md[m].bd.gd,sizeof(int),mo->md[m].bd.Ne+1,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.ren[i][j],sizeof(double),3,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.wen[i][j],sizeof(double),3,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fread(mo->md[m].bd. Pi[i],sizeof(double complex),4,fp);
    for(i=0;i<=mo->md[m].bd.Ne;i++) fread(mo->md[m].bd.dPi[i],sizeof(double complex),4,fp);
    init_elem_const(&(mo->md[m].bd)); // setup
    // sub domain data
    malloc_sub_domain(&(mo->md[m])); // malloc
    for(d=0;d<=mo->md[m].MN;d++){
      fread(&tmp,sizeof(int),1,fp);
      fread(mo->md[m].bd.sb[d].sid,sizeof(int),mo->md[m].bd.sb[d].Ne+1,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) fread(mo->md[m].bd.sb[d]. P[i],sizeof(double complex),4,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) fread(mo->md[m].bd.sb[d].dP[i],sizeof(double complex),4,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.sb[d]. pv[i][j],sizeof(double complex),3,fp);
      for(i=0;i<=mo->md[m].bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(mo->md[m].bd.sb[d].dpv[i][j],sizeof(double complex),3,fp);
    }

    mo->md[m].cm.MN=mo->md[m].MN;
    mo->md[m].cm.na=2*mo->md[m].bd.NN;
    initalize_cmd2(&(mo->md[m].cm));
  }

  fclose(fp); 
}

void mo_output_node_particles(char *fname,MOBJ *mo)
{
  FILE *fp;
  int s1,s2,oid,i,j;
  char *sd,fo[128]="";

  sd=strrchr(fname,'.');
  if(sd==NULL){ // no file extension
    sprintf(fo,"%s.particles",fname);
  }
  else {
    s1=strlen(fname);
    s2=strlen(sd);
    strncpy(fo,fname,s1-s2);
    sprintf(fo,"%s.particles",fo);
  }
  
  if((fp=fopen(fo,"wt"))==NULL){    printf("Can not open the %s file.\n",fo);    exit(1);  }
  fprintf(fp,"# x y z object_id\n");
  
  for(oid=0;oid<mo->N;oid++){
    for(i=1;i<=mo->md[oid].bd.Ne;i++){
      for(j=0;j<4;j++){
        fprintf(fp,"%15.14e %15.14e %15.14e %d\n",mo->md[oid].bd.ren[i][j][0],mo->md[oid].bd.ren[i][j][1],mo->md[oid].bd.ren[i][j][2],oid);
      }
    }
  }

  fclose(fp);
}

////////////////////////////////////////////////////////////////////////
void init_boundary_data2(DMDA *ad)
{
  double cr[3][4],cw[3][3],r[3],w[3];
  int i,j,l,m;
  
  // element node and incident field data
  for(i=1;i<=ad->bd.Ne;i++){
    for(l=0;l<3;l++){
      for(m=0;m<4;m++) cr[l][m]=ad->bd.cr[i][l][m];
      for(m=0;m<3;m++) cw[l][m]=ad->bd.cw[i][l][m];
    }

    if(ad->bd.ed[i][3]==0){ // linear triangular element
      for(j=0;j<4;j++){
        lit_rw_zeta_eta(r,w,ad->bd.zt_34[j],ad->bd.et_34[j],cr,cw);
        for(l=0;l<3;l++){
          ad->bd.ren[i][j][l]=r[l];
          ad->bd.wen[i][j][l]=w[l];
        }
      }
    }
    else { // bi-linear element
      for(j=0;j<4;j++){
        bil_rw_zeta_eta(r,w,ad->bd.zt_44[j],ad->bd.et_44[j],cr,cw);
        for(l=0;l<3;l++){
          ad->bd.ren[i][j][l]=r[l];
          ad->bd.wen[i][j][l]=w[l];
        }
      }
    }
  }
}

void mo_malloc(MOBJ *mo)
{
  int i,N;

  N=mo->N;
  mo->md=(DMDA *)m_alloc2(N,sizeof(DMDA),"bem3_aw_b2.c, mo_malloc(), mo->md");

  mo->mod_fn=(char **)m_alloc2(N,sizeof(char *),"bem3_aw_b2.c, mo_malloc(), mo->mod_fn");
  mo->rs=(double **)m_alloc2(N,sizeof(double *),"bem3_aw_b2.c, mo_malloc(), mo->rs");
  for(i=0;i<N;i++){
    mo->mod_fn[i]=(char *)m_alloc2(256,sizeof(char),"bem3_aw_b2.c, mo_malloc(), mo->mod_fn[i]");
    mo->rs[i]=(double *)m_alloc2(3,sizeof(double),"bem3_aw_b2.c, mo_malloc(), mo->rs[i]");
  }
}

void dat_read2(char *fname,DMDA *md)
{
  void malloc_node(BOUD *bd); // bem3_aw_b1.c
  void malloc_elem(BOUD *bd); // bem3_aw_b1.c
  void malloc_sub_domain(DMDA *ad); // bem3_aw_b1.c
  void initalize_cmd2(CMD *cm); // bem3_aw_b2_solve_bieq.c

  FILE *fp;
  int i,j,d,tmp;

  if((fp=fopen(fname,"rb"))==NULL){    printf("bem3_aw_b2.c, dat_read2(), Failed to open the %s file.\n",fname);    exit(1);  }

  fread(md,sizeof(DMDA),1,fp);
  // medium data
  md->rho0=(double *)m_alloc2(md->MN+1,sizeof(double),"bem3_aw_b2.c, dat_read2(),md->rho0");
  md->c0  =(double *)m_alloc2(md->MN+1,sizeof(double),"bem3_aw_b2.c, dat_read2(),md->c0");
  md->k0  =(double *)m_alloc2(md->MN+1,sizeof(double),"bem3_aw_b2.c, dat_read2(),md->k0");
  md->K0  =(double *)m_alloc2(md->MN+1,sizeof(double),"bem3_aw_b2.c, dat_read2(),md->K0");
  md->k1=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"bem3_aw_b2.c, dat_read2(),md->k1");
  md->k2=(double complex *)m_alloc2(md->MN+1,sizeof(double complex),"bem3_aw_b2.c, dat_read2(),md->k2");
  fread(md->rho0,sizeof(double),md->MN+1,fp);
  fread(md->c0,sizeof(double),md->MN+1,fp);
  fread(md->k0,sizeof(double),md->MN+1,fp);
  fread(md->K0,sizeof(double),md->MN+1,fp);
  fread(md->k1,sizeof(double complex),md->MN+1,fp);
  fread(md->k2,sizeof(double complex),md->MN+1,fp);
  // beam data
  md->aw.bd.pw=(Apw *)m_alloc2(md->aw.n_pw,sizeof(Apw),"bem3_aw_b2.c,dat_read2(),md->aw.bd.pw");
  fread(md->aw.bd.pw,sizeof(Apw),md->aw.n_pw,fp);
  md->aw.bd.bb=(Abb *)m_alloc2(md->aw.n_bb,sizeof(Abb),"bem3_aw_b2.c,dat_read2(),md->aw.bd.bb");
  fread(md->aw.bd.bb,sizeof(Abb),md->aw.n_bb,fp);
  md->aw.bd.fb=(Afb *)m_alloc2(md->aw.n_fb,sizeof(Afb),"bem3_aw_b2.c,dat_read2(),md->aw.bd.fb");
  fread(md->aw.bd.fb,sizeof(Afb),md->aw.n_fb,fp);
  // reset beam data
  init_maw(&(md->aw));
  read_data_maw(&(md->aw));   
  setup_maw(&(md->aw));
  // check parameter of incident field
  if(md->rho0[0]!=md->aw.rho0 || md->c0[0]!=md->aw.c0 || md->k0[0]!=md->aw.k0){
    printf("incident field data mismatched. check incident field datafile.\n");
    printf("model rho0[0]=% 8.6g, incident field rho0=% 8.6g\n",md->rho0[0],md->aw.rho0);
    printf("model c0[0]  =% 8.6g, incident field c0  =% 8.6g\n",md->c0[0],md->aw.c0);
    printf("model k0[0]  =% 8.6g, incddent field k0  =% 8.6g\n",md->k0[0],md->aw.k0);
    printf("Exit...\n");
    exit(0);
  }
  // BOUD
  malloc_node(&(md->bd)); // malloc
  malloc_elem(&(md->bd)); // malloc
  for(i=0;i<=md->bd.Nn;i++) fread(md->bd.rn[i],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) fread(md->bd.ed[i],sizeof(int),4,fp);
  for(i=0;i<=md->bd.Ne;i++) fread(md->bd.eni[i],sizeof(int),4,fp);
  fread(md->bd.md,sizeof(int),md->bd.Ne+1,fp);
  fread(md->bd.sd,sizeof(int),md->bd.Ne+1,fp);
  fread(md->bd.gd,sizeof(int),md->bd.Ne+1,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fread(md->bd.ren[i][j],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fread(md->bd.wen[i][j],sizeof(double),3,fp);
  // sub domain data
  malloc_sub_domain(md); // malloc
  for(d=0;d<=md->MN;d++){
    fread(&tmp,sizeof(int),1,fp);
    fread(md->bd.sb[d].sid,sizeof(int),md->bd.sb[d].Ne+1,fp);
  }
  // coefficient matrix data
  md->cm.MN=md->MN;
  md->cm.na=2*md->bd.NN;
  initalize_cmd2(&(md->cm));
  for(i=0;i<=md->MN;i++){
    fread(md->cm.tgfn[i],sizeof(char),128,fp); // G
    fread(md->cm.thfn[i],sizeof(char),128,fp); // H
    fread(md->cm.tdgfn[i],sizeof(char),128,fp); // dG
    fread(md->cm.tdhfn[i],sizeof(char),128,fp); // dH
    fread(md->cm.tdffn[i],sizeof(char),128,fp); // dF
  }
  fread(md->cm.lupfn,sizeof(char),128,fp); // LU + pivot

  fclose(fp);
}

void init_boundary_data_u2(DMDA *md)
{
  double complex p,dp,cp;
  double r[3],w[3];
  int i,j,l;

  cp=-1.0/md->aw.k2;

  // element node and incident field data
  for(i=1;i<=md->bd.Ne;i++){
    if(md->bd.ed[i][3]==0){ // linear triangular element
      for(j=0;j<4;j++){
        for(l=0;l<3;l++){
          r[l]=md->bd.ren[i][j][l];
          w[l]=md->bd.wen[i][j][l];
        }
        vuni_d(w);
        calc_maw_dpdn(&p,&dp,r,w,&(md->aw)); // sound pressure 
        p*=cp; // velocity potential
        dp*=cp;
        md->bd.Pi [i][j]=p;
        md->bd.dPi[i][j]=dp;
        md->bd.Pi0 [i][j]=p;
        md->bd.dPi0[i][j]=dp;
        md->bd.Pip [i][j]=p;
      }
    }
    else { // bi-linear element
      for(j=0;j<4;j++){
        for(l=0;l<3;l++){
          r[l]=md->bd.ren[i][j][l];
          w[l]=md->bd.wen[i][j][l];
        }
        vuni_d(w);
        calc_maw_dpdn(&p,&dp,r,w,&(md->aw)); // sound pressure
        p*=cp; // velocity potential
        dp*=cp;
        md->bd.Pi [i][j]=p;
        md->bd.dPi[i][j]=dp;
        md->bd.Pi0 [i][j]=p;
        md->bd.dPi0[i][j]=dp;
        md->bd.Pip [i][j]=p;
      }
    }
  }
}
