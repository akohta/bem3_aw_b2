#include "bem3_aw_b2.h"

void solve_bieq_dmda2(DMDA *ad,char *ofn)
{
  void create_cmatrix(CMD *cm,DMDA *md); // bem3_aw_b1_solve_bieq.c
  void initalize_cmd2(CMD *cm);
  void create_Amatrix_csr(DMDA *md);
  void lu_dec_A(CMD *cm);
  void dat_write2(char *fname,DMDA *md);
  void finalize_cmd2(CMD *cm);
  
  time_t start,end,ms,me;

  printf("\ncreate matrix data \n");
  time(&start);
  
  printf("  coefficient matrix          "); fflush(stdout);
  time(&ms);
  ad->cm.type=2; // precision setting 0:4p GL,1:9p or 7p(triangular) GL, 2:GLN p GL, 3: GHN p GL
  ad->cm.MN=ad->MN;
  initalize_cmd2(&(ad->cm));
  create_cmatrix(&(ad->cm),ad);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  printf("  inverse matrix              "); fflush(stdout);
  time(&ms);
  ad->cm.nn=(size_t)ad->bd.NN;
  ad->cm.na=ad->cm.nn*2;
  create_Amatrix_csr(ad);
  lu_dec_A(&(ad->cm));
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));
  dat_write2(ofn,ad);

  finalize_cmd2(&(ad->cm));
  time(&end);
  printf("Total elapsed time : %g (sec)\n",difftime(end,start));
}

void mo_solve_iter(MOBJ *mo)
{
  void solve_coefficient(int flg,DMDA *md);
  void solve_coefficient_vp(int flg,DMDA *md);
  void solve_coefficient_pv(int flg,DMDA *md);
  void reset_incident_bd(DMDA *md);
  void add_scattered_vp(DMDA *dst,DMDA *src);
  double ccd_f(MOBJ *mo);

  time_t start,end,ms,me;
  double cc;
  int iter,m,d,itermax;

  itermax=ITER_MAX;

  if(mo->N==1){ // single object
    printf("\nsingle object\n");
    solve_coefficient(1,&(mo->md[0]));
  }
  else { // multi object
    printf("\ninitialize boundary value\n");
    for(m=0;m<mo->N;m++){ // calc boundary value each object
      printf("object id = %2d\n",m);
      solve_coefficient_vp(1,&(mo->md[m]));
    }

    printf("\nunit of iterative operation\n");
    time(&start);
    for(m=0;m<mo->N;m++){ // reset incidient field
      reset_incident_bd(&(mo->md[m]));
    }
    for(m=0;m<mo->N;m++){ // add scattered field to initial incidient field
      printf("object id = %2d\n",m);
      printf("  add scattered field       "); fflush(stdout);
      time(&ms);
      for(d=0;d<mo->N;d++) {
        if(d==m) continue;
        add_scattered_vp(&(mo->md[d]),&(mo->md[m]));
      }
      time(&me);
      printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));
    }
    for(m=0;m<mo->N;m++){
      printf("object id = %2d\n",m);
      solve_coefficient_vp(1,&(mo->md[m]));
    }
    time(&end);
    printf("Estimated time per iteration : %5g (sec)\n\n",difftime(end,start));

    printf("iterative operation start. convergence criterion cc < %g, iter>%d\n",CCD,ITER_MIN);
    time(&start);
    for(iter=0;iter<itermax;iter++){
      printf("iter = %2d",iter);
      for(m=0;m<mo->N;m++){ // renew incidient field
        reset_incident_bd(&(mo->md[m]));
      }
      for(m=0;m<mo->N;m++){ // add scattered field to incidient field
        for(d=0;d<mo->N;d++) {
          if(d==m) continue;
          add_scattered_vp(&(mo->md[d]),&(mo->md[m]));
        }
      }
      for(m=0;m<mo->N;m++){
        solve_coefficient_vp(0,&(mo->md[m]));
      }

      cc=ccd_f(mo);
      printf(", cc = %g\n",cc);
      if(cc<CCD && iter>ITER_MIN) break;
    }
    for(m=0;m<mo->N;m++) solve_coefficient_pv(0,&(mo->md[m]));

    time(&end);
    printf("finished. Elapsed time : %5g (sec)\n",difftime(end,start));
    if(iter==itermax){
      printf("The number of iterations has reached the upper limit. The solution has not converged.\n");
    } 
  }
}


///////////////////////////////////////////////////////////////////////
void initalize_cmd2(CMD *cm)
{
  int i,MN;

  MN=cm->MN;
  cm->tgfn=(char **)m_alloc2(MN+1,sizeof(char *),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->tgfn");
  cm->thfn=(char **)m_alloc2(MN+1,sizeof(char *),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->thfn");
  cm->tdgfn=(char **)m_alloc2(MN+1,sizeof(char *),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->tdgfn");
  cm->tdhfn=(char **)m_alloc2(MN+1,sizeof(char *),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->tdhfn");
  cm->tdffn=(char **)m_alloc2(MN+1,sizeof(char *),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->tdffn");
  for(i=0;i<=MN;i++){
    cm->tgfn[i]=(char *)m_alloc2(128,sizeof(char ),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->tgfn[i]");
    cm->thfn[i]=(char *)m_alloc2(128,sizeof(char ),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->thfn[i]");
    sprintf(cm->tgfn[i],"tmpG_%05d.cmat",i);
    sprintf(cm->thfn[i],"tmpH_%05d.cmat",i);
    cm->tdgfn[i]=(char *)m_alloc2(128,sizeof(char ),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->tdgfn[i]");
    cm->tdhfn[i]=(char *)m_alloc2(128,sizeof(char ),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->tdhfn[i]");
    cm->tdffn[i]=(char *)m_alloc2(128,sizeof(char ),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->tdffn[i]");
    sprintf(cm->tdgfn[i],"tmpdG_%05d.cmat",i);
    sprintf(cm->tdhfn[i],"tmpdH_%05d.cmat",i);
    sprintf(cm->tdffn[i],"tmpdF_%05d.cmat",i);
  }
  cm->lupfn=(char *)m_alloc2(128,sizeof(char ),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->tgfn[i]");
  sprintf(cm->lupfn,"tmpinvA.cmat");

  cm->aval=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->aval");
  cm->aptr=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->aptr");
  cm->aidx=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_b2_solve_bieq.c, initialize_cmd2(), cm->aidx");
  sprintf(cm->aval,"tmpAval.dat");
  sprintf(cm->aptr,"tmpAprt.dat");
  sprintf(cm->aidx,"tmpAidx.dat");
  cm->b=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_b2_solve_bieq.c, initialize_cmd(), cm->b");
  sprintf(cm->b,"tmpB.dat");
}

void finalize_cmd2(CMD *cm)
{
  int i;

  remove(cm->aval);
  remove(cm->aptr);
  remove(cm->aidx);
  remove(cm->b);

  // free memory
  for(i=0;i<=cm->MN;i++){
    free(cm->tgfn[i]);    free(cm->thfn[i]);
    free(cm->tdgfn[i]);    free(cm->tdhfn[i]);  free(cm->tdffn[i]);
  }
  free(cm->tgfn);  free(cm->thfn);
  free(cm->tdgfn);  free(cm->tdhfn); free(cm->tdffn);
  free(cm->lupfn);

  cm->MN=0;

  free(cm->aval);
  free(cm->aptr);
  free(cm->aidx);
  cm->nn=0;
  cm->na=0;
  cm->nnz=0;
  free(cm->b);
}

void create_Amatrix_csr(DMDA *md)
{
  void create_Amatrix_csr_dac(int did,FILE *av,FILE *ap,FILE *ai,CMD *cm,DMDA *md);

  FILE *av,*ai,*ap;

  size_t did;

  if((av=fopen(md->cm.aval,"wb"))==NULL){    printf("bem3_aw_b2_solve_bieq.c, create_Amatrix_csr(),*av. Failed to create %s file.\n",md->cm.aval);    exit(1);  }
  if((ai=fopen(md->cm.aidx,"wb"))==NULL){    printf("bem3_aw_b2_solve_bieq.c, create_Amatrix_csr(),*ai. Failed to create %s file.\n",md->cm.aidx);    exit(1);  }
  if((ap=fopen(md->cm.aptr,"wb"))==NULL){    printf("bem3_aw_b2_solve_bieq.c, create_Amatrix_csr(),*ap. Failed to create %s file.\n",md->cm.aidx);    exit(1);  }

  md->cm.nnz=0; // initialize nnz
  for(did=0;did<=md->cm.MN;did++) create_Amatrix_csr_dac(did,av,ap,ai,&(md->cm),md);

  fwrite(&(md->cm.nnz),sizeof(size_t),1,ap);

  fclose(av);
  fclose(ai);
  fclose(ap);
}

void create_Amatrix_csr_dac(int did,FILE *av,FILE *ap,FILE *ai,CMD *cm,DMDA *md)
{
  FILE *fg,*fh;

  double complex *tG,*tH,*tA,k2m,k2s,hk2;
  size_t Ne,t,tn,s,asd,mdid,sdid,etype,l,cc,*ti;
  int td,sd;

  if((fg=fopen(cm->tgfn[did],"rb"))==NULL){    printf("bem3_aw_b2_solve_bieq.c, create_Amatrix_csr_dac(),*fg. Failed to open %s file.\n",cm->tgfn[did]);    exit(1);  }
  if((fh=fopen(cm->thfn[did],"rb"))==NULL){    printf("bem3_aw_b2_solve_bieq.c, create_Amatrix_csr_dac(),*fh. Failed to open %s file.\n",cm->thfn[did]);    exit(1);  }

  Ne=(size_t)md->bd.sb[did].Ne;
  tG=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"bem3_aw_b2_solve_bieq.c, create_Amatrix_csr_dac(),tG"); // malloc
  tH=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"bem3_aw_b2_solve_bieq.c, create_Amatrix_csr_dac(),tH"); // malloc

  tA=(double complex *)m_alloc2(cm->na,sizeof(double complex),"bem3_aw_b2_solve_bieq.c, create_Amatrix_csr_dac(),tA"); // malloc
  ti=(size_t *)m_alloc2(cm->na,sizeof(size_t),"bem3_aw_b2_solve_bieq.c, create_Amatrix_csr_dac(),ti"); // malloc

  for(t=1;t<=Ne;t++){
    td=md->bd.sb[did].sid[t];

    for(tn=0;tn<4;tn++){
      if(fread(tG,sizeof(double complex),Ne*4,fg)!=Ne*4){
        printf("bem3_aw_b2_solve_bieq.c, create_Amatrix_csr_dac(), failed to read the tG. exit...\n");
        exit(1);
      }
      if(fread(tH,sizeof(double complex),Ne*4,fh)!=Ne*4){
        printf("bem3_aw_b2_solve_bieq.c, create_Amatrix_csr_dac(), failed to read the tG. exit...\n");
        exit(1);
      }
      if( tn==3 && ELT3==check_element_type(td,&(md->bd)) )  continue;

      fwrite(&(cm->nnz),sizeof(size_t),1,ap); // write A pointer
      for(l=0;l<cm->na;l++) tA[l]=0.0;

      for(s=1;s<=Ne;s++){
        sd=md->bd.sb[did].sid[s]; // signed element id
        asd=(size_t)abs(sd);
        mdid=(size_t)md->bd.md[asd]; // main domain id
        sdid=(size_t)md->bd.sd[asd]; // sub domain id
        etype=check_element_type(sd,&(md->bd));        //printf("sd = %d, mdid=%d, sdid=%d\n",sd,mdid,sdid); // test

        if(did==mdid){ // main domain
          if(etype==ELT3){ // linear-triangular element
            for(l=0;l<3;l++) {
              tA[ cm->nn*0 + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
              tA[ cm->nn*1 + md->bd.eni[asd][l] ]=-tG[(s-1)*4+l];
            }
          }
          else { // bi-linear element
            for(l=0;l<4;l++) {
              tA[ cm->nn*0 + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
              tA[ cm->nn*1 + md->bd.eni[asd][l] ]=-tG[(s-1)*4+l];
            }
          }
        } // end main domain
        else { // sub domain
          k2m=md->k2[mdid];
          k2s=md->k2[sdid];
          hk2=k2m/k2s;
          if(etype==ELT3){ // linear-triangular element
            for(l=0;l<3;l++){
              tA[ cm->nn*0 + md->bd.eni[asd][l] ]=hk2*tH[(s-1)*4+l];
              tA[ cm->nn*1 + md->bd.eni[asd][l] ]=    tG[(s-1)*4+l];
            }
          }
          else { // bi-linear element
            for(l=0;l<4;l++){
              tA[ cm->nn*0 + md->bd.eni[asd][l] ]=hk2*tH[(s-1)*4+l];
              tA[ cm->nn*1 + md->bd.eni[asd][l] ]=    tG[(s-1)*4+l];
            }
          }
        } // end sub domain
      } // end for s

      // compress and store data
      cc=0;
      for(l=0;l<cm->na;l++){
        if( creal(tA[l])==0.0 && cimag(tA[l])==0.0) continue;
        tA[cc]=tA[l];
        ti[cc]=l;
        cc+=1;
      }
      fwrite(tA,sizeof(double complex),cc,av);
      fwrite(ti,sizeof(size_t),cc,ai);
      cm->nnz+=cc;

    } // end for tn
  } // end for t

  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);
  free(tA);
  free(ti);
}

void lu_dec_A(CMD *cm)
{
  FILE *fa,*fxa,*fas;
  MKL_Complex16 *A,tc;
  MKL_INT is,ie,i,j,p;

  A=(MKL_Complex16 *)m_alloc2(cm->na*cm->na,sizeof(MKL_Complex16),"bem3_aw_b2_solve_bieq.c, lu_dec_A(),A");
  // read matrix A
  if((fa=fopen(cm->aval,"rb"))==NULL){     printf("bem3_aw_b2_solve_bieq.c,lu_dec_A(), Failed to open the %s file.\n",cm->aval);    exit(1); }
  if((fxa=fopen(cm->aptr,"rb"))==NULL){     printf("bem3_aw_b2_solve_bieq.c, lu_dec_A(), Failed to open the %s file.\n",cm->aptr);    exit(1); }
  if((fas=fopen(cm->aidx,"rb"))==NULL){     printf("bem3_aw_b2_solve_bieq.c, lu_dec_A(), Failed to open the %s file.\n",cm->aidx);    exit(1); }

  if(fread(&is,sizeof(MKL_INT),1,fxa)!=1){
    printf("bem3_aw_b2_solve_bieq.c, lu_dec_A(), failed to read the is. exit...\n");
    exit(1);
  }    
  
  for(j=0;j<cm->na;j++){
    if(fread(&ie,sizeof(MKL_INT),1,fxa)!=1){
      printf("bem3_aw_b2_solve_bieq.c, lu_dec_A(), failed to read the ie. exit...\n");
      exit(1);
    }
    for(p=is;p<ie;p++){
      if(fread(&i,sizeof(MKL_INT),1,fas)!=1){
        printf("bem3_aw_b2_solve_bieq.c, lu_dec_A(), failed to read the i. exit...\n");
        exit(1);
      }
      if(fread(&tc,sizeof(MKL_Complex16),1,fa)!=1){
        printf("bem3_aw_b2_solve_bieq.c, lu_dec_A(), failed to read the A. exit...\n");
        exit(1);
      }
      A[j*cm->na+i]=tc;
    }
    is=ie;
  }
  fclose(fa);
  fclose(fxa);
  fclose(fas);

  MKL_INT N,lda,info,*ipiv;

  N=(MKL_INT)cm->na;
  lda=N;
  ipiv=(MKL_INT*)m_alloc2(cm->na,sizeof(MKL_INT),"bem3_aw_b2_solve_bieq.c, lu_dec_A(),ipiv");

  // LU factorization
  info=LAPACKE_zgetrf(LAPACK_ROW_MAJOR , N , N , A , lda , ipiv );
  if(info!=0){
    printf("bem3_aw_b2_solve_bieq.c, lu_dec_A(), LAPCKE_zgetrf() error! info=%d. Exit...\n",(int)info);
    exit(1);
  }
  // inverse matrix
  info=LAPACKE_zgetri(LAPACK_ROW_MAJOR , N , A , lda , ipiv );
  if(info!=0){
    printf("bem3_aw_b2_solve_bieq.c, lu_dec_A(), LAPCKE_zgetri() error! info=%d. Exit...\n",(int)info);
    exit(1);
  }

  // output
  if((fa=fopen(cm->lupfn,"wb"))==NULL){     printf("bem3_aw_b2_solve_bieq.c, lu_dec_A(), Failed to open the %s file.\n",cm->lupfn);    exit(1); }
  fwrite(A,sizeof(MKL_Complex16),cm->na*cm->na,fa); // matrix A
  fclose(fa);

  free(A);
  free(ipiv);
}

void dat_write2(char *fname,DMDA *md)
{
  FILE *fp;
  char tfn[128],pfn[125];
  int i,j,d;

  if((fp=fopen(fname,"wb"))==NULL){    printf("solve_bieq2.c, dat_write2(), Failed to create the %s file.\n",fname);    exit(1);  }
  
  fwrite(md,sizeof(DMDA),1,fp);
  // medium data
  fwrite(md->rho0,sizeof(double),md->MN+1,fp);
  fwrite(md->c0,sizeof(double),md->MN+1,fp);
  fwrite(md->k0,sizeof(double),md->MN+1,fp);
  fwrite(md->K0,sizeof(double),md->MN+1,fp);
  fwrite(md->k1,sizeof(double complex),md->MN+1,fp);
  fwrite(md->k2,sizeof(double complex),md->MN+1,fp);
  // beam data
  fwrite(md->aw.bd.pw,sizeof(Apw),md->aw.n_pw,fp);
  fwrite(md->aw.bd.bb,sizeof(Abb),md->aw.n_bb,fp);
  fwrite(md->aw.bd.fb,sizeof(Afb),md->aw.n_fb,fp);
  // BOUD
  for(i=0;i<=md->bd.Nn;i++) fwrite(md->bd.rn[i],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.ed[i],sizeof(int),4,fp);
  for(i=0;i<=md->bd.Ne;i++) fwrite(md->bd.eni[i],sizeof(int),4,fp);
  fwrite(md->bd.md,sizeof(int),md->bd.Ne+1,fp);
  fwrite(md->bd.sd,sizeof(int),md->bd.Ne+1,fp);
  fwrite(md->bd.gd,sizeof(int),md->bd.Ne+1,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.ren[i][j],sizeof(double),3,fp);
  for(i=0;i<=md->bd.Ne;i++) for(j=0;j<4;j++) fwrite(md->bd.wen[i][j],sizeof(double),3,fp);
  // sub domain data
  for(d=0;d<=md->MN;d++){
    fwrite(&(md->bd.sb[d].Ne),sizeof(int),1,fp);
    fwrite(md->bd.sb[d].sid,sizeof(int),md->bd.sb[d].Ne+1,fp);
  }

  // matrix data
  for(i=0;i<128;i++){
    if(fname[i]=='.'){
      pfn[i]='\0';
      break;
    }
    pfn[i]=fname[i];
  }
  for(i=0;i<=md->MN;i++){
    sprintf(tfn,"%s_%s",pfn,md->cm.tgfn[i]+3);  rename(md->cm.tgfn[i],tfn);     fwrite(tfn,sizeof(char),128,fp); // G
    sprintf(tfn,"%s_%s",pfn,md->cm.thfn[i]+3);  rename(md->cm.thfn[i],tfn);     fwrite(tfn,sizeof(char),128,fp); // H
    sprintf(tfn,"%s_%s",pfn,md->cm.tdgfn[i]+3); rename(md->cm.tdgfn[i],tfn);    fwrite(tfn,sizeof(char),128,fp); // dG
    sprintf(tfn,"%s_%s",pfn,md->cm.tdhfn[i]+3); rename(md->cm.tdhfn[i],tfn);    fwrite(tfn,sizeof(char),128,fp); // dH
    sprintf(tfn,"%s_%s",pfn,md->cm.tdffn[i]+3); rename(md->cm.tdffn[i],tfn);    fwrite(tfn,sizeof(char),128,fp); // dF
  }
  sprintf(tfn,"%s_%s",pfn,md->cm.lupfn+3); rename(md->cm.lupfn,tfn);    fwrite(tfn,sizeof(char),128,fp); // LU + pivot

  fclose(fp);
}



void solve_coefficient(int flg,DMDA *md)
{
  void tmatrix_bd_store(double complex *X,DMDA *md); // bem3_aw_b1_solve_bieq.c
  void solve_pv_bv(CMD *cm,DMDA *ad); // bem3_aw_b1_solve_bieq.c
  void solve_dpv_bv(CMD *cm,DMDA *ad); // bem3_aw_b1_solve_bieq.c
  void create_Bmatrix(double complex *B,DMDA *md);
  void solve_ABmatrix(double complex *Bc,DMDA *md);

  time_t ms,me;
  double complex *B;

  if(flg==1){ printf("  solve VP boundary value   "); fflush(stdout); }
  time(&ms);
  B=(double complex *)m_alloc2(md->cm.na,sizeof(double complex),"bem3_aw_b2_solve_bieq.c, solve_coefficient(),B");
  create_Bmatrix(B,md);
  solve_ABmatrix(B,md);
  tmatrix_bd_store(B,md);
  free(B);

  time(&me);
  if(flg==1){ printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms)); }

  if(flg==1){ printf("  solve PV boundary value   "); fflush(stdout); }
  time(&ms);
  solve_pv_bv(&(md->cm),md);
  solve_dpv_bv(&(md->cm),md);
  time(&me);
  if(flg==1){ printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms)); }
}

void create_Bmatrix(double complex *B,DMDA *md)
{
  void create_Bmatrix_dac(size_t *j,size_t did,double complex *B,DMDA *md);

  size_t did,j;

  j=0;
  for(did=0;did<=md->MN;did++) create_Bmatrix_dac(&j,did,B,md);
}

void create_Bmatrix_dac(size_t *j,size_t did,double complex *B,DMDA *md)
{
  FILE *fg,*fh;

  double complex *tG,*tH,k2m,k2s,hk2;
  size_t Ne,t,tn,s,asd,mdid,sdid,etype,l;
  int td,sd;

  if((fg=fopen(md->cm.tgfn[did],"rb"))==NULL){    printf("bem3_aw_b2_solve_bieq.c, create_Bmatrix_dac(),*fg. Failed to open %s file.\n",md->cm.tgfn[did]);    exit(1);  }
  if((fh=fopen(md->cm.thfn[did],"rb"))==NULL){    printf("bem3_aw_b2_solve_bieq.c, create_Bmatrix_dac(),*fh. Failed to open %s file.\n",md->cm.thfn[did]);    exit(1);  }

  Ne=(size_t)md->bd.sb[did].Ne;
  tG=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"bem3_aw_b2_solve_bieq.c, create_Bmatrix_dac(),tG"); // malloc
  tH=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"bem3_aw_b2_solve_bieq.c, create_Bmatrix_dac(),tH"); // malloc

  for(t=1;t<=Ne;t++){
    td=md->bd.sb[did].sid[t];

    for(tn=0;tn<4;tn++){
      if(fread(tG,sizeof(double complex),Ne*4,fg)!=Ne*4){
        printf("bem3_aw_b2_solve_bieq.c, create_Bmatrix_dac(), failed to read the tG. exit...\n");
        exit(1);
      }
      if(fread(tH,sizeof(double complex),Ne*4,fh)!=Ne*4){
        printf("bem3_aw_b2_solve_bieq.c, create_Bmatrix_dac(), failed to read the tH. exit...\n");
        exit(1);
      }
      if( tn==3 && ELT3==check_element_type(td,&(md->bd)) )  continue;

      B[*j]=0.0;

      for(s=1;s<=Ne;s++){
        sd=md->bd.sb[did].sid[s]; // signed element id
        asd=(size_t)abs(sd);
        mdid=(size_t)md->bd.md[asd]; // main domain id
        sdid=(size_t)md->bd.sd[asd]; // sub domain id
        etype=check_element_type(sd,&(md->bd));      

        if(did!=mdid){ // sub domain
          k2m=md->k2[mdid];
          k2s=md->k2[sdid];
          hk2=k2m/k2s;
          if(etype==ELT3){ // linear-triangular element
            for(l=0;l<3;l++){
              if(mdid==0) B[*j]+=-hk2*tH[(s-1)*4+l]*md->bd.Pi[asd][l]-tG[(s-1)*4+l]*md->bd.dPi[asd][l];
            }
          }
          else { // bi-linear element
            for(l=0;l<4;l++){
              if(mdid==0) B[*j]+=-hk2*tH[(s-1)*4+l]*md->bd.Pi[asd][l]-tG[(s-1)*4+l]*md->bd.dPi[asd][l];
            }
          }
        } // end sub domain
      } // end for s
      *j+=1;
    } // end for tn
  } // end for t

  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);
}

void solve_ABmatrix(double complex *Bc,DMDA *md)
{
  FILE *fa;
  double complex *tA,*X;
  size_t N,i,j;

  N=md->cm.na;

  tA=(double complex *)m_alloc2(N,sizeof(double complex),"bem3_aw_b2_solve_bieq.c, solve_ABmatrix(),tA");
  X=(double complex *)m_alloc2(N,sizeof(double complex),"bem3_aw_b2_solve_bieq.c, solve_ABmatrix(),X");

  if((fa=fopen(md->cm.lupfn,"rb"))==NULL){     printf("bem3_aw_b2_solve_bieq.c, solve_ABmatrix(), Failed to open the %s file.\n",md->cm.lupfn);    exit(1); }

  for(i=0;i<N;i++){
    if(fread(tA,sizeof(double complex),N,fa)!=N){
      printf("bem3_aw_b2_solve_bieq.c, solve_ABmatrix(), failed to read the tA. exit...\n");
      exit(1);
    }
    X[i]=0.0;
    for(j=0;j<N;j++) X[i]+=tA[j]*Bc[j];
  }

  for(i=0;i<N;i++){
    Bc[i]=X[i];
  }

  free(tA);
  free(X);
}

void solve_coefficient_vp(int flg,DMDA *md)
{
  void tmatrix_bd_store(double complex *X,DMDA *md); // bem3_aw_b1_solve_bieq.c
  void create_Bmatrix(double complex *B,DMDA *md);
  void solve_ABmatrix(double complex *Bc,DMDA *md);

  time_t ms,me;
  double complex *B;

  if(flg==1){ printf("  solve VP boundary value   "); fflush(stdout); }
  time(&ms);
  B=(double complex *)m_alloc2(md->cm.na,sizeof(double complex),"bem3_aw_b2_solve_bieq.c, solve_coefficient(),B");
  create_Bmatrix(B,md);
  solve_ABmatrix(B,md);
  tmatrix_bd_store(B,md);
  free(B);

  time(&me);
  if(flg==1) printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));
}

void reset_incident_bd(DMDA *md)
{
  int i,j;

  for(i=1;i<=md->bd.Ne;i++){
    for(j=0;j<4;j++){
      md->bd.Pip[i][j]=md->bd.Pi[i][j];
      md->bd. Pi[i][j]=md->bd. Pi0[i][j];
      md->bd.dPi[i][j]=md->bd.dPi0[i][j];
    }
  }
}

void add_scattered_vp(DMDA *dst,DMDA *src)
{
  int itr_vp_s(double complex *p,double complex *dp,double *rt,double *n,int type,DMDA *md);

  double complex p,dp;
  double r[3],w[3];
  int i,j,l,did;

  int type=PREC_DEF;

  #pragma omp parallel for schedule(dynamic) private(j,l,r,w,p,dp,did)
  for(i=1;i<=dst->bd.Ne;i++){
    if(dst->bd.ed[i][3]==0){ // linear triangular element
      for(j=0;j<4;j++){
        if(dst->bd.md[i]==0){
          for(l=0;l<3;l++){
            r[l]=dst->bd.ren[i][j][l];
            w[l]=dst->bd.wen[i][j][l];
          }
          vuni_d(w);
          did=itr_vp_s(&p,&dp,r,w,type,src);
          if(did!=0){
            printf("bem3_aw_b2_solve_bieq.c, add_scattered_vp(), scattered field domain id error!(overlapped). domain id=%d. Exit...\n",did);
            exit(1);
          }
          dst->bd. Pi[i][j]+= p;
          dst->bd.dPi[i][j]+=dp;
        }
      }
    }
    else { // bi-linear element
      for(j=0;j<4;j++){
        if(dst->bd.md[i]==0){
          for(l=0;l<3;l++){
            r[l]=dst->bd.ren[i][j][l];
            w[l]=dst->bd.wen[i][j][l];
          }
          vuni_d(w);
          did=itr_vp_s(&p,&dp,r,w,type,src);
          if(did!=0){
            printf("bem3_aw_b2_solve_bieq.c, add_scattered_vp(), scattered field domain id error!(overlapped). domain id=%d. Exit...\n",did);
            exit(1);
          }
          dst->bd. Pi[i][j]+= p;
          dst->bd.dPi[i][j]+=dp;
        }
      }
    }
  }
}

int itr_vp_s(double complex *p,double complex *dp,double *rt,double *nv,int type,DMDA *md)
{
  int domain_id_m_dmda(double *rt,DMDA *ad);
  
  double complex CC[9],dCC[9],kc;
  double F,dF;
  int did,s,sd,l,n;

  did=domain_id_m_dmda(rt,md);

  *p=0.0;
  *dp=0.0;

  F=0.0;
  dF=0.0;
  kc=md->k0[did];

  for(s=1;s<=md->bd.sb[did].Ne;s++){
    sd=md->bd.sb[did].sid[s];
    dcoef_rt(CC,dCC,rt,nv,sd,kc,type,&(md->bd));

    for(n=0;n<4;n++){
      *p += CC[n+0]*md->bd.sb[did].dP[s][n]- CC[n+4]*md->bd.sb[did].P[s][n];
      *dp+=dCC[n+0]*md->bd.sb[did].dP[s][n]-dCC[n+4]*md->bd.sb[did].P[s][n];
    }

    F+=creal(CC[8]);
    dF+=creal(dCC[8]);
  }

  if(did==0){
    for(l=0;l<4;l++){
      *p/=1.0+F;
      *dp=((*dp)-(*p)*dF)/(1.0+F);
    }
  }
  else{
    for(l=0;l<4;l++){
      *p/=F;
      *dp=((*dp)-(*p)*dF)/F;
    }
  }

  return did;
}

double ccd_f(MOBJ *mo)
{
  double complex tc,tp,dt;
  double sum,bc,bp,db,ds;
  int i,s,p,sc;

  sum=0.0;
  sc=0;

  for(i=0;i<mo->N;i++){
    for(s=1;s<=mo->md[i].bd.Ne;s++){
      if(mo->md[i].bd.md[s]==0){
        for(p=0;p<4;p++){
          tc=mo->md[i].bd. Pi[s][p];
          tp=mo->md[i].bd.Pip[s][p];
          dt=tc-tp;

          bc=creal(tc*conj(tc));
          bp=creal(tp*conj(tp));
          db=creal(dt*conj(dt));
          ds=2.0*db/(bc+bp);
          sum+=ds;

          sc+=1;
        }
      }
    }
  }

  return sum/(double)sc;
}

void solve_coefficient_pv(int flg,DMDA *md)
{
  void solve_pv_bv (CMD *cm,DMDA *ad); // bem3_aw_b1_solve_bieq.c
  void solve_dpv_bv(CMD *cm,DMDA *ad); // bem3_aw_b1_solve_bieq.c

  time_t ms,me;

  if(flg==1){ printf("  solve PV boundary value   "); fflush(stdout); }
  time(&ms);
  solve_pv_bv(&(md->cm),md);
  solve_dpv_bv(&(md->cm),md);
  time(&me);
  if(flg==1) printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));
}
