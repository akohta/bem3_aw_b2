#include "bem3_aw_b2.h"


void solve_bieq_dmda(DMDA *ad)
{
  void initalize_cmd(CMD *cm);
  void finalize_cmd(CMD *cm);
  void create_cmatrix(CMD *cm,DMDA *ad); // coefficient matrix
  void create_tmatrix_csr(CMD *cm,DMDA *ad); // create total matrix, Compressed Sparse Row(CSR) data format
  void solve_tmatrix_csr(CMD *cm,DMDA *ad);
  void solve_pv_bv(CMD *cm,DMDA *ad);
  void solve_dpv_bv(CMD *cm,DMDA *ad);

  time_t start,end,ms,me;
  CMD cm;

  printf("\nsolve acoustic wave boundary value \n");
  time(&start);

  // all coefficient matrix
  printf("  coefficient matrix          "); fflush(stdout);
  time(&ms);
  cm.type=2; // precision setting 0:4p GL,1:9p or 7p(triangular) GL, 2:GLN p GL, 3: GHN p GL
  cm.MN=ad->MN;
  initalize_cmd(&cm);
  create_cmatrix(&cm,ad);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  // velocity potential boundary value
  printf("  solve VP boundary value     "); fflush(stdout);
  time(&ms);
  cm.nn=(size_t)ad->bd.NN;
  cm.na=cm.nn*2;
  create_tmatrix_csr(&cm,ad);
  solve_tmatrix_csr(&cm,ad);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  // particle velocity boundary value
  printf("  solve PV boundary value     "); fflush(stdout);
  time(&ms);
  solve_pv_bv(&cm,ad);
  solve_dpv_bv(&cm,ad);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  finalize_cmd(&cm);
  time(&end);
  printf("Total elapsed time : %g (sec)\n",difftime(end,start));
}

////////////////////////////////////////////////////
void initalize_cmd(CMD *cm)
{
  int i,MN;

  MN=cm->MN;
  cm->tgfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq.c, initialize_cmd(), cm->tgfn");
  cm->thfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq.c, initialize_cmd(), cm->thfn");
  cm->tdgfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq.c, initialize_cmd(), cm->tdgfn");
  cm->tdhfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq.c, initialize_cmd(), cm->tdhfn");
  cm->tdffn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq.c, initialize_cmd(), cm->tdffn");
  for(i=0;i<=MN;i++){
    cm->tgfn[i]=(char *)m_alloc2(16,sizeof(char ),"solve_bieq.c, initialize_cmd(), cm->tgfn[i]");
    cm->thfn[i]=(char *)m_alloc2(16,sizeof(char ),"solve_bieq.c, initialize_cmd(), cm->thfn[i]");
    sprintf(cm->tgfn[i],"tmpG_%05d.dat",i);
    sprintf(cm->thfn[i],"tmpH_%05d.dat",i);
    cm->tdgfn[i]=(char *)m_alloc2(16,sizeof(char ),"solve_bieq.c, initialize_cmd(), cm->tdgfn[i]");
    cm->tdhfn[i]=(char *)m_alloc2(16,sizeof(char ),"solve_bieq.c, initialize_cmd(), cm->tdhfn[i]");
    cm->tdffn[i]=(char *)m_alloc2(16,sizeof(char ),"solve_bieq.c, initialize_cmd(), cm->tdffn[i]");
    sprintf(cm->tdgfn[i],"tmpdG_%05d.dat",i);
    sprintf(cm->tdhfn[i],"tmpdH_%05d.dat",i);
    sprintf(cm->tdffn[i],"tmpdF_%05d.dat",i);
  }

  cm->aval=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq.c, initialize_cmd(), cm->aval");
  cm->aptr=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq.c, initialize_cmd(), cm->aptr");
  cm->aidx=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq.c, initialize_cmd(), cm->aidx");
  sprintf(cm->aval,"tmpAval.dat");
  sprintf(cm->aptr,"tmpAprt.dat");
  sprintf(cm->aidx,"tmpAidx.dat");
  cm->b=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq.c, initialize_cmd(), cm->b");
  sprintf(cm->b,"tmpB.dat");
}

void finalize_cmd(CMD *cm)
{
  int i;

  // delete temporary file
  for(i=0;i<=cm->MN;i++){
    remove(cm->tgfn[i]);
    remove(cm->thfn[i]);
    remove(cm->tdgfn[i]);
    remove(cm->tdhfn[i]);
    remove(cm->tdffn[i]);
  }
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

  cm->MN=0;

  free(cm->aval);
  free(cm->aptr);
  free(cm->aidx);
  cm->nn=0;
  cm->na=0;
  cm->nnz=0;
  free(cm->b);
}

void create_cmatrix(CMD *cm,DMDA *ad)
{
  void create_cmatrix_domain(int did,CMD *cm,DMDA *ad);

  int i;

  for(i=0;i<=ad->MN;i++) create_cmatrix_domain(i,cm,ad);
}

void create_cmatrix_domain(int did,CMD *cm,DMDA *ad)
{
  FILE *fg,*fh,*fdg,*fdh,*fdf;
  double complex *tG,*tH,*tdG,*tdH,CC[9],dCz[9],dCe[9],kc;
  double F,vtz[3],vte[3],dFz,dFe,*tdF;
  size_t Ne,N,t,s,tn,tl,i;
  int td,sd;

  Ne=(size_t)ad->bd.sb[did].Ne;
  N=4*Ne;

  tG=(double complex *)m_alloc2(N*N,sizeof(double complex),"solve_bieq.c, create_matrix_domain(),tG"); // malloc
  tH=(double complex *)m_alloc2(N*N,sizeof(double complex),"solve_bieq.c, create_matrix_domain(),tH"); // malloc
  if((fg=fopen(cm->tgfn[did],"wb"))==NULL){    printf("solve_bieq.c, create_matrix_domain(),*fg. Failed to create %s file.\n",cm->tgfn[did]);    exit(1);  }
  if((fh=fopen(cm->thfn[did],"wb"))==NULL){    printf("solve_bieq.c, create_matrix_domain(),*fh. Failed to create %s file.\n",cm->thfn[did]);    exit(1);  }

  tdG=(double complex *)m_alloc2(2*N*N,sizeof(double complex),"solve_bieq.c, create_matrix_domain(),tdG"); // malloc
  tdH=(double complex *)m_alloc2(2*N*N,sizeof(double complex),"solve_bieq.c, create_matrix_domain(),tdH"); // malloc
  tdF=(double *)m_alloc2(3*N,sizeof(double),"solve_bieq.c, create_matrix_domain(),tdF"); // malloc
  if((fdg=fopen(cm->tdgfn[did],"wb"))==NULL){    printf("solve_bieq.c, create_matrix_domain(),*fdg. Failed to create %s file.\n",cm->tdgfn[did]);    exit(1);  }
  if((fdh=fopen(cm->tdhfn[did],"wb"))==NULL){    printf("solve_bieq.c, create_matrix_domain(),*fdh. Failed to create %s file.\n",cm->tdhfn[did]);    exit(1);  }
  if((fdf=fopen(cm->tdffn[did],"wb"))==NULL){    printf("solve_bieq.c, create_matrix_domain(),*fdf. Failed to create %s file.\n",cm->tdffn[did]);    exit(1);  }


  kc=(double complex)ad->k0[did];
#pragma omp parallel for schedule(dynamic) private(td,tl,tn,vtz,vte,F,dFz,dFe,s,sd,CC,dCz,dCe,i)
  for(t=1;t<=Ne;t++){
    td=ad->bd.sb[did].sid[t];
    if( ELT3==check_element_type(td,&(ad->bd)) ) tl=3;
    else tl=4;

    for(tn=0;tn<tl;tn++){
      tz_te_bd_node(vtz,vte,td,tn,&(ad->bd));

      F=0.0;
      dFz=0.0;
      dFe=0.0;
      for(s=1;s<=Ne;s++){
        sd=ad->bd.sb[did].sid[s];

        dcoef_bd_node_t2(CC,dCz,dCe,td,tn,vtz,vte,sd,kc,cm->type,&(ad->bd)); // precision check with derivative coefficient

        for(i=0;i<4;i++){
          tG[(t-1)*4*N+tn*N+(s-1)*4+i]=CC[i+0];
          tH[(t-1)*4*N+tn*N+(s-1)*4+i]=CC[i+4];

          tdG[(t-1)*2*4*N+tn*2*N+N*0+(s-1)*4+i]=dCz[i+0];
          tdH[(t-1)*2*4*N+tn*2*N+N*0+(s-1)*4+i]=dCz[i+4];
          tdG[(t-1)*2*4*N+tn*2*N+N*1+(s-1)*4+i]=dCe[i+0];
          tdH[(t-1)*2*4*N+tn*2*N+N*1+(s-1)*4+i]=dCe[i+4];
        }
        F+=creal(CC[8]);
        dFz+=creal(dCz[8]);
        dFe+=creal(dCe[8]);
      }

      if(did==0){
        tH[(t-1)*4*N+tn*N+(t-1)*4+tn]+=1.0+F;
        tdF[(t-1)*4*3+tn*3+0]=1.0+F;
      }
      else {
        tH[(t-1)*4*N+tn*N+(t-1)*4+tn]+=F;
        tdF[(t-1)*4*3+tn*3+0]=F;
      }
      tdF[(t-1)*4*3+tn*3+1]=dFz;
      tdF[(t-1)*4*3+tn*3+2]=dFe;
    }
  }

  fwrite(tG,sizeof(double complex),N*N,fg);
  fwrite(tH,sizeof(double complex),N*N,fh);
  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);

  fwrite(tdG,sizeof(double complex),2*N*N,fdg);
  fwrite(tdH,sizeof(double complex),2*N*N,fdh);
  fwrite(tdF,sizeof(double),3*N,fdf);
  fclose(fdg);
  fclose(fdh);
  fclose(fdf);
  free(tdG);
  free(tdH);
  free(tdF);
}

void create_tmatrix_csr(CMD *cm,DMDA *ad)
{
  void create_matrix_csr_dac(int did,FILE *av,FILE *ap,FILE *ai,FILE *b,CMD *cm,DMDA *ad);

  FILE *av,*ai,*ap,*b;

  size_t did;

  if((av=fopen(cm->aval,"wb"))==NULL){    printf("solve_bieq.c, create_tmatrix_csr(),*av. Failed to create %s file.\n",cm->aval);    exit(1);  }
  if((ai=fopen(cm->aidx,"wb"))==NULL){    printf("solve_bieq.c, create_tmatrix_csr(),*ai. Failed to create %s file.\n",cm->aidx);    exit(1);  }
  if((ap=fopen(cm->aptr,"wb"))==NULL){    printf("solve_bieq.c, create_tmatrix_csr(),*ap. Failed to create %s file.\n",cm->aidx);    exit(1);  }
  if((b=fopen(cm->b,"wb"))==NULL){    printf("solve_bieq.c, create_tmatrix_csr(),*b. Failed to create %s file.\n",cm->b);    exit(1);  }

  cm->nnz=0; // initialize nnz
  for(did=0;did<=cm->MN;did++) create_matrix_csr_dac(did,av,ap,ai,b,cm,ad);

  fwrite(&(cm->nnz),sizeof(size_t),1,ap);

  fclose(av);
  fclose(ai);
  fclose(ap);
  fclose(b);
}

void create_matrix_csr_dac(int did,FILE *av,FILE *ap,FILE *ai,FILE *b,CMD *cm,DMDA *ad)
{
  FILE *fg,*fh;

  double complex *tG,*tH,*tA,tB,k2m,k2s,hk2;
  size_t Ne,t,tn,s,asd,mdid,sdid,etype,l,cc,*ti;
  int td,sd;

  if((fg=fopen(cm->tgfn[did],"rb"))==NULL){    printf("solve_bieq.c, create_tmatrix_csr_dac(),*fg. Failed to open %s file.\n",cm->tgfn[did]);    exit(1);  }
  if((fh=fopen(cm->thfn[did],"rb"))==NULL){    printf("solve_bieq.c, create_tmatrix_csr_dac(),*fh. Failed to open %s file.\n",cm->thfn[did]);    exit(1);  }

  Ne=(size_t)ad->bd.sb[did].Ne;
  tG=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"solve_bieq.c, create_tmatrix_csr_dac(),tG"); // malloc
  tH=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"solve_bieq.c, create_tmatrix_csr_dac(),tH"); // malloc

  tA=(double complex *)m_alloc2(cm->na,sizeof(double complex),"solve_bieq.c, create_tmatrix_csr_dac(),tA"); // malloc
  ti=(size_t *)m_alloc2(cm->na,sizeof(size_t),"solve_bieq.c, create_tmatrix_csr_dac(),ti"); // malloc

  for(t=1;t<=Ne;t++){
    td=ad->bd.sb[did].sid[t];

    for(tn=0;tn<4;tn++){
      fread(tG,sizeof(double complex),Ne*4,fg);
      fread(tH,sizeof(double complex),Ne*4,fh);
      if( tn==3 && ELT3==check_element_type(td,&(ad->bd)) )  continue;

      fwrite(&(cm->nnz),sizeof(size_t),1,ap); // write A pointer
      for(l=0;l<cm->na;l++) tA[l]=0.0;
      tB=0.0;

      for(s=1;s<=Ne;s++){
        sd=ad->bd.sb[did].sid[s]; // signed element id
        asd=(size_t)abs(sd);
        mdid=(size_t)ad->bd.md[asd]; // main domain id
        sdid=(size_t)ad->bd.sd[asd]; // sub domain id
        etype=check_element_type(sd,&(ad->bd));        

        if(did==mdid){ // main domain
          if(etype==ELT3){ // linear-triangular element
            for(l=0;l<3;l++) {
              tA[ cm->nn*0 + ad->bd.eni[asd][l] ]= tH[(s-1)*4+l];
              tA[ cm->nn*1 + ad->bd.eni[asd][l] ]=-tG[(s-1)*4+l];
            }
          }
          else { // bi-linear element
            for(l=0;l<4;l++) {
              tA[ cm->nn*0 + ad->bd.eni[asd][l] ]= tH[(s-1)*4+l];
              tA[ cm->nn*1 + ad->bd.eni[asd][l] ]=-tG[(s-1)*4+l];
            }
          }
        } // end main domain
        else { // sub domain
          k2m=ad->k2[mdid];
          k2s=ad->k2[sdid];
          hk2=k2m/k2s;
          if(etype==ELT3){ // linear-triangular element
            for(l=0;l<3;l++){
              tA[ cm->nn*0 + ad->bd.eni[asd][l] ]=hk2*tH[(s-1)*4+l];
              tA[ cm->nn*1 + ad->bd.eni[asd][l] ]=    tG[(s-1)*4+l];
              if(mdid==0) tB+=-hk2*tH[(s-1)*4+l]*ad->bd.Pi[asd][l]-tG[(s-1)*4+l]*ad->bd.dPi[asd][l];
            }
          }
          else { // bi-linear element
            for(l=0;l<4;l++){
              tA[ cm->nn*0 + ad->bd.eni[asd][l] ]=hk2*tH[(s-1)*4+l];
              tA[ cm->nn*1 + ad->bd.eni[asd][l] ]=    tG[(s-1)*4+l];
              if(mdid==0) tB+=-hk2*tH[(s-1)*4+l]*ad->bd.Pi[asd][l]-tG[(s-1)*4+l]*ad->bd.dPi[asd][l];
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

      fwrite(&tB,sizeof(double complex),1,b);

    } // end for tn
  } // end for t

  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);
  free(tA);
  free(ti);
}

void solve_tmatrix_csr(CMD *cm,DMDA *ad)
{
  int d3b1_mkl_solver_pardiso_d(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x); // b_mkl_solver.c
  int d3b1_mkl_solver_pardiso_s(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x); // b_mkl_solver.c
  void tmatrix_bd_store(double complex *X,DMDA *ad);

  int err;
  double complex *x; // results
  x=(double complex *)m_alloc2(cm->na,sizeof(double complex),"solve_bieq.c, solve_tmatrix_csr(),x");

  if(PARDISO_PREC==0) err=d3b1_mkl_solver_pardiso_d(cm->na,cm->nnz,cm->aval,cm->aptr,cm->aidx,cm->b,x);
  else err=d3b1_mkl_solver_pardiso_s(cm->na,cm->nnz,cm->aval,cm->aptr,cm->aidx,cm->b,x);
  if(err!=0){
    printf("solve_bieq.c, solve_tmatrix_csr(), b_mkl_solver_pardiso() returned error. err=%d. Exit...\n",err);
    exit(1);
  }

/*
  // for test using lapacke
  int i;
  int b_mkl_solver_lapacke(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x); // b_mkl_solfer.c
  err=b_mkl_solver_lapacke(cm->na,cm->nnz,cm->aval,cm->aptr,cm->aidx,cm->b,x);
  for(i=10;i<15;i++){
    printf("x[%d]=%15.14e + %15.14e\n",i,creal(x[i]),cimag(x[i]));
  }
*/

  // store data
  tmatrix_bd_store(x,ad);

  free(x);
}

void tmatrix_bd_store(double complex *X,DMDA *ad)
{
  double complex k2m,k2s,hk2;
  size_t d,s,l,nn;
  int sd,asd,etype,mdid,sdid;

  nn=ad->bd.NN;

  for(d=0;d<=ad->MN;d++){
    for(s=1;s<=ad->bd.sb[d].Ne;s++){
      sd=ad->bd.sb[d].sid[s];
      asd=abs(sd);
      mdid=ad->bd.md[asd];
      sdid=ad->bd.sd[asd];
      etype=check_element_type(sd,&(ad->bd));

      if(d==mdid){ // main domain
        if(etype==ELT3){ // linear triangular element
          for(l=0;l<3;l++){ // node id
            ad->bd.sb[d]. P[s][l]=X[nn*0+ad->bd.eni[asd][l]];
            ad->bd.sb[d].dP[s][l]=X[nn*1+ad->bd.eni[asd][l]];
          }
        }
        else { // bi-linear element
          for(l=0;l<4;l++){ // node id
            ad->bd.sb[d]. P[s][l]=X[nn*0+ad->bd.eni[asd][l]];
            ad->bd.sb[d].dP[s][l]=X[nn*1+ad->bd.eni[asd][l]];
          }
        }
      }
      else { // subdomain
        k2m=ad->k2[mdid];
        k2s=ad->k2[sdid];
        hk2=k2m/k2s;
        if(etype==ELT3){ // linear triangular element
          for(l=0;l<3;l++){ // node id
            ad->bd.sb[d]. P[s][l]= hk2*X[nn*0+ad->bd.eni[asd][l]];
            ad->bd.sb[d].dP[s][l]=-1.0*X[nn*1+ad->bd.eni[asd][l]];
            if(mdid==0){
              ad->bd.sb[d]. P[s][l]+= hk2*ad->bd. Pi[asd][l];
              ad->bd.sb[d].dP[s][l]+=-1.0*ad->bd.dPi[asd][l];
            }
          }
        }
        else { // bi-linear element
          for(l=0;l<4;l++){ // node id
            ad->bd.sb[d]. P[s][l]= hk2*X[nn*0+ad->bd.eni[asd][l]];
            ad->bd.sb[d].dP[s][l]=-1.0*X[nn*1+ad->bd.eni[asd][l]];
            if(mdid==0){
              ad->bd.sb[d]. P[s][l]+= hk2*ad->bd. Pi[asd][l];
              ad->bd.sb[d].dP[s][l]+=-1.0*ad->bd.dPi[asd][l];
            }
          }
        }
      } // end subdomain
    }  // end s
  } // end d
}

void solve_pv_bv(CMD *cm,DMDA *ad)
{
  void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte);
  void dpdt_fcm(double complex *dPz,double complex *dPe,int t,int tn,int Ne,
      double complex **P,double complex **dP,double complex *tdG,double complex *tdH,double *tdF);

  FILE *fdg,*fdh,*fdf;
  double complex *tdG,*tdH,dPz,dPe;
  double sig,tdF[3],vtz[3],vte[3],vn[3],a[9];
  int d,Ne,t,tn,at,atd,i;

  for(d=0;d<=ad->MN;d++){
    if((fdg=fopen(cm->tdgfn[d],"rb"))==NULL){    printf("solve_bieq.c, solve_pv_bv(),*fdg. Failed to open %s file.\n",cm->tdgfn[d]);    exit(1);  }
    if((fdh=fopen(cm->tdhfn[d],"rb"))==NULL){    printf("solve_bieq.c, solve_pv_bv(),*fdh. Failed to open %s file.\n",cm->tdhfn[d]);    exit(1);  }
    if((fdf=fopen(cm->tdffn[d],"rb"))==NULL){    printf("solve_bieq.c, solve_pv_bv(),*fdf. Failed to open %s file.\n",cm->tdffn[d]);    exit(1);  }

    Ne=ad->bd.sb[d].Ne;
    tdG=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"solve_bieq.c, solve_pv_bv(),tdG"); // malloc
    tdH=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"solve_bieq.c, solve_pv_bv(),tdH"); // malloc

    for(t=1;t<=Ne;t++){
      at=ad->bd.sb[d].sid[t];
      if(at<0) sig=-1.0;
      else sig=1.0;
      atd=abs(at);

      for(tn=0;tn<4;tn++){
        fread(tdG,sizeof(double complex),2*Ne*4,fdg);
        fread(tdH,sizeof(double complex),2*Ne*4,fdh);
        fread(tdF,sizeof(double),3,fdf);
        if( tn==3 && ELT3==check_element_type(atd,&(ad->bd)) )  continue;

        // tangential vector
        tz_te_bd_node(vtz,vte,atd,tn,&(ad->bd));
        // normal vector
        for(i=0;i<3;i++) vn[i]=sig*ad->bd.wen[atd][tn][i];
        vuni_d(vn);
        // inv-matrix
        inv_matrix_33(a,vn,vtz,vte);
        // dPdtz,dPdte
        dpdt_fcm(&dPz,&dPe,t,tn,Ne,ad->bd.sb[d].P,ad->bd.sb[d].dP,tdG,tdH,tdF);
        // particle velocity
        for(i=0;i<3;i++) ad->bd.sb[d].pv[t][tn][i]=-(ad->bd.sb[d].dP[t][tn]*a[i*3+0]+dPz*a[i*3+1]+dPe*a[i*3+2]);

      } // end for tn
    } // end for t
    fclose(fdg);
    fclose(fdh);
    fclose(fdf);
    free(tdG);
    free(tdH);
  } // end for d
}

void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte)
{
  double i_det;
  int i;

  i_det=1.0/(vn[0]*vtz[1]*vte[2]+vn[1]*vtz[2]*vte[0]+vn[2]*vtz[0]*vte[1]
          -( vn[2]*vtz[1]*vte[0]+vn[1]*vtz[0]*vte[2]+vn[0]*vtz[2]*vte[1]));

  ma[0*3+0]= vtz[1]*vte[2]-vtz[2]*vte[1];
  ma[0*3+1]=- vn[1]*vte[2]+ vn[2]*vte[1];
  ma[0*3+2]=  vn[1]*vtz[2]- vn[2]*vtz[1];

  ma[1*3+0]=-vtz[0]*vte[2]+vtz[2]*vte[0];
  ma[1*3+1]=  vn[0]*vte[2]- vn[2]*vte[0];
  ma[1*3+2]=- vn[0]*vtz[2]+ vn[2]*vtz[0];

  ma[2*3+0]= vtz[0]*vte[1]-vtz[1]*vte[0];
  ma[2*3+1]=- vn[0]*vte[1]+ vn[1]*vte[0];
  ma[2*3+2]=  vn[0]*vtz[1]- vn[1]*vtz[0];

  for(i=0;i<9;i++) ma[i]*=i_det;
}

void dpdt_fcm(double complex *dPz,double complex *dPe,int t,int tn,int Ne,double complex **P,double complex **dP,double complex *tdG,double complex *tdH,double *tdF)
{
  int s,sn;

  *dPz=0.0;
  *dPe=0.0;
  for(s=1;s<=Ne;s++){
    for(sn=0;sn<4;sn++){
      *dPz+=tdG[Ne*4*0+4*(s-1)+sn]*dP[s][sn]-tdH[Ne*4*0+4*(s-1)+sn]*P[s][sn];
      *dPe+=tdG[Ne*4*1+4*(s-1)+sn]*dP[s][sn]-tdH[Ne*4*1+4*(s-1)+sn]*P[s][sn];
    }
  }
  *dPz=(*dPz-P[t][tn]*tdF[1])/tdF[0];
  *dPe=(*dPe-P[t][tn]*tdF[2])/tdF[0];
}


void solve_dpv_bv(CMD *cm,DMDA *ad)
{
  FILE *fg,*fh;
  MKL_Complex16 *A,*B;
  double complex *tG,*tH;
  int Ne,d,t,atd,tn,nn,nc,nr,s,ns,asd,i;
  MKL_INT *ipiv,nrhs,lda,ldb,info;

  nrhs=3;

  for(d=0;d<=ad->MN;d++){
    Ne=ad->bd.sb[d].Ne;
    if(Ne==0) continue;

    if((fg=fopen(cm->tgfn[d],"rb"))==NULL){    printf("solve_bieq.c, solve_dpv_bv(),*fg. Failed to open %s file.\n",cm->tgfn[d]);    exit(1);  }
    if((fh=fopen(cm->thfn[d],"rb"))==NULL){    printf("solve_bieq.c, solve_dpv_bv(),*fh. Failed to open %s file.\n",cm->thfn[d]);    exit(1);  }

    // matrix size
    nn=0;
    for(s=1;s<=Ne;s++){
      asd=abs(ad->bd.sb[d].sid[s]);
      if( ELT3==check_element_type(asd,&(ad->bd)) ) nn+=3;
      else nn+=4;
    }

    lda=nn;
    ldb=nrhs;

    A =(MKL_Complex16 *)m_alloc2(nn*lda,sizeof(MKL_Complex16),"solve_bieq.c, solve_dpv_bv(),A"); // malloc
    B =(MKL_Complex16 *)m_alloc2(nn*ldb,sizeof(MKL_Complex16),"solve_bieq.c, solve_dpv_bv(),B"); // malloc
    ipiv=(MKL_INT *)m_alloc2(nn,sizeof(MKL_INT),"solve_bieq.c, solve_dpv_bv(),ipiv"); // malloc
    tG=(double complex *)m_alloc2(4*Ne,sizeof(double complex),"solve_bieq.c, solve_dpv_bv(),tG"); // malloc
    tH=(double complex *)m_alloc2(4*Ne,sizeof(double complex),"solve_bieq.c, solve_dpv_bv(),tH"); // malloc

    nc=0;
    for(t=1;t<=Ne;t++){
      atd=abs(ad->bd.sb[d].sid[t]);

      for(tn=0;tn<4;tn++){
        fread(tG,sizeof(double complex),Ne*4,fg);
        fread(tH,sizeof(double complex),Ne*4,fh);
        if( tn==3 && ELT3==check_element_type(atd,&(ad->bd)) )  continue;

        for(i=0;i<ldb;i++){
          B[ldb*nc+i].real=0.0;
          B[ldb*nc+i].imag=0.0;
        }
        nr=0;
        for(s=1;s<=Ne;s++){
          asd=abs(ad->bd.sb[d].sid[s]);
          if( ELT4==check_element_type(asd,&(ad->bd)) ){
            for(ns=0;ns<4;ns++) {
              A[nn*nc+nr].real=creal(tG[4*(s-1)+ns]);
              A[nn*nc+nr].imag=cimag(tG[4*(s-1)+ns]);
              nr+=1;
              for(i=0;i<3;i++){
                B[ldb*nc+i].real+=creal(tH[4*(s-1)+ns]*ad->bd.sb[d].pv[s][ns][i]);
                B[ldb*nc+i].imag+=cimag(tH[4*(s-1)+ns]*ad->bd.sb[d].pv[s][ns][i]);
              }
            }
          }
          else {
            for(ns=0;ns<3;ns++) {
              A[nn*nc+nr].real=creal(tG[4*(s-1)+ns]);
              A[nn*nc+nr].imag=cimag(tG[4*(s-1)+ns]);
              nr+=1;
              for(i=0;i<3;i++){
                B[ldb*nc+i].real+=creal(tH[4*(s-1)+ns]*ad->bd.sb[d].pv[s][ns][i]);
                B[ldb*nc+i].imag+=cimag(tH[4*(s-1)+ns]*ad->bd.sb[d].pv[s][ns][i]);
              }
            }
          }
        } // end for s
        nc+=1;
      } // end for tn

    } // end for t

    // solve mkl lapack
    info = LAPACKE_zgesv( LAPACK_ROW_MAJOR, nn, nrhs, A, lda, ipiv, B, ldb );
    if(info!=0){
      printf("solve_bieq.c, solve_dpv_bv(), LAPACKE_zgesv(). info=%lld error. Exit...\n",info);
      exit(1);
    }

    // store data
    nc=0;
    for(s=1;s<=Ne;s++){
      asd=abs(ad->bd.sb[d].sid[s]);
      if( ELT4==check_element_type(asd,&(ad->bd)) ){
        for(ns=0;ns<4;ns++){
          for(i=0;i<3;i++) {
            ad->bd.sb[d].dpv[s][ns][i]=B[ldb*nc+i].real+B[ldb*nc+i].imag*I;
          }
          nc+=1;
        }
      }
      else {
        for(ns=0;ns<3;ns++){
          for(i=0;i<3;i++) {
            ad->bd.sb[d].dpv[s][ns][i]=B[ldb*nc+i].real+B[ldb*nc+i].imag*I;
          }
          nc+=1;
        }
      }
    }

    fclose(fg);
    fclose(fh);
    free(A);
    free(B);
    free(ipiv);
    free(tG);
    free(tH);
  } // end for d
}
