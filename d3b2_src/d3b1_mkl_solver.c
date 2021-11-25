/*
 * b_mkl_solver.c
 *
 *  Created on: Jan 31, 2019
 *      Author: ohta
 */

#include "d3b1_const.h"


void d3b1_mkl_pardiso_create_matrixd(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,
        MKL_Complex16 *A,MKL_INT *xa,MKL_INT *asub,MKL_Complex16 *B); // double precision
void d3b1_mkl_pardiso_create_matrixs(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,
        MKL_Complex8  *A,MKL_INT *xa,MKL_INT *asub,MKL_Complex8  *B); // single precision

void d3b1_mkl_lapacke_create_matrix(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,
        MKL_Complex16 *A,MKL_Complex16 *B);


int d3b1_mkl_solver_pardiso_d(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x)
{
  MKL_Complex16 *A,*B,*X,cdum;
  MKL_INT N,mtype,nrhs,i,j,idum,*xa,*asub;
  void *pt[64]; // internal solver memory pointer
  // pardiso control param
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;

  // matrix setting
  N=(MKL_INT)n; // matrix size
  mtype=13; // matrix type, complex unsymmetric matrix

  // matrix a
  A=(MKL_Complex16 *)m_alloc2(nnz,sizeof(MKL_Complex16),"b_mkl_solver(),A");
  xa=(MKL_INT *)m_alloc2(n+1,sizeof(MKL_INT),"b_mkl_solver(),xa");
  asub=(MKL_INT *)m_alloc2(nnz,sizeof(MKL_INT),"b_mkl_solver(),asub");

  // rhs
  nrhs=1;
  B=(MKL_Complex16 *)m_alloc2(n,sizeof(MKL_Complex16),"b_mkl_solver(),B");
  X=(MKL_Complex16 *)m_alloc2(n,sizeof(MKL_Complex16),"b_mkl_solver(),X");

  // read matrix data from file
  d3b1_mkl_pardiso_create_matrixd(n,nnz,fn_a,fn_xa,fn_asub,fn_b,A,xa,asub,B);

  // Setup Pardiso control parameters
  for(i=0;i<64;i++) iparm[i]=0;

  iparm[0] = 1;         /* 1:No solver default */
  iparm[1] = 0;         /* 0:fill-in reordering from the minimum degree algorithm, 2: METIS, 3:parallel of 2 option */
  iparm[3] = 0;         /* No iterative-direct algorithm */
  iparm[4] = 0;         /* No user fill-in reducing permutation */
  iparm[5] = 0;         /* Write solution into x */
  iparm[6] = 0;         /* Not in use */
  iparm[7] = 2;         /* Max numbers of iterative refinement steps */
  iparm[8] = 0;         /* Not in use */
  iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;        /* Conjugate transposed/transpose solve */
  iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
  iparm[13] = 0;        /* Output: Number of perturbed pivots */
  iparm[14] = 0;        /* Not in use */
  iparm[15] = 0;        /* Not in use */
  iparm[16] = 0;        /* Not in use */
  iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1;       /* Output: Mflops for LU factorization */
  iparm[19] = 0;        /* Output: Numbers of CG Iterations */
  iparm[59] =0;         /* 0:in-core mode,1:variable, 2:out of core mode, required config file. */

  maxfct = 1;           /* Maximum number of numerical factorizations.  */
  mnum = 1;         /* Which factorization to use. */

  msglvl = PRDISO_STAT;    /* 1: Print statistical information  */
  error = 0;            /* Initialize error flag */

  //Initialize the internal solver memory pointer
  for(i=0;i<64;i++) pt[i]=0;

  //Reordering and Symbolic Factorization
  phase = 11;

  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      &N, A, xa, asub, &idum, &nrhs, iparm, &msglvl, &cdum, &cdum, &error);

  if (error != 0){
    printf ("\nERROR during symbolic factorization: %lld", error);
    exit (1);
  }
  if(msglvl==1){
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors  = %lld", iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %lld", iparm[18]);
    printf("\n");
  }

  // numerical factorization
  phase = 22;

  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
       &N, A, xa, asub, &idum, &nrhs, iparm, &msglvl, &cdum, &cdum, &error);

  if (error != 0){
    printf ("\nERROR during numerical factorization: %lld", error);
    exit (2);
  }
  if(msglvl==1) printf ("\nFactorization completed ...\n ");

  // Back substitution and iterative refinement
  phase = 33;
  // init X
  for (j = 0; j < N; j++){
    X[j].real = 0.0;
    X[j].imag = 0.0;
  }

  if(msglvl==1) printf ("\n\nSolving the system...\n");
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      &N, A, xa, asub, &idum, &nrhs, iparm, &msglvl, B, X, &error);

  if (error != 0){
    printf ("\nERROR during solution: %lld", error);
    exit (3);
  }

  if(msglvl==1){
    printf ("\nThe solution of the system is: \n");
    for (j = 10; j < 15; j++){
      printf ("X [%lld] = % 15.14e % 15.14e\n", j, X[j].real, X[j].imag);
    }
    printf ("\n");
  }

  // termination
  phase = -1;           /* Release internal memory. */

  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
       &N, &cdum, xa, asub, &idum, &nrhs,
       iparm, &msglvl, &cdum, &cdum, &error);


  for(i=0;i<n;i++) x[i]=X[i].real+X[i].imag*I;

  free(A);
  free(xa);
  free(asub);

  free(B);
  free(X);
  return (int)error;
}

int d3b1_mkl_solver_pardiso_s(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x)
{
  MKL_Complex8 *A,*B,*X,cdum;
  MKL_INT N,mtype,nrhs,i,j,idum,*xa,*asub;
  void *pt[64]; // internal solver memory pointer
  // pardiso control param
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;

  // matrix setting
  N=(MKL_INT)n; // matrix size
  mtype=13; // matrix type, complex unsymmetric matrix

  // matrix a
  A=(MKL_Complex8 *)m_alloc2(nnz,sizeof(MKL_Complex8),"b_mkl_solver(),A");
  xa=(MKL_INT *)m_alloc2(n+1,sizeof(MKL_INT),"b_mkl_solver(),xa");
  asub=(MKL_INT *)m_alloc2(nnz,sizeof(MKL_INT),"b_mkl_solver(),asub");

  // rhs
  nrhs=1;
  B=(MKL_Complex8 *)m_alloc2(n,sizeof(MKL_Complex8),"b_mkl_solver(),B");
  X=(MKL_Complex8 *)m_alloc2(n,sizeof(MKL_Complex8),"b_mkl_solver(),X");

  // read matrix data from file
  d3b1_mkl_pardiso_create_matrixs(n,nnz,fn_a,fn_xa,fn_asub,fn_b,A,xa,asub,B);

  // Setup Pardiso control parameters
  for(i=0;i<64;i++) iparm[i]=0;

  iparm[0] = 1;         /* 1:No solver default */
  iparm[1] = 0;         /* 0:fill-in reordering from the minimum degree algorithm, 2: METIS, 3:parallel of 2 option */
  iparm[3] = 0;         /* No iterative-direct algorithm */
  iparm[4] = 0;         /* No user fill-in reducing permutation */
  iparm[5] = 0;         /* Write solution into x */
  iparm[6] = 0;         /* Not in use */
  iparm[7] = 2;         /* Max numbers of iterative refinement steps */
  iparm[8] = 0;         /* Not in use */
  iparm[9] = 13;        /* Perturb the pivot elements  */
  iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;        /* Conjugate transposed/transpose solve */
  iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
  iparm[13] = 0;        /* Output: Number of perturbed pivots */
  iparm[14] = 0;        /* Not in use */
  iparm[15] = 0;        /* Not in use */
  iparm[16] = 0;        /* Not in use */
  iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1;       /* Output: Mflops for LU factorization */
  iparm[19] = 0;        /* Output: Numbers of CG Iterations */
  iparm[27] = 1;        /* 0:double precision, 1:single precision. */
  iparm[59] =0;         /* 0:in-core mode,1:variable, 2:out of core mode, required config file. */

  maxfct = 1;           /* Maximum number of numerical factorizations.  */
  mnum = 1;         /* Which factorization to use. */

  msglvl = PRDISO_STAT;     /* 1: Print statistical information  */
  error = 0;            /* Initialize error flag */

  //Initialize the internal solver memory pointer
  for(i=0;i<64;i++) pt[i]=0;

  //Reordering and Symbolic Factorization
  phase = 11;

  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      &N, A, xa, asub, &idum, &nrhs, iparm, &msglvl, &cdum, &cdum, &error);

  if (error != 0){
    printf ("\nERROR during symbolic factorization: %lld", error);
    exit (1);
  }
  if(msglvl==1){
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors  = %lld", iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %lld", iparm[18]);
    printf("\n");
  }

  // numerical factorization
  phase = 22;

  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
       &N, A, xa, asub, &idum, &nrhs, iparm, &msglvl, &cdum, &cdum, &error);

  if (error != 0){
    printf ("\nERROR during numerical factorization: %lld", error);
    exit (2);
  }
  if(msglvl==1) printf ("\nFactorization completed ...\n ");

  // Back substitution and iterative refinement
  phase = 33;
  // init X
  for (j = 0; j < N; j++){
    X[j].real = 0.0;
    X[j].imag = 0.0;
  }

  if(msglvl==1) printf ("\n\nSolving the system...\n");
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      &N, A, xa, asub, &idum, &nrhs, iparm, &msglvl, B, X, &error);

  if (error != 0){
    printf ("\nERROR during solution: %lld", error);
    exit (3);
  }

  if(msglvl==1){
    printf ("\nThe solution of the system is: \n");
    for (j = 10; j < 15; j++){
      printf ("X [%lld] = % 15.14e % 15.14e\n", j, X[j].real, X[j].imag);
    }
    printf ("\n");
  }

  // termination
  phase = -1;           /* Release internal memory. */

  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
       &N, &cdum, xa, asub, &idum, &nrhs,
       iparm, &msglvl, &cdum, &cdum, &error);


  for(i=0;i<n;i++) x[i]=X[i].real+X[i].imag*I;

  free(A);
  free(xa);
  free(asub);

  free(B);
  free(X);
  return (int)error;
}


void d3b1_mkl_pardiso_create_matrixd(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,
        MKL_Complex16 *A,MKL_INT *xa,MKL_INT *asub,MKL_Complex16 *B)
{
  FILE *fp;
  size_t i;

  // matrix A
  if((fp=fopen(fn_a,"rb"))==NULL){     printf("b_mkl_pardiso_create_matrix(),Failed to open the %s file.\n",fn_a);    exit(1); }
  fread(A,sizeof(MKL_Complex16),nnz,fp);
  fclose(fp);

  // vector xa
  if((fp=fopen(fn_xa,"rb"))==NULL){     printf("b_mkl_pardiso_create_matrix(),Failed to open the %s file.\n",fn_xa);    exit(1); }
  fread(xa,sizeof(MKL_INT),n+1,fp);
  fclose(fp);
  // offset 1
  for(i=0;i<=n;i++) xa[i]+=1;

  // vector asub
  if((fp=fopen(fn_asub,"rb"))==NULL){     printf("b_mkl_pardiso_create_matrix(),Failed to open the %s file.\n",fn_asub);    exit(1); }
  fread(asub,sizeof(MKL_INT),nnz,fp);
  fclose(fp);
  // offset 1
  for(i=0;i<nnz;i++) asub[i]+=1;

  // vector b
  if((fp=fopen(fn_b,"rb"))==NULL){     printf("b_mkl_pardiso_create_matrix(),Failed to open the %s file.\n",fn_b);    exit(1); }
  fread(B,sizeof(MKL_Complex16),n,fp);
  fclose(fp);
}

void d3b1_mkl_pardiso_create_matrixs(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,
        MKL_Complex8 *A,MKL_INT *xa,MKL_INT *asub,MKL_Complex8 *B)
{
  FILE *fp;
  double complex tc;
  size_t i;

  // matrix A
  if((fp=fopen(fn_a,"rb"))==NULL){     printf("b_mkl_pardiso_create_matrix(),Failed to open the %s file.\n",fn_a);    exit(1); }
  for(i=0;i<nnz;i++){
    fread(&tc,sizeof(double complex),1,fp);
    A[i].real=creal(tc);
    A[i].imag=cimag(tc);
  }
  fclose(fp);

  // vector xa
  if((fp=fopen(fn_xa,"rb"))==NULL){     printf("b_mkl_pardiso_create_matrix(),Failed to open the %s file.\n",fn_xa);    exit(1); }
  fread(xa,sizeof(MKL_INT),n+1,fp);
  fclose(fp);
  // offset 1
  for(i=0;i<=n;i++) xa[i]+=1;

  // vector asub
  if((fp=fopen(fn_asub,"rb"))==NULL){     printf("b_mkl_pardiso_create_matrix(),Failed to open the %s file.\n",fn_asub);    exit(1); }
  fread(asub,sizeof(MKL_INT),nnz,fp);
  fclose(fp);
  // offset 1
  for(i=0;i<nnz;i++) asub[i]+=1;

  // vector b
  if((fp=fopen(fn_b,"rb"))==NULL){     printf("b_mkl_pardiso_create_matrix(),Failed to open the %s file.\n",fn_b);    exit(1); }
  for(i=0;i<n;i++){
    fread(&tc,sizeof(double complex),1,fp);
    B[i].real=creal(tc);
    B[i].imag=cimag(tc);
  }
  fclose(fp);
}

int d3b1_mkl_solver_lapacke(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x)
{
  MKL_Complex16 *A,*B;
  MKL_INT N,nrhs,lda,ldb,info,*ipiv,i;

  N=(MKL_INT)n;
  nrhs=1;
  lda=N;
  ldb=1;
  ipiv=(MKL_INT*)malloc(sizeof(MKL_INT)*n);

  A=(MKL_Complex16 *)m_alloc2(n*n,sizeof(MKL_Complex16),"b_mkl_solver_lapacke(),A");
  B=(MKL_Complex16 *)m_alloc2(n,sizeof(MKL_Complex16),"b_mkl_solver_lapacke(),B");

  d3b1_mkl_lapacke_create_matrix(n,nnz,fn_a,fn_xa,fn_asub,fn_b,A,B);

  info = LAPACKE_zgesv( LAPACK_ROW_MAJOR, n, nrhs, A, lda, ipiv, B, ldb );

  for(i=0;i<n;i++) x[i]=B[i].real+B[i].imag*I;

  free(ipiv);
  free(A);
  free(B);

  return (int)info;
}

void d3b1_mkl_lapacke_create_matrix(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,
        MKL_Complex16 *A,MKL_Complex16 *B)
{
  FILE *fa,*fxa,*fas;
  MKL_Complex16 tc;
  MKL_INT is,ie,i,j,p;

  // matrix A
  if((fa=fopen(fn_a,"rb"))==NULL){     printf("b_mkl_lapacke_create_matrix(), Failed to open the %s file.\n",fn_a);    exit(1); }
  if((fxa=fopen(fn_xa,"rb"))==NULL){     printf("b_mkl_lapacke_create_matrix(), Failed to open the %s file.\n",fn_xa);    exit(1); }
  if((fas=fopen(fn_asub,"rb"))==NULL){     printf("b_mkl_lapacke_create_matrix(), Failed to open the %s file.\n",fn_asub);    exit(1); }

  fread(&is,sizeof(MKL_INT),1,fxa);
  fread(&ie,sizeof(MKL_INT),1,fxa);
  for(j=0;j<n;j++){
    for(p=is;p<ie;p++){
      fread(&i,sizeof(MKL_INT),1,fas);
      fread(&tc,sizeof(MKL_Complex16),1,fa);
      A[j*n+i]=tc;
    }
    is=ie;
    fread(&ie,sizeof(MKL_INT),1,fxa);
  }
  fclose(fa);
  fclose(fxa);
  fclose(fas);

  // vector B
  if((fa=fopen(fn_b,"rb"))==NULL){     printf("b_mkl_lapacke_create_matrix(), Failed to open the %s file.\n",fn_b);    exit(1); }
  fread(B,sizeof(MKL_Complex16),n,fa);
  fclose(fa);
}


