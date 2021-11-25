/*
 * bem3_aw_b2.h
 *
 *  Created on: Nov 23, 2021
 *      Author: ohta
 */

#ifndef BEM3_AW_B2_H_
#define BEM3_AW_B2_H_

#include "d3b1_const.h"
#include "d3b1_elem.h"

#define ITER_MIN 10 // minimum number of iterative operations
#define ITER_MAX 50 // maximum number of iterative operations
#define PREC_DEF 2  // scattered field precision settings for iterative operation, 0:4p GL, 1:9p GL, 2:GLNp GL,3:GHNp GL
#define CCD 1.0e-7  // convergence condition determination parameter

typedef struct object_data{
  int N; // number of objects

  char **mod_fn; // object data file name
  double **rs;   // translation vector of each object

  DMDA *md;      // object data
}MOBJ;

// ----- create_matrix -----
// bem3_aw_b2.c
void read_dmda2(int argc,char **argv,DMDA *ad); // read datafile for single object.
void print_dmda2(DMDA *ad);                     // print data.
void initialize_dmda2(DMDA *ad);                // memory allocation and initialize coefficients.
void finalize_dmda2(DMDA *ad);                  // memory free.
// bem3_aw_b2_solve_bieq.c
void solve_bieq_dmda2(DMDA *ad,char *ofn);      // create coefficient matrix and its inverse matrix. 


// ----- solver -----
void mo_initialize(int argc,char **argv,MOBJ *mo);  // initialize multi-object data.
void mo_print_data(MOBJ *mo);                       // print data.
void mo_finalize(MOBJ *mo);                        // memory free.
void mo_dat_write(char *fname,MOBJ *mo);            // write datafile. 
void mo_dat_read(char *fname,MOBJ *mo);             // read datafile. 
void mo_output_node_particles(char *fname,MOBJ *mo);// outputs the nodes as point cloud data ( .particles file ) 
// bem3_aw_b2_solve_bieq.c
void mo_solve_iter(MOBJ *mo);                       // solve boundary integral equations.


// -- bem3_aw_b2_force.c --
int mo_force_FN(double *F,double *N,double *rc,int oid,int type,MOBJ *mo);


// -- bem3_aw_b2_field.c --
void mo_object_domain_id(int *oid,int *did,double *rt,MOBJ *mo);
int mo_phi_s(double complex *phi,double *rt,int type,MOBJ *mo); // scattered or internal field
int mo_phi_t(double complex *phi,double *rt,int type,MOBJ *mo); // total(incident + scattered) field
int mo_phi_i(double complex *phi,double *rt,int type,MOBJ *mo); // incident field
// return object id.
int mo_p_s(double complex *p,double *rt,int type,MOBJ *mo); // scattered or internal field
int mo_p_t(double complex *p,double *rt,int type,MOBJ *mo); // total field
int mo_p_i(double complex *p,double *rt,int type,MOBJ *mo); // incident field

int mo_pv_s(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo); // scattered or internal field
int mo_pv_t(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo); // total field
int mo_pv_i(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo); // incident field

int mo_pv_s_dbieq(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo); // for far field
int mo_pv_t_dbieq(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo); // 
int mo_pv_i_dbieq(double complex *p,double complex *pv,double *rt,int type,MOBJ *mo); //

void mo_pv_s_bd(double complex *p,double complex *pv,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo);
void mo_pv_t_bd(double complex *p,double complex *pv,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo);
void mo_pv_i_bd(double complex *p,double complex *pv,int oid,int did,int t,double zeta_t,double eta_t,int type,MOBJ *mo);

#endif
