/* focused beam
 * aw_fb.h
 *
 *  Created on: Jan 16, 2020
 *      Author: ohta
 */
#ifndef AW_FB_H_
#define AW_FB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "my_utils.h"
#include "gauleg.h"

// Transducer thickness to avoid singularities 
#define TD_TH 5.0e-3     // [m]

// trapezoidal rule settings
#define TRAP_MNUM 12     // maximum step of composite trapezoidal rule
#define TRAP_LNUM 2      // minimum step of composite trapezoidal rule
#define TRAP_EPS 1.0e-15 // for convergence criterion

typedef struct aw_focusedbeam_adata{
  double omega;         // angular frequency [rad/s], omega=2*M_PI*f.
  double K0;            // bulk modulus [N/m^2], K0 = c0^2*rho0.
  double complex k1,k2; // k1=I*omega/K_0, k2=I*omega*rho0, k0^2=-k1*k2.
  double k0;            // wavenumber, k0=omega/c0.
  double lambda0;       // wavelength.
  double alpha;         // alpha=sin^{-1}(B/F).
  double R[9];          // rotation matrix.
  
  double *wt,*ct,*st;   // gauss-legendre quadrature nodes and weights.
  double *cp,*sp;       // nodes for trapezoidal rule. 
}AfbD;

typedef struct aw_focusedbeam{
  double rho0;          // medium density [kg/m^3].
  double c0;            // speed of sound in the medium [m/s].
  double f;             // frequency [Hz].
  double complex p;     // amplitude including initial phase.
  double B;             // radius of focused transducer.
  double F;             // curvature radius (focal length) of the transducer. 
  double tv[3];         // translation vector, tv[0]=x, tv[1]=y, tv[2]=z.
  double theta,psi;     // parameter for rotation of acoustic axis.
  int nn;               // sampling number for Gauss-Legendre quadrature.

  AfbD wd;
}Afb;

void read_data_afb(char *fname, Afb *pw); // read datafile
void print_data_afb(Afb *pw);             // print aw_focusedbeam data
void setup_afb(Afb *pw);                  // memory allocation and initialize Afb
void free_afb(Afb *pw);                   // memory free

int calc_afb_pv(double complex *p,double complex *v,double *r,Afb *pw);
// input 
// *r  : coordinate of calculation point, r[0]=x, r[1]=y, r[2]=z.
// *pw : pointer of the Apw object.
// output
// *p : sound pressure.
// *v : particle velocity, v[0]=v_x, v[1]=v_y, v[2]=v_z. 
// return 
// 1 : inside the transducer, 0 : outside the transducer.

int calc_afb_dpdn(double complex *p,double complex *dpdn,double *r,double *n,Afb *pw);
// input 
// *r  : coordinate of calculation point, r[0]=x, r[1]=y, r[2]=z.
// *n  : direction for directional derivative, n[0]=n_x, n[1]=n_y, n[2]=n_z, |n|=1.
// *pw : pointer of the Apw object.
// output
// *p  : sound pressure.
// *dpdn : n-directional derivative of p, dpdn = nabla p * n.
// return 
// 1 : inside the transducer, 0 : outside the transducer.

#endif
