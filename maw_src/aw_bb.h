/* Bessel beam
 * aw_bb.h
 *
 *  Created on: Jan 16, 2020
 *      Author: ohta
 */
#ifndef AW_BB_H_
#define AW_BB_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_specfunc.h>


typedef struct aw_besselbeam_adata{
  double omega;         // angular frequency [rad/s], omega=2*M_PI*f.
  double K0;            // bulk modulus [N/m^2], K0 = c0^2*rho0.
  double complex k1,k2; // k1=I*omega/K_0, k2=I*omega*rho0, k0^2=-k1*k2.
  double k0;            // wavenumber, k0=omega/c0.
  double lambda0;       // wavelength.
  double kr;            // r-component of wave vector without rotation.
  double kz;            // z-component of wave vector without rotation. 
  double R[9];          // rotation matrix.
}AbbD;

typedef struct aw_besselbeam{
  double rho0;          // medium density [kg/m^3].
  double c0;            // speed of sound in the medium [m/s].
  double f;             // frequency [Hz].
  double complex p;     // amplitude including initial phase.
  double d_angle;       // deflection angle
  double tv[3];         // translation vector, tv[0]=x, tv[1]=y, tv[2]=z.
  double theta,psi;     // parameter for rotation of acoustic axis.

  AbbD wd;
}Abb;

void read_data_abb(char *fname, Abb *pw); // read datafile 
void print_data_abb(Abb *pw);             // print aw_besselbeam data
void setup_abb(Abb *pw);                  // setup Abb

void calc_abb_pv(double complex *p,double complex *v,double *r,Abb *pw);
// input 
// *r  : coordinate of calculation point, r[0]=x, r[1]=y, r[2]=z.
// *pw : pointer of the Apw object.
// output
// *p : sound pressure.
// *v : particle velocity, v[0]=v_x, v[1]=v_y, v[2]=v_z. 

void calc_abb_dpdn(double complex *p,double complex *dpdv, double *r ,double *n, Abb *pw);
// input 
// *r  : coordinate of calculation point, r[0]=x, r[1]=y, r[2]=z.
// *n  : direction for directional derivative, n[0]=n_x, n[1]=n_y, n[2]=n_z, |n|=1.
// *pw : pointer of the Apw object.
// output
// *p  : sound pressure.
// *dpdn : n-directional derivative of p, dpdn = nabla p * n.

#endif
