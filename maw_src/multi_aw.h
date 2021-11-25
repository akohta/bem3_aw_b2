/* Analysis of sound pressure and particle velocity of plane waves,
 * Bessel beams, and focused beams. 
 * 
 * multi_aw.h
 *
 *  Created on: Nov 01, 2021
 *      Author: ohta
 */
#ifndef MULTI_AW_H_
#define MULTI_AW_H_

#include "aw_pw.h"
#include "aw_bb.h"
#include "aw_fb.h"

// default datafile name
#define fn_pw "mpw.txt" // plane wave.
#define fn_bb "mbb.txt" // bessel beam.
#define fn_fb "mfb.txt" // focused beam.

typedef struct beam_data{
  Apw *pw;
  Abb *bb;
  Afb *fb;
}Bdata;

typedef struct beam_object{
  int n_pw;  char fname_pw[28]; // number of defined beams, readed beam datafile name.
  int n_bb;  char fname_bb[28]; 
  int n_fb;  char fname_fb[28];

  // common data
  double rho0;          // medium density [kg/m^3].
  double c0;            // speed of sound in the medium [m/s].
  double f;             // frequency [Hz].
  double lambda0;       // wavelength
  double k0;            // wavenumber
  double complex k1,k2; // k1=I*omega/K_0, k2=I*omega*rho0, k0^2=-k1*k2.
  
  Bdata bd;
}Maw;


void init_maw(Maw *aw);        // initialize the struct data.
void read_data_maw(Maw *aw);   // beam datafile is automatically searched and readed. search path is current directry.
void print_data_maw(Maw *aw);  // print defined beam data.
void setup_maw(Maw *aw);       // memory allocation and calculation of coefficients.
void  free_maw(Maw *aw);       // free allocated memory.


int calc_maw_pv(double complex *p,double complex *v,double *r,Maw *aw);
// input 
// *r  : coordinate of calculation point, r[0]=x, r[1]=y, r[2]=z.
// *aw : pointer of the Maw object.
// output
// *p : sound pressure.
// *v : particle velocity, v[0]=v_x, v[1]=v_y, v[2]=v_z. 
// return (only avairable for focused beam)
// 1< : inside the transducer, the number returned is the focused beam (transducer) identification number. 
// 0  : outside the transducer. 
// -1 : it is not a focused beam.

int calc_maw_dpdn(double complex *p,double complex *dpdn,double *r,double *n,Maw *aw);
// input 
// *r  : coordinate of calculation point, r[0]=x, r[1]=y, r[2]=z.
// *n  : direction for directional derivative, n[0]=n_x, n[1]=n_y, n[2]=n_z, |n|=1.
// *aw : pointer of the Maw object.
// output
// *p  : sound pressure.
// *dpdn : n-directional derivative of p, dpdn = nabla p * n.
// return (only avairable for focused beam)
// 1< : inside the transducer, the number returned is the focused beam (transducer) identification number. 
// 0  : outside the transducer.
// -1 : it is not a focused beam.


// manually set beam datafile. select proper function. returnd value is the number of defined beams.
int read_data_maw_pw(char *fname,Maw *aw); // for plane wave.
int read_data_maw_bb(char *fname,Maw *aw); // for Besse beam.
int read_data_maw_fb(char *fname,Maw *aw); // for focused beam.

#endif
