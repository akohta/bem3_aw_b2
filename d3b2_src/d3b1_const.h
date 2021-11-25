/*
 * d3b1_const.h
 *
 *  Created on: Dec 14, 2018
 *      Author: ohta
 */

#ifndef BEM1_CONST_H_
#define BEM1_CONST_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "multi_aw.h"
#include "de_int.h"
#include "gauleg.h"
// mkl
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include "mkl_lapacke.h"


// epsilon
#define MEPS   1.0e-15 // machine epsilon
#define GEPS   1.0e-15 // geometrical precision,must be same as gmsh -tol (tolerance) option
#define IEPS   1.0e-12 // epsilon for DE integration
#define DEPS   0.1     // D_epsilon ( domain division coefficient )
#define CDFH   1.0e-3  // numerical difference H

// constant
#define I_TH    3.333333333333333333333e-1 // 1.0/3.0
#define I_FP    7.957747154594766788444e-2 // 1.0/(4.0*M_PI)
#define I_EP    3.978873577297383394222e-2 // 1.0/(8.0*M_PI)
#define SQ3     1.732050807568877293528e+0 // sqrt(3.0)
#define I_SPSQ3 3.062938307898844719507e-2 // 1.0/(6.0*M_PI*sqrt(3.0))

// gmsh setting
#define MSHVER  2.2 // version number : 2.2 later
#define MSHASCI 0   // file type is ascii
#define MSHPREC 8   // data size(precision), 8 = sizeof(double)
#define ELT3    2   // element type : 3-node triangle ( linear triangular element )
#define ELT4    3   // element type : 4-node quadrangle ( bi-linear element )
#define OPENDID 99  // open domain id 

// gaussian quadrature node and weight data
// 4node and 9node rule for quadrangular
#define P44_N  0.5773502691896257645091 // sqrt(1/3)
#define P44_W  1.0
#define P49_N  0.7745966692414833770359 // sqrt(3/5)
#define P49_W0 0.3086419753086419753086 // 25/81
#define P49_W1 0.4938271604938271604938 // 40/81
#define P49_W2 0.7901234567901234567901 // 64/81
// 4node and 7node rule for triangular
#define P34_N0 0.2                      // 1/5
#define P34_N1 0.3464101615137754587055 // sqrt(3)/5
#define P34_W0 0.5208333333333333333333 // 25/48
#define P34_W1 0.5625                   // 27/48
#define P37_N0 0.3480702390148154917985 // (sqrt(15)+1)/14
#define P37_N1 0.2052130961576726346557 // (sqrt(15)-1)/14
#define P37_N2 0.6028753385763033130539 // sqrt(3)*(sqrt(15)+1)/14
#define P37_N3 0.3554395089236065568357 // sqrt(3)*(sqrt(15)-1)/14
#define P37_W0 0.1259391805448271525957 // (155-sqrt(15))/1200
#define P37_W1 0.1323941527885061807376 // (155+sqrt(15))/1200
#define P37_W2 0.225                    // 270/1200

// gauss-legendre quadrature node(order) number
#define GLN 12
#define GHN 64

// for coefficient of boundary integral equation(BIEQ)
#define SW_BDR_BIEQ 0 // r-directional integration parameter.  0: GLN order GL, 1: GHN order GL
#define SW_BDT_BIEQ 0 // theta-directional integration parameter, 0: GHN order GL, 1: DE

// for coefficient of derivative boundary integral equation(DBIEQ)
#define CBD_CDF   100.0 // cut off parameter for near-boundary
#define CBD_SDF     0.5 // threshold for DBIEQ coefficient. if dF/dv > CBD_SDF, switch to GHN order GL
#define SW_DBIEQ      0 // 0:principal value method, 1:differential method
#define SW_DBDR_BDE   0 // r-directional integration parameter for domain bar{D_eps}. 0:GHN, 1:DE
#define SW_DBDR_DE    0 // r-directional integration parameter for domain D_eps. 0:GLN, 1:GHN
#define SW_DBDT_BDE   0 // theta-directional integration parameter for domain bar{D_eps}. 0: GHN, 1: DE
#define SW_DBDT_DE    0 // theta-directional integration (0~2pi) parameter for domain D_eps. 0: GHN, 1:DE

// for PARDISO
#define PARDISO_PREC  0 // 0:double precision complex, 1:single precision complex.
#define PRDISO_STAT   0 // 0:no print, 1:print statistical information


#endif /* BEM1_CONST_H_ */
