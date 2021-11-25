#if !defined DE_INT_H
#define DE_INT_H

#include <math.h>
#include <complex.h>

#define mmax 256
#define efs  0.1
#define hoff 10.0

double complex deintz(double complex (*f)(double,void *),double a,double b,void *pa, double eps,double *err); 
double         deintd(double         (*f)(double,void *),double a,double b,void *pa, double eps,double *err); 
/* integral of f(x) over (a,b) 
 f         : function pointer of integrand 
 a         : lower limit of integration
 b         : upper limit of integration
 pa        : pointer of parameter list 
 eps       : relative error requested 
 *err      : estimation of the absolute error 
 
 f(x) needs to be analytic over (a,b).
 error message
   err >= 0 : normal termination. returned value is estimated abosolute error. 
   err < 0  : abnormal termination (exceed mmax). 
*/

double complex deintiz(double complex (*f)(double,void *),double a,void *pa,double eps,double *err);
double         deintid(double         (*f)(double,void *),double a,void *pa,double eps,double *err);
/* integral of f(x) over (a,infinity) */

#endif
