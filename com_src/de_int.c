#include "de_int.h"
#include <math.h>

double complex deintz(double complex (*f)(double,void *),double a,double b,void *pa,double eps,double *err)
{
  int m;
  double complex fa,fb,Ih,ir,irback,iback;
  double pi2, epsln, epsh, h0, ehp, ehm, epst, ba, h, t, ep, em,  errt, errh, errd;
  double xw,xa,wg;

  if(a==b){
    *err=0.0;
    return 0.0;
  }

  pi2 = 2.0*M_PI;
  epsln = 1.0 - log(efs * eps);
  epsh = sqrt(efs * eps);
  h0 = hoff / epsln;
  ehp = exp(h0);
  ehm = 1.0 / ehp;
  epst = exp(-ehm * epsln);
  ba = b - a;
  ir = (*f)((a + b) * 0.5,pa) * (ba * 0.25);
  Ih = ir * (2.0 * pi2);
  *err = cabs(Ih) * epst;
  h = 2.0 * h0;
  m = 1;
  do {
    iback = Ih;
    irback = ir;
    t = h * 0.5;
    do {
      em = exp(t);
      ep = pi2 * em;
      em = pi2 / em;
      do {
        xw = 1.0 / (1.0 + exp(ep - em));
        xa = ba * xw;
	wg = xa * (1.0 - xw);
        fa = (*f)(a + xa,pa) * wg;
        fb = (*f)(b - xa,pa) * wg;
        ir += fa + fb;
        Ih += (fa + fb) * (ep + em);
        errt = (cabs(fa) + cabs(fb)) * (ep + em);
        if (m == 1) *err += errt * epst;
        ep *= ehp;
        em *= ehm;
      } while (errt > *err || xw > epsh);
      t += h;
    } while (t < h0);
    if (m == 1) {
      errh = (*err / epst) * epsh * h0;
      errd = 1.0 + 2.0 * errh;
    } else {
      errd = h * (cabs(Ih - 2.0 * iback) + 4.0 * cabs(ir - 2.0 * irback));
    }
    h *= 0.5;
    m *= 2;
  } while (errd > errh && m < mmax);
  Ih *= h;
  if (errd > errh) {
    *err = -errd * (double)m;
  } else {
    *err = errh * epsh * (double)m / (2.0 * efs);
  }
  return Ih;
}

double deintd(double (*f)(double,void *),double a,double b,void *pa,double eps,double *err)
{
  int m;
  double Ih;
  double pi2, epsln, epsh, h0, ehp, ehm, epst, ba, ir, h, iback,
    irback, t, ep, em, xw, xa, wg, fa, fb, errt, errh, errd;

  if(a==b){
    *err=0.0;
    return 0.0;
  }

  pi2 = 2.0 * M_PI;
  epsln = 1.0 - log(efs * eps);
  epsh = sqrt(efs * eps);
  h0 = hoff / epsln;
  ehp = exp(h0);
  ehm = 1.0 / ehp;
  epst = exp(-ehm * epsln);
  ba = b - a;
  ir = (*f)((a + b) * 0.5,pa) * (ba * 0.25);
  Ih = ir * (2.0 * pi2);
  *err = fabs(Ih) * epst;
  h = 2 * h0;
  m = 1;
  do {
    iback = Ih;
    irback = ir;
    t = h * 0.5;
    do {
      em = exp(t);
      ep = pi2 * em;
      em = pi2 / em;
      do {
	xw = 1.0 / (1.0 + exp(ep - em));
	xa = ba * xw;
	wg = xa * (1.0 - xw);
	fa = (*f)(a + xa,pa) * wg;
	fb = (*f)(b - xa,pa) * wg;
	ir += fa + fb;
	Ih += (fa + fb) * (ep + em);
	errt = (fabs(fa) + fabs(fb)) * (ep + em);
	if (m == 1) *err += errt * epst;
	ep *= ehp;
	em *= ehm;
      } while (errt > *err || xw > epsh);
      t += h;
    } while (t < h0);
    if (m == 1) {
      errh = (*err / epst) * epsh * h0;
      errd = 1.0 + 2.0 * errh;
    } else {
      errd = h * (fabs(Ih - 2.0 * iback) + 4.0 * fabs(ir - 2.0 * irback));
    }
    h *= 0.5;
    m *= 2;
  } while (errd > errh && m < mmax);
  Ih *= h;
  if (errd > errh) {
    *err = -errd * m;
  } else {
    *err = errh * epsh * m / (2.0 * efs);
  }
  return Ih;
}

double complex deintiz(double complex (*f)(double,void *),double a,void *pa,double eps,double *err)
{
  int m;
  double complex ir,Ih,iback,irback,fp,fm;
  double pi4, epsln, epsh, h0, ehp, ehm, epst, h, 
    t, ep, em, xp, xm, errt, errh, errd;
    
  pi4 = atan(1.0);
  epsln = 1.0 - log(efs * eps);
  epsh = sqrt(efs * eps);
  h0 = hoff / epsln;
  ehp = exp(h0);
  ehm = 1.0 / ehp;
  epst = exp(-ehm * epsln);
  ir = (*f)(a + 1.0,pa);
  Ih = ir * (2.0 * pi4);
  *err = cabs(Ih) * epst;
  h = 2.0 * h0;
  m = 1;
  do {
    iback = Ih;
    irback = ir;
    t = h * 0.5;
    do {
      em = exp(t);
      ep = pi4 * em;
      em = pi4 / em;
      do {
        xp = exp(ep - em);
        xm = 1.0 / xp;
        fp = (*f)(a + xp,pa) * xp;
        fm = (*f)(a + xm,pa) * xm;
        ir += fp + fm;
        Ih += (fp + fm) * (ep + em);
        errt = (cabs(fp) + cabs(fm)) * (ep + em);
        if (m == 1) *err += errt * epst;
        ep *= ehp;
        em *= ehm;
      } while (errt > *err || xm > epsh);
      t += h;
    } while (t < h0);
    if (m == 1) {
      errh = (*err / epst) * epsh * h0;
      errd = 1.0 + 2.0 * errh;
    } else {
      errd = h * (cabs(Ih - 2.0 * iback) + 4.0 * cabs(ir - 2.0 * irback));
    }
    h *= 0.5;
    m *= 2;
  } while (errd > errh && m < mmax);
  Ih *= h;
  if (errd > errh) {
    *err = -errd * (double)m;
  } else {
    *err = errh * epsh * (double)m / (2.0 * efs);
  }
  return Ih;
}

double deintid(double(*f)(double,void *),double a,void*pa,double eps,double *err)
{
  int m;
  double Ih;
  double pi4, epsln, epsh, h0, ehp, ehm, epst, ir, h, iback, irback, 
    t, ep, em, xp, xm, fp, fm, errt, errh, errd;
    
  pi4 = atan(1.0);
  epsln = 1.0 - log(efs * eps);
  epsh = sqrt(efs * eps);
  h0 = hoff / epsln;
  ehp = exp(h0);
  ehm = 1.0 / ehp;
  epst = exp(-ehm * epsln);
  ir = (*f)(a + 1.0,pa);
  Ih = ir * (2.0 * pi4);
  *err = fabs(Ih) * epst;
  h = 2.0 * h0;
  m = 1;
  do {
    iback = Ih;
    irback = ir;
    t = h * 0.5;
    do {
      em = exp(t);
      ep = pi4 * em;
      em = pi4 / em;
      do {
	xp = exp(ep - em);
	xm = 1.0 / xp;
	fp = (*f)(a + xp,pa) * xp;
	fm = (*f)(a + xm,pa) * xm;
	ir += fp + fm;
	Ih += (fp + fm) * (ep + em);
	errt = (fabs(fp) + fabs(fm)) * (ep + em);
	if (m == 1) *err += errt * epst;
	ep *= ehp;
	em *= ehm;
      } while (errt > *err || xm > epsh);
      t += h;
    } while (t < h0);
    if (m == 1) {
      errh = (*err / epst) * epsh * h0;
      errd = 1.0 + 2.0 * errh;
    } else {
      errd = h * (fabs(Ih - 2.0 * iback) + 4.0 * fabs(ir - 2.0 * irback));
    }
    h *= 0.5;
    m *= 2;
  } while (errd > errh && m < mmax);
  Ih *= h;
  if (errd > errh) {
    *err = -errd * (double)m;
  } else {
    *err = errh * epsh * (double)m / (2.0 * efs);
  }
  return Ih;
}
