/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nr.h.  Do not confuse this file with the same-named
   file nr.h that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#ifndef _NR_H_
#define _NR_H_

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {float r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

#include <stdio.h>


void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi);
void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp);

void bsstep(double y[], double dydx[], int nv, double *xx, double htry,
	double eps, double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));

float chebev(float a, float b, float c[], int m, float x);

void derivs(double x, double y[], double dydx[]);
void jacobn(double x, double y[], double dfdx[], double **dfdy, int n);

void lubksb(double **a, int n, int *indx, double b[]);
void ludcmp(double **a, int n, int *indx, double *d);

double midpnt(double (*func)(double), double a, double b, int n);
void mmid(double y[], double dydx[], int nvar, double xs, double htot,
	int nstep, double yout[], void (*derivs)(double, double[], double[]));

void odeint(double ystart[], int nvar, double x1, double x2,
	double eps, double h1, double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double,
	double [], double *, double *, void (*)(double, double [], double [])));
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
	int nv);
double qromo(double (*func)(double), double a, double b,
	     double (*choose)(double(*)(double), double, double, int));
void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
	void (*derivs)(double, double [], double []));
void rkck(double y[], double dydx[], int n, double x, double h,
	double yout[], double yerr[], void (*derivs)(double, double [], double []));
void rkdumb(double vstart[], int nvar, double x1, double x2, int nstep,
	    void (*derivs)(double, double [], double []));
void rkqs(double y[], double dydx[], int n, double *x,
	  double htry, double eps, double yscal[], double *hdid, double *hnext,
	  void (*derivs)(double, double [], double []));

void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[],
	   void (*derivs)(double, double [], double []));

void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	    double yscal[], double *hdid, double *hnext,
	    void (*derivs)(double, double [], double []));

#endif /* _NR_H_ */
