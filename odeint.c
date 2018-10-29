#include <math.h>
#include "defs.h"
#define NRANSI
#include "nrutil.h"
#define MAXSTP 1000000000
#define TINY 1.0e-30

extern int kmax,kount;
extern double *xp,**yp,dxsav;

void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])))
{
	int nstp,i;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;

	yscal=dvector(1,nvar);
	y=dvector(1,nvar);
	dydx=dvector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
	        (*derivs)(x,y,dydx);
#if PROPERTIME	
        	for (i=1;i<=4;i++){
		  //yscal[i]=10.;
		  //yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		  // yscal[i]=fabs(y[i])+TINY;
		  yscal[i]=FMAX(0.01*fabs(y[i]), 0.01);
		}
		for (i=5;i<=N;i++){
		  //yscal[i]=10.;
		  //yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		  // yscal[i]=fabs(y[i])+TINY;
		  yscal[i]=FMAX(0.1*fabs(y[i]), 0.1);
		}
	
#else
        	for (i=1;i<=3;i++){
		  //yscal[i]=10.;
		  yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		  // yscal[i]=fabs(y[i])+TINY;
		}
		yscal[1]=0.1;
		yscal[4]=0.01;
		yscal[5]=y[5];
		yscal[6]=1.;
		/*
		for (i=4;i<=nvar;i++){
		  //yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		  yscal[i]=FMAX(fabs(0.001*h*dydx[i]),0.00001);
	
		  // yscal[i]=FMAX(FMIN(0.1*fabs(y[i]),0.01),1.0);
		}
		*/
#endif
	
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
		
		        xp[++kount]=x;
			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax) {
				xp[++kount]=x;
				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			}
			free_dvector(dydx,1,nvar);
			free_dvector(y,1,nvar);
			free_dvector(yscal,1,nvar);
			return;
		}
		if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
		h=hnext;
	
		/*	if(hnext>1.0)
			h=1.0;
			else 
			h=hnext;*/
	}
	nrerror("Too many steps in routine odeint");
}
#undef MAXSTP
#undef TINY
#undef NRANSI
