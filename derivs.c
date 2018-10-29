#include <math.h>
#include "functions.h"
#include "defs.h"

//void derivs(double x, double y[], double dydx[], double B[], double E[])
void derivs(double x, double y[], double dydx[])
{
  /* y[0-2] are positions, y[3-5] are velocities */
  double invgam;
  double B[3], E[3];

  invgam=1./sqrt(1+y[3]*y[3]+y[4]*y[4]+y[5]*y[5]);
  
  E[0]=100.*cos(20.*(x-y[2]));
  E[1]=0.;
  E[2]=0.;
  B[0]=0.;
  B[1]=E[0];
  B[2]=0.;


  dydx[0] = y[3]*invgam;
  dydx[1] = y[4]*invgam;
  dydx[2] = y[5]*invgam;
  dydx[3] = E[0] + invgam*(B[2]*y[4]-B[1]*y[5]);
  dydx[4] = E[1] + invgam*(B[0]*y[5]-B[2]*y[3]);
  dydx[5] = E[2] + invgam*(B[1]*y[3]-B[0]*y[4]);
}
