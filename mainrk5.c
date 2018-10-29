/************************************************************
 *   Calculate synchrotron/jitter spectrum from ensemble
 *   of particles. 
 *
 *
 *   Time/length units are in terms of gyrofrequency/gyroradius 
 *   set to nonrelativistic value omg=eB/mc, rg=c/omg
 *
 ************************************************************/


#include "defs.h"
#include "functions.h"

int main(int argc, char **argv){


  /* Particle variables - change to arrays later
     u,v,w are spatial components of 4-momentum u=p/mc
   */
  
  double 
    x0, y0, z0,
    gamma, gamma0,
    betax, betay, betaz;

  double 
    *y, *dydx,
    *yscal;

 
  /* useful parameters */
  double
    omg, omgstar,  //nonrel + rel gyrofrequencies
    rg, rgstar,    // and gyro-radii
    freq;          // frequency of wave in units of gyrofrequency

  /* Particle parameters and physical constants */
  double
    q0, m0, c;

  /* index integers */

  int
    i, j, k;

  /* update parameters */

  double
    tstep, htry, tnew,
    time, x,
    eps, hdid, hnext;


  FILE
    *trace;

 

  /*********************************************************************/
  

  y=(double *) calloc((unsigned) (N),sizeof(double));
  dydx=(double *) calloc((unsigned) (N),sizeof(double));
  yscal=(double *) calloc((unsigned) (N),sizeof(double));


  c=m0=q0=1.;


  /* initialise particle position and velocity */
  x0=y0=z0=0.;

  gamma0=20.;
  betax=0.0;
  betay=0.0;
  betaz=sqrt(1.-1./(gamma0*gamma0));
  printf("initial Lorentz factor of particle: %f\n", gamma0);
  printf("initial velocity of particle: %16f\n", betaz);
 

  gamma=gamma0;

  y[0]=x0;
  y[1]=y0;
  y[2]=z0;
  y[3]=betax*gamma;
  y[4]=betay*gamma;
  y[5]=betaz*gamma;



  freq=20.;


  trace=fopen("trace.dat","w");

#if RK5

  eps=1.0e-6;

  for (i=0;i<N;i++) yscal[i]=1.0;
  htry=0.9;
  
  
  i=0;
  x=0.0;
  while(x<200){
    bsstep(y,dydx,N,&x,htry,eps,yscal,&hdid,&hnext,derivs);
    if(i%50==0){
      gamma=sqrt(1.0+y[3]*y[3]+y[4]*y[4]+y[5]*y[5]);
      fprintf(trace,"%f %.32f %.32f %.32f %.32f %.32f %.32f %.32f\n",	\
	      x, y[0], y[1], y[2], gamma, y[3], y[4], y[5]);
    }
    htry=hnext;
    i++;
  }
 
  printf("finished RK routine\n");
 
#endif /* end RK5 */
 
  
  fclose(trace);

  free(y);
  free(dydx);
  free(yscal);


  return(0);
 


} /* code exits */


