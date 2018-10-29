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
    u0, v0, w0,
    u1, v1, w1,
    xp, yp, zp, 
    up, vp, wp,
    gamma, invgamma,
    gamma0,
    betax, betay, betaz;

  double 
    *y, *dydx,
    *yscal;

  /* Field parameters - will call from function later */

  double
    ex0, ey0, ez0,
    ex, ey, ez,
    bx, by, bz;

  /* useful parameter */
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
    time, x, f,
    eps, hdid, hnext;


  FILE
    *trace;

  c=m0=q0=1.;

  /*********************************************************************/
  
  trace=fopen("trace.dat","w");

  /* initialise particle position and velocity */
  x0=y0=z0=0.;
  xp=x0;
  yp=y0;
  zp=z0;

  y=(double *) calloc((unsigned) (N),sizeof(double));
  //  if (!y) nrerror("allocation error 1 in dmatrix()");
  dydx=(double *) calloc((unsigned) (N),sizeof(double));
  yscal=(double *) calloc((unsigned) (N),sizeof(double));

  betax=0.0;
  betay=0.0;
  betaz=0.999;
  gamma0=1.0/sqrt(1.-(betax*betax+betay*betay+betaz*betaz));
  printf("initial Lorentz factor of particle: %f\n", gamma0);

  gamma=gamma0;
 

  up=betax*gamma;
  vp=betay*gamma;
  wp=betaz*gamma;

  y[0]=x0;
  y[1]=y0;
  y[2]=z0;
  y[3]=up;
  y[4]=vp;
  y[5]=wp;

  /* calculate field quantities */

  ex0=100.0;
  ey0=0.0;
  ez0=0.0;

  freq=20.;

  time=0.0;
  htry=0.001;


#if BORIS
 
  /*                     Boris pusher                      */

  /* Particle velocity is defined at t=n-1/2               */
  /* Particle position is defined at t=n                   */
  /* Fields must be determined at t=n                      */
  /* Update gives velocity at n+1/2 and position at n+1    */
  
  for(i=0;i<5000000;i++){
   
    ex=ex0*cos(freq*(time-zp));
    ey=ey0*sin(freq*(time-zp));
    bx=-ey0*sin(freq*(time-zp));
    by=ex0*cos(freq*(time-zp));;
    ez=0.0;
    bz=0.0;
   
    tstep=0.5*htry;

    /* First half of acceleration by electric impulses */
    u0=up+ex*tstep;
    v0=vp+ey*tstep;
    w0=wp+ez*tstep;
    
    /* First half of rotation by B */
    gamma=sqrt(1.0+u0*u0+v0*v0+w0*w0);
    invgamma=1.0/gamma;
    bz=invgamma*tstep*bz;
    by=invgamma*tstep*by;
    bx=invgamma*tstep*bx;
    
    f=2./(1.+bx*bx+by*by+bz*bz);
    
    u1=(u0+v0*bz-w0*by)*f;
    v1=(v0+w0*bx-u0*bz)*f;
    w1=(w0+u0*by-v0*bx)*f;
    
    /* Second half of rotation and acceleration */
    
    u0=u0+v1*bz-w1*by+ex*tstep; 
    v0=v0+w1*bx-u1*bz+ey*tstep;
    w0=w0+u1*by-v1*bx+ez*tstep;
    
    /* use updated particle velocities */
    up=u0;
    vp=v0;
    wp=w0;
    
    /* update position using new velocities*/
    gamma=sqrt(1.0+u0*u0+v0*v0+w0*w0);
    invgamma=1.0/gamma;

    xp=x0+htry*up*invgamma;
    yp=y0+htry*vp*invgamma;
    zp=z0+htry*wp*invgamma;
    
    if(i%50==0){
      fprintf(trace,"%f %.32f %.32f %.32f %.32f %.32f %.32f %.32f\n",	\
	      time, xp, yp, zp, gamma, up, vp, wp);
    }
    time+=htry;

    x0=xp;
    y0=yp;
    z0=zp;
    
  }

#endif /* end BORIS */
 

#if RK5

  eps=1.0e-6;

  for (i=0;i<N;i++) yscal[i]=1.0;
  htry=0.6;
  
  
  
  x=0.0;
  for(i=0;i<1000000;i++){
    rkqs(y,dydx,N,&x,htry,eps,yscal,&hdid,&hnext,derivs);
    if(i%50==0){
      gamma=sqrt(1.0+y[3]*y[3]+y[4]*y[4]+y[5]*y[5]);
      fprintf(trace,"%f %.32f %.32f %.32f %.32f %.32f %.32f %.32f\n",	\
	      x, y[0], y[1], y[2], gamma, y[3], y[4], y[5]);
    }
  }
 
  printf("finished RK routine\n");
 
#endif /* end RK5 */
 
  
  fclose(trace);

  free(y);
  free(dydx);
  free(yscal);


  return(0);
 


} /* code exits */


