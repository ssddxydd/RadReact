/*********************************************************************

  Code for pushing particles using classical radiation reaction 
	
*************************************************************************/

#include <stdio.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"
#include "functions.h"
#include "defs.h"
#include <math.h>

#define RADCOEFF 0.01    /* 2/3 alpha_f B_0/B_crit */
#define CHARGE -1
/* GLOBALS */
double dxsav,*xp,**yp;   /* defining declarations */
int kmax,kount;

double xlog;

int nrhs;   /* counts function evaluations */


            
int main(void)
{


 /* Particle variables - change to arrays later
     u,v,w are spatial components of 4-momentum u=p/mc
   */
  
  int nbad,nok;

  double
    eps=1.0e-6,h1=0.1,hmin=0.0,
    x1,x2,*yvec;

  double 
    x_0, y_0, z_0,
    gamma, gamma0,
    invgam,
    betax, betay, betaz, beta,
    utheta, uphi;

  double 
    *y, *dydx,
    *yscal, 
    *B, *E, 
    e_x, e_y, e_z, b_x, b_y, b_z;

  double 
    u0,v0,w0, 
    u1,v1,w1,
    up,vp,wp,
    x_p,y_p,z_p;

  /* useful parameters */
  double
    omg, omgstar,  //nonrel + rel gyrofrequencies
    rg, rgstar,    // and gyro-radii
    freq,          // frequency of wave in units of gyrofrequency
    omega, om_crit,
    alpha, dalpha,
    Bmag, Emag, umag,
    a,b,c,mod, F, rad_force;

  /* Particle parameters and physical constants */
  double
    q0=CHARGE, m0, clight;

  /* index integers */

  int
    i, j, k;

  /* update parameters */

  double
    tstep, htry, tnew,
    time, x, dummy,
    hdid, hnext;

  FILE
    *trace, *trace2,
    *spectrum, *spectot,
    *radangle;


  /****************************************************************/
  yvec=dvector(1,N);

  B=dvector(1,3);
  E=dvector(1,3);


  /* particle begin and end times (in proper time this may be very long) */
  x1=0.0;
  x2=2.*PI*1.e3;

  /* maximum number of output points from odeint */
  kmax=1.e4;
  dxsav=(x2-x1)/(kmax+1);

  xp=dvector(1,kmax);
  yp=dmatrix(1,N,1,kmax);


  /**************************************************************/

  /* initialise particle position and velocity */
  x_0=y_0=z_0=0.;

  gamma0=100.; //1.973200071;
  betax=0.0;
  betaz=0.0;
  betay=sqrt(1.-1./(gamma0*gamma0));

  beta=sqrt(betax*betax + betay*betay + betaz*betaz);
 
  
  printf("initial Lorentz factor of particle: %f\n", gamma0);
  printf("initial velocity of particle: %16f  %16f  %16f\n", betax, betay, betaz);



#if ODE_INTEGRATOR

  gamma=gamma0;
#if CARTESIAN
  yvec[1]=x_0;
  yvec[2]=y_0;
  yvec[3]=z_0;
  yvec[4]=betax*gamma;
  yvec[5]=betay*gamma;
  yvec[6]=betaz*gamma;
#endif
#if SPHERICAL
  yvec[1]=x_0;
  yvec[2]=y_0;
  yvec[3]=z_0;
  yvec[4]=beta*gamma;
  yvec[5]=betax/beta;
  yvec[6]=atan2(betay,betaz);
#endif


 

  nrhs=0;

  printf("starting integrator\n");
  
#if BS5
  odeint(yvec,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);
#endif
 
#if RK5
  odeint(yvec,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs);
#endif
 
  printf("\n%s %13s %3d\n","successful steps:"," ",nok);
  printf("%s %20s %3d\n","bad steps:"," ",nbad);
  printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
  printf("\n%s %3d\n","stored intermediate values:    ",kount);
  


  /**************************************************************/


  
  trace=fopen("trace.dat","w");
  

  for (i=1;i<=kount;i++){
    
#if CARTESIAN
    umag=sqrt( yp[4][i]* yp[4][i]+ yp[5][i]* yp[5][i]+ yp[6][i]* yp[6][i]);
    gamma=sqrt(1.+umag*umag);

    fprintf(trace,"%f %.32f %.32f %.32f %.32f %.32f %.32f %.32f\n", 
	    xp[i], yp[1][i], yp[2][i], yp[3][i], yp[4][i], yp[5][i], yp[6][i], gamma);

#endif
#if SPHERICAL

    gamma=sqrt(1.+yp[4][i]* yp[4][i]);

    fprintf(trace,"%f %.32f %.32f %.32f %.32f %.32f %.32f %.32f\n",	\
	    xp[i], yp[1][i], yp[2][i], yp[3][i], yp[4][i], yp[5][i],\
	    yp[6][i], gamma);
#endif


  }

  fclose(trace);


#endif  /* ODE_INTEGRATOR */


#if BORIS
 
  /*                     Boris pusher                      */

  /* Particle velocity is defined at t=n-1/2               */
  /* Particle position is defined at t=n                   */
  /* Fields must be determined at t=n                      */
  /* Update gives velocity at n+1/2 and position at n+1    */
  

  yvec[1]=x_0;
  yvec[2]=y_0;
  yvec[3]=z_0;
  yvec[4]=betax*gamma0;
  yvec[5]=betay*gamma0;
  yvec[6]=betaz*gamma0;
  /* note - these quantities are staggered */
  up=yvec[4];vp=yvec[5];wp=yvec[6]; 


  time=x1;
  htry=.1;

  trace=fopen("trace.dat","w");
  i=0;
  while(time<x2){

    field_eval(time, yvec, B, E);

    e_x=E[1]; e_y=E[2]; e_z=E[3];
    b_x=B[1]; b_y=B[2]; b_z=B[3];
    // up=yvec[4];vp=yvec[5];wp=yvec[6]; 

    tstep=0.5*q0*htry;

    /* First half of acceleration by electric impulses */
    u0=up+e_x*tstep;
    v0=vp+e_y*tstep;
    w0=wp+e_z*tstep;
    
    /* First half of rotation by B */
    gamma=sqrt(1.0+u0*u0+v0*v0+w0*w0);
    invgam=1.0/gamma;
    b_z*=invgam*tstep;
    b_y*=invgam*tstep;
    b_x*=invgam*tstep;
    
    mod=2./(1.+b_x*b_x+b_y*b_y+b_z*b_z);
    
    u1=(u0+v0*b_z-w0*b_y)*mod;
    v1=(v0+w0*b_x-u0*b_z)*mod;
    w1=(w0+u0*b_y-v0*b_x)*mod;
    
    /* Second half of rotation and acceleration */
    
    up=u0+v1*b_z-w1*b_y+e_x*tstep; 
    vp=v0+w1*b_x-u1*b_z+e_y*tstep;
    wp=w0+u1*b_y-v1*b_x+e_z*tstep;
    
    
   




#if RADREACT
    u1=0.5*(yvec[4]+up);
    v1=0.5*(yvec[5]+vp);
    w1=0.5*(yvec[6]+wp);
    gamma=sqrt(1.0+u1*u1+v1*v1+w1*w1);
    invgam=1.0/gamma;

    a=E[1]+(v1*B[3]-w1*B[2])*invgam;
    b=E[2]+(w1*B[1]-u1*B[3])*invgam;
    c=E[3]+(u1*B[2]-v1*B[1])*invgam;


    mod=a*a+b*b+c*c;
    dummy=(E[1]*u1+E[2]*v1+E[3]*w1)*invgam;
    mod-=dummy*dummy;

    rad_force=gamma*mod*RADCOEFF*htry;
  
  // dummy=cos(x-y[3])*cos(x-y[3])*(gam-y[6])*(gam-y[6])/(gam*gam);
  // printf("Rad force: mod= %e f= %e\n", mod, rad_force );
  
     up -= rad_force*u1;
     vp -= rad_force*v1;
     wp -= rad_force*w1;

#endif

 /* update position using new velocities*/
    gamma=sqrt(1.0+up*up+vp*vp+wp*wp);
    invgam=1.0/gamma;

    x_p=x_0+htry*up*invgam;
    y_p=y_0+htry*vp*invgam;
    z_p=z_0+htry*wp*invgam;
    
    if(i%20==0){
      fprintf(trace,"%f %.32f %.32f %.32f %.32f %.32f %.32f %.32f\n",	\
	      time, x_p, y_p, z_p, gamma, up, vp, wp);
    }
    time+=htry;
    i++;

    x_0=x_p;
    y_0=y_p;
    z_0=z_p;


    yvec[1]=x_p;yvec[2]=y_p;yvec[3]=z_p; 
    yvec[4]=up;yvec[5]=vp;yvec[6]=wp; 
 
    
  }

#endif /* end BORIS */
 




  /****************** Free vectors ******************************/

  free_dvector(B,1,3);
  free_dvector(E,1,3);
  free_dvector(yvec,1,N);
  free_dvector(xp,1,kmax);
  free_dmatrix(yp,1,N,1,kmax);
  return 0;
}



/***************************************************************/





/*  ---------------        Subroutines     ------------------  */




/***************************************************************/


void field_eval(double x, double y[], double B[], double E[]){

  int i;

  E[1]=0.;
  E[2]=0.;
  E[3]=0.;
  B[1]=1.;
  B[2]=0.;
  B[3]=0.;


  
}



void derivs(double x,double y[],double dydx[])   {

/* y[1-3] are positions, y[4-6] are momenta */
  double umag2, gam, invgam, invu,
    gam2, bet, bet2 , invgam2,
    dummy, sin2;

  double *B, *E;

  double q=CHARGE;

  double mod, rad_force;

  double ex,ey,ez;
  double a, b, c;

  extern int nrhs;


  B=dvector(1,3);
  E=dvector(1,3);


  nrhs++;

  field_eval(x, y, B, E);

#if CARTESIAN
  
  /*    1 2 3 4   5     6       */
  /* y=[x,y,z,ux,uy,uz] */
  umag2=y[4]*y[4]+y[5]*y[5]+y[6]*y[6];
  gam2=1.+umag2;
  gam=sqrt(gam2);
  invgam = 1./gam;
 
  


#endif

#if SPHERICAL

  /*    1 2 3 4   5     6        */
  /* y=[x,y,z,u,mu,phi] */

  /*
  if(y[6]==0.)
    nrerror("alligned with z-axis");
  */

  gam2=1.+y[4]*y[4];
  gam=sqrt(gam2);

  invgam = 1./gam;
  invu = 1./y[4];
  bet=y[4]*invgam;

  ex = y[5];
  sin2=1.-ex*ex;
  dummy=sqrt(sin2);
  ey = dummy*cos(y[6]);
  ez = dummy*sin(y[6]);

  if(sin2==0)
    nrerror("particle alligned with x-axis");

#endif

  /*****************************************************************/


#if CARTESIAN


  /* (F_ij u_j)^2 */
  a=E[1]+(y[5]*B[3]-y[6]*B[2])*invgam;
  b=E[2]+(y[6]*B[1]-y[4]*B[3])*invgam;
  c=E[3]+(y[4]*B[2]-y[5]*B[1])*invgam;
  
  dydx[1] = y[4]*invgam;
  dydx[2] = y[5]*invgam;
  dydx[3] = y[6]*invgam;
  dydx[4] = q*a;
  dydx[5] = q*b;
  dydx[6] = q*c;
  

#if RADREACT
  mod=a*a+b*b+c*c;
  dummy=(E[1]*y[4]+E[2]*y[5]+E[3]*y[6])*invgam;
  mod-=dummy*dummy;

  rad_force=gam*mod*RADCOEFF;
  
  // dummy=cos(x-y[3])*cos(x-y[3])*(gam-y[6])*(gam-y[6])/(gam*gam);
  // printf("Rad force: mod= %e f= %e\n", mod, rad_force );
  
  dydx[4] -= rad_force*y[4];
  dydx[5] -= rad_force*y[5];
  dydx[6] -= rad_force*y[6];

#endif
 
#endif

#if SPHERICAL


  dydx[1] = bet*ex;
  dydx[2] = bet*ey;
  dydx[3] = bet*ez;
  dydx[4] = q*(E[1]*ex+E[2]*ey+E[3]*ez);
  dydx[5] = (E[1]*sin2-ex*(E[2]*ey+E[3]*ez))*invu;
  dydx[5] += (B[3]*ey-B[2]*ez)*invgam;
  dydx[5] *= q;
  
  dydx[6] = ey*(E[3]*invu+invgam*(ex*B[2]-ey*B[1]))-ez*(E[2]*invu+invgam*(ez*B[1]-ex*B[3]));
  dydx[6] *= q/sin2;
  
  /*
    dydx[5] =  invgam*((cos(y[6])*B[1]+sin(y[6])*B[2])/tan(y[7])-B[3]);
    dydx[6] =  invgam*(sin(y[6])*B[1]-cos(y[6])*B[2]);
  */
  

#if RADREACT
  ex*=y[4]; ey*=y[4]; ez*=y[4];
  a=E[1]+(ey*B[3]-ez*B[2])*invgam;
  b=E[2]+(ez*B[1]-ex*B[3])*invgam;
  c=E[3]+(ex*B[2]-ey*B[1])*invgam;
  mod=a*a+b*b+c*c;
  dummy=(E[1]*ex+E[2]*ey+E[3]*ez)*invgam;
  mod-=dummy*dummy;
  
  rad_force=gam*mod*RADCOEFF;
  
  // dummy=cos(x-y[3])*cos(x-y[3])*(gam-y[6])*(gam-y[6])/(gam*gam);
  // printf("Rad force: mod= %e f= %e\n", mod, rad_force );
  
  dydx[4] -= rad_force*y[4];


#endif

#endif

   
 




  free_dvector(B,1,3);
  free_dvector(E,1,3);

}








#undef NRANSI
