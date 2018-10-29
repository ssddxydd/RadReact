/* Driver for routine bsstep */

#include <stdio.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"
#include "functions.h"
#include "defs.h"
#include <math.h>


/* GLOBALS */
double dxsav,*xp,**yp;   /* defining declarations */
int kmax,kount;

double xlog;

int nrhs;   /* counts function evaluations */

double Bunits;   /* magnetic field strength in units of m^2c^4/e^3 */
                       /* Which corresponds to ~ 1.6e17 G  */
            
int main(void)
{


 /* Particle variables - change to arrays later
     u,v,w are spatial components of 4-momentum u=p/mc
   */
  
  int nbad,nok;

  double
    eps=1.0e-6,h1=0.1,hmin=0.0,
    x1,x2,*ystart;

  double 
    x0, y0, z0,
    gamma, gamma0,
    invgam,
    betax, betay, betaz;

  double 
    *y, *dydx,
    *yscal, 
    *westfold,
    *B, *E;

  double
    logxmin, logxmax, dlogx;

  /* useful parameters */
  double
    omg, omgstar,  //nonrel + rel gyrofrequencies
    rg, rgstar,    // and gyro-radii
    freq,          // frequency of wave in units of gyrofrequency
    omega, om_crit,
    alpha, dalpha,
    Bmag, Emag, umag,
    a,b,c,mod;

  /* Particle parameters and physical constants */
  double
    q0, m0, clight;

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
    *spectrum, *spectot;


  /****************************************************************/
  ystart=dvector(1,N);

  B=dvector(1,3);
  E=dvector(1,3);
  westfold=dvector(1,NSPEC);

   
  Bunits=6.e-7; /* Roughly 1G field */

  for(j=1;j<=NSPEC;j++){
    westfold[j]=0.;
  }

  /* for radiation spectrum */
  logxmin=-2.;
  logxmax=4.;
  dummy=(double)NSPEC - 1.;
  dlogx=(logxmax-logxmin)/dummy;

  /* particle begin and end times */
  x1=0.0;
  x2=1000.;

  /* maximum number of output points from odeint */
  kmax=100;
  dxsav=(x2-x1)/(kmax+1);

  xp=dvector(1,kmax);
  yp=dmatrix(1,N,1,kmax);


  /**************************************************************/

  /* initialise particle position and velocity */
  x0=y0=z0=0.;

  gamma0=10.; //1.973200071;
  betax=0.0;
  betay=0.0;
  betaz=sqrt(1.-1./(gamma0*gamma0));
 

  printf("initial Lorentz factor of particle: %f\n", gamma0);
  printf("initial velocity of particle: %16f\n", betaz);
 
  gamma=gamma0;

  ystart[1]=x0;
  ystart[2]=y0;
  ystart[3]=z0;
  ystart[4]=betax*gamma;
  ystart[5]=betay*gamma;
  ystart[6]=betaz*gamma;



 

  nrhs=0;

  printf("starting integrator\n");
  
#if BS5
  odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep);
#endif
 
#if RK4
  odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs);
#endif
 
  printf("\n%s %13s %3d\n","successful steps:"," ",nok);
  printf("%s %20s %3d\n","bad steps:"," ",nbad);
  printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
  printf("\n%s %3d\n","stored intermediate values:    ",kount);
  


  /**************************************************************/


  
  trace=fopen("trace.dat","w");
  spectrum=fopen("spectrum.dat","w");
  
  for (i=1;i<=kount;i++){
    

    /* calculate critical frequency */
    umag=sqrt( yp[4][i]* yp[4][i]+ yp[5][i]* yp[5][i]+ yp[6][i]* yp[6][i]);
    gamma=sqrt(1.+umag*umag);
    invgam=1./gamma;

    ystart[3]=yp[3][i];
    field_eval(xp[i],ystart,B,E);
    
    Bmag=sqrt(B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);

    /* find E + vxB */
    a=E[1]+(yp[5][i]*B[3]-yp[6][i]*B[2])*invgam;
    b=E[2]+(yp[6][i]*B[1]-yp[4][i]*B[3])*invgam;
    c=E[3]+(yp[4][i]*B[2]-yp[5][i]*B[1])*invgam;

    mod=a*a+b*b+c*c;
    dummy=(E[1]*yp[4][i]+E[2]*yp[5][i]+E[3]*yp[6][i])*invgam;
    mod-=dummy*dummy;
    mod=sqrt(mod);

    om_crit=1.5*gamma*gamma*mod;

    fprintf(trace,"%f %.32f %.32f %.32f %.32f %.32f %.32f %.32f %.32f\n",	\
	    xp[i], yp[1][i], yp[2][i], yp[3][i], yp[4][i], yp[5][i],\
	    yp[6][i], gamma, om_crit);
    omega=logxmin;

#if IOSPECTRUM
  
    for(j=1;j<=NSPEC;j++){
        
      x=omega-log10(om_crit);
      dummy = fairy(pow(10.,x));
      dummy*=mod;
 
      westfold[j]+=dummy;

      fprintf(spectrum,"%f %f %f %f \n",xp[i], pow(10.,omega), om_crit, dummy);

      omega+=dlogx;
    }
    fprintf(spectrum,"\n");
#endif
  }

  fclose(trace);
  fclose(spectrum);

  /*output summed spectrum */
#if IOSPECTRUM  
  spectot=fopen("tot_spectrum.dat","w");
  omega=logxmin;

  for(j=1;j<=NSPEC;j++){
    fprintf(spectot, "%.32f %.32f\n",  pow(10.,omega), westfold[j]);
    omega+=dlogx;
  }
  
  fclose(spectot);
#endif



  free_dvector(B,1,3);
  free_dvector(E,1,3);
  free_dvector(ystart,1,N);
  free_dvector(xp,1,kmax);
  free_dmatrix(yp,1,N,1,kmax);
  free_dvector(westfold,1,NSPEC);
  return 0;
}



/***************************************************************/

void derivs(double x,double y[],double dydx[])
{

  /* y[1-3] are positions, y[4-6] are velocities */
  double gam, invgam, betay, dummy;
  double *B, *E;
  double a,b,c,mod;
  double rad_force;

  B=dvector(1,3);
  E=dvector(1,3);

  nrhs++;
 
  
  field_eval(x, y, B, E);
  gam=sqrt(1.+y[4]*y[4]+y[5]*y[5]+y[6]*y[6]);
  invgam=1./gam;

#if RADREACT

  /* (F_ij u_j)^2 */
  a=E[1]+(y[5]*B[3]-y[6]*B[2])*invgam;
  b=E[2]+(y[6]*B[1]-y[4]*B[3])*invgam;
  c=E[3]+(y[4]*B[2]-y[5]*B[1])*invgam;
  
  mod=a*a+b*b+c*c;
  dummy=(E[1]*y[4]+E[2]*y[5]+E[3]*y[6])*invgam;
  mod-=dummy*dummy;
  
  rad_force=gam*mod*Bunits;

  //  printf("Rad force %f\n", rad_force);

  dydx[1] = y[4]*invgam;
  dydx[2] = y[5]*invgam;
  dydx[3] = y[6]*invgam;
  dydx[4] = E[1] + invgam*(B[3]*y[5]-B[2]*y[6])-rad_force*y[4];
  dydx[5] = E[2] + invgam*(B[1]*y[6]-B[3]*y[4])-rad_force*y[5];
  dydx[6] = E[3] + invgam*(B[2]*y[4]-B[1]*y[5])-rad_force*y[6];

#endif

#if LORENTZ

  dydx[1] = y[4]*invgam;
  dydx[2] = y[5]*invgam;
  dydx[3] = y[6]*invgam;
  dydx[4] = E[1] + invgam*(B[3]*y[5]-B[2]*y[6]);
  dydx[5] = E[2] + invgam*(B[1]*y[6]-B[3]*y[4]);
  dydx[6] = E[3] + invgam*(B[2]*y[4]-B[1]*y[5]);

#endif
 
  free_dvector(B,1,3);
  free_dvector(E,1,3);
}

void field_eval(double x, double y[], double B[], double E[]){

  E[1]=100.*cos(20.*(x-y[3]));
  E[2]=0.;
  E[3]=0.;
  B[1]=0.;
  B[2]=E[1];
  B[3]=0.;
  
 

}


#undef NRANSI
