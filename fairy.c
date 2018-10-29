#include "nr.h"
#include "functions.h"
#include <math.h>



double fairy(double x){

  extern double xlog;

  double 
    fairy, 
    zlow, zhigh, 
    answer;

  xlog=log10(x);
  if(x<=1.){
    zlow=0.;
    zhigh=1.;
    answer=qromo(s1func,zlow,zhigh,midpnt);
    fairy=3./2.*x*x*answer;
  }
  else{
    zlow=0.;
    zhigh=1.;
    answer=qromo(s2func,zlow,zhigh,midpnt);
    fairy=x*answer;
  }

  return(fairy);

}

/*****************************************************************/



/*****************************************************************/

double s1func(double z){
  
  extern double xlog;

  double
    mys1func,
    x,y,xarg,xnu,
    rk,ri,rip,rkp;

  y=pow(z,(-3./2.));
  x=pow(10., xlog);
  xarg=x*y;
  xnu=5./3.;
  


  if(xarg < 0.001)
    rk=0.902745/2.*pow((xarg/2.),(-5./3.));
  else if(xarg<=10.)
    bessik(xarg,xnu,&ri,&rk,&rip,&rkp);
  else
    rk=sqrt(2.*atan(1.)/xarg)*exp(-xarg);
  
  mys1func=pow(z,(-5./2.))*rk;
    
  return(mys1func);
}

double s2func(double z){
  
  extern double xlog;
 
  double
    mys2func,x,
    y,xarg,xnu,rk,ri,rip,
    rkp;
   

  y=-log(z);
  x=pow(10., xlog);
  xarg=x+y;


  xnu=5./3.;

  if(xarg < 0.001)
    rk=0.902745/2.*pow((xarg/2.),(-5./3.));
  else if(xarg <= 10.)
    bessik(xarg,xnu,&ri,&rk,&rip,&rkp);
  else
    rk=sqrt(2.*atan(1.)/xarg)*exp(-xarg);
     
  mys2func=rk/z;
  return(mys2func);
}


