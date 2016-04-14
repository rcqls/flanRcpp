
#include "MATH_Function.h"

#include <math.h>
#include "integration.h"
#include "stdafx.h"

MATH_Integration::MATH_Integration(){
      mFunctionCalledNumber=0;
      mIntervalsNumber=0;
    }


void MATH_Integration::integralFunction(double x, double xminusa, double bminusx, double &y, void *ptr) {
      // this callback calculates f(x)=exp(x)
      //y=(dynamic_cast<const MATH_Integration*>(ptr))->computeFunction(x);
      y=static_cast<FLAN_Function*>(ptr)->computeFunction(x);
  }

bool MATH_Integration::integrate(const double& a,const double& b,double& v) {

      bool succeeds=true;
      
      alglib::autogkstate s;
      
      alglib::autogkreport rep;
      double x;
      // Integration of a smooth function F(x) on a finite interval [a,b].
      //Fast-convergent algorithm based on a Gauss-Kronrod formula is used. Result
      //is calculated with accuracy close to the machine precision.
      //Algorithm works well only with smooth integrands.  It  may  be  used  with
      //continuous non-smooth integrands, but with  less  performance.
      //It should never be used with integrands which have integrable singularities
      //at lower or upper limits - algorithm may crash. Use AutoGKSingular in such cases.
      alglib::autogksmooth((double) a, (double) b, s);
      alglib::autogkintegrate(s,integralFunction,getFunction());
      alglib::autogkresults(s, x, rep);
      // analyse the rep structure
      // terminationtype=-5:non-convergence of Gauss-Kronrod nodes calculation subroutine.
      // terminationtype=-1:incorrect parameters specified
      // terminationtype=1: ok
      mFunctionCalledNumber=rep.nfev;
      mIntervalsNumber=rep.nintervals;
      succeeds=succeeds && (rep.terminationtype==1);
      
      v=x;
      return succeeds;
}
