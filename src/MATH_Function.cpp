
#include "MATH_Function.h"

#include <math.h>

using namespace Rcpp;

MATH_Integration::MATH_Integration(List fns,double reltol_) {
      reltol=reltol_;
      integrate=new Function("integrate");
      fcts=fns;
      integrand=new Function("identity");
}


void MATH_Integration::setFunctionName(std::string name) {
  *integrand = fcts[name];
}


double MATH_Integration::integralFunction(double a, double b) {
  setFunctionName("test");
  List res=(*integrate)(*integrand, a, b,_["rel.tol"] = reltol);
  double integ=as<double>(res["value"]);
  return integ;
}
