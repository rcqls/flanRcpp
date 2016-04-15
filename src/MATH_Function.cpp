
#include "MATH_Function.h"

#include <math.h>

using namespace Rcpp;

MATH_Integration::MATH_Integration(List fns,double reltol_) {
      reltol=reltol_;
      integrate=new Function("integrate");
      fcts=fns;
      integrand=new Function("identity");  // allouer
}


void MATH_Integration::setFunctionName(std::string name) {
  *integrand = fcts[name];
}


double MATH_Integration::integralFunction(double a, double b,double rho,double delta) {
  setFunctionName("CF_GY_WD");
  List res=(*integrate)(*integrand, a, b, _["rel.tol"] = reltol,_["rho"] =rho,_["delta"] =delta);
  double integ=as<double>(res["value"]);
  return integ;
}
