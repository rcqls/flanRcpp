
#include "MATH_Function.h"

#include <math.h>

using namespace Rcpp;

MATH_Integration::MATH_Integration(double reltol_) {
      mReltol=reltol_;
      mIntegrate=new Function("integrate");
      List fns=Environment::global_env().get(".integrands");
      fcts=List::create(
        _["CLONE_P0_WD"]=fns["CLONE_P0_WD"],
        _["CLONE_PK_WD"]=fns["CLONE_PK_WD"],
        _["CLONE_dP0_dr_WD"]=fns["CLONE_dP0_dr_WD"],
        _["CLONE_dPK_dr_WD"]=fns["CLONE_dPK_dr_WD"]
      );
//       fcts=ListOf<Function>(Environment::global_env()[".integrands"]);

//       Rcpp_PreserveObject(as<SEXP>(fcts));
//       fcts=fns;
//       std::cout<<"Size de fns="<<fns.size()<<std::endl;
//       std::cout<<"Size de fcts="<<fcts.size()<<std::endl;
      std::cout<<"toto"<<fcts.size()<<std::endl;
      mIntegrand=new Function("identity");  // allouer
      std::cout<<"toto2"<<fcts.size()<<std::endl;

}



void MATH_Integration::setFunctionName(std::string name) {
  std::cout<<"reltol=<"<<mReltol<<">"<<std::endl;
      std::cout<<"Size de fcts="<<fcts.size()<<std::endl;
//   Function func(fcts[name]);
  *mIntegrand = fcts[name];
}

// List MATH_Integration::getFns(){
//  return fcts;
// }


double MATH_Integration::testintegralFunction(double a, double b,double rho,double delta) {
  std::string name="CLONE_P0_WD";
  setFunctionName(name);
  List res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["rho"] =rho,_["delta"] =delta);
  double integ=as<double>(res["value"]);
  std::cout<<integ<<std::endl;

  name="CLONE_dP0_dr_WD";
  setFunctionName(name);
  res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["rho"] =rho,_["delta"] =delta);
  integ=as<double>(res["value"]);

  return integ;
}


double MATH_Integration::integralFunction(double a, double b,double rho,double delta,int k) {

  List res;
  if(k > 0) res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["rho"] =rho,_["delta"] =delta,_["k"]=k);
  else res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["rho"] =rho,_["delta"] =delta);
  double integ=as<double>(res["value"]);

  return integ;
}
