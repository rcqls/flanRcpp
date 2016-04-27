
#include "MATH_Function.h"

#include <math.h>

using namespace Rcpp;


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


double MATH_Integration::integralFunction(double a, double b,double rho,double delta,double k) {

  List res;
  if(k > 0) res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["rho"] =rho,_["delta"] =delta,_["k"]=k);   
  else res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["rho"] =rho,_["delta"] =delta);
  double integ=as<double>(res["value"]);
//   double integ=res["value"];
  return integ;
}


void MATH_Polynom::square(){
  
    int n=mCoef.size();
    std::vector <double> Sq(2*n-1);
 
   // Take ever term of first polynomial
   int i=0;
//    for (int i=0; i<m; i++)
   for(std::vector<double>::iterator it1=mCoef.begin() ; it1!=mCoef.end() ; it1++,i++) {
     // Multiply the current term of first polynomial
     // with every term of second polynomial.
     int j=0;
     for(std::vector<double>::iterator it2=mCoef.begin() ; it2!=mCoef.end() ; it2++,j++) {
         Sq[i+j] += (*it1)*(*it2);
      }
   }
   setCoef(Sq);
}






