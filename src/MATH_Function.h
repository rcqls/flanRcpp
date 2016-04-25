#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

// #include "FLAN_Function.h"
#include <Rcpp.h>
using namespace Rcpp;

class MATH_Integration {

private:

  double mReltol;

  Function* mIntegrate;

  Function* mIntegrand;

  List fcts;


protected:

    // inline void setIntervalsNumber(const int& n) {
    //   mIntervalsNumber=n;
    // };
    //
    //
    // inline void setFunctionCalledNumber(const int& n) {
    //     mFunctionCalledNumber=n;
    // };
    //
    // inline double computeFunction(const double& x) {
    //     return getFunction()->computeFunction(x);
    // };



public:

    MATH_Integration() {};
    MATH_Integration(double reltol_);

    ~MATH_Integration(){};

List getFns(){return fcts;};

    void setFunctionName(std::string name);

  /*
   * Integrals functions
   */

  double testintegralFunction(double a, double b,double rho, double delta);
  double integralFunction(double a, double b,double rho, double delta,int k);

};

#endif
