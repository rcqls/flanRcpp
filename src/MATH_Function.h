#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

// #include "FLAN_Function.h"
#include <Rcpp.h>
using namespace Rcpp;

class MATH_Integration {

private:

  double reltol;

  Function* integrate;

  Function* integrand;

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
    MATH_Integration(List fns,double reltol_);

    ~MATH_Integration(){};


    void setFunctionName(std::string name);

  /*
   * Integrals functions
   */

  double integralFunction(double a, double b,double rho, double delta);


};

#endif
