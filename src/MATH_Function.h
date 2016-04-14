#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

#include "FLAN_Function.h"


class MATH_Integration {

private:
  
  int mFunctionCalledNumber;
  int mIntervalsNumber;
  
  FLAN_Function *mFunction;
  
  
protected:
  
    inline void setIntervalsNumber(const int& n) {
      mIntervalsNumber=n;
    };
    
    
    inline void setFunctionCalledNumber(const int& n) {
        mFunctionCalledNumber=n;
    };
    
    inline double computeFunction(const double& x) {
        return getFunction()->computeFunction(x);
    };
  
    
  
public:
  
    MATH_Integration();
    
    ~MATH_Integration(){};
  
    inline void setFunction(FLAN_Function *f) {
        mFunction=f;
    }; 
    
    inline FLAN_Function* getFunction() const {
        return mFunction;
    };
    
    inline int getFunctionCalledNumber() const {
        return  mFunctionCalledNumber;
    };
    /*! \brief get intervals number
     */
    inline int getIntervalsNumber() const {
        return  mIntervalsNumber;
    };
    
    
  
  /*
   * Integrals functions
   */
  
  static void integralFunction(double x, double xminusa, double bminusx, double &y, void *ptr);

  bool integrate(const double& a,const double& b,double& v);
};

#endif