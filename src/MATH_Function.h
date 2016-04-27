#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

#include <Rcpp.h>
using namespace Rcpp;

class MATH_Integration {

private:

  double mReltol;

  Function* mIntegrate;

  Function* mIntegrand;

  List fcts;


protected:


public:

    MATH_Integration() {};
    MATH_Integration(List fns, double reltol_){
      mReltol=reltol_;
      mIntegrate=new Function("integrate");
//       List fns=Environment::global_env().get(".integrands");
      fcts=fns;

      mIntegrand=new Function("identity");  // allouer


    };

    ~MATH_Integration(){};

// List getFns(){return fcts;};

    void setFunctionName(std::string name){
    //   std::cout<<"reltol=<"<<mReltol<<">"<<std::endl;
    //       std::cout<<"Size de fcts="<<fcts.size()<<std::endl;
    //   Function func(fcts[name]);
      *mIntegrand = fcts[name];
    };

  /*
   * Integrals functions
   */

  double testintegralFunction(double a, double b,double rho, double delta);
  double integralFunction(double a, double b,double rho, double delta,double k);

};



class MATH_Polynom {
  
  
protected:
   std::vector<double> mCoef;


public:

    MATH_Polynom() {
      std::vector<double> C(1);
      C[0]=0;
      mCoef=C;
    };
    
    MATH_Polynom(int deg) {
      std::vector<double> C(deg);
      std::vector<double>::iterator itC;
      for(itC=C.begin() ; itC != C.end()-1 ; ++itC) *itC=0;
	++itC;
	*itC=1;
	mCoef=C;
    };
    
    MATH_Polynom(std::vector<double> C) {
      mCoef=C;
    };
    
    ~MATH_Polynom(){};
    
    
    void setCoef(std::vector<double> C){
      mCoef=C;
    }
    
    int getDegree(){
      return mCoef.size()-1;
    };
    
    
    double& operator[](int i) {
      return mCoef[i];
    };
    
    MATH_Polynom& operator*=(double f) {
      int i=0;
      for(std::vector<double>::iterator it=mCoef.begin() ; it!=mCoef.end() ; ++it,i++) {
	(*it)*=f;
      }
    };
    
    MATH_Polynom& operator+=(double f) {
      mCoef[0]+=f;
    };
    
    void reduce(double eps){
      int i=0;
      int dmax=0;
      for(std::vector<double>::iterator it=mCoef.begin() ; it!=mCoef.end() ; ++it,i++) {
	if(*it <= eps) *it=0;
	if(*it > 0) dmax=i;
      }
      mCoef.resize(dmax+1);
    };
    
   

    void square();


};


#endif
