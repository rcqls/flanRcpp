/**
  * FLAN Software
  *
  * @author 2015-2020 Adrien Mazoyer  <adrien.mazoyer@imag.fr> 
  * @see The GNU Public License (GPL)
  */
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
// #include <Rcpp.h>

#ifndef FLAN_MUTATION_MODEL_H
#define FLAN_MUTATION_MODEL_H

#include "FLAN_Clone.h"
using namespace Rcpp ;


class FLAN_MutationModel {
  
private:
    
    FLAN_Clone* mClone;    // Clone object
    
protected:
    
    double mMutNumber;     // Mean number of mutations
    double mFitness;       // Relative fitness
    double mDeath;         // Death probability
    
    bool mLT;              // logical: if TRUE probabilities are P[X <= x]
			   //   otherwise, P[X > x]
    
    double mZ4;            // tuning parameter to unbias estimation
    
public: 

    //  Create an object 
    FLAN_MutationModel();
    
    // Create object for GF method
    FLAN_MutationModel(double death,std::string model);
    
    // Create object for GF method
    FLAN_MutationModel(double rho,double death,std::string model);
    
    // Create object for distribution
    FLAN_MutationModel(List args);
    
    // Destructor
    ~FLAN_MutationModel();   
    
    
    // Set attributes
    inline void setMutNumber(double alpha) {
      mMutNumber=alpha;
    };
    inline void setFitness(double rho){
      mFitness=rho;
    };
    inline void setDeath(double death){
      mDeath=death;
    };
    inline void setClone(FLAN_Clone* clone){
      mClone=clone;
    };
    
    List getFns(){
      return mClone->get();
}
    
    // Get attributes
    
    inline double getMutNumber(){
      return mMutNumber;
    };
    inline double getFitness(){
      return mFitness;
    };
    inline double getDeath(){
      return mDeath;    
    };
    inline FLAN_Clone* getClone(){
      return mClone;
    };
    
    
    
    // --------------------------
    // Probability methods
    //---------------------------
    
    
    std::vector<double> computeProbability(int m) ;

    std::vector<double> deduceProbability(int m,std::vector<double>& pClone) ;
			
			    
    bool computeProbability1DerivativeAlpha(int m,
					    std::vector<double>& Q,
					    std::vector<double>& dQ_da)  ;    
    bool deduceProbability1DerivativeAlpha(int m,
					    std::vector<double>& pClone,
					    std::vector<double>& Q,
					    std::vector<double>& dQ_da)  ;
    
					      
    bool computeProbability1DerivativeRho(int m,
					  std::vector<double>& Q,
					  std::vector<double>& dQ_dr)  ;
    bool deduceProbability1DerivativeRho(int m,
					  std::vector<double>& pClone,
					  std::vector<double>& dpClone_r,
					  std::vector<double>& Q,
					  std::vector<double>& dQ_dr) ;
    
					    
    bool computeProbability1DerivativesAlphaRho(int m,
						std::vector<double>& Q,
						std::vector<double>& dQ_da,
						std::vector<double>& dQ_dr)  ;
    bool deduceProbability1DerivativesAlphaRho(int m,
						std::vector<double>& pClone,
						std::vector<double>& dpClone_r,
						std::vector<double>& Q,
						std::vector<double>& dQ_da,
						std::vector<double>& dQ_dr)  ;
    
//     void computeProbability1Derivatives(int m,double alpha,double rho,double death) ;
    
    std::vector<double> computeCumulativeFunction(int m) ;
    
    
    
    //  compute the Generating function and its derivative
//     double computeGeneratingFunction(double z)  ;
//     double computeGeneratingFunctionDerivative(double z)  ;
//    
//     // --------------------
//     // GF ESTIMATION covariance methods
//     // -------------------
// 
//     double covariance(double z1, double z2) ;
//     
//     void covariance(double z1,double z2,double z3,
//                     double cov[2][2]) ;
//     
// 		    
//     // Unbiased estimation of pi and its standart deviation if fluctuation of final counts
//     void unbiasPiEstimation(double z,double cvfn,
// 			      double pi,double sd_pi);

};
#endif
