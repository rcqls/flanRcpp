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

#include "FLAN_Sim.h"

using namespace Rcpp;

FLAN_Sim::FLAN_Sim() {}

FLAN_Sim::FLAN_Sim(List args){
      mMut=as<double>(args["mutations"]);
      mFitness=as<double>(args["fitness"]);
      mDeath=as<double>(args["death"]);
      
      List dist=args["dist"];
      FLAN_Dist* mDist=new FLAN_Dist(dist);          // Lifetime Distribution
      mDist->adjustGrowthRate(mDeath);   // Rescales parameter to unit growth rate
      
      mMfn=as<double>(args["mfn"]);
      mCvfn=as<double>(args["cvfn"]);

      mClone=new FLAN_SimClone(mFitness,mDeath,mDist);

}


FLAN_Sim::~FLAN_Sim() {}


    // --------------------------------------------------------------------------------------------------------------
    //////////////////////////////////////////// Samples computation ////////////////////////////////////////////////
    // --------------------------------------------------------------------------------------------------------------

    /* Generates a sample of size n of couples (mutant count ; final count)
     * The final counts are sampled with the log-normal
     * distribution with coefficient of variation cvfn and mean mfn.
     *
     * n: number of experiments
     * mfn: mean of final counts
     * Cfn: variation number of final counts
     * alphapi: mean number of mutations probability of mutation (depends of cvfn)
     * rho: fitness parameter
     * death: probability of death
     */


Rcpp::List FLAN_Sim::computeSamplesMutantsFinalsNumber(int n)  {

  RNGScope rngScope;
  
  
  std::vector<long int> mutantCount(n);
  std::vector<double> finalCount(n);

  if (mCvfn>0) {
      double sdLog2=log(1+mCvfn*mCvfn);
      double sdLog=sqrt(sdLog2);
      double meanLog=log(mMfn)-sdLog2/2;


      for(std::vector<double>::iterator it = finalCount.begin();
	  it != finalCount.end(); ++it) {


	    *it=R::rlnorm(meanLog,sdLog);
      }

      mutantCount=computeSampleMutantsNumber(n,finalCount);

  } else {
      for(std::vector<double>::iterator it = finalCount.begin();
	  it != finalCount.end(); ++it) {

	    *it=mMfn;
      }
      mutantCount=computeSampleMutantsNumber(n);
  }

  //export in R
  return Rcpp::List::create(
    _["mc"]=NumericVector(mutantCount.begin(),mutantCount.end()),
    _["fn"]=NumericVector(finalCount.begin(),finalCount.end())
  );


}



    /* Generates a sample of size n of mutant counts
     * n: sample size
     * alpha: mean number of mutations
     * rho: fitness parameter
     * death: probability of death
     */



std::vector<long int> FLAN_Sim::computeSampleMutantsNumber(int n)  {

    
    std::vector<long int> mutantCount(n);
//     FLAN_Clone *clone=new FLAN_Clone(mFitness,mDeath,mDist);

    // poisson number of mutation

    std::vector<double> sample;
    double s;
    // set the size of the mutants
    double mcs;
    int mc;
    for (std::vector<long int>::iterator it = mutantCount.begin();
	  it != mutantCount.end(); ++it) {
//       int i=0;i<n;i++) {
        // simulate a poisson number

	  mcs=R::rpois(mMut);
	  mc=(int)(mcs);

        if (mc>0) {
            //Poisson sum of Yule variables
            sample=mClone->computeSample(mc);

            //s = sum_j sample[j]
//             int m=sample.getSize();
            s=0;
	    int j=0;
	    bool testneg=false;
//             for (int j=0;j<m;j++){
	    while(j<mc && !testneg){
	      if(sample[j] < 0){
		testneg=true;
		s=-10^5;
	      } else {
		s+=sample[j];
		j++;
	      }
	    }
	    *it=s;
        } else
            *it=0;
    }
    
    return mutantCount;
}


    /* Generates a sample of size n of mutant counts, given a sample of final counts.
     *
     * n: number of experiments
     * pi: mutation probability
     * rho: fitness parameter
     * death: probability of death
     * finalCount:
     */

std::vector<long int> FLAN_Sim::computeSampleMutantsNumber(int n,
					  const std::vector<double>& finalCount)  {

    std::vector<long int> mutantCount(n);
//     FLAN_SimClone* clone=new FLAN_SimClone(mFitness,mDeath,mDist);

    // poisson number of mutation

    std::vector<double> sample;
    double s;
    // set the size of the mutants
//     mutantCount.resize(n);

    double mcs;
    int mc;
    int i=0;
//     RNGScope rngScope;

    for (std::vector<long int>::iterator it = mutantCount.begin();
	  it != mutantCount.end(); ++it, i++) {
        // simulate a poisson number

      mcs=R::rpois(mMut*finalCount[i]);
      mc=(int)(mcs);
// 	  std::cout<<"mc ="<<mc<<std::endl<<std::endl;

        if (mc>0) {
            //Poisson sum of Yule variables
            sample=mClone->computeSample(mc);

            //s = sum_j sample[j]
//             int m=sample.getSize();
            s=0;
	    int j=0;
	    bool testneg=false;
//             for (int j=0;j<m;j++){
	    while(j<mc && !testneg){
	      if(sample[j] < 0){
		testneg=true;
		s=-10^5;
	      } else {
		s+=sample[j];
		j++;
	      }
	    }
	    *it=s;
        } else
            *it=0;
    }
    
    return mutantCount;

}
