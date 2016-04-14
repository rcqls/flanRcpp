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

#ifndef FLAN_SIM_H
#define FLAN_SIM_H


#include "FLAN_Clone.h"
using namespace Rcpp ;


class FLAN_Sim {
public:

    /*! \brief  create an object */

    FLAN_Sim();

    FLAN_Sim(List args);

    // DESTRUCTORS

    /*! \brief  destroy an object.
     */
    ~FLAN_Sim();



private:

  double mMut;
  double mFitness;
  double mDeath;
  FLAN_SimClone *mClone;        
  double mMfn;
  double mCvfn;

    // --------------------------------------------------------------------------------------------------------------
    //////////////////////////////////////////// Samples computation ////////////////////////////////////////////////
    // --------------------------------------------------------------------------------------------------------------
    /* Generates a sample of size n of mutant counts
     * n: sample size
     * alpha: mean number of mutations
     * rho: fitness parameter
     * death: probability of death
     */

    std::vector<long int> computeSampleMutantsNumber(int n) ;

    /* Generates a sample of size n of mutant counts, given a sample of final counts.
     *
     * n: number of experiments
     * pi: mutation probability
     * rho: fitness parameter
     * death: probability of death
     * finalCount:
     */
    std::vector<long int> computeSampleMutantsNumber(int n,
				    const std::vector<double>& finalCount) ;
public:
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
    List computeSamplesMutantsFinalsNumber(int n)  ;

};

#endif
