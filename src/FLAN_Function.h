/**
  * FLAN Software
  *
  * @author 2015-2020 Stephane Despreaux <stephane.despreaux@imag.fr>
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
#ifndef FLAN_FUNCTION_H
#define FLAN_FUNCTION_H

#include <Rcpp.h>

// DEFINE_SPTR(MRE_Function);

class FLAN_Function {

//   SP_OBJECT(FLAN_Function);

    // ATTRIBUTES

public:
    static const unsigned char CF_GY;
    static const unsigned char CFD_GY;
    static const unsigned char LOG_NORMAL;
    static const unsigned char CF_GY_WD;
    static const unsigned char CF_GY_WD_1;
    static const unsigned char CF_GY_WD_K;
    static const unsigned char CF_GY_WD_DR;
    static const unsigned char CF_GY_WD_1_DR;
    static const unsigned char CF_GY_WD_K_DR;
    static const unsigned char CF_GY_WD_DD;
    static const unsigned char CF_GY_WD_1_DD;
    static const unsigned char CF_GY_WD_K_DD;

private:

    unsigned char mFunctionName;
    double mZ;
    double mRho;
    double mDeath;
    int mK;

    // ASSOCIATIONS


    // CONSTRUCTORS
protected:

public:
  
    FLAN_Function();

    // DESTRUCTORS


    /* destroy a function.
     */
    ~FLAN_Function();
    
    /* set the type of the function
     */
    inline void setFunctionName(const unsigned char name) {
        mFunctionName=name;
    }

    inline void setParameter(double z,double rho,double death) {
        mZ=z;
        mRho=rho;
        mDeath=death;
    };
    inline void setParameter(double z,double rho) {
        mZ=z;
        mRho=rho;
        mDeath=0;
    };
    inline void setParameter(int k) {
        mK=k;
    };

    /* compute the function to optimize at parameters value
     */
    double computeFunction(const double& x);


};
#endif
