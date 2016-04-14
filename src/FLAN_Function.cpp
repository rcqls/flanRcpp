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
#include "FLAN_Function.h"

const unsigned char FLAN_Function::CF_GY=1;
const unsigned char FLAN_Function::CF_GY_WD=2;
const unsigned char FLAN_Function::CF_GY_WD_1=3;
const unsigned char FLAN_Function::CF_GY_WD_K=4;
const unsigned char FLAN_Function::CF_GY_WD_DR=5;
const unsigned char FLAN_Function::CF_GY_WD_1_DR=6;
const unsigned char FLAN_Function::CF_GY_WD_K_DR=7;
const unsigned char FLAN_Function::CF_GY_WD_DD=8;
const unsigned char FLAN_Function::CF_GY_WD_1_DD=9;
const unsigned char FLAN_Function::CF_GY_WD_K_DD=10;
const unsigned char FLAN_Function::CFD_GY=11;
const unsigned char FLAN_Function::LOG_NORMAL=12;


FLAN_Function::FLAN_Function() {
  mRho=0;
  mZ=0;
  mK=0;
  mDeath=0;
}

 
FLAN_Function::~FLAN_Function() {
}


double FLAN_Function::computeFunction(const double& x) {
    double a=1;
    double umz=1-mZ;
    double zstar=mZ/umz;
    switch(mFunctionName) {
    case CF_GY:
	zstar=(mZ-mDeath/(1-mDeath))/umz;
	return pow(x,mRho)/(1+zstar*x);
// 	a=(1.-mDeath)*mZ-mDeath;
// 	return pow(x,mRho-1.)*((mDeath*umz+a*x)/((1.-mDeath)*umz+a*x));
        break;
    case CFD_GY:
      zstar=(mZ-mDeath/(1-mDeath))/umz;
      return pow(x,mRho)*log(x)/(1+zstar*x);
//         a=(1-mDeath)*mZ-mDeath;
//         return (1+mRho*log(x))*pow(x,mRho-1)*((mDeath*umz+a*x)/((1-mDeath)*umz+a*x));     
        break;
    case LOG_NORMAL:
        return exp(mRho*x+mZ)*exp(-x*x/2)*(1./sqrt(2*M_PI));
    case CF_GY_WD:
        return (1.-x)*pow(x,mRho-1)/(1.-mZ*x);
    case CF_GY_WD_DR:
        return log(x)*(1.-x)*pow(x,mRho-1)/(1.-mZ*x);
    case CF_GY_WD_DD:
        a=(1.-mZ*x);
        return (1.-x)*pow(x,mRho)/(a*a);
    case CF_GY_WD_1:
        a=(1.-mZ*x);
        return pow(x,mRho)/(a*a);
    case CF_GY_WD_1_DR:
        a=(1.-mZ*x);
        return log(x)*pow(x,mRho)/(a*a);
    case CF_GY_WD_1_DD:
        a=(1.-mZ*x);
        return pow(x,mRho+1)/(a*a*a);    
    case CF_GY_WD_K:
        return pow(x,mRho)*pow(1-x,mK-1)/pow(1.-mZ*x,mK+1);
    case CF_GY_WD_K_DR:
        return log(x)*pow(x,mRho)*pow(1-x,mK-1)/pow(1.-mZ*x,mK+1);
    case CF_GY_WD_K_DD:
        return pow(x,mRho+1)*pow(1-x,mK-1)/pow(1.-mZ*x,mK+2);
    default:
        return 0;
    }
}
