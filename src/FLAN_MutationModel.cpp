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
#include "FLAN_MutationModel.h"
    // --------------------------
    // Probability methods
    //---------------------------
    

std::vector<double> FLAN_MutationModel::computeProbability(int m) {
    
    //vector of clone probabilities
    std::vector<double> P(m+1);
    
    P=mClone->computeProbability(m);
    
    return deduceProbability(m,P);
}


std::vector<double> FLAN_MutationModel::deduceProbability(int m,std::vector<double>& P) {
    
    // initialize Q
    std::vector<double> Q(m+1);
    
//     for(std::vector<double>::iterator it=Q.begin(); it != Q.end(); ++it) *it=0;
//     for (int k=0;k<=m;k++) Q[k]=0;
    
    //initial probability q0
    Q[0]=exp(-mMutNumber*(1.-P[0]));
    
    if(m == 0) return Q;
    
    double s=0;
    for (int k=1;k<=m;k++) {
        // convolution
        s=0;
        for (int i=1;i<=k;i++) {
            s+=i*P[i]*Q[k-i];
         
        }
        Q[k]=(mMutNumber/k)*s;
    }
    
    return Q;
}


List FLAN_MutationModel::computeProbability1DerivativeAlpha(int m
// 							    std::vector<double>& Q,
// 							    std::vector<double>& dQ_da
									  ) {
    
                                                 
    
    std::vector<double> P(m+1);
    // compute the probabilities of mClone
    P=mClone->computeProbability(m);
    
    return deduceProbability1DerivativeAlpha(m,P);
}


List FLAN_MutationModel::deduceProbability1DerivativeAlpha(int m,
							    std::vector<double>& P
// 							    std::vector<double>& Q,
// 							    std::vector<double>& dQ_da
							  ) {
    
                                                    
    //initial probability Q0
    std::vector<double> Q(m+1);
    std::vector<double> dQ_da(m+1);
//     Q.resize(m+1);
//     dQ_da.resize(m+1);
    
//     for(std::vector<double>::iterator it=Q.begin(); it != Q.end(); ++it) *it=0;
    
    //initial probability q0
    Q[0]=exp(-mMutNumber*(1.-P[0]));
        // first derivatives of Q with respect to mMutNumber
    dQ_da[0]=-(1-P[0])*Q[0];
   
//     if (m == 0) return true;
    if (m == 0) return List::create(_["Q"]=Q[0],
				    _["dQ_da"]=dQ_da[0]
				   );
    
    double s,ds_da; 
    for (int k=1;k<=m;k++) {
        s=0;ds_da=0;
        for (int i=1;i<=k;i++) {
            s+=i*P[i]*Q[k-i];
	    ds_da+=P[i]*Q[k-i];
        }
	Q[k]=(mMutNumber/k)*s;
        dQ_da[k]=ds_da-Q[k];
    }
    return List::create(_["Q"]=NumericVector(Q.begin(),Q.end()),
			_["dQ_da"]=NumericVector(dQ_da.begin(),dQ_da.end())
		       );
//     return true;
}


List FLAN_MutationModel::computeProbability1DerivativeRho(int m
// 							  std::vector<double>& Q,
// 							  std::vector<double>& dQ_dr
									) {
    
                                                 
    
    std::vector<double> P(m+1);
    std::vector<double> dP_dr(m+1);
    
    // compute the probabilities and their derivative with respect to mFitness of mClone
    P=mClone->computeProbability(m);
    dP_dr=mClone->computeProbability1DerivativeRho(m,P);
    
    // compute the probabilities Q

    return deduceProbability1DerivativeRho(m,P,dP_dr);
}




List FLAN_MutationModel::deduceProbability1DerivativeRho(int m,
							  std::vector<double>& P,
							  std::vector<double>& dP_dr
// 							  std::vector<double>& Q,
// 							  std::vector<double>& dQ_dr
								       ) {
    
                                                    
    //initial probability Q0

    std::vector<double> Q(m+1);
    std::vector<double> dQ_dr(m+1);
//     Q.resize(m+1);
//     dQ_dr.resize(m+1);
    
    // first derivatives of Q with respect to mFitness
    Q[0]=exp(-mMutNumber*(1.-P[0]));
    dQ_dr[0]=mMutNumber*(dP_dr[0])*Q[0];
    
//     if (m==0) return true;
    if (m==0)  return List::create(_["Q"]=Q[0],
				   _["dQ_dr"]=dQ_dr[0]
				  );
    
    double s=0,ds_dr=0;
    
    for (int k=1;k<=m;k++) {
        s=0;ds_dr=0;
        for (int i=1;i<=k;i++) {
            s +=i*P[i]*Q[k-i];
	    ds_dr+=dP_dr[i]*Q[k-i];
//             ds_dr+=i*dP_dr[i]*Q[k-i]+i*P[i]*dQ_dr[k-i];
        }
        Q[k] =(mMutNumber/k)*s;
	dQ_dr[k]=mMutNumber*ds_dr;
//         dQ_dr[k]=(mMutNumber/k)*ds_dr;

    }
    return List::create(_["Q"]=NumericVector(Q.begin(),Q.end()),
			_["dQ_dr"]=NumericVector(dQ_dr.begin(),dQ_dr.end())
			);
//     return true;
}



List FLAN_MutationModel::computeProbability1DerivativesAlphaRho(int m
// 								std::vector<double>& Q,
// 								std::vector<double>& dQ_da,
// 								std::vector<double>& dQ_dr
							       ) {
    
                                                 
    // compute the derivative with respect to mFitness of GY
    
    std::vector<double> P;
    std::vector<double> dP_dr;
    P=mClone->computeProbability(m);
    dP_dr=mClone->computeProbability1DerivativeRho(m,P);
    
    
    return deduceProbability1DerivativesAlphaRho(m,P,dP_dr);
}


List FLAN_MutationModel::deduceProbability1DerivativesAlphaRho(int m,
								std::vector<double>& P,
								std::vector<double>& dP_dr
// 								std::vector<double>& Q,
// 								std::vector<double>& dQ_da,
// 								std::vector<double>& dQ_dr
							      ) {
    
                                                    
    //initial probability Q0
    std::vector<double> Q(m+1);
    std::vector<double> dQ_da(m+1);
    std::vector<double> dQ_dr(m+1);
//     Q.resize(m+1);
//     dQ_da.resize(m+1);
//     dQ_dr.resize(m+1);
    
    
    Q[0]=exp(-mMutNumber*(1-P[0]));
    // first derivatives of Q with respect to mMutNumber
    dQ_da[0]=-(1-P[0])*Q[0];
    
    // first derivatives of Q with respect to mFitness
    dQ_dr[0]=mMutNumber*(dP_dr[0])*Q[0];

//     if (m==0) return true;
    if (m==0) return List::create(_["Q"]=Q[0],
				  _["dQ_da"]=dQ_da[0],
				  _["dQ_dr"]=dQ_dr[0]
				 );
    
    double s,ds_da,ds_dr;
    
    for (int k=1;k<=m;k++) {
        s=0;ds_da=0;ds_dr=0;
        for (int i=1;i<=k;i++) {
            s +=i*P[i]*Q[k-i];
	    ds_da+=P[i]*Q[k-i];
	    ds_dr+=dP_dr[i]*Q[k-i];
        }
        Q[k]=(mMutNumber/k)*s;
	dQ_da[k]=ds_da-Q[k];
	dQ_dr[k]=mMutNumber*ds_dr;
    }
    return List::create(_["Q"]=NumericVector(Q.begin(),Q.end()),
			_["dQ_da"]=NumericVector(dQ_da.begin(),dQ_da.end()),
			_["dQ_dr"]=NumericVector(dQ_dr.begin(),dQ_dr.end())
	    );
}


  

std::vector<double> FLAN_MutationModel::computeCumulativeFunction(int m) {
		
    std::vector<double> cumsum(m+1);
    
    std::vector<double> Q(m+1);
    
    Q=computeProbability(m);
    
    std::partial_sum(Q.begin(),Q.end(),cumsum.begin(),std::plus<double>());
    
    //initial probability q0
//     cumsum[0]=Q[0];
// //     sum_pk+=pk[0];
//     
//     for (int k=1;k<=m;k++) {
//       cumsum[k]=cumsum[k-1]+Q[k];
//     }
//     
    if(mLT) {
      for(std::vector<double>::iterator it=cumsum.begin();it!=cumsum.end();++it) {
	(*it)*=-1;
	(*it)+=1;
      }
    }
    
    return cumsum;
}


// GF method


// List FLAN_MutationModel::MutationGFEstimation(std::vector& mutantsCount){
//   
//   
//   
//   
// }




// double FLAN_MutationModel::computeGeneratingFunction(double z) {
//     return exp(mMutNumber*(mClone->computeGeneratingFunction(z)-1));
// }
// 
// 
// double FLAN_MutationModel::covariance(double z1, double z2) {
//     return computeGeneratingFunction(z1*z2)-
//         computeGeneratingFunction(z1)*
//         computeGeneratingFunction(z2);
// }
// 
// void FLAN_MutationModel::covariance(double z1, double z2,double z3,
// 				    double cov[2][2]) {
// 
// //     double cvfn2=cvfn*cvfn;
//     
//     double M[3][3];
//     M[0][0]=covariance(z1,z1);
//     M[0][1]=covariance(z1,z2);
//     M[0][2]=covariance(z1,z3);
// 
//   
//     M[1][1]=covariance(z2,z2);
//     M[1][2]=covariance(z2,z3);
// 
//     M[2][2]=covariance(z3,z3);
//     for (int i=0;i<3;i++) {
//         for (int j=0;j<i;j++) {
//             M[i][j]=M[j][i];
//         }
//     }
// 
//     double ccf_z1=mClone->computeGeneratingFunction(z1);
//     double ccf_z2=mClone->computeGeneratingFunction(z2);
//     double ccf_z3=mClone->computeGeneratingFunction(z3);
// 
//     double dccf_z1=mClone->computeGeneratingFunction1DerivativeRho(z1);
//     double dccf_z2=mClone->computeGeneratingFunction1DerivativeRho(z2);
//     double dccf_z3=mClone->computeGeneratingFunction1DerivativeRho(z3);
//     
// 
//     double d01=(ccf_z2-1)*dccf_z1-(ccf_z1-1)*dccf_z2;
//     double d11=(ccf_z1-1)*dccf_z2-(ccf_z2-1)*dccf_z1;
//     
//     
//     double CO[3][2];
//     CO[0][1]= (ccf_z2-1)/(mMutNumber*d01*computeGeneratingFunction(z1));
//     
//     CO[1][1]=(ccf_z1-1)/(mMutNumber*d11*computeGeneratingFunction(z2));
//     CO[2][1]= 0;
//     
//     CO[0][0]=(mMutNumber*dccf_z3*CO[0][1])/(1-ccf_z3);
//     CO[1][0]=(mMutNumber*dccf_z3*CO[1][1])/(1-ccf_z3);
//     CO[2][0]=1./(computeGeneratingFunction(z3)*(ccf_z3-1));
//     
// //     cov.setSize(2,2);
//     double d, c;
//     for (int i=0;i<2;i++) {
//         for (int j=0;j<2;j++) {
//             d=0;
//             for (int k=0;k<3;k++) {
//                 c=0;
//                 for (int s=0;s<3;s++) {
//                     c+=M[k][s]*CO[s][j];
//                 }
//                 d+=CO[k][i]*c;
//             }
//             cov[i][j]=d;
//         }
//     }
//     
// }
// 
// 
//  // Unbiased estimation of pi and its standart deviation if fluctuation of final counts
// void FLAN_MutationModel::unbiasPiEstimation(double z,double cvfn,
//                                               double pi,double sd_pi)  {
//     
//     double umds=1-mDeath/(1-mDeath);
// //     double f=mMutNumber*umds*cvfn*cvfn;
//     double f=mMutNumber*(1-mClone->computeGeneratingFunction(z))*cvfn*cvfn;
//     pi*=1+f/2;
//     sd_pi*=1+f;
// }
