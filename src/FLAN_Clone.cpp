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


#include "FLAN_Clone.h"


FLAN_SimClone::FLAN_SimClone(){}

FLAN_SimClone::FLAN_SimClone(double rho,double death, FLAN_Dist *dist){

  mFitness=rho;
  mDeath=death;
  mDist=dist;

}

FLAN_SimClone::~FLAN_SimClone(){}

const double FLAN_SimClone::DEATH_EPS_SIM=1.e-4;

std::vector<double> FLAN_SimClone::computeSample(int n) {

  std::vector<double> sample(n);

  std::string name=mDist->getDistName();

  for(std::vector<double>::iterator it = sample.begin();
	  it != sample.end(); ++it) {
      //mExponentialD->setProperties(mFitness);
      //mExponentialD->computeSample(n,sample);
      *it=R::rexp(mFitness);
  }

  if(name.compare("dirac") == 0){

    // specialized method
    if (mDeath<DEATH_EPS_SIM) {
        // case without mDeath
//         for (int k=0 ; k<n ; k++) {
	for(std::vector<double>::iterator it = sample.begin();
	  it != sample.end(); ++it) {
// 	  sample[k] = pow(2,floor(sample[k]/log(2)));
	  *it = pow(2.,floor(*it/log(2.)));
        }
    } else {
        // case with mDeath
        double t=0,a=log(2*(1-mDeath));
//         for (int k=0 ; k<n ; k++) {
	for(std::vector<double>::iterator it = sample.begin();
	  it != sample.end(); ++it) {

	  double cg=0;
	  int m=1;
	  int nAlives=0,nDeaths=0;

// 	  t=sample[k];
	  t=*it;
            while ((cg<t-a) && (m>0)) {
//                 mBinomialD->setProperties(m,mDeath);
                // number of mDeaths in [0,m]
//                 nDeaths=mBinomialD->random();
		nDeaths=R::rbinom(m,mDeath);
                // the number of alives among m
                nAlives=m-nDeaths;
                // the new alive cells
                m+=nAlives-nDeaths;
                // current generation
                cg+=a;
            }
//             std::cout<<"cg final ="<<cg<<" | m final ="<<m<<" | Tps de mort (t-a) = "<<t-a<<std::endl;
//             sample[k] =m;
	      *it=m;
        }
    }
   // end if dirac distribution
  } else if(name.compare("exp") == 0) {

      double p,up;
      double geoD,beD;

      if (mDeath<DEATH_EPS_SIM) {
  //         for (int k=0;k<n;k++) {
	for(std::vector<double>::iterator it = sample.begin();
	    it != sample.end(); ++it) {
	      p=exp(-(*it));

	      geoD=R::rgeom(p);

	      if(geoD>=0) *it=geoD+1;
  // 	    else {
  // // 	      std::cout<<"Geom négative avec param p="<<p<<endl;
  // // 	      std::cout<<"Rappel p=exp(-sample) et sample="<<sample[k]<<endl;
  // 	      sample[k]="NaN";
  // 	    }

  // 	    if(sample[k]<0) std::cout<<sample[k]<<endl;

	  }
      } else {
	for(std::vector<double>::iterator it = sample.begin();
	    it != sample.end(); ++it) {

	      p=exp(-(*it));
	      up=(1-2*mDeath)/(1-mDeath*(1+p));

	      beD=R::rbinom(1,up);

  //             mBernoulliD->setProperties(up);
  //             mBernoulliD->computeSample(1,beD);

	      if (beD==1){

		  geoD=R::rgeom(p*up);

  // 	      if(geoD[0]<0) {
  // 		std::cout<<"Geom négative avec param p="<<p<<endl;
  // 		std::cout<<"Rappel p=exp(-sample) et sample="<<sample[k]<<endl;
    // 	      std::cout<<"tel quel"<<geoD[0]<<endl;
    // 	      std::cout<<"après castage"<<(int)(geoD[0])<<endl;
  // 	      }
		if(geoD>=0) *it=geoD+1;

	      } else *it=0;
  // 	  if(sample[k]<0) std::cout<<sample[k]<<endl;
	  }
      }
    // end if exponential distribution
  } else {

    // transform the sample

    std::vector<double> splitTimesList;
    int nc=0;
    for(std::vector<double>::iterator it = sample.begin();
	  it != sample.end(); ++it) {
	     nc=splitTimes(*it,splitTimesList);
        //sample[k]=splitTimesList.getSize()-1;
        *it=nc;
    }

  }

  return sample;
}


int FLAN_SimClone::splitTimes(double t,std::vector<double>& splitTimesList) {

    std::string name=mDist->getDistName();
    std::vector<double> params=mDist->getDistParams();

    // return the clone size
    int cloneSize=0;

    // initialize return vector
    splitTimesList.resize(1);

    //current generation
    int g=0;

    // make a sample of size 1
    //division times of current generation
    std::vector<double> st(1);

    if(name.compare("lnorm") == 0) st[0]=R::rlnorm(params[0],params[1]);
    if(name.compare("gamma") == 0) st[0]=R::rgamma(params[0],params[1]);



    // add initial time
//     splitTimesList.add(0);
    splitTimesList[0]=0;

    // take only the times if before t
    // ng is the number of simulated value less than t
    // number of alive cells
    int ng=0;

    std::vector<double> u(1);
//     if (mDeath>0) runif->computeSample(1,u);
    if (mDeath>0) u[0]=R::runif(0.,1.);

    if ((st[0]<t)&& ((mDeath==0) || ((mDeath>0) && (u[0]>mDeath)))) {
        // division at  t of an alive cell
//         splitTimesList.add(st[0]);
      splitTimesList.push_back(st[0]);
        ng=1;
    } else if (st[0]>=t) {
        // no division before t
        ng=0;
        cloneSize++;

    } else if ((mDeath>0) && (u[0]<=mDeath)) {
        // the cell dies
        ng=0;
    }



    std::vector<double> gtf;
    int nActiveAlives=ng;

    // main loop
    while ((ng>0) && (ng<1.e+6)) {
        //next generation
        g++;

        // compute a sample of size 2*ng
        //division times of daughters
//         distribution.computeSample(2*ng,gtf);
	gtf.resize(2*ng);

	for(std::vector<double>::iterator it=gtf.begin() ; it != gtf.end() ;
	    ++it) {
	  if(name.compare("lnorm") == 0) *it=R::rlnorm(params[0],params[1]);
	  if(name.compare("gamma") == 0) *it=R::rgamma(params[0],params[1]);
	}

        //  split times of daughters with duplicating split times
        for (int i=0;i<ng;i++) {
            gtf[2*i]+=st[i];
            gtf[2*i+1]+=st[i];
        }

        ng*=2;
        // keep only values less than t
        // daughters that still divide
        // keep their split times
//         if (mDeath>0) runif->computeSample(ng,u);
	if (mDeath>0) {
	  u.resize(ng);
	  for(std::vector<double>::iterator it=u.begin() ; it != u.end() ;
	    ++it) {
	    *it=R::runif(0.,1.);
	  }
	}

        st.resize(ng);
        nActiveAlives=0;
        for (int i=0;i<ng;i++) {

            if ((gtf[i]<t)&&
                ((mDeath==0) || ((mDeath>0) && u[i]>mDeath))) {
                // division at  t of an alive cell
                st[nActiveAlives]=gtf[i];
                nActiveAlives++;
            } else if (gtf[i]>=t) {
                // no division before t
                cloneSize++;

            }
        }
        // resize st to usefull length
        st.resize(nActiveAlives);


        // append st to spliTimes list
        //stack split times
//         splitTimesList.append(st);
	splitTimesList.insert(splitTimesList.end(), st.begin(), st.end());

        // ng is the size of st
        //dividing cells in next generation
        ng=nActiveAlives;
    }
    //reorder split times
//     splitTimesList.sort();

    std::sort(splitTimesList.begin(),splitTimesList.end());

    //return split times before t
//     splitTimesList.add(t);
    splitTimesList.push_back(t);

    return cloneSize;
}



  /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class for the distribution of a clone size
  */


const double FLAN_Clone::DEATH_EPS_DIST=1.e-4;

FLAN_Clone::FLAN_Clone() {
}


FLAN_Clone::FLAN_Clone(double death) {
  mDeath=death;
}

// create object for distribution
FLAN_Clone::FLAN_Clone(double rho, double death) {
  mFitness=rho;
  mDeath=death;


}


FLAN_Clone::~FLAN_Clone() {}



 /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class when the lifetime model is supposed to be exponential
  */


// create object for GF method



std::vector<double> FLAN_ExponentialClone::computeProbability(int m){

  std::vector<double> P(m+1);

  if(mDeath<DEATH_EPS_DIST){

    P[0]=0;
    if(m == 0) return P;

    int k=1;
    for(std::vector<double>::iterator it=P.begin()+1 ; it!= P.end() ; ++it,k++) *it=R::beta(mFitness+1,k);

  }
//    else {
//     double d1=mDeath/(1-mDeath);
//     int m_max=1000;
//     // set the function type
//     mFunction->setFunctionName(FLAN_Function::CF_GY_WD);
//     mFunction->setParameter(d1,mFitness);
//
//     //integrate the function in [0,1]
//     double i;
//     // m=0 probability
//     mIntegrator->integrate(0,1,i);
//     P[0]=i*d1*mFitness;
//
// //     sum_pk+=p_k;
//     if (m==0) return P;
//
//     //m=1 probability
//     double d2=(1.-2.*mDeath)/(1-mDeath);
//     d2*=d2;
//     mFunction->setFunctionName(FLAN_Function::CF_GY_WD_1);
//     mIntegrator->integrate(0.,1.,i);
//     P[1]=i*d2*mFitness;
// //     sum_pk+=p_k;
//     if (m==1) return P;
//
//     //m>1 probability
//     int m1=m;
//     if (m1>=m_max) m1=m_max;
//     int k=2;
// //     for (int k=2;k<=m1;k++) {
//     for(std::vector<double>::iterator it=P.begin()+2 ; it!= P.end() ; ++it, k++){
//         mFunction->setParameter(k);
//         mFunction->setFunctionName(FLAN_Function::CF_GY_WD_K);
//         mIntegrator->integrate(0.,1.,i);
//         *it=i*d2*mFitness;
// //         sum_pk+=p_k;
//     }
//
//     // equivalent computation
//     double a=pow(d2,(1.-mFitness)/2.)*mFitness*R::gammafn(mFitness+1);
//     k=m1+1;
//     for(std::vector<double>::iterator it=P.begin()+m1+1 ; it!= P.end() ; ++it, k++){
//     //     for (int k=m1+1;k<=m;k++) {
//         *it=a*pow(k,-mFitness-1);
// //         sum_pk+=p_k;
//     }
//
//  }

  return P;

}


std::vector<double> FLAN_ExponentialClone::computeProbability1DerivativeRho(int m,const std::vector<double> P){

  std::vector<double>dP_dr(m+1);

  if(mDeath<DEATH_EPS_DIST){

    dP_dr[0]=0;
    if(m == 0) return dP_dr;

    int k=1;
    double dg=R::digamma(mFitness+1);
    for(std::vector<double>::iterator it=dP_dr.begin()+1 ; it!= dP_dr.end(); ++it,k++)
      *it=P[k]*(1/mFitness+dg-R::digamma(mFitness+k+1));

  } else {
    double dstar=mDeath/(1-mDeath);

    //                  [.....]
  }

  return dP_dr;

}




/*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class when the lifetime model is suppsoed to be constant
  */

std::vector<double> FLAN_DiracClone::computeProbability(int m){

  std::vector<double> P(m+1);

  if(mDeath<DEATH_EPS_DIST){

    P[0]=0;

    if(m == 0) return P;

    for(std::vector<double>::iterator it=P.begin()+1;it!=P.end(); ++it) *it=0;
    int n=(int) (log(m)/log(2));
        //take only the indices index=2^ind for ind in [0,n] where index < m
    for (int k=0;k<=n;k++) P[pow(2,k)]=(1-pow(2,-mFitness))*pow(2,-k*mFitness);

  } else {
    double dstar=mDeath/(1-mDeath);

    //                  [.....]
  }

  return P;

}


std::vector<double> FLAN_DiracClone::computeProbability1DerivativeRho(int m,const std::vector<double> P){

  std::vector<double>dP_dr(m+1);

  if(mDeath<DEATH_EPS_DIST){

    dP_dr[0]=0;
    if(m == 0) return dP_dr;

    for(std::vector<double>::iterator it=dP_dr.begin()+1;it!=dP_dr.end(); ++it) *it=0;
    int n=(int) (log(m)/log(2));
    double index;

    for(int k=0;k<=n;k++) {
      index=pow(2,k);
      dP_dr[index]=log(2)*(pow(2,-mFitness*(k+1))-k*P[index]);
    }

  } else {
    double dstar=mDeath/(1-mDeath);

    //                  [.....]
  }
  return dP_dr;

}


// double FLAN_DiracClone::computeGeneratingFunction(double z){
//
//     double eps=1.e-8;
//
//     // z=0 return 0
//     if (fabs(z)<eps) return 0;
//
//     // z=1 return 1
//     if (fabs(1-z)<eps) return 1;
//
//     // otherwize
//     double s=0;
//     if (mDeath<eps) {
// //         double a=pow((int)2,-mFitness);
// 	double a=pow(2,-mFitness);
//         int n=(int) (4-log(fabs(log(z)))/log(2));
//
//         for (int k=0;k<=n;k++) {
//             s+=pow(z,pow(2,k))*pow(a,k);
//         }
//         s*=(1-a);
//     } else {
//         double a=log(2*(1-mDeath));
//         double dstar=mDeath/(1-mDeath);
//         int n=(int)(-log(eps)/(mFitness*a));
//         double tp=exp(-mFitness*a);
//         double tpi=1;
//         double bi=z;
//         s=z;
//         for (int i=1;i<=n;i++) {
//             bi=mDeath+(1-mDeath)*bi*bi;
//             tpi*=tp;
//             s+=tpi*bi;
//         }
//         s*=(1-tp);
//         s+=dstar*tpi*tp;
//     }
//
//     return s;
//
// }
//
//
// double FLAN_DiracClone::computeGeneratingFunction1DerivativeRho(double z) {
//
//     double eps=1.e-8;
//     // z=0 return 0
//     if (fabs(z)<eps) return 0;
//     // z=1 return 1
//     if (fabs(1-z)<eps) return 0;
//     // otherwise
//
//     if (mDeath<eps) {
//         double a=pow((int)2,-mFitness);
//         int n=(int) (4-log(fabs(log(z)))/log(2));
//         double s1=0,s2=0,t=0;
//         for (int k=0;k<=n;k++) {
//             t=pow(z,pow(2,k))*pow(a,k);
//             s1+=t;
//             s2+=k*t;
//
//         }
//         return log(2)*(a*s1-(1-a)*s2);
//     } else {
//         double a=log(2*(1-mDeath));
//         double dstar=mDeath/(1-mDeath);
//         int n=(int) (-log(eps)/(mFitness*a));
//         double tp=exp(-mFitness*a);
//         double tpi=1;
//         double bi=z;
//         double s=z;
//         double ds=0;
//         for (int i=1;i<=n;i++) {
//             bi=mDeath+(1-mDeath)*bi*bi;
//             tpi*=tp;
//             s+=tpi*bi;
//             ds-=i*tpi*bi;
//         }
//         s=a*(ds*(1-tp)+s*tp);
//         // no usefull when n is big
// //         s-=dstar*a*tpi*tp*(n+1);
// 	s+=dstar*tpi*tp;
//         return s;
//     }
// }
//
//
//
