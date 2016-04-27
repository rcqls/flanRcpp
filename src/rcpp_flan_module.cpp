#include "rcpp_flan_module.h"

// static void finalizer_of_flan( Flan* ptr ){
//       printf("finalizer sim_vam has been called\n");
// }

RCPP_MODULE(flan_module) {
  class_<FLAN_Sim>("FlanSim")
	.constructor<List>()
	.method("rflan",&FLAN_Sim::computeSamplesMutantsFinalsNumber,"compute sample mutants")
  ;

  class_<FLAN_MutationModel>("FlanMutMod")
	.constructor<List>()
	.constructor<double,double,std::string>()
	.method("pflan",&FLAN_MutationModel::computeCumulativeFunction,"compute cumulative function")
	.method("dflan",&FLAN_MutationModel::computeProbability,"compute probability")
	.method("dflanda",&FLAN_MutationModel::computeProbability1DerivativeAlpha,"compute probability derivative wtt alpha")
	.method("dflandr",&FLAN_MutationModel::computeProbability1DerivativeRho,"compute probability derivative wtt rho")
	.method("dflangrad",&FLAN_MutationModel::computeProbability1DerivativesAlphaRho,"compute probability derivative wtt alpha and rho")
// 	.method("getfcts",&FLAN_MutationModel::getFns,"get function")
// 	.method("unbias_mutprob_estimation",&FLAN_MutationModel::unbiasPiEstimation,"unbias mutprob estimation")
  ;
  
  class_<FLAN_ExponentialClone>("FlanExpClone")
	.constructor<double,double>()
	.method("dclone",&FLAN_ExponentialClone::computeProbability,"compute probability")
	.method("dclonedr",&FLAN_ExponentialClone::computeProbability1DerivativeRho,"compute probability")
  ;
  class_<FLAN_SimClone>("FlanSimClone")
	.constructor<double,double,List>()
	.method("rclone",&FLAN_SimClone::computeSample,"compute clone sample")
  ;
  
  
  
  
  class_<MATH_Polynom>("MathPol")
  .constructor< std::vector<double> >()
  .method("square.pol",&MATH_Polynom::square,"compute prod of a polynom with itself")
  ;

  //function( "newMaintenancePolicy", &newMaintenancePolicy );

}
