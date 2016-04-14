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
	.method("pflan",&FLAN_MutationModel::computeCumulativeFunction,"compute cumulative function")
	.method("dflan",&FLAN_MutationModel::computeProbability,"compute probability")
// 	.method("unbias_mutprob_estimation",&FLAN_MutationModel::unbiasPiEstimation,"unbias mutprob estimation")
  ;

  //function( "newMaintenancePolicy", &newMaintenancePolicy );

}
