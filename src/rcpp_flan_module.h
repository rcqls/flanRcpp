#ifndef RCPP_FLAN_MODULE_H
#define RCPP_FLAN_MODULE_H

#include "FLAN_Sim.h"
#include "FLAN_MutationModel.h"
#include "MATH_Function.h"
// #include "rcpp_mle_vam.h"
// #include "rcpp_vam_model.h"
//
RCPP_EXPOSED_AS(FLAN_Sim);
RCPP_EXPOSED_WRAP(FLAN_Sim);
//
RCPP_EXPOSED_AS(FLAN_MutationModel);
RCPP_EXPOSED_WRAP(FLAN_MutationModel);
//
// RCPP_EXPOSED_AS(VamModel);
// RCPP_EXPOSED_WRAP(VamModel);

#endif //RCPP_FLAN_MODULE_H
