.onLoad <- function(libname, pkgname) {
	loadRcppModules()
# 	local({
# 		.flantol=.Machine$double.eps^0.5
# 		.integrands=list(
# 				CLONE_P0_WD=function(x,rho,delta) {(1-x)*x^(rho-1)/(1-delta*x)},
# 				CLONE_PK_WD=function(x,rho,delta,k) {x^rho*(1-x/k)^(k-1)/(1-x/k*delta)^(k+1)},
# # 				(1-x)^(k-1)*x^rho/(1-delta*x/k)^(k+1)
# 				CLONE_dP0_dr_WD=function(x,rho,delta) {(1-x)*x^(rho-1)/(1-delta*x)*log(x)},
# 				CLONE_dPK_dr_WD=function(x,rho,delta,k) {x^rho*(1-x/k)^(k-1)/(1-x/k*delta)^(k+1)*log(x/k)}
# 				)
# 				},globalenv())
}

.onAttach <- function(libname, pkgname) {
	local({
		.flantol=.Machine$double.eps^0.5
		.integrands=list(
				CLONE_P0_WD=function(x,rho,delta) {(1-x)*x^(rho-1)/(1-delta*x)},
				CLONE_PK_WD=function(x,rho,delta,k) {x^rho*(1-x/k)^(k-1)/(1-x/k*delta)^(k+1)},
# 				(1-x)^(k-1)*x^rho/(1-delta*x/k)^(k+1)
				CLONE_dP0_dr_WD=function(x,rho,delta) {(1-x)*x^(rho-1)/(1-delta*x)*log(x)},
				CLONE_dPK_dr_WD=function(x,rho,delta,k) {x^rho*(1-x/k)^(k-1)/(1-x/k*delta)^(k+1)*log(x/k)}
				)
				},globalenv())
}
