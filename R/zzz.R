.onLoad <- function(libname, pkgname) {
	loadRcppModules()
	local({
		.flantol=.Machine$double.eps^0.5
		.integrands<-list(
				CF_GY_WD=function(x,rho,delta){(1-x)*x^(rho-1)/(1-delta/(1-delta)*x)}
				)
				},globalenv())
}
