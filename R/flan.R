  ## --------------------------------------------------------------------------------------------------------------
  ############################################## flantest object class ##############################################
  ## --------------------------------------------------------------------------------------------------------------

# Constructor of the class flantest
flantest <- function(Tstat,parameter,estimate,p.value,conf.int,estimates,null.value,alternative,model,method,data.name,sample.name,nsamples){
  obj <- list(Tstat=Tstat,estimate=estimate,parameter=parameter,conf.int=conf.int,p.value=p.value,
  null.value=null.value,alternative=alternative,data.name=data.name,model=model,method=method,nsamples=nsamples)
  class(obj) <- "flantest"
  return(obj)
}

# Print function for flantest objects, inspired from print function for htest objects
print.flantest <- function(object){
  digits <- getOption("digits")
  prefix <- "\t"
  cat("\n")
  if(object$nsamples == 2){
    cat(strwrap(paste("Paired ",object$method,"-test"," (",object$model," model)",sep=""), prefix=prefix), sep="\n")
  } else {
    cat(strwrap(paste("One sample ",object$method,"-test"," (",object$model," model)",sep=""), prefix=prefix), sep="\n")
  }
  cat("\n")
  cat("------------------------------------- Data -----------------------------------\n")

  if(object$nsamples == 2){
    if(length(object$data.name) == 1){
      cat("data:  ",object$data.name[[1]][1]," and ",object$data.name[[1]][2],"\n", sep="")
    }
    if(length(object$data.name) == 2){
      if(length(object$data.name[[2]]) == 2){
	cat("data:  ",object$data.name[[1]][1]," with ",object$data.name[[2]][1],
	" and ",object$data.name[[1]][2]," with ",object$data.name[[2]][2],"\n", sep="")
      } else {
	cat("data:  ",object$data.name[[1]][1]," with ",object$data.name[[2]][1],
	" and ",object$data.name[[1]][2],"\n", sep="")
      }
    }
  }
  if(object$nsamples == 1){
    if(length(object$data.name) == 1){
      cat("data:  ",object$data.name[[1]][1],"\n", sep="")
    }
    if(length(object$data.name) == 2){
      cat("data:  ",object$data.name[[1]][1]," with ",object$data.name[[2]][1],"\n", sep="")
    }
  }

  out <- character()

  if(!is.null(object$null.value)){
    if(length(object$null.value) == 1){
      if(!is.null(object$parameter)){
	if(object$nsamples == 2){
	  Nparam <- dim(object$param)[2]
	  cat("Sample 1 parameters : ")
	  for(i in 1:Nparam){
	    out <- c(out,paste(colnames(object$parameter)[i], "=", format(signif(object$parameter[1,i],
	      max(1, digits - 2)))))
	  }
	  cat(strwrap(paste(out, collapse=", ")), sep="\n")
	  out <- character()
	  cat("Sample 2 parameters : ")
	  for(i in 1:Nparam){
	    out <- c(out,paste(colnames(object$parameter)[i], "=", format(signif(object$parameter[2,i],
	      max(1, digits - 2)))))
	  }
	} else {
	  Nparam <- length(object$parameter)
	  for(i in 1:Nparam){
	    if(names(object$parameter[i])=="mfn"){
	      out <- c(out,paste(names(object$parameter[i]), "=", format(object$parameter[i],
	      scientific=TRUE,digit=5)))
	    } else {
	      out <- c(out,paste(names(object$parameter[i]), "=", format(signif(object$parameter[i],
	      max(1, digits - 2)))))
	    }
	  }
	}
      }
      cat(strwrap(paste(out, collapse=", ")), sep="\n")
      cat("---------------------------------- Statistics --------------------------------\n")
      if(!is.null(object$Tstat)){
	cat("Tstat =", format(signif(object$Tstat,
	  max(1, digits - 2))),"\n")
      }
      if (!is.null(object$p.value)) {
	fp <- format.pval(object$p.value, digits= max(1, digits -
	    3))
	cat("p-value for",names(object$null.value[1]), if (substr(fp, 1, 1) ==
	  "<") fp else paste("=", fp),"\n")
      }

      if (!is.null(object$alternative)) {
	cat("Alternative hypothesis: ")
	alt.char <- switch(object$alternative, two.sided="not equal to",
		less="less than", greater="greater than")
	cat("true ", names(object$null.value), " is ", alt.char,
	  " ", object$null.value, "\n", sep="")
      }

      if (!is.null(object$conf.int)) {
	cat(format(100 * attr(object$conf.int, "conf.level")), " percent confidence interval:\n",
	paste(format(object$conf.int),collapse="  "),
	"\n")
      }
      if (!is.null(object$estimate)) {
	if(object$nsamples == 1){
	  cat("Sample estimates : \n")
	  print(object$estimate)
	}
	if(object$nsamples == 2){
	  cat("Sample 1 estimates : \n")
	  print(object$estimate[1,])
	  cat("Sample 2 estimates : \n")
	  print(object$estimate[2,])
	}
      }
    } else if(length(object$null.value)==2){

      if(!is.null(object$parameter)){
	if(object$nsamples == 2){
	  cat("Sample 1 parameters : ")
	  Nparam <- dim(object$param)[2]
	  for(i in 1:Nparam){
	    out <- c(out,paste(colnames(object$parameter)[i], "=", format(signif(object$parameter[1,i],
	      max(1, digits - 2)))))
	  }
	  cat(strwrap(paste(out, collapse=", ")), sep="\n")
	  out <- character()
	  cat("Sample 2 parameters : ")
	  for(i in 1:Nparam){
	    if(colnames(object$parameter)[i] == "mfn"){
	      out <- c(out,paste(colnames(object$parameter)[i], "=", format(object$parameter[2,i],
	      scientific=TRUE,digit=5)))
	    } else {
	      out <- c(out,paste(colnames(object$parameter)[i], "=", format(signif(object$parameter[2,i],
	      max(1, digits - 2)))))
	    }
	  }
	}
	if(object$nsamples == 1){
	  cat("Sample parameters : ")
	  Nparam <- length(object$parameter)

	  for(i in 1:Nparam){
	    out <- c(out,paste(names(object$parameter[i]), "=", format(signif(object$parameter[i],
	      max(1, digits - 2)))))
	  }
	}
      }
      cat(strwrap(paste(out, collapse=", ")), sep="\n")
      cat("---------------------------------- Statistics --------------------------------\n")
      if(!is.null(object$Tstat)){
	cat("Tstat = (",format(signif(object$Tstat[1],
	  max(1, digits - 2)))," , ",format(signif(object$Tstat[2],
	  max(1, digits - 2))),")","\n",sep="")
      }
      if (!is.null(object$p.value)) {
	fp1 <- format.pval(object$p.value[1], digits=max(1, digits -
	    3))
	fp2 <- format.pval(object$p.value[2], digits=max(1, digits -
	    3))
	cat("p-value for",names(object$null.value[1]), if (substr(fp1, 1, 1) ==
	  "<") fp1 else paste("=", fp1),"\np-value for",names(object$null.value[2]),if (substr(fp2, 1, 1) ==
	  "<") fp2 else paste("=", fp2),"\n")
      }
      if (!is.null(object$alternative)) {
	cat("Alternative hypotheses: ")
	alt.char1 <- switch(object$alternative[1], two.sided="not equal to",
		less="less than", greater="greater than")
	alt.char2 <- switch(object$alternative[2], two.sided="not equal to",
		less="less than", greater="greater than")
	cat("true ", names(object$null.value[1]), " is ", alt.char1,
	  " ", object$null.value[1], "\n", sep="")
	cat("                        ")
	cat("true ", names(object$null.value[2]), " is ", alt.char2,
	  " ", object$null.value[2], "\n", sep="")
      }
      if (!is.null(object$conf.int)) {
	cat(format(100 * attr(object$conf.int, "conf.level")), " percent confidence interval for ",names(object$null.value[1]),": \n",
	  paste(format(object$conf.int[,1]),collapse="   "),
	  "\n",
	  format(100 * attr(object$conf.int, "conf.level")), " percent confidence interval for ",names(object$null.value[2]),": \n",
	  paste(format(object$conf.int[,2]),collapse="   "),
	  "\n",sep="")
      }
      if (!is.null(object$estimate)) {
	if(object$nsamples == 1){
	  cat("Sample estimates : \n")
	  print(object$estimate[1,])
	}
	if(object$nsamples == 2){
	  cat("Sample 1 estimates : \n")
	  print(object$estimate[1,])
	  cat("Sample 2 estimates : \n")
	  print(object$estimate[2,])
	}
      }
    }
    cat("\n")
  }
}


MutationsNumberP0Estimation <- function(mc,death){

# Return the estimate of alpha for a sample mc
# by p0-method under cell deaths with probability delta
# # when delta is known.
  if(death > 0){
    dstar <- death/(1-death)			# extinction probability
    epgf <- function(z) mean(z^mc)		# empirical PGF

    a <- -log(epgf(dstar))/(1-dstar)	#  estimate alpha

    sda <- (1-dstar)^(-2)*(epgf(dstar^2)/(epgf(dstar)^2)-1)	# variance of alpha's estimate

    sda <- sqrt(sda/length(mc))

  } else {
    p0 <- mean(mc==0)
    a <- -log(p0)
    sda <- sqrt((1/p0-1)/length(mc))
  }

  list(mutations=a,sd.mutations=sda)

}


    ## --------------------------------------------------------------------------------------------------------------------------
    ############################################## Estimations and tests functions ##############################################
    ## --------------------------------------------------------------------------------------------------------------------------

  # Estimation function of the mean number of mutations. Compute also the probability of mutation if with.prob is TRUE and the fitness parameter if fitness is empty.
    # Arguments:
    # mc:  a (non-empty) numeric vector of mutants counts.
    # fn:  an optional (non-empty) numeric vector with same length as mc of final numbers of cells.
    # mfn:  mean final number of cells which. Ignored if fn is non-missing.
    # cvfn:  coefficient of variation of final number of cells. Ignored if fn is non-missing.
    # fitness:  fitness parameter: ratio of growth rates of normal and mutant cells. If fitness is NULL (default) then the fitness will be estimated. Otherwise, the given value will be used to estimate the mean mutation number mutations
    # death:  death probability. Must be smaller than 0.5.
    # method:  estimation method as a character string: one of ML (default), P0, or GF.
    # winsor:  winsorization parameter: positive integer. Only used when method is ML or when method is P0 and fitness is NULL.
    # model:  statistical lifetime model as a character string: one of LD (default) for Luria-Delbrück model with exponential lifetimes, or H for Haldane model with constant lifetimes.
    # Assumptions for the estimation
    #   - fitness (if given)
    #   - death probability
    #   - lifetimes distribution model : exponentialy distributed lifetime ("LD") or constant lifetime "H"
    #   - sample of mutants counts
    #   - sample of finals counts (or a mean and coefficient of variation for finals counts)
    # Three possible methods :
    #   - p0-method ("P0")   (winsorize the mc sample if rho is estimated with the parameter winsor)
    #   - Generating function method ("GF")
    #   - Maximum of Likelihood method ("ML")  (winsorize the mc sample with the parameter winsor)
    #
    # Returns :
    #   - Estimate of alpha (or pi) and rho if non given
    #   - Standard deviation of alpha, (or pi) and rho if non given

mutestim <- function(mc,fn=NULL,mfn=NULL,cvfn=NULL,                 # user's data
                  fitness=NULL,death=0,                # user's parameters
                  method=c("ML","GF","P0"),winsor=512, # estimation method
                  model=c("LD","H"))                   # clone growth model
                  {

    if(missing(method)){method <- "ML"}
    if(missing(model)){model <- "LD"}
    method <- match.arg(method)
    model <- match.arg(model)
    with.prob.fn <- FALSE             # if fn is given, deduce estimate of mutprob from estimate of mutations. Moreover, if method is "ML", compute directly estimate of mutprob
    with.prob.stat <- FALSE           # if mfn or cvfn is given, deduce estimate of mutprob from estimate of mutations

    if(is.null(mc)){
      stop("'mc' is empty...")
    }
    if(sum(floor(mc)==mc) != length(mc)) stop("'mc' must be a vector of integer.")
    if(!is.null(fn)){
      if(length(fn) != length(mc)){
	stop("'fn must have the same length as 'mc'.")
      }
      if(sd(fn) == 0) {
	mfn <- mean(fn)
	cvfn <- 0
	fn <- NULL
	with.prob.stat <- TRUE
      } else with.prob.fn <- TRUE
    }
    if(!is.null(mfn)){
      if(mfn < 0 | length(mfn) > 1){
	stop("'mfn' must be a single positive number.")
      }
      if(is.null(cvfn)) cvfn <- 0
      with.prob.stat <- TRUE
    }
    if(!is.null(cvfn)){
      if(cvfn < 0 | length(cvfn) > 1){
	stop("'cvfn' must be a single positive number.")
      }
      if(is.null(mfn)) mfn <- 1e9
      with.prob.stat <- TRUE
    }
    if(!is.null(fitness)){
      if(fitness < 0 | length(fitness) > 1){
	stop("'fitness' must be empty or a single positive number.")
      }
    }
    if(death < 0 | death >= 0.5 | length(death) > 1){
      stop("'death must be a single positive and < 0.5 number.")
    }
    if(winsor < 0 | trunc(winsor) != winsor | length(winsor) > 1){
      stop("'winsor' must be a single positive integer.")
    }
    if(method=="P0" & sum(mc==0) == 0 & death == 0){
      stop(paste("P0 method can not be used if 'mc' does not contain any null counts and if 'death' is null.",sep=" "))
    }

    if(with.prob.fn) {
      with.prob.stat <- FALSE      # If fn is non empty: ignore mfn and cvfn even if they are non empty
    }


    if (is.null(fn)) fn <- c()


    optimizer="bfgs"    # Optimization method for maximum of likelihood

    if(is.null(fitness)){   # If fitness is empty : compiute estimates of mean number of mutations (mutation probaility), fitness

      # Maximum of Likelihood estimators
      if(method == "ML"){
	    output <- .Call("MutationFitnessMLOptimization",
		      list(mc=mc,
			  fn=as.double(fn),
			  optimizer=as.character(optimizer),
			  model=as.character(model),
			  death=as.double(death),
			  winsor=as.integer(winsor),
			  with.prob.fn=as.logical(with.prob.fn),
			  with.prob.stat=as.logical(with.prob.stat),
			  mfn=as.double(mfn),
			  cvfn=as.double(cvfn)
			)
		    )
	if(!output$succeeds) warning("Initialization of 'fitness' with 'GF'-method may have failed.")
	output$succeeds <- NULL
      }
      # P0 method
      if(method == "P0"){

	a.est <- MutationsNumberP0Estimation(mc,death)	 # P0 estimator of mutations

	mut <- a.est$mutations
	sdmut <- a.est$sd.mutations

	if(with.prob.fn | with.prob.stat){     # P0 estimator of mutprob
	  if(with.prob.fn){
	    mfn <- mean(fn)
	    cvfn <- sd(fn)/mfn
	  }
	  mut <- mut/mfn
	  sdmut <- sdmut/mfn

	  if(cvfn != 0){
	    umds <- 1-death/(1-death)			# extinction probability
	    temp <- a.est$mutations*(cvfn^2)*umds         # relative bias on mutprob
	    mut <- mut*(1+temp/2)		# unbiased estimator of mutrate
	    sdmut <- sdmut*(1+temp)		# standard deviation of unbiased estimator
	  }
	}
	# Maximum of likelihood estimator of fitness
	r.est <- .Call("FitnessP0Optimization",
		      list(mc=mc,
			   fn=as.double(fn),
			   optimizer=as.character(optimizer),
			   model=as.character(model),
			   mut=as.double(mut),
			   death=as.double(death),
			   winsor=as.integer(winsor)
			)
		    )
	if(!r.est$succeeds) warning("Initialization of 'fitness' with 'GF'-method may have failed.")

	if(with.prob.fn | with.prob.stat) output <- list(mutprob=mut,sd.mutprob=sdmut,fitness=r.est$fitness,sd.fitness=r.est$sd.fitness)
	else output <- list(mutations=mut,sd.mutations=sdmut,fitness=r.est$fitness,sd.fitness=r.est$sd.fitness)
      }

      # GF method
      if(method == "GF"){
	output <- .Call("MutationFitnessGFEstimation",
		      list(mc=mc,
			    fn=as.double(fn),
			    model=as.character(model),
			    death=as.double(death),
			    with.prob.fn=as.logical(with.prob.fn),
			    with.prob.stat=as.logical(with.prob.stat),
			    mfn=as.double(mfn),
			    cvfn=as.double(cvfn)
			)
		    )

	# Sometime the resolution of fixed point equation does not converge
	if(!output$succeeds) warning(paste("Impossible to estimate 'fitness' of",deparse(substitute(mc)),"with 'GF'-method : 'fitness' is set to default value 1 and only",if(with.prob.fn | with.prob.stat){"mutation probability and mutation number are"}else{"mutation number is"},"estimated.",sep=" "))
	output$succeeds <- NULL
      }

    } else {    # Else : compute estimate(s) of mean number of mutations, or mutation probability

      # Maximum of Likelihood estimator of mutations or mutprob
      if(method == "ML"){
	    output <- .Call("MutationMLOptimization",
		      list(mc=mc,
			    fn=as.double(fn),
			    optimizer=as.character(optimizer),
			    model=as.character(model),
			    fitness=as.double(fitness),
			    death=as.double(death),
			    winsor=as.integer(winsor),
			    with.prob.fn=as.logical(with.prob.fn),
			    with.prob.stat=as.logical(with.prob.stat),
			    mfn=as.double(mfn),
			    cvfn=as.double(cvfn)
			)
		    )


      }
      # P0 method
      if(method == "P0"){

	a.est <- MutationsNumberP0Estimation(mc,death)     # P0 estimator of mutations

	if(with.prob.fn | with.prob.stat){     # P0 estimator of mutprob
	  if(with.prob.fn){
	    mfn <- mean(fn)
	    cvfn <- sd(fn)/mfn
	  }
	  mut <- a.est$mutations/mfn
	  sdmut <- a.est$sd.mutations/mfn

	  if(cvfn != 0){
	    umds <- 1-death/(1-death)			# extinction probability
	    temp <- a.est$mutations*(cvfn^2)*umds         # relative bias on mutprob
	    mut <- mut*(1+temp/2)		# unbiased estimator of mutrate
	    sdmut <- sdmut*(1+temp)		# standard deviation of unbiased estimator
	  }
	} else {
	  mut <- a.est$mutations
	  sdmut <- a.est$sd.mutations
	}

	if(with.prob.fn | with.prob.stat) output <- list(mutprob=mut,sd.mutprob=sdmut)
	else output <- list(mutations=mut,sd.mutations=sdmut)
      }

      # GF method
      if(method == "GF"){
	output <- .Call("MutationGFEstimation",
		      list(mc=mc,
			    fn=as.double(fn),
			    model=as.character(model),
			    fitness=as.double(fitness),
			    death=as.double(death),
			    with.prob.fn=as.logical(with.prob.fn),
			    with.prob.stat=as.logical(with.prob.stat),
			    mfn=as.double(mfn),
			    cvfn=as.double(cvfn)
			)
		    )
      }

    }

    output
}


    # One-sample or two-sample tests for mean number of mutations, fitness and mutation probability using mutestim estimates
    # The possible assumptions are the same as for the mutestim function.
    # Arguments :
    # Returns a "flantest" object, which contains:
    # Tstat:  the value of the computed statistic.
    # parameter:  the values of the parameter of the model: fitness(if not tested), death, mfn (if needed) and cvfn
    # p.value:  the p-value(s) of the test.
    # conf.int:  a confidence interval for the parameter(s) of interest appropriate to the specified alternative hypothesis.
    # estimates:  the estimate(s) of interest.
    # null.value:  the specified hypothesized value(s).
    # alternative:  a character string describing the alternative hypothesis.
    # model:  the statistical lifetime model.
    # method:  method used to compute the estimation(s) of the parameter(s) of interest.
    # data.name:  a character string giving the name of the complete data.


flan.test <- function(mc,fn=NULL,mfn=NULL,cvfn=NULL,                      # user's data
               fitness=NULL,death=0,                         # user's parameters
               mutations0=1,mutprob0=NULL,fitness0=1,       # null hypotheses
               conf.level=0.95,                              # confidence level
               alternative=c("two.sided","less","greater"),  # alternative
               method=c("ML","GF","P0"),winsor=512,          # estimation method
               model=c("LD","H"))                            # clone growth model
               {

  with.prob <- FALSE               # Boolean: if TRUE (if fn, mfn, or cvfn are given), mutprob is tested instead of mutations

  if(is.null(mc)){
    stop("'mc' is empty...")
  } else {
    if(is.list(mc)){
      if(length(mc) == 1){
	nsamples <- 1
	dname <- list(deparse(substitute(mc)))
	if(is.null(mc[[1]])) stop("'mc[[1]]' is empty...")
	if(is.list(fn)){
	  if(length(fn) != nsamples) stop("'mc' and 'fn' must have the same length.")
	  if(!is.null(fn[[1]])){
	    dname <- c(dname,deparse(substitute(fn)))
	    with.prob <- TRUE
	    mfn <- mean(fn)
	    cvfn <- sd(fn)/mfn
	  }
	} else {
	  if(!is.null(fn)) stop("'fn' must have the same type as 'mc'.")
	  fn <- list(fn)
	}
      }
      if(length(mc) == 2){
	nsamples <- 2
	dname <- list(c(paste(deparse(substitute(mc)),"1",sep=""),paste(deparse(substitute(mc)),"2",sep="")))
	if(is.null(mc[[1]])) stop("'mc[[1]]' is empty...")
	if(is.null(mc[[2]])) stop("'mc[[2]]' is empty...")

	if(is.list(fn)){
	  if(length(fn) != nsamples) stop("'fn' must have the same length as 'mc'.")

	  if(!is.null(fn[[1]])) {
	    if(length(fn[[1]]) != length(mc[[1]])) stop("'fn[[1]]' must have the same length as 'mc[[1]]'.")

	    if(is.null(fn[[2]])){
	      warning("'fn[[2]]' is empty : 'fn[[1]]' is ignored.")
	      fn <- list(NULL,NULL)
	    } else {
	      if(length(fn[[2]]) != length(mc[[2]])) stop("'fn[[2]]' must have the same length as 'mc[[2]]'.")

	      dname <- c(dname,list(c(paste(deparse(substitute(fn)),"1",sep=""),paste(deparse(substitute(fn)),"2",sep=""))))
	      with.prob <- TRUE
	      mfn <- unlist(lapply(fn,mean))
	      cvfn <- unlist(lapply(fn,sd))/mfn
	    }
	  } else {
	    if(!is.null(fn[[2]])) {
	      warning("'fn[[1]]' is empty : 'fn[[2]]' is ignored.")
	      fn <- list(NULL,NULL)
	    }
	  }
	} else {
	  if(!is.null(fn)) stop("'fn' must have the same type as 'mc'.")
	  fn <- list(fn,fn)
	}
      }
    } else {
      nsamples <- 1
      dname <- list(deparse(substitute(mc)))
      if(!is.null(fn)){
	if(is.list(fn)) stop("'fn' must have the same type as 'mc'.")
	if(length(fn) != length(mc)) stop("'fn' must have the same length as 'mc'.")
	dname <- c(dname,deparse(substitute(fn)))
	with.prob <- TRUE
	mfn <- mean(fn)
	cvfn <- sd(fn)/mfn
      }
      mc <- list(mc)
      fn <- list(fn)
    }
  }

  if(!is.null(mfn)){
    if(sum(mfn < 0) != 0 | length(mfn) > 2) stop("if given, 'mfn' must be a vector with size <= 2 of positive numbers.")

    if(is.null(cvfn)) cvfn <- rep(0,length(mfn))
    with.prob <- TRUE
  }
  if(!is.null(cvfn)){
    if(sum(cvfn < 0) != 0 | length(cvfn) > 2){
      stop("if given, 'cvfn' must be a vector with size <= 2 of positive numbers.")
    }
    if(is.null(mfn)) mfn <- rep(1e9,length(cvfn))
    with.prob <- TRUE
  }
  if(!is.null(fitness)){
    if(sum(fitness < 0) != 0 | length(fitness) > 2){
      stop("if given, 'fitness' must be a vector with size <= 2 of positive numbers.")
    }
    fitness0 <- NULL
  }
  if(length(mutations0) > 1 | mutations0 < 0){
    stop("'mutations0' must be a single positive number.")
  }
  if(!is.null(fitness0)){
    if(length(fitness0) > 1 | fitness0 < 0){
      stop("'fitness0' must be a single positive number.")
    }
  }
  if(!is.null(mutprob0)){
    if(length(mutprob0) > 1 | mutprob0 < 0 | mutprob0 >= 1){
      stop("if given, 'mutprob0' must be a positive and <= 1 number.")
    }
    with.prob <- TRUE
    if(is.null(mfn)) mfn <- 1e9
    if(is.null(cvfn)) cvfn <- 0
  } else {
    if(with.prob){
      if(is.null(mfn)) mfn <- 1e9
      if(is.null(cvfn)) cvfn <- 0
      if(nsamples == 1) mutprob0 <- mutations0/mfn
    }
  }


  # Default values if two-samples test
  if(nsamples == 2){
    if(missing(mutations0)) mutations0 <- 0
    if(is.null(fitness)){
      if(missing(fitness0)) fitness0 <- 0
    }
    if(with.prob) {
      if(missing(mutprob0)) mutprob0 <- 0
    }
  }

  if((length(conf.level) > 1) | conf.level > 1 | conf.level < 0){
    stop("'conf.level' must be a single positive and <= 1 numbers.")
  }
  if(sum(death < 0 | death >= 0.5) != 0 | length(death) > 2){
    stop("'death' must be a vector with size <= 2 of positive and < 0.5 numbers.")
  }
  if(winsor < 0 | trunc(winsor) != winsor | length(winsor) > 1){
    stop("'winsor' must be a single positive integer.")
  }

  H0 <- c(if(with.prob) mutprob0 else mutations0,fitness0)     # Vector of null hypothesises

  np <- length(H0)                                             # Number of tested values

  if(missing(alternative)) {alternative <- rep("two.sided",np)}
  if(missing(method)) {method <- "ML"}
  if(missing(model)) {model <- "LD"}

  alternative <- match.arg(alternative,several.ok=TRUE)
  method <- match.arg(method)
  model <- match.arg(model)


  if(nsamples == 1){

    parameter <- NULL
    names <- character()

    if(np == 1){
      names(H0) <- if(with.prob) "mutation probability" else "mutation number"
      parameter <- c(fitness,death)
      names <- c("fitness","death")
      if(with.prob){
	parameter <- c(parameter,mfn,cvfn)
	names <- c(names,"mfn","cvfn")
      }
    }
    if(np == 2) {
      names(H0) <- c(if(with.prob) "mutation probability" else "mutation number", "fitness")
      parameter <- death
      names <- "death"
      if(with.prob){
	parameter <- c(parameter,mfn,cvfn)
	names <- c(names,"mfn","cvfn")
      }
    }
    names(parameter) <- names
  }
  if(nsamples == 2){
    if(np == 1 & length(fitness) == 1){fitness <- c(fitness,fitness)}
    if(length(death) == 1){death <- c(death,death)}
    if(length(mfn) == 1){mfn <- c(mfn,mfn)}
    if(length(cvfn) == 1){cvfn <- c(cvfn,cvfn)}

    parameter <- NULL
    names <- character()


    if(np == 1){
      names(H0) <- if(with.prob) "mutprob difference" else "mutations difference"
      parameter <- cbind(fitness,death)
      names <- c("fitness","death")
      if(with.prob) {
	parameter <- cbind(parameter,mfn,cvfn)
	names <- c(names,"mfn","cvfn")
      }
    }
    if(np == 2){
      names(H0) <- c(if(with.prob) "mutprob difference" else "mutations difference","fitness difference")
      parameter <- cbind(death)
      names <- "death"
      if(with.prob){
	parameter <- cbind(parameter,mfn,cvfn)
	names <- c(names,"mfn","cvfn")
      }
    }
    colnames(parameter) <- names
  }

  if(is.null(mfn)) mfn <- list(mfn)
  if(is.null(cvfn)) cvfn <- list(cvfn)
  if(is.null(fitness)) fitness <- list(fitness)

   # Estimates mean number of mutations alpha, given fitness parameter rho
  ests <- mapply(function(x,y,m,c,f,d){
    mutestim(mc=x,fn=y,mfn=m,cvfn=c,method=method,model=model,fitness=f,death=d,winsor=winsor)
  },mc,fn,mfn,cvfn,fitness,death)


  if(np == 1){   # If only mutations (or mutprob) is tested
    if(length(alternative) > 1) alternative <- alternative[1]

    Tstat <- unlist(ests[1,])    # Extract the estimate to compute statistic of the test
    sds <- unlist(ests[2,])      # Standard deviation of the estimate
    names(Tstat) <- NULL
    ests <- Tstat                # Keep the estimate

    if(nsamples == 2){    # If Two-sample test
      ests <- cbind(ests)
      Tstat <- -diff(Tstat)      # Statistic of the test is build with difference of estimates
      sds <- sqrt(sum(sds^2))    # Standard deviation of the difference between estimates
      colnames(ests) <- if(with.prob) "mutprob" else "mutations"
    }
    names(ests) <- if(with.prob) "mutprob" else "mutations"


  }
  if(np == 2){    # If fitness is also tested
    if(length(alternative) == 1) alternative <- rep(alternative,2)
    if(length(alternative) > 2) alternative <- alternative[c(1,2)]

    Tstat <- rbind(unlist(ests[1,]),unlist(ests[3,]))   # Extract estimates to compute statistic of the test
    sds <- rbind(unlist(ests[2,]),unlist(ests[4,]))     # Standard deviation of the estimates
    colnames(Tstat) <- NULL
    ests <- cbind(Tstat[1,],Tstat[2,])                  # Keep the estimates

    if(nsamples == 2){               # If Two-sample test
      Tstat <-- apply(Tstat,1,diff)                     # Statistic of the test is build with difference of estimates
      sds <- apply(sds,1,function(s){sqrt(sum(s^2))})   # Standard deviation of the difference between estimates
    }
    colnames(ests) <- c(if(with.prob) "mutprob" else "mutations","fitness")
  }

  cint <- mapply(function(e,s,alt){
    if(alt == "less"){                                        # Confidence interval(s) of the test
      c(-Inf,e+s*qnorm(conf.level))                           # with respect to the confidence level
    } else if(alt == "greater"){                              # and the alternative
      c(e-s*qnorm(conf.level),Inf)
    } else if(alt == "two.sided"){
      e+s*c(-1,1)*qnorm((1+conf.level)/2)
    }
  },Tstat,sds,alternative)



  Tstat <- (Tstat-H0)/sds                                      # Statistic(s) of the test

  pval <- mapply(function(alt,tstat){
      if(alt == "less"){                                        # p-value(s) of the test
	pnorm(tstat)                                            # with respect to the alternative
      } else if(alt == "greater"){
	pnorm(tstat,lower.tail=FALSE)
      } else if(alt == "two.sided"){
	2*pnorm(-abs(tstat))
      }
    },alternative,Tstat)

    if(nsamples == 2){
      names(Tstat) <- sapply(colnames(ests),function(noun){paste(noun,"difference")})
    } else {
      names(Tstat) <- names(ests)
    }

    colnames(cint) <- names(Tstat)
    names(pval) <- names(Tstat)
    names(alternative) <- names(Tstat)

  rownames(cint) <- c("bInf","bSup")


  attr(cint,"conf.level") <- conf.level

  # Build object with 'flantest' class

  flantest(Tstat=Tstat,estimate=ests,parameter=parameter,conf.int=cint,p.value=pval,
  null.value=H0,alternative=alternative,data.name=dname,model=model,method=method,nsamples=nsamples)

}

## --------------------------------------------------------------------------------------------------------------
    ############################ Density, distribution function, quantile function and random #######################
    ################################## generation for the mutants counts ############################################
    ## --------------------------------------------------------------------------------------------------------------

    # There are two cases :
    #  - If alpha is given, then returns a sample mcs of mutants counts with parameters :
    #   * mutations : mean number of mutation of mutation
    #   * fitness : fitness parameter
    #   * death : probability of death
    #   * dist : lifetimes distribution
    #
    #  - If alpha is not given, then compute first a sample fn of final counts with a Log-Normal distribution with mean mfn and
    #    coefficient of variation cvfn. Then compute a sample mcs of mutants counts, where the ith count has the following parameters :
    #   * mutprob*fni : mean number of mutation
    #                      (where fni is the corresponding of final count, and pi the probability of mutation)
    #   * fitness : fitness parameter
    #   * death : probability of death
    #   *dist : lifetimes distribution
    #
    #   Returns the two samples mcs and fn
    #
rflan <- function(n,mutations=1,mutprob=NULL,fitness=1,death=0,
        dist=list(name="lnorm",meanlog=-0.3795851,sdlog=0.3016223),
        mfn=1e9,cvfn=0) {

    if(n <= 0) stop("'n' must be a positive integer.")

    if(length(n) > 1) n <- length(n)

    if(mutations < 0 | length(mutations) > 1) stop("'mutations' must be a single positive number")

    if(fitness < 0 | length(fitness) > 1) stop("'fitness' must be a single positive number")

    if(mfn < 0 | length(mfn) > 1) stop("'mfn' must be a single positive number")

    if(cvfn < 0 | length(cvfn) > 1) stop("'cvfn' must be a single positive number")

    if(!is.list(dist)) stop("'dist' must be a list of a character chain followed by its arguments")


    names(dist)[1] <- "name"
    if(dist$name == "exp"){
      names(dist)[2] <- "rate"
    } else if(dist$name == "dirac"){
      names(dist)[2] <- "location"
    } else if(dist$name == "lnorm"){
      names(dist)[2] <- "meanlog"
      names(dist)[3] <- "sdlog"
    } else if(dist$name == "gamma"){
      names(dist)[2] <- "shape"
      names(dist)[3] <- "scale"
    } else stop("'dist[[1]]' must be a character chain among 'exp', 'dirac', 'lnorm', 'gamma'")


    if(death < 0 | death >= 0.5 | length(death) > 1) stop("'death' must be a single positive and < 0.5 number ")


    if(!is.null(mutprob)){
      if(mutprob < 0 | mutprob > 1 | length(mutprob) > 1) stop("'mutprob' must be a single positive and <= 1 number")
      mutations <- mutprob*(if(cvfn == 0) mfn else 1)              # Sample with mutprob instead of mutations if cvfn > 0
    } else {
      if(cvfn > 0) mutations <- mutations/mfn    # Default value of mutprob if missing and cvfn > 0
    }



      # output <- .Call("rflan",list(
			#   n=as.integer(n),
			#   mutations=as.double(mutations),
			#   fitness=as.double(fitness),
			#   death=as.double(death),
			#   distribution=dist,
			#   mfn=as.double(mfn),
			#   cvfn=as.double(cvfn)
			#   )
		  # )

    flan.sim <- new(FlanSim,list(
      mutations=mutations,
      fitness=fitness,
      death=death,
      dist=dist,
      mfn=mfn,
      cvfn=cvfn
      ))
      
    output <-flan.sim$rflan(n)


    mc <- output$mc
    indNA <- which(mc < 0)
    if(length(indNA) > 0){
      warning("Production of NaNs in geometric sampling.")
      output$mc[indNA] <- NaN
    }

    output
}


  # Quantile function of the mutants count with parameter
  #   mutations : mean number of mutations
  #   fitness : fitness parameter
  #   death : death probability
  #   model : lifetimes distribution model : exponentialy distributed lifetimes ("LD") or constant lifetimes ("H")

qflan <- function(p,mutations=1,fitness=1,death=0,model=c("LD","H"),lower.tail=TRUE){

  if((sum(p < 0)+sum(p > 1)) != 0){
    stop("'p' must be a vector of positive and <= 1 numbers")
  }
  if(mutations < 0 | length(mutations) > 1){
    stop("'mutations' must be a single positive number.")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }
  if(missing(model)){model="LD"}
  model <- match.arg(model)
#   if(!lower.tail) p <- 1-p


  m <- 100
  P <- pflan(0:m,mutations,fitness,death,model,lower.tail)
  sapply(p,function(pp){
    if(pp == 1) return(Inf)
    if (lower.tail & pp <= P[1]) return(0)
    if (!lower.tail & pp >= P[1]) return(0)
    else {
      k <- if(lower.tail) max(which(P<pp)) else max(which(P>pp))          # quantile if in the actual table
      while (k>=m){                   # if p not yet in the table
	m2 <- 2*m                       # double the table
	P2 <- pflan(m:m2,mutations,fitness,death,model,lower.tail)      # truncated distribution function
	if (lower.tail & pp <= P2[1]) return(k)
	if (!lower.tail & pp >= P2[1]) return(k)
	else {
	  k <- k-1+if(lower.tail) max(which(P2<pp)) else max(which(P2>pp))           # quantile if in the table
	  m <- m2
	}
      }                              # end while
    }                                  # end if
    return(k)
  })

}

  # Distribution function of the mutants count with parameter
  #   mutations : mean number of mutations
  #   fitness : fitness parameter
  #   death : death probability
  #   model : lifetimes distribution model : exponentialy distributed lifetimes ("LD") or constant lifetimes ("H")


pflan <- function(m,mutations=1,fitness=1,death=0,model=c("LD","H"),lower.tail=TRUE){

  if(sum(m < 0) > 0 | sum(trunc(m) != m) > 0){
    stop("'m' must be a vector of positive integers")
  }
  if(mutations < 0 | length(mutations) > 1){
    stop("'mutations' must be a single positive number.")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }
  if(missing(model)){model="LD"}
  model <- match.arg(model)

#   output=.Call("pflan",
#     list(m=as.integer(max(m)),mutations=as.double(mutations),fitness=as.double(fitness),death=as.double(death),model=as.character(model),lowerTail=as.logical(lower.tail),bool=as.logical(length(m)>1)))

  flan.mutmodel <- new(FlanMutMod,list(
    mutations=mutations,
    fitness=fitness,
    death=death,
    model=model,
    lt=lower.tail
  ))
  
  output <- flan.mutmodel$pflan(max(m))

  output[m+1]
}


  # Density function of the mutants count with parameter
  #   mutations : mean number of mutations
  #   fitness : fitness parameter
  #   death : death probability
  #   model : lifetimes distribution model : exponentialy distributed lifetimes ("LD") or constant lifetimes ("H")


dflan <- function(m,mutations=1,fitness=1,death=0,model=c("LD","H")){

  if(sum(m < 0) > 0 | sum(trunc(m) != m) > 0){
    stop("'m' must be a vector of positive integers")
  }
  if(mutations < 0 | length(mutations) > 1){
    stop("'mutations' must be a single positive number.")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }
  if(missing(model)){model="LD"}
  model <- match.arg(model)

#   output=.Call("dflan",
#   list(m=as.integer(max(m)),mutations=as.double(mutations),fitness=as.double(fitness),death=as.double(death),model=as.character(model),bool=as.logical(length(m)>1)))
  
  flan.mutmodel <- new(FlanMutMod,list(
    mutations=mutations,
    fitness=fitness,
    death=death,
    model=model,
    lt=TRUE
  ))
  
  
  output <- flan.mutmodel$dflan(max(m))

  output[m+1]
}
