\name{flan-package}
\alias{flan-package}
\alias{flan}
\docType{package}
\encoding{UTF-8}



\title{Fluctuation analysis of mutant cell counts}

\description{
 Statistical tools for fluctuation analysis of mutant cells counts.
}

\details{
  \tabular{ll}{
    Package: \tab flan\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2016-02\cr
    License: \tab GPL\cr
  }
  
  This package contains functions dedicated to :
  
  - Distribution (dflan, pflan, qflan, and rflan): built as the classic distribution functions, like dnorm, pnorm, qnorm, and rnorm.
  
  - Parametric estimation (mutestim): estimates of the parameters of interest, which are the mean number of mutations (or the mutation probability) and the fitness parameter. 
  
  - Parametric testing (flan.test): built as the classic test functions, like t.test.
}

\author{
  Adrien Mazoyer, Stéphane Despreaux, and Bernard Ycart
  
  Maintainer: Adrien Mazoyer  <adrien.mazoyer@imag.fr>
}

\references{
  B. Ycart and N. Veziris: Unbiased estimates of mutation rates under fluctuating final counts. 
PLoS one 9(7) e101434 (2014) 
  
  B. Ycart: Fluctuation analysis with cell deaths.
 J. Applied Probab. Statist, 9(1):12-28 (2014)
 
  B. Ycart: Fluctuation analysis: can estimates be trusted? 
One PLoS one 8(12) e80958 (2013)
 
  A. Hamon and B. Ycart: Statistics for the Luria-Delbrück distribution. 
Elect. J. Statist., 6:1251-1272 (2012)
}


\keyword{package}

