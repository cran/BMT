#'@title Maximum Likelihood Fit of the BMT Distribution to Non-censored Data.
#'  
#'@description Fit of the BMT distribution to non-censored data by maximum 
#'  likelihood estimation (mle).
#'  
#'@name BMTfit.mle
#'  
#'@details This function is not intended to be called directly but is internally
#'  called in \code{\link{BMTfit}} when used with the maximum likelihood method.
#'  
#'  \code{BMTfit.mle} is based on the function \code{\link[fitdistrplus]{mledist}} from the 
#'  package \pkg{fitdistrplus} but it focuses on the maximum likelihood 
#'  parameter estimation for the BMT distribution (see \pkg{BMT} for
#'  details about the BMT distribution and \code{\link[fitdistrplus]{mledist}} for details
#'  about maximum likelihood fit of univariate distributions).
#'  
#'@param data A numeric vector with the observed values for non-censored data.
#'@param start A named list giving the initial values of parameters of the BMT 
#'  distribution or a function of data computing initial values and returning a 
#'  named list. (see the 'details' section of 
#'  \code{\link[fitdistrplus]{mledist}}).
#'@param fix.arg An optional named list giving the values of fixed parameters of
#'  the BMT distribution or a function of data computing (fixed) parameter 
#'  values and returning a named list. Parameters with fixed value are thus NOT 
#'  estimated. (see the 'details' section of 
#'  \code{\link[fitdistrplus]{mledist}}).
#'@param type.p.3.4 Type of parametrization associated to p3 and p4. "t w" means 
#'  tails weights parametrization (default) and "a-s" means asymmetry-steepness 
#'  parametrization.
#'@param type.p.1.2 Type of parametrization associated to p1 and p2. "c-d" means 
#'  domain parametrization (default) and "l-s" means location-scale 
#'  parametrization.
#'@param optim.method \code{"default"} (see the 'details' section of 
#'  \code{\link[fitdistrplus]{mledist}}) or optimization method to pass to 
#'  \code{\link{optim}}.
#'@param custom.optim A function carrying the optimization (see the 'details' 
#'  section of \code{\link[fitdistrplus]{mledist}}).
#'@param silent A logical to remove or show warnings when bootstrapping.
#'@param \dots Further arguments to be passed to generic functions or to the 
#'  function \code{"mledist"}. See \code{\link[fitdistrplus]{mledist}} for details.
#'  
#'@return \code{BMTfit.mle} returns a list with following components,
#'  
#'  \item{estimate}{ the parameter estimates.}
#'  
#'  \item{convergence}{ an integer code for the convergence of 
#'  \code{\link{optim}}/\code{\link{constrOptim}} defined as below or defined by
#'  the user in the user-supplied optimization function.
#'  
#'  \code{0} indicates successful convergence.
#'  
#'  \code{1} indicates that the iteration limit of \code{\link{optim}} has been 
#'  reached.
#'  
#'  \code{10} indicates degeneracy of the Nealder-Mead simplex.
#'  
#'  \code{100} indicates that \code{\link{optim}} encountered an internal error.
#'  }
#'  
#'  \item{loglik}{the log-likelihood value.}
#'  
#'  \item{hessian}{a symmetric matrix computed by \code{\link{optim}} as an 
#'  estimate of the Hessian at the solution found or computed in the 
#'  user-supplied optimization function. It is used in \code{\link{BMTfit}} to estimate
#'  standard errors. }
#'  
#'  \item{optim.function}{the name of the optimization function used for maximum
#'  likelihood.}
#'  
#'  \item{optim.method}{when \code{\link{optim}} is used, the name of the 
#'  algorithm used, \code{NULL} otherwise.}
#'  
#'  \item{fix.arg}{the named list giving the values of parameters of the named 
#'  distribution that must kept fixed rather than estimated or \code{NULL} if there are no such parameters. }
#'  
#'  \item{fix.arg.fun}{the function used to set the value of \code{fix.arg} or 
#'  \code{NULL}.}
#'  
#'  \item{weights}{the vector of weights used in the estimation process or 
#'  \code{NULL}.}
#'  
#'  \item{counts}{A two-element integer vector giving the number of calls to the
#'  log-likelihood function and its gradient respectively. This excludes those 
#'  calls needed to compute the Hessian, if requested, and any calls to 
#'  log-likelihood function to compute a finite-difference approximation to the 
#'  gradient. \code{counts} is returned by \code{\link{optim}} or the 
#'  user-supplied function or set to \code{NULL}.}
#'  
#'  \item{optim.message}{A character string giving any additional information 
#'  returned by the optimizer, or \code{NULL}. To understand exactly the 
#'  message, see the source code.}
#'  
#'@references Torres-Jimenez, C. J. (2017, September), \emph{Comparison of estimation
#'  methods for the BMT distribution}. ArXiv e-prints.
#'  
#'  Torres-Jimenez, C. J. (2018), \emph{The BMT Item Response Theory model: A 
#'  new skewed distribution family with bounded domain and an IRT model based on
#'  it}, PhD thesis, Doctorado en ciencias - Estadistica, Universidad Nacional 
#'  de Colombia, Sede Bogota.
#'  
#'@seealso See \code{\link{BMT}} for the BMT density, distribution, quantile 
#'  function and random deviates. See \code{\link{BMTfit.mme}}, 
#'  \code{\link{BMTfit.qme}}, \code{\link{BMTfit.mge}}, 
#'  \code{\link{BMTfit.mpse}} and \code{\link{BMTfit.mqde}} for other estimation
#'  methods. See \code{\link{optim}} and \code{\link{constrOptim}} for 
#'  optimization routines. See \code{\link{BMTfit}} and \code{\link[fitdistrplus]{fitdist}} 
#'  for functions that return an object of class \code{"fitdist"}.
#'  
#'@author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co}
#'  
#'@source Based on the function \code{\link[fitdistrplus]{mledist}} of the R package: 
#'  \pkg{fitdistrplus}
#'  
#'  Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package 
#'  for Fitting Distributions}. Journal of Statistical Software, 64(4), 1-34.
#'  
#' @examples
#' # (1) basic fit by maximum likelihood estimation
#' set.seed(1234)
#' x1 <- rBMT(n=100, p3 = 0.25, p4 = 0.75)
#' BMTfit.mle(x1)
#' 
#' # (2) how to change the optimisation method?
#' BMTfit.mle(x1, optim.method="L-BFGS-B") 
#' BMTfit.mle(x1, custom.optim="nlminb")
#' 
#' # (3) estimation of the tails weights parameters of the BMT 
#' # distribution with domain fixed at [0,1]
#' BMTfit.mle(x1, start=list(p3=0.5, p4=0.5), fix.arg=list(p1=0, p2=1))
#' 
#' # (4) estimation of the asymmetry-steepness parameters of the BMT 
#' # distribution with domain fixed at [0,1]
#' BMTfit.mle(x1, start=list(p3=0, p4=0.5), type.p.3.4 = "a-s", 
#'            fix.arg=list(p1=0, p2=1))
#' 
#'@keywords distribution

#####################
#' @rdname BMTfit.mle
#' @export
BMTfit.mle <- function(data, 
                       start = list(p3 = 0.5, p4 = 0.5, p1 = min(data) - 0.1, p2 = max(data) + 0.1),
                       fix.arg = NULL, type.p.3.4 = "t w", type.p.1.2 = "c-d", 
                       optim.method = "Nelder-Mead", custom.optim = NULL, silent = TRUE, ...){
  # Control data
  if (!(is.vector(data) & is.numeric(data) & length(data) > 1))
    stop("data must be a numeric vector of length greater than 1")
  # Further arguments to be passed
  my3dots <- list(...)    
  if (length(my3dots) == 0) 
    my3dots <- NULL
  # Control weights
  if(!is.null(my3dots$weights))
    stop("Estimation with weights is not considered yet")
  # Control type.p.3.4. It allows partial match.
  TYPE.P.3.4 <- c("t w", "a-s") # tail weights or asymmetry-steepness
  int.type.p.3.4 <- pmatch(type.p.3.4, TYPE.P.3.4)
  if (is.na(int.type.p.3.4))
    stop("invalid type of parametrization for parameters 3 and 4")
  if (int.type.p.3.4 == -1)
    stop("ambiguous type of parametrization for parameters 3 and 4")
  # mle only allows parametrization "c-d" 
  # because all data have to be inside the estimated domain.
  if(type.p.1.2 != "c-d")
    stop("maximum likelihood estimation only allows parametrization \"c-d\"")
  # Type of parametrizations are fixed parameters
  fix.arg$type.p.3.4 <- type.p.3.4
  fix.arg$type.p.1.2 <- "c-d"
  # Establish box constraints according to parameters in start
  stnames <- names(start)
  m <- length(stnames)
  # Initialize all box constraints: (0, 1)
  lower <- rep(0 + .epsilon, m)
  upper <- rep(1 - .epsilon, m)
  # domain parametrization
  # c has to be inside (-Inf, min(data)) 
  lower[stnames == "p1"] <- -Inf
  upper[stnames == "p1"] <- min(data) - .epsilon
  # d has to be inside (max(data), Inf)
  lower[stnames == "p2"] <- max(data) + .epsilon
  upper[stnames == "p2"] <- Inf
  # asymmetry-steepness parametrization
  if(int.type.p.3.4 == 2) {
    # asymmetry has to be inside (-1, 1) 
    lower[stnames == "p3"] <- -1 + .epsilon
  }
  # nlminb optimization method
  if(!is.null(custom.optim))
    if(custom.optim=="nlminb")
      custom.optim <- .m.nlminb
  # mledist function of fitdistplus
  mle <- fitdistrplus::mledist(data, "BMT", start = start, fix.arg = fix.arg, 
                 optim.method = optim.method, lower = lower, upper = upper,
                 custom.optim = custom.optim, silent = silent, ...)
  return(mle)
}
