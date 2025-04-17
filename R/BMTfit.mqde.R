#'@title Minimum Quantile Distance Fit of the BMT Distribution to Non-censored 
#'  Data.
#'  
#'@description Fit of the BMT distribution to non-censored data by minimum 
#'  quantile distance (mqde), which can also be called maximum quantile 
#'  goodness-of-fit.
#'  
#'@name BMTfit.mqde
#'  
#'@details This function is not intended to be called directly but is internally
#'  called in \code{\link{BMTfit}} when used with the minimum quantile distance 
#'  method.
#'  
#'  \code{BMTfit.mqde} is based on the function \code{\link{mqdedist}} but it 
#'  focuses on the minimum quantile distance parameter estimation for the BMT 
#'  distribution (see \code{\link{BMT}} for details about the BMT distribution 
#'  and \code{\link{mqdedist}} for details about minimum quantile distance fit 
#'  of univariate distributions).
#'  
#'  Given the close-form expression of the quantile 
#'  function, two optimization methods were added when the euclidean distance is
#'  selected: Coordinate descend (\code{"CD"}) and Newton-Rhapson (\code{"NR"}).
#'  
#'@param data A numeric vector with the observed values for non-censored data.
#'@param probs A numeric vector of the probabilities for which the minimum 
#'  quantile distance estimation is done. \eqn{p[k] = (k - 0.5) / n} (default).
#'@param qtype The quantile type used by the R \code{\link{quantile}} function 
#'  to compute the empirical quantiles. Type 5 (default), i.e. \eqn{x[k]} is 
#'  both the \eqn{k}th order statistic and the type 5 sample quantile of 
#'  \eqn{p[k] = (k - 0.5) / n}.
#'@param dist The distance measure between observed and theoretical quantiles to
#'  be used. This must be one of "euclidean" (default), "maximum", or 
#'  "manhattan". Any unambiguous substring can be given.
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
#'  \code{\link{optim}}. Given the close-form expression of the quantile 
#'  function, two optimization methods were added when the euclidean distance is
#'  selected: Coordinate descend (\code{"CD"}) and Newton-Rhapson (\code{"NR"}).
#'@param custom.optim A function carrying the optimization (see the 'details' 
#'  section of \code{\link[fitdistrplus]{mledist}}).
#'@param weights an optional vector of weights to be used in the fitting process. 
#'  Should be \code{NULL} or a numeric vector with strictly positive numbers. 
#'  If non-\code{NULL}, weighted mqde is used, otherwise ordinary mqde.
#'@param silent A logical to remove or show warnings when bootstrapping.
#'@param \dots Further arguments to be passed to generic functions or to the 
#'  function \code{"mqdedist"}. See \code{\link{mqdedist}} for details.
#'  
#'@return \code{BMTfit.mqde} returns a list with following components,
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
#'  \item{value}{the value of the corresponding objective function of the 
#'  estimation method at the estimate.}
#'  
#'  \item{hessian}{a symmetric matrix computed by \code{\link{optim}} as an 
#'  estimate of the Hessian at the solution found or computed in the 
#'  user-supplied optimization function.}
#'  
#'  \item{loglik}{the log-likelihood value.}
#'  
#'  \item{probs}{ the probability vector on which observed and theoretical 
#'  quantiles were calculated. }
#'  
#'  \item{dist}{ the name of the distance between observed and theoretical 
#'  quantiles used. }
#'  
#'  \item{optim.function}{the name of the optimization function used for maximum
#'  product of spacing.}
#'  
#'  \item{optim.method}{when \code{\link{optim}} is used, the name of the 
#'  algorithm used, \code{NULL} otherwise.}
#'  
#'  \item{fix.arg}{the named list giving the values of parameters of the named 
#'  distribution that must kept fixed rather than estimated or \code{NULL} if 
#'  there are no such parameters. }
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
#'  \code{\link{BMTfit.mle}}, \code{\link{BMTfit.mge}}, 
#'  \code{\link{BMTfit.mpse}} and \code{\link{BMTfit.qme}} for other estimation 
#'  methods. See \code{\link{optim}} and \code{\link{constrOptim}} for 
#'  optimization routines. See \code{\link{BMTfit}} and \code{\link[fitdistrplus]{fitdist}} 
#'  for functions that return an object of class \code{"fitdist"}.
#'  
#'@author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co}
#'  
#'@source Based on the function \code{\link{mqdedist}} which in turn is based on
#'  the function \code{\link[fitdistrplus]{mledist}} of the R package: 
#'  \pkg{fitdistrplus}
#'  
#'  Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package 
#'  for Fitting Distributions}. Journal of Statistical Software, 64(4), 1-34.
#'  
#' @examples
#' # (1) basic fit by minimum quantile distance estimation
#' set.seed(1234)
#' x1 <- rBMT(n=100, p3=0.25, p4=0.75)
#' BMTfit.mqde(x1)
#' 
#' # (2) quantile matching is a particular case of minimum quantile distance
#' BMTfit.mqde(x1, probs=c(0.2,0.4,0.6,0.8), qtype=7)
#' 
#' # (3) maximum or manhattan instead of euclidean distance
#' BMTfit.mqde(x1, dist="maximum")
#' BMTfit.mqde(x1, dist="manhattan")
#' 
#' # (4) how to change the optimisation method?
#' BMTfit.mqde(x1, optim.method="L-BFGS-B") 
#' BMTfit.mqde(x1, custom.optim="nlminb")
#' 
#' # (5) estimation of the tails weights parameters of the BMT 
#' # distribution with domain fixed at [0,1]
#' BMTfit.mqde(x1, start=list(p3=0.5, p4=0.5), fix.arg=list(p1=0, p2=1))
#' 
#' # (6) estimation of the asymmetry-steepness parameters of the BMT 
#' # distribution with domain fixed at [0,1]
#' BMTfit.mqde(x1, start=list(p3=0, p4=0.5), type.p.3.4 = "a-s", 
#'             fix.arg=list(p1=0, p2=1))
#' 
#'@keywords distribution

######################
#' @rdname BMTfit.mqde
#' @export
BMTfit.mqde <- function(data, probs = (1:length(data)-0.5)/length(data), qtype = 5, dist = "euclidean", 
                        start = list(p3 = 0.5, p4 = 0.5, p1 = min(data) - 0.1, p2 = max(data) + 0.1),
                        fix.arg = NULL, type.p.3.4 = "t w", type.p.1.2 = "c-d", 
                        optim.method = "Nelder-Mead", custom.optim = NULL, weights = NULL, silent = TRUE, ...){
  # Control data
  if (!(is.vector(data) & is.numeric(data) & length(data) > 1))
    stop("data must be a numeric vector of length greater than 1")
  # Control probs, qtype
  if (!(is.vector(probs) & is.numeric(probs)) | anyNA(probs) | any(probs < 0 | probs > 1))
    stop("probs must be a numeric vector with all elements greater than zero and less than one")
  probs <- unique(sort(probs))
  if (qtype < 0 || qtype > 9) 
    stop("wrong type for the R quantile function")
  # Control dist
  int.dist <- pmatch(dist, c("euclidean", "maximum", "manhattan"))
  if (is.na(int.dist))
    stop("invalid distance measure to be used")
  if (int.dist == -1)
    stop("ambiguous distance measure to be used")
  # Control optim.method
  if (is.null(custom.optim)) 
    optim.method <- match.arg(optim.method, c("default", "Nelder-Mead", "BFGS", "CG", 
                                              "L-BFGS-B", "SANN", "Brent", "CD","NR"))
  # Control start and fix.arg
  start.arg <- start
  if (is.vector(start.arg)) 
    start.arg <- as.list(start.arg)
  stnames <- names(start.arg)
  fixnames <- names(fix.arg)
  
  # Further arguments to be passed
  my3dots <- list(...)    
  if (length(my3dots) == 0) 
    my3dots <- NULL
  if (is.vector(data)) {
    n <- length(data)
    if (!(is.numeric(data) & n > 1)) 
      stop("data must be a numeric vector of length greater than 1")
  }
  else
    stop("Minimum quantile distance estimation is not yet available for censored data.")
  # Control weights
  if (!is.null(weights)) {
    if (any(weights <= 0)) 
      stop("weights should be a vector of numbers greater than 0")
    if (length(weights) != n) 
      stop("weights should be a vector with a length equal to the observation number")
    w <- sum(weights)
  }
  # Control maximum number of iterations
  maxit <- ifelse(is.null(my3dots$control$maxit),3000,my3dots$control$maxit)
  # Control type.p.3.4. It allows partial match.
  TYPE.P.3.4 <- c("t w", "a-s") # tail weights or asymmetry-steepness
  int.type.p.3.4 <- pmatch(type.p.3.4, TYPE.P.3.4)
  if (is.na(int.type.p.3.4))
    stop("invalid type of parametrization for parameters 3 and 4")
  if (int.type.p.3.4 == -1)
    stop("ambiguous type of parametrization for parameters 3 and 4")
  # Control type.p.1.2. It allows partial match.
  TYPE.P.1.2 <- c("c-d", "l-s") # domain or location-scale
  int.type.p.1.2 <- pmatch(type.p.1.2, TYPE.P.1.2)
  if (is.na(int.type.p.1.2))
    stop("invalid type of parametrization for parameters 1 and 2")
  if (int.type.p.1.2 == -1)
    stop("ambiguous type of parametrization for parameters 1 and 2")
  # Type of parametrizations are fixed parameters
  fix.arg$type.p.3.4 <- type.p.3.4
  fix.arg$type.p.1.2 <- type.p.1.2
  # Establish box constraints according to parameters in start
  npar <- length(stnames)
  # Initialize all box constraints: (-Inf, Inf)
  if(is.null(my3dots$lower)){
    lower <- rep(-Inf, npar)
    # domain parametrization
    if (int.type.p.1.2 == 1) {
      # c has to be inside (-Inf, min(data)) 
      # d has to be inside (max(data), Inf)
      lower[stnames == "p2"] <- max(data) + .epsilon
    }
    # location-scale parametrization
    else{
      # sigma has to be inside (0, Inf) 
      lower[stnames == "p2"] <- 0 + .epsilon
    }
    # tail weights parametrization
    if (int.type.p.3.4 == 1) {
      # Both tail weights have to be inside (0,1) 
      lower[stnames == "p3" | stnames == "p4"] <- 0 + .epsilon
    }
    # asymmetry-steepness parametrization
    else{
      # asymmetric has to be inside (-1, 1) 
      # steepness has to be inside (0, 1)
      lower[stnames == "p3"] <- -1 + .epsilon
      lower[stnames == "p4"] <- 0 + .epsilon
    }
  }
  else{
    lower <- my3dots$lower
    my3dots$lower <- NULL
  }
  if(is.null(my3dots$upper)){
    upper <- rep(Inf, npar)
    # domain parametrization
    if (int.type.p.1.2 == 1) {
      # c has to be inside (-Inf, min(data)) 
      upper[stnames == "p1"] <- min(data) - .epsilon
      # d has to be inside (max(data), Inf)
    }
    # location-scale parametrization
    # sigma has to be inside (0, Inf) 
    # tail weights parametrization
    if (int.type.p.3.4 == 1) {
      # Both tail weights have to be inside (0,1) 
      upper[stnames == "p3" | stnames == "p4"] <- 1 - .epsilon
    }
    # asymmetry-steepness parametrization
    else{
      # asymmetric has to be inside (-1, 1) 
      # steepness has to be inside (0, 1)
      upper[stnames == "p3" | stnames == "p4"] <- 1 - .epsilon
    }
  }
  else{
    upper <- my3dots$upper
    my3dots$upper <- NULL
  }
  names(upper) <- names(lower) <- stnames
  # nlminb optimization method
  if(!is.null(custom.optim)){
    if(custom.optim=="nlminb")
      custom.optim <- .m.nlminb
    # mqdedist function
    mqde <- do.call(mqdedist, append(list(data, "BMT", probs = probs, qtype = qtype, dist = dist, 
                     start = start, fix.arg = fix.arg, optim.method = optim.method, 
                     lower = lower, upper = upper, custom.optim = custom.optim, 
                     weights = weights, silent = silent), my3dots))
  }
  else{
    if(optim.method == "CD"){
      if(int.dist != 1)
        stop("Coordinate descend (CD) optimization metod is only considered with euclidean distance")
      # Change to tails curvature and domain parameterization
      par <- append(start.arg, fix.arg)
      if(int.type.p.3.4 != 1){
        c.par <- BMTchangepars(par$p3, par$p4, par$type.p.3.4, par$p1, par$p2, par$type.p.1.2)
        par$p3 <- c.par$p3
        par$p4 <- c.par$p4
        par$type.p.3.4 <- c.par$type.p.3.4
      }
      # Sample quantiles
      if(qtype == 0)
        sq <- data
      else
        sq <- quantile(data, probs = probs, type = qtype, names = FALSE)
      t <- 0.5-cos((acos(2*probs-1)-2*pi)/3)
      a <- 3 * t * (t - 1)^2
      b <- 3 * t * t * (t - 1) 
      c <- t * t * (3 - 2 * t)
      # 
      flag <- FALSE
      iter <- 0
      conv <- 0
      message <- NULL
      if(is.null(weights)){
        # theorical quantiles
        theoq <- do.call("qBMT", c(list(p = probs), as.list(par[stnames]), as.list(par[names(fix.arg)])))
        # objective function
        value <- sum((sq - theoq)^2)
        repeat{
          value.old <- value
          iter <- iter + 1
          if(!("p3" %in% fixnames) | !("p4" %in% fixnames)){
            sx <- (sq - par$p1)/(par$p2 - par$p1)
            if(!("p4" %in% fixnames)){
              a.b <- sum(a*b)
              b.b <- sum(b*b)
              b.c_sx <- sum(b*(c-sx))
              if(!("p3" %in% fixnames)){
                a.a <- sum(a*a)
                a.c_sx <- sum(a*(c-sx)) 
                par$p3 <- (a.b * b.c_sx - b.b * a.c_sx ) / (a.a * b.b - a.b * a.b)
                par$p3 <- unname(ifelse(par$p3 > upper["p3"], upper["p3"], par$p3))
                par$p3 <- unname(ifelse(par$p3 < lower["p3"], lower["p3"], par$p3))
              }
              par$p4 <- (- par$p3 * a.b - b.c_sx) / b.b
              par$p4 <- unname(ifelse(par$p4 > upper["p4"], upper["p4"], par$p4))
              par$p4 <- unname(ifelse(par$p4 < lower["p4"], lower["p4"], par$p4))
            }
            else{
              if(!("p3" %in% fixnames)){
                par$p3 <- (- par$p4 * sum(a*b) - sum(a*(c-sx))) / sum(a*a)
                par$p3 <- unname(ifelse(par$p3 > upper["p3"], upper["p3"], par$p3))
                par$p3 <- unname(ifelse(par$p3 < lower["p3"], lower["p3"], par$p3))
              }
            }
          }
          else
            flag <- TRUE
          if(!("p1" %in% fixnames) | !("p2" %in% fixnames)){
            tx <- a * par$p3 + b * par$p4 + c
            mean.sq <- mean(sq)
            mean.tx <- mean(tx)
            if(!("p2" %in% fixnames)){
              mean.tx2 <- mean(tx^2)
              mean.sq.tx <- mean(sq * tx)
              if(!("p1" %in% fixnames)){
                par$p1 <- (mean.sq * mean.tx2 - mean.sq.tx * mean.tx) / (mean.tx2 - mean.tx^2)
                if(int.type.p.1.2 == 1)
                  par$p1 <- unname(ifelse(par$p1 > upper["p1"], upper["p1"], par$p1))
              }
              par$p2 <- par$p1 + (mean.sq.tx - par$p1 * mean.tx) / mean.tx2
              if(int.type.p.1.2 == 1)
                par$p2 <- unname(ifelse(par$p2 < lower["p2"], lower["p2"], par$p2))
            }
            else{
              if(!("p1" %in% fixnames)){
                par$p1 <- (mean.sq - par$p2 * mean.tx) / (1 - mean.tx)
                if(int.type.p.1.2 == 1)
                  par$p1 <- unname(ifelse(par$p1 > upper["p1"], upper["p1"], par$p1))
              }
            }
          }
          else
            flag <- TRUE
          # theorical quantiles
          theoq <- do.call("qBMT", c(list(p = probs), as.list(par[stnames]), as.list(par[names(fix.arg)])))
          # objective function
          value <- sum((sq - theoq)^2)
          if(abs(value - value.old) < 1e-12 | iter == maxit | flag == TRUE)
            break
        }
      }
      else{
        # theorical quantiles
        theoq <- do.call("qBMT", c(list(p = probs), as.list(par[stnames]), as.list(par[names(fix.arg)])))
        # objective function
        value <- sum(weights*(sq - theoq)^2)
        repeat{
          value.old <- value
          iter <- iter + 1
          if(!("p3" %in% fixnames) | !("p4" %in% fixnames)){
            sx <- (sq - par$p1)/(par$p2 - par$p1)
            if(!("p4" %in% fixnames)){
              a.b <- sum(weights*a*b)
              b.b <- sum(weights*b*b)
              b.c_sx <- sum(weights*b*(c-sx))
              if(!("p3" %in% fixnames)){
                a.a <- sum(weights*a*a)
                a.c_sx <- sum(weights*a*(c-sx)) 
                par$p3 <- (a.b * b.c_sx - b.b * a.c_sx ) / (a.a * b.b - a.b * a.b)
                par$p3 <- unname(ifelse(par$p3 > upper["p3"], upper["p3"], par$p3))
                par$p3 <- unname(ifelse(par$p3 < lower["p3"], lower["p3"], par$p3))
              }
              par$p4 <- (- par$p3 * a.b - b.c_sx) / b.b
              par$p4 <- unname(ifelse(par$p4 > upper["p4"], upper["p4"], par$p4))
              par$p4 <- unname(ifelse(par$p4 < lower["p4"], lower["p4"], par$p4))
            }
            else{
              if(!("p3" %in% fixnames)){
                par$p3 <- (- par$p4 * sum(weights*a*b) - sum(weights*a*(c-sx))) / sum(weights*a*a)
                par$p3 <- unname(ifelse(par$p3 > upper["p3"], upper["p3"], par$p3))
                par$p3 <- unname(ifelse(par$p3 < lower["p3"], lower["p3"], par$p3))
              }
            }
          }
          else
            flag <- TRUE
          if(!("p1" %in% fixnames) | !("p2" %in% fixnames)){
            tx <- a * par$p3 + b * par$p4 + c
            mean.sq <- sum(weights * sq) / w
            mean.tx <- sum(weights * tx) / w
            if(!("p2" %in% fixnames)){
              mean.tx2 <- sum(weights * tx^2) / w
              mean.sq.tx <- sum(weights * sq * tx) / w
              if(!("p1" %in% fixnames)){
                par$p1 <- (mean.sq * mean.tx2 - mean.sq.tx * mean.tx) / (mean.tx2 - mean.tx^2)
                if(int.type.p.1.2 == 1)
                  par$p1 <- unname(ifelse(par$p1 > upper["p1"], upper["p1"], par$p1))
              }
              par$p2 <- par$p1 + (mean.sq.tx - par$p1 * mean.tx) / mean.tx2
              if(int.type.p.1.2 == 1)
                par$p2 <- unname(ifelse(par$p2 < lower["p2"], lower["p2"], par$p2))
            }
            else{
              if(!("p1" %in% fixnames)){
                par$p1 <- (mean.sq - par$p2 * mean.tx) / (1 - mean.tx)
                if(int.type.p.1.2 == 1)
                  par$p1 <- unname(ifelse(par$p1 > upper["p1"], upper["p1"], par$p1))
              }
            }
          }
          else
            flag <- TRUE
          # theorical quantiles
          theoq <- do.call("qBMT", c(list(p = probs), as.list(par[stnames]), as.list(par[names(fix.arg)])))
          # objective function
          value <- sum(weights*(sq - theoq)^2)
          if(abs(value - value.old) < 1e-12 | iter == maxit | flag == TRUE)
            break
        }
      }
      if(iter == maxit){
        conv <- 1
        message <- "Maximum number of iterations"
      }
      # Restore to the given parameterization
      if(int.type.p.3.4 != 1){
        c.par <- BMTchangepars(par$p3, par$p4, par$type.p.3.4, par$p1, par$p2, par$type.p.1.2)
        par$p3 <- c.par$p3
        par$p4 <- c.par$p4
        par$type.p.3.4 <- c.par$type.p.3.4
      }
      mqde <- list(estimate = unlist(par[stnames]), convergence = conv, value = value, 
                   hessian = NULL, probs = probs, dist = dist, 
                   optim.function = NULL, fix.arg = fix.arg, 
                   loglik = .loglik(par[stnames], fix.arg, data, "dBMT"), 
                   optim.method = "CD", fix.arg.fun = NULL, 
                   counts = c(iter, NA), optim.message = message)
    }
    else if(optim.method == "NR"){
      if(int.dist != 1)
        stop("Newton-Raphson (NR) optimization method is only considered with euclidean distance")
      # Change to tails curvature and domain parameterization
      par <- append(start.arg, fix.arg)
      if(int.type.p.3.4 != 1 | int.type.p.1.2 != 1){
        c.par <- BMTchangepars(par$p3,par$p4,par$type.p.3.4,par$p1,par$p2,par$type.p.1.2)
        if(int.type.p.3.4 != 1){
          par$p3 <- c.par$p3
          par$p4 <- c.par$p4
          par$type.p.3.4 <- c.par$type.p.3.4
        }
        if(int.type.p.1.2 != 1){
          par$p1 <- c.par$p1
          par$p2 <- c.par$p2
          par$type.p.1.2 <- c.par$type.p.1.2
        }
      }
      # Sample quantiles
      if(qtype == 0)
        sq <- data
      else
        sq <- quantile(data, probs = probs, type = qtype, names = FALSE)
      t <- 0.5-cos((acos(2*probs-1)-2*pi)/3)
      a <- 3 * t * (t - 1)^2
      b <- 3 * t * t * (t - 1) 
      c <- t * t * (3 - 2 * t)
      # 
      iter <- 0
      conv <- 0
      message <- NULL
      if(is.null(weights)){
        repeat{
          iter <- iter + 1
          R <- par$p2 - par$p1
          sx <- (sq - par$p1) / R
          tx <- a * par$p3 + b * par$p4 + c
          sx_tx <- sx - tx
          sxsx_tx <- 2*sx - tx
          a.sx_tx <- sum(a * sx_tx)
          b.sx_tx <- sum(b * sx_tx)
          sx.sx_tx <- sum(sx * sx_tx)
          kl.1_kl <- par$p3 * (1 - par$p3)
          kr.1_kr <- par$p4 * (1 - par$p4)
          gradient <- -2 * c(kl.1_kl * a.sx_tx,
                             kr.1_kr * b.sx_tx,
                             sx.sx_tx,
                             sum(sx_tx) / R)
          hessian <- diag(c(kl.1_kl * ((2 * par$p3 - 1) * a.sx_tx + kl.1_kl * sum(a*a)),
                            kr.1_kr * ((2 * par$p4 - 1) * b.sx_tx + kr.1_kr * sum(b*b)),
                            sum(sx * sxsx_tx),
                            n / R^2))
          hessian[lower.tri(hessian)] <- 2*c(kl.1_kl * kr.1_kr * sum(a * b),
                                             kl.1_kl * sum(a * sx),
                                             kl.1_kl * sum(a) / R,
                                             kr.1_kr * sum(b * sx),
                                             kr.1_kr * sum(b) / R,
                                             sum(sxsx_tx) / R)
          hessian <- hessian + t(hessian)
          delta <- solve(hessian,gradient)
          theta <- c(log(par$p3/(1-par$p3)), log(par$p4/(1-par$p4)), log(R), par$p1)
          theta <- theta - delta
          par$p3 <- plogis(theta[1])
          par$p4 <- plogis(theta[2])
          par$p1 <- theta[4]
          par$p2 <- par$p1 + exp(theta[3])
          if(all(abs(delta) < 1e-12) | iter == maxit)
            break
        }
      }
      else{
        #
      }
      # theorical quantiles
      theoq <- do.call("qBMT", c(list(p = probs), as.list(par[stnames]), as.list(fix.arg)))
      # objective function
      value <- sum((sq - theoq)^2)
      if(iter == maxit){
        conv <- 1
        message <- "Maximum number of iterations"
      }
      # Restore to the given parameterization
      if(int.type.p.3.4 != 1 | int.type.p.1.2 != 1){
        c.par <- BMTchangepars(par$p3,par$p4,par$type.p.3.4,par$p1,par$p2,par$type.p.1.2)
        if(int.type.p.3.4 != 1){
          par$p3 <- c.par$p3
          par$p4 <- c.par$p4
          par$type.p.3.4 <- c.par$type.p.3.4
        }
        if(int.type.p.1.2 != 1){
          par$p1 <- c.par$p1
          par$p2 <- c.par$p2
          par$type.p.1.2 <- c.par$type.p.1.2
        }
      }
      mqde <- list(estimate = unlist(par[stnames]), convergence = conv, value = value, 
                   hessian = NULL, probs = probs, dist = dist, 
                   optim.function = NULL, fix.arg = fix.arg, 
                   loglik = .loglik(par[stnames], fix.arg, data, "dBMT"), 
                   optim.method = "NR", fix.arg.fun = NULL, 
                   counts = c(iter, iter), optim.message = message)
    }
    else{
      # mqdedist function
      mqde <- mqdedist(data, "BMT", probs = probs, qtype = qtype, dist = dist, 
                       start = start, fix.arg = fix.arg, optim.method = optim.method, 
                       lower = lower, upper = upper, custom.optim = custom.optim, 
                       weights = weights, silent = silent, ...)
    }
  }
  # Estimation with location-scale parametrization might allow data outside the estimated domain
  par <- append(mqde$estimate,fix.arg)
  if (int.type.p.1.2 == 2)
    par <- BMTchangepars(par$p3, par$p4, par$type.p.3.4, par$p1, par$p2, par$type.p.1.2)
  n.obs <- sum(data < par$p1 | data > par$p2)
  if(n.obs > 0){
    text <- paste("The resultant estimated domain is [", round(par$p1, 4), ",", round(par$p2, 4),
                  "] and there are ", n.obs, " observations out of it.", sep="")
    warning(text)
  }
  return(mqde)
}
