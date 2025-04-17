#' @title Fit of the BMT Distribution to Non-censored Data.
#'   
#' @description Fit of the BMT distribution to non-censored data by maximum 
#'   likelihood (mle), moment matching (mme), quantile matching (qme), maximum 
#'   goodness-of-fit (mge), also known as minimum distance, maximum product of 
#'   spacing (mpse), also called maximum spacing, and minimum quantile distance 
#'   (mqde), which can also be called maximum quantile goodness-of-fit.
#'   
#' @name BMTfit
#'   
#' @details This function is based on the function \code{\link[fitdistrplus]{fitdist}} 
#'   but it focuses on the parameter 
#'   estimation for the BMT distribution (see \code{\link{BMT}} for details). It
#'   has six possible fitting methods: maximum likelihood (mle), moment matching
#'   (mme), quantile matching (qme), maximum goodness-of-fit (mge), also known 
#'   as minimum distance, maximum product of spacing (mpse), also called maximum
#'   spacing, and minimum quantile distance (mqde), which can also be called 
#'   maximum quantile goodness-of-fit. These fitting methods are carried out in 
#'   \code{\link{BMTfit.mle}}, \code{\link{BMTfit.mme}}, 
#'   \code{\link{BMTfit.qme}}, \code{\link{BMTfit.mge}}, 
#'   \code{\link{BMTfit.mpse}}, and \code{\link{BMTfit.mqde}}, respectively (see
#'   each function for details). \code{BMTfit} returns an object of class 
#'   \code{"fitdist"} (see \code{\link[fitdistrplus]{fitdist}} for details). Therefore, it 
#'   benefits of all the developed functions and methods for that class (see 
#'   \code{\link[fitdistrplus]{fitdistrplus}} for details).
#'   
#'   Generic methods of a \code{\link[fitdistrplus]{fitdist}} object are \code{print}, 
#'   \code{plot}, \code{summary}, \code{quantile}, \code{logLik}, \code{vcov} 
#'   and \code{coef}.
#'   
#' @param data A numeric vector with the observed values for non-censored data.
#' @param method A character string coding for the fitting method: \code{"mle"} 
#'   for 'maximum likelihood estimation', \code{"mme"} for 'moment matching 
#'   estimation', \code{"qme"} for 'quantile matching estimation', \code{"mge"} 
#'   for 'maximum goodness-of-fit estimation', \code{"mpse"} for 'maximum 
#'   product of spacing estimation', and \code{"mqde"} for 'minimum quantile 
#'   estimation'.
#' @param start A named list giving the initial values of parameters of the BMT 
#'   distribution or a function of data computing initial values and returning a
#'   named list. (see the 'details' section of 
#'   \code{\link[fitdistrplus]{mledist}}).
#' @param fix.arg An optional named list giving the values of fixed parameters 
#'   of the BMT distribution or a function of data computing (fixed) parameter 
#'   values and returning a named list. Parameters with fixed value are thus NOT
#'   estimated. (see the 'details' section of 
#'   \code{\link[fitdistrplus]{mledist}}).
#' @param type.p.3.4 Type of parametrization associated to p3 and p4. "t w" means
#'   tails weights parametrization (default) and "a-s" means asymmetry-steepness
#'   parametrization.
#' @param type.p.1.2 Type of parametrization associated to p1 and p2. "c-d" means
#'   domain parametrization (default) and "l-s" means location-scale 
#'   parametrization.
#' @param optim.method \code{"default"} (see the 'details' section of 
#'   \code{\link[fitdistrplus]{mledist}}) or optimization method to pass to 
#'   \code{\link{optim}}.
#' @param custom.optim A function carrying the optimization (see the 'details' 
#'   section of \code{\link[fitdistrplus]{mledist}}).
#' @param keepdata A logical. If \code{TRUE}, dataset is returned, otherwise 
#'   only a sample subset is returned.
#' @param keepdata.nb When \code{keepdata=FALSE}, the length (>1) of the subset 
#'   returned.
#' @param \dots Further arguments to be passed to generic functions, or to one 
#'   of the functions \code{"BMTfit.mle"},  \code{"BMTfit.mme"}, 
#'   \code{"BMTfit.qme"}, \code{"BMTfit.mge"}, \code{"BMTfit.mpse"}, or 
#'   \code{"BMTfit.mqde"} depending of the chosen method. See 
#'   \code{\link{BMTfit.mle}}, \code{\link{BMTfit.mme}}, 
#'   \code{\link{BMTfit.qme}}, \code{\link{BMTfit.mge}}, 
#'   \code{\link{BMTfit.mpse}}, \code{\link{BMTfit.mqde}} for details on 
#'   parameter estimation.
#'   
#' @return \code{fitdist} returns an object of class \code{"fitdist"} with the 
#'   following components:
#'   
#'   \item{estimate }{ the parameter estimates.}
#'   
#'   \item{method }{ the character string coding for the fitting method : 
#'   \code{"mle"} for 'maximum likelihood estimation', \code{"mme"} for 'moment 
#'   matching estimation', \code{"qme"} for 'quantile matching estimation', 
#'   \code{"mge"} for 'maximum goodness-of-fit estimation', \code{"mpse"} for 
#'   'maximum product of spacing estimation', and \code{"mqde"} for 'minimum 
#'   quantile estimation'.}
#'   
#'   \item{sd}{ the estimated standard errors, \code{NA} if numerically not 
#'   computable or \code{NULL} if not available.}
#'   
#'   \item{cor}{ the estimated correlation matrix, \code{NA} if numerically not 
#'   computable or \code{NULL} if not available.}
#'   
#'   \item{vcov}{ the estimated variance-covariance matrix, \code{NULL} if not 
#'   available.}
#'   
#'   \item{loglik}{ the log-likelihood.}
#'   
#'   \item{aic}{ the Akaike information criterion.}
#'   
#'   \item{bic}{ the the so-called BIC or SBC (Schwarz Bayesian criterion).}
#'   
#'   \item{n}{ the length of the data set.}
#'   
#'   \item{data}{ the data set.}
#'   
#'   \item{distname}{ the name of the distribution (BMT).}
#'   
#'   \item{fix.arg}{ the named list giving the values of parameters of the named
#'   distribution that must be kept fixed rather than estimated or \code{NULL}
#'   if there are no such parameters. }
#'   
#'   \item{fix.arg.fun}{the function used to set the value of \code{fix.arg} or 
#'   \code{NULL}.}
#'   
#'   \item{discrete}{ the input argument or the automatic definition by the 
#'   function to be passed to functions \code{\link[fitdistrplus]{gofstat}}, 
#'   \code{\link[fitdistrplus]{plotdist}} and \code{\link[fitdistrplus]{cdfcomp}}. }
#'   
#'   \item{dots}{ the list of  further arguments passed in \dots to be used in 
#'   \code{\link[fitdistrplus]{bootdist}} in iterative calls to \code{\link[fitdistrplus]{mledist}}, 
#'   \code{\link[fitdistrplus]{mmedist}}, \code{\link[fitdistrplus]{qmedist}}, \code{\link[fitdistrplus]{mgedist}},
#'   \code{\link{mpsedist}}, \code{\link{mqdedist}} or \code{NULL} if no such
#'   arguments.}
#'   
#'   \item{weights}{the vector of weights used in the estimation process or 
#'   \code{NULL}.}
#'   
#' @references Torres-Jimenez, C. J. (2017, September), \emph{Comparison of 
#'   estimation methods for the BMT distribution}. ArXiv e-prints.
#'   
#'   Torres-Jimenez, C. J. (2018), \emph{The BMT Item Response Theory model: A 
#'   new skewed distribution family with bounded domain and an IRT model based 
#'   on it}, PhD thesis, Doctorado en ciencias - Estadistica, Universidad 
#'   Nacional de Colombia, Sede Bogota.
#'   
#' @seealso See \code{\link{BMT}} for the BMT density, distribution, quantile 
#'   function and random deviates. See \code{\link{BMTfit.mle}}, 
#'   \code{\link{BMTfit.mme}}, \code{\link{BMTfit.qme}}, 
#'   \code{\link{BMTfit.mge}}, \code{\link{BMTfit.mpse}} and 
#'   \code{\link{BMTfit.mqde}} for details on parameter estimation. See 
#'   \code{\link[fitdistrplus]{fitdist}} for details on the object fitdist and its methods 
#'   \code{print}, \code{plot}, \code{summary}, \code{quantile}, \code{logLik}, 
#'   \code{vcov} and \code{coef}, and \pkg{fitdistrplus} for an overview
#'   of the package to which that object belongs to.
#'   
#' @author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co}
#'   
#' @source Based on the function \code{\link[fitdistrplus]{fitdist}} of the R package: 
#'   \pkg{fitdistrplus}
#'   
#'   Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package 
#'   for Fitting Distributions}. Journal of Statistical Software, 64(4), 1-34.
#'   
#' @examples
#' # (1) fit of the BMT distribution by maximum likelihood estimation
#' data(groundbeef)
#' serving <- groundbeef$serving
#' fit.mle <- BMTfit(serving)
#' summary(fit.mle)
#' plot(fit.mle)
#' plot(fit.mle, demp = TRUE)
#' plot(fit.mle, histo = FALSE, demp = TRUE)
#' cdfcomp(fit.mle, addlegend=FALSE)
#' denscomp(fit.mle, addlegend=FALSE)
#' ppcomp(fit.mle, addlegend=FALSE)
#' qqcomp(fit.mle, addlegend=FALSE)
#' 
#' # (2) Comparison of various estimation methods
#' fit.mme <- BMTfit(serving, method="mme")
#' fit.mpse <- BMTfit(serving, method="mpse")
#' fit.mqde <- BMTfit(serving, method="mqde")
#' summary(fit.mme)
#' summary(fit.mpse)
#' summary(fit.mqde)
#' cdfcomp(list(fit.mle, fit.mme, fit.mpse, fit.mqde), 
#'         legendtext=c("mle", "mme", "mpse", "mqde"))
#' denscomp(list(fit.mle, fit.mme, fit.mpse, fit.mqde), 
#'          legendtext=c("mle", "mme", "mpse", "mqde"))
#' qqcomp(list(fit.mle, fit.mme, fit.mpse, fit.mqde), 
#'        legendtext=c("mle", "mme", "mpse", "mqde"))
#' ppcomp(list(fit.mle, fit.mme, fit.mpse, fit.mqde), 
#'        legendtext=c("mle", "mme", "mpse", "mqde"))
#' gofstat(list(fit.mle, fit.mme, fit.mpse, fit.mqde), 
#'         fitnames=c("mle", "mme", "mpse", "mqde"))
#' 
#' # (3) how to change the optimisation method?
#' BMTfit(serving, optim.method="Nelder-Mead")
#' BMTfit(serving, optim.method="L-BFGS-B") 
#' BMTfit(serving, custom.optim="nlminb")
#' 
#' # (4) estimation of the tails weights parameters of the BMT distribution 
#' # with domain fixed at [9,201] using Kolmogorov-Smirnov
#' fit.KS <- BMTfit(serving, method="mge", gof="KS", 
#'                  start=list(p3=0.5, p4=0.5), fix.arg=list(p1=9, p2=201))
#' summary(fit.KS)
#' plot(fit.KS)
#' 
#' # (5) estimation of the asymmetry-steepness parameters of the BMT 
#' # distribution with domain fixed at [9,201] using minimum quantile distance 
#' # with a closed formula (optim.method="CD")
#' fit.mqde.CD <- BMTfit(serving, method="mqde", optim.method="CD", 
#'                       start=list(p3=0.5, p4=0.5), type.p.3.4 = "a-s", 
#'                       fix.arg=list(p1=9, p2=201))
#' summary(fit.mqde.CD)
#' plot(fit.mqde.CD)
#' 
#' @keywords distribution

#################
#' @import stats
#' @import utils
#' @import partitions
#' @import fitdistrplus
#' @rdname BMTfit
#' @export 
BMTfit <- function(data, method = c("mle","mme","qme","mge","mpse","mqde"),
                   start = list(p3 = 0.5, p4 = 0.5, p1 = min(data) - 0.1, p2 = max(data) + 0.1),
                   fix.arg = NULL, type.p.3.4 = "t w", type.p.1.2 = "c-d",
                   optim.method = "Nelder-Mead", custom.optim = NULL,
                   keepdata = TRUE, keepdata.nb = 100, ...) {
  # Control keepdata and keepdata.nb
  if (!is.logical(keepdata) || !is.numeric(keepdata.nb) || keepdata.nb < 2)
    stop("wrong arguments 'keepdata' and 'keepdata.nb'")
  # Further arguments to be passed
  my3dots <- list(...)    
  if (length(my3dots) == 0) 
    my3dots <- NULL
  # Lenght of data
  n <- length(data)
  # Control method
  method <- match.arg(method, c("mle","mme","qme","mge","mpse","mqde"))
  # Separation for each estimation method
  res <- switch (method, 
                 mle = BMTfit.mle(data, start=start, fix.arg=fix.arg, 
                                  type.p.3.4=type.p.3.4, type.p.1.2=type.p.1.2, 
                                  optim.method=optim.method, custom.optim=custom.optim, ...),
                 mme = BMTfit.mme(data, start=start, fix.arg=fix.arg, 
                                  type.p.3.4=type.p.3.4, type.p.1.2=type.p.1.2, 
                                  optim.method=optim.method, custom.optim=custom.optim, ...),
                 qme = BMTfit.qme(data, start=start, fix.arg=fix.arg, 
                                  type.p.3.4=type.p.3.4, type.p.1.2=type.p.1.2, 
                                  optim.method=optim.method, custom.optim=custom.optim, ...),
                 mge = BMTfit.mge(data, start=start, fix.arg=fix.arg, 
                                  type.p.3.4=type.p.3.4, type.p.1.2=type.p.1.2, 
                                  optim.method=optim.method, custom.optim=custom.optim, ...),
                 mpse = BMTfit.mpse(data, start=start, fix.arg=fix.arg, 
                                    type.p.3.4=type.p.3.4, type.p.1.2=type.p.1.2, 
                                    optim.method=optim.method, custom.optim=custom.optim, ...),
                 mqde = BMTfit.mqde(data, start=start, fix.arg=fix.arg, 
                                    type.p.3.4=type.p.3.4, type.p.1.2=type.p.1.2, 
                                    optim.method=optim.method, custom.optim=custom.optim, ...))
  # Optimization method message
  if(!is.null(res$optim.message))
    cat("\noptim.message:",res$optim.message,"\n\n")
  # Unsuccesful convergence
  if (res$convergence > 0)
    stop("Unsuccesful convergence with the error code ", res$convergence, 
         ".\nAnother optimization method could succeed.\n")
  # Parameters covariance, correlation and standard error
  sd <- correl <- varcovar <- NULL
  # Maximum likelihood and maximum product of spacing
  if (method == "mle" || method == "mpse"){
    if (!is.null(res$hessian)){
      # check for NA values and invertible Hessian
      if (all(!is.na(res$hessian)) && qr(res$hessian)$rank == NCOL(res$hessian)) {
        # Parameters covariance, correlation and standard error
        varcovar <- solve(res$hessian)
        sd <- sqrt(diag(varcovar))
        correl <- cov2cor(varcovar)
      }
      else
        sd <- correl <- varcovar <- NA
    }
    else
      sd <- correl <- varcovar <- NA
  }
  # Object fitdist
  if (!keepdata){
    n2keep <- min(keepdata.nb, n) - 2
    imin <- which.min(data)
    imax <- which.max(data)
    subdata <- data[sample((1:n)[-c(imin, imax)], size = n2keep, replace = FALSE)]
    data <- c(subdata, data[c(imin, imax)])
  }
  npar <- length(res$estimate)
  aic <- -2 * res$loglik + 2 * npar
  bic <- -2 * res$loglik + log(n) * npar
  # Optimization method goes to dots
  if(is.null(custom.optim)){
    my3dots$optim.method <- optim.method
  }
  else
    my3dots$custom.optim <- custom.optim
  reslist <- list(estimate = res$estimate, method = method, sd = sd, cor = correl, 
                  vcov = varcovar, loglik = res$loglik, aic = aic, bic = bic, 
                  n = n, data = data, distname = "BMT", fix.arg = res$fix.arg, 
                  fix.arg.fun = res$fix.arg.fun, dots = my3dots, 
                  convergence = res$convergence, discrete = FALSE, weights = res$weights)
  return(structure(reslist, class = "fitdist"))
}

############################
##### Hidden functions #####

# Wrapper for function nlminb in order to work as a custom.optim
.m.nlminb <- function(fn , par, ...) {
  opt <- nlminb(start = par, objective = fn, ...)
  return(list(par = opt$par, convergence = opt$convergence, value = opt$objective,
              hessian = NULL, counts = as.vector(opt$evaluations), 
              message = opt$message))
}

# Compute log-likelihood
.loglik <- function(par, fix.arg, obs, ddistnam) {
  sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)))))
}
