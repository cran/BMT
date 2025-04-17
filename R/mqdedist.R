#' @title Minimum Quantile Distance Fit of Univariate Distributions.
#' @description Fit of univariate distributions for non-censored data using 
#'   minimum quantile distance estimation (mqde), which can also be called 
#'   maximum quantile goodness-of-fit estimation.
#' 
#' @name mqdedist
#'   
#' @details The \code{mqdedist} function carries out the minimum quantile 
#'   distance estimation numerically, by minimization of a distance between 
#'   observed and theoretical quantiles.
#'   
#'   The optimization process is the same as 
#'   \code{\link[fitdistrplus]{mledist}}, see the 'details' section of that 
#'   function.
#'   
#'   Optionally, a vector of \code{weights} can be used in the fitting process. 
#'   By default (when \code{weights=NULL}), ordinary mqde is carried out, 
#'   otherwise the specified weights are used to compute a weighted distance.
#'   
#'   We believe this function should be added to the package 
#'   \pkg{fitdistrplus}. Until it is accepted and incorporated into that
#'   package, it will remain in the package \pkg{BMT}. This function is 
#'   internally called in \code{\link{BMTfit.mqde}}.
#'   
#' @param data A numeric vector with the observed values for non-censored data.
#' @param distr A character string \code{"name"} naming a distribution for which
#'   the corresponding quantile function \code{qname} and the corresponding 
#'   density distribution \code{dname} must be classically defined.
#' @param probs A numeric vector of the probabilities for which the minimum 
#'   quantile distance estimation is done. \eqn{p[k] = (k - 0.5) / n} (default).
#' @param qtype The quantile type used by the R \code{\link{quantile}} function 
#'   to compute the empirical quantiles. Type 5 (default), i.e. \eqn{x[k]} is 
#'   both the \eqn{k}th order statistic and the type 5 sample quantile of 
#'   \eqn{p[k] = (k - 0.5) / n}.
#' @param dist The distance measure between observed and theoretical quantiles 
#'   to be used. This must be one of "euclidean" (default), "maximum", or 
#'   "manhattan". Any unambiguous substring can be given.
#' @param start A named list giving the initial values of parameters of the 
#'   named distribution or a function of data computing initial values and 
#'   returning a named list. This argument may be omitted (default) for some 
#'   distributions for which reasonable starting values are computed (see the 
#'   'details' section of \code{\link[fitdistrplus]{mledist}}).
#' @param fix.arg An optional named list giving the values of fixed parameters 
#'   of the named distribution or a function of data computing (fixed) parameter
#'   values and returning a named list. Parameters with fixed value are thus NOT
#'   estimated.
#' @param optim.method \code{"default"} (see details) or optimization method to pass to 
#'   \code{\link{optim}}.
#' @param lower Left bounds on the parameters for the \code{"L-BFGS-B"} method 
#'   (see \code{\link{optim}}) or the \code{\link{constrOptim}} function (as an 
#'   equivalent linear constraint).
#' @param upper Right bounds on the parameters for the \code{"L-BFGS-B"} method 
#'   (see \code{\link{optim}}) or the \code{\link{constrOptim}} function (as an 
#'   equivalent linear constraint).
#' @param custom.optim A function carrying the optimization (see details).
#' @param weights An optional vector of weights to be used in the fitting 
#'   process. Should be \code{NULL} or a numeric vector with strictly positive 
#'   numbers. If non-\code{NULL}, weighted mqde is used, otherwise ordinary 
#'   mqde.
#' @param silent A logical to remove or show warnings when bootstrapping.
#' @param gradient A function to return the gradient of the optimization 
#'   objective function for the \code{"BFGS"}, \code{"CG"} and \code{"L-BFGS-B"}
#'   methods. If it is \code{NULL}, a finite-difference approximation will be 
#'   used, see \code{\link{optim}}.
#' @param \dots Further arguments passed to the \code{\link{optim}}, 
#'   \code{\link{constrOptim}} or \code{custom.optim} function.
#'   
#' @return \code{mqdedist} returns a list with following components,
#'   
#'   \item{estimate}{ the parameter estimates.}
#'   
#'   \item{convergence}{ an integer code for the convergence of 
#'   \code{\link{optim}} defined as below or defined by the user in the 
#'   user-supplied optimization function.
#'   
#'   \code{0} indicates successful convergence.
#'   
#'   \code{1} indicates that the iteration limit of \code{\link{optim}} has been
#'   reached.
#'   
#'   \code{10} indicates degeneracy of the Nealder-Mead simplex.
#'   
#'   \code{100} indicates that \code{\link{optim}} encountered an internal 
#'   error. }
#'   
#'   \item{value}{the value of the optimization objective function at the solution found.}
#'   
#'   \item{hessian}{ a symmetric matrix computed by \code{\link{optim}} as an
#'   estimate of the Hessian at the solution found or computed in the
#'   user-supplied optimization function. }
#'   
#'   \item{probs}{ the probability vector on which observed and theoretical quantiles were calculated. }
#'   
#'   \item{dist}{ the name of the distance between observed and theoretical quantiles used. }
#'   
#'   \item{optim.function }{ the name of the optimization function used.  }
#'   
#'   \item{fix.arg}{ the named list giving the values of parameters of the named
#'   distribution that must kept fixed rather than estimated by maximum
#'   likelihood or \code{NULL} if there are no such parameters. }
#'   
#'   \item{loglik}{ the log-likelihood. }
#'   
#'   \item{optim.method}{when \code{\link{optim}} is used, the name of the
#'   algorithm used, \code{NULL} otherwise.}
#'   
#'   \item{fix.arg.fun}{the function used to set the value of \code{fix.arg} or
#'   \code{NULL}.}
#'   
#'   \item{weights}{the vector of weights used in the estimation process or 
#'   \code{NULL}.}
#'   
#'   \item{counts}{A two-element integer vector giving the number of calls to
#'   the log-likelihood function and its gradient respectively. This excludes
#'   those calls needed to compute the Hessian, if requested, and any calls to
#'   log-likelihood function to compute a finite-difference approximation to the
#'   gradient. \code{counts} is returned by \code{\link{optim}} or the
#'   user-supplied optimization function, or set to \code{NULL}.}
#'   
#'   \item{optim.message}{A character string giving any additional information 
#'   returned by the optimizer, or \code{NULL}. To understand exactly the 
#'   message, see the source code.}
#'   
#' @references LaRiccia, V. N. (1982). \emph{Asymptotic Properties of Weighted 
#'   $L^2$ Quantile Distance Estimators}. The Annals of Statistics, 10(2), 
#'   621-624.
#'   
#'   Torres-Jimenez, C. J. (2017, September), \emph{Comparison of estimation methods 
#'   for the BMT distribution}. ArXiv e-prints.
#'   
#' @seealso \code{\link{mpsedist}}, \code{\link[fitdistrplus]{mledist}}, 
#'   \code{\link[fitdistrplus]{mmedist}}, \code{\link[fitdistrplus]{qmedist}}, 
#'   \code{\link[fitdistrplus]{mgedist}}, \code{\link{optim}}, 
#'   \code{\link{constrOptim}}, and \code{\link{quantile}}.
#'   
#' @author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co}
#'   
#' @source Based on the function mledist of the R package: \pkg{fitdistrplus}
#'   
#'   Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package 
#'   for Fitting Distributions}. Journal of Statistical Software, 64(4), 1-34.
#'   
#'   Functions \code{checkparam} and \code{startargdefault} are needed and 
#'   were copied from the same package (fitdistrplus version: 1.0-9).
#'   
#' @examples
#' # (1) basic fit of a normal distribution 
#' set.seed(1234)
#' x1 <- rnorm(n = 100)
#' mean(x1); sd(x1)
#' mqde1 <- mqdedist(x1, "norm")
#' mqde1$estimate
#' 
#' # (2) defining your own distribution functions, here for the Gumbel 
#' # distribution for other distributions, see the CRAN task view dedicated 
#' # to probability distributions
#' dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
#' pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
#' qgumbel <- function(p, a, b) a-b*log(-log(p))
#' mqde1 <- mqdedist(x1, "gumbel", start = list(a = 10, b = 5))
#' mqde1$estimate
#' 
#' # (3) fit a discrete distribution (Poisson)
#' set.seed(1234)
#' x2 <- rpois(n = 30, lambda = 2)
#' mqde2 <- mqdedist(x2, "pois")
#' mqde2$estimate
#' 
#' # (4) fit a finite-support distribution (beta)
#' set.seed(1234)
#' x3 <- rbeta(n = 100, shape1 = 5, shape2 = 10)
#' mqde3 <- mqdedist(x3, "beta")
#' mqde3$estimate
#' 
#' # (5) fit frequency distributions on USArrests dataset.
#' x4 <- USArrests$Assault
#' mqde4pois <- mqdedist(x4, "pois")
#' mqde4pois$estimate
#' mqde4nbinom <- mqdedist(x4, "nbinom")
#' mqde4nbinom$estimate
#' 
#' # (6) weighted fit of a normal distribution 
#' set.seed(1234)
#' w1 <- runif(100)
#' weighted.mean(x1, w1)
#' mqde1 <- mqdedist(x1, "norm", weights = w1)
#' mqde1$estimate
#' 
#' @keywords distribution


###################
#' @rdname mqdedist
#' @export
mqdedist <- function (data, distr, probs = (1:length(data)-0.5)/length(data), qtype = 5, 
                      dist = "euclidean", start = NULL, fix.arg = NULL, optim.method = "default", 
                      lower = -Inf, upper = Inf, custom.optim = NULL, weights = NULL, 
                      silent = TRUE, gradient = NULL, ...) 
{
  if (!is.character(distr)) 
    stop("distr must be a character string naming a distribution")
  else distname <- distr
  qdistname <- paste("q", distname, sep = "")
  if (!exists(qdistname, mode = "function")) 
    stop(paste("The ", qdistname, " function must be defined"))
  ddistname <- paste("d", distname, sep = "")
  if (!exists(ddistname, mode = "function")) 
    stop(paste("The ", ddistname, " function must be defined"))
  if (!(is.vector(probs) & is.numeric(probs)) | anyNA(probs) | any(probs <= 0 | probs >= 1))
    stop("probs must be a numeric vector with all elements greater than zero and less than one")
  probs <- unique(sort(probs))
  if (qtype < 0 || qtype > 9) 
    stop("wrong type for the R quantile function")
  if (is.null(custom.optim)) 
    optim.method <- match.arg(optim.method, c("default", "Nelder-Mead", "BFGS", "CG", 
                                              "L-BFGS-B", "SANN", "Brent"))
  int.dist <- pmatch(dist, c("euclidean", "maximum", "manhattan"))
  if (is.na(int.dist))
    stop("invalid distance measure to be used")
  if (int.dist == -1)
    stop("ambiguous distance measure to be used")
  start.arg <- start
  if (is.vector(start.arg)) 
    start.arg <- as.list(start.arg)
  my3dots <- list(...)
  if (!is.null(weights)) {
    if (any(weights <= 0)) 
      stop("weights should be a vector of numbers greater than 0")
    if (length(weights) != NROW(probs)) 
      stop("weights should be a vector with a length equal to the the probabilities probs")
    warning("weights are not taken into account in the default initial values")
  }
  if (is.vector(data)) {
    cens <- FALSE
    if (!(is.numeric(data) & length(data) > 1)) 
      stop("data must be a numeric vector of length greater than 1 for non-censored data")
#            \n            or a dataframe with two columns named left and right and more than one line for censored data")
  }
  else {
    stop("Minimum quantile distance estimation is not yet available for censored data.")
#     cens <- TRUE
#     censdata <- data
#     if (!(is.vector(censdata$left) & is.vector(censdata$right) & 
#           length(censdata[, 1]) > 1)) 
#       stop("data must be a numeric vector of length greater than 1 for non censored data\n        or a dataframe with two columns named left and right and more than one line for censored data")
#     pdistname <- paste("p", distname, sep = "")
#     if (!exists(pdistname, mode = "function")) 
#       stop(paste("The ", pdistname, " function must be defined to apply maximum likelihood to censored data"))
  }
#   if (cens) {
#     lcens <- censdata[is.na(censdata$left), ]$right
#     if (any(is.na(lcens))) 
#       stop("An observation cannot be both right and left censored, coded with two NA values")
#     rcens <- censdata[is.na(censdata$right), ]$left
#     ncens <- censdata[censdata$left == censdata$right & 
#                         !is.na(censdata$left) & !is.na(censdata$right), 
#                       ]$left
#     icens <- censdata[censdata$left != censdata$right & 
#                         !is.na(censdata$left) & !is.na(censdata$right), 
#                       ]
#     data <- c(rcens, lcens, ncens, (icens$left + icens$right)/2)
#   }
  argqdistname <- names(formals(qdistname))
  chfixstt <- checkparam(start.arg = start.arg, fix.arg = fix.arg, 
                         argdistname = argqdistname, errtxt = NULL, 
                         data10 = head(data, 10), distname = distname)
  if (!chfixstt$ok) 
    stop(chfixstt$txt)
  if (is.function(chfixstt$start.arg)) 
    vstart <- unlist(chfixstt$start.arg(data))
  else vstart <- unlist(chfixstt$start.arg)
  if (is.function(fix.arg)) {
    fix.arg.fun <- fix.arg
    fix.arg <- fix.arg(data)
  }
  else fix.arg.fun <- NULL
  if(qtype == 0)
    s <- data
  else
    s <- quantile(data, probs = probs, type = qtype, names = FALSE)
#   if (!cens) {
  if(is.null(weights)){
    if (int.dist == 1) # euclidean
      fnobj <- function(par, fix.arg, obs, qdistnam, probs, weights) {
        theoq <- do.call(qdistnam, c(list(p = probs), as.list(par), as.list(fix.arg)))
        sum((obs - theoq)^2)
      }
    else if (int.dist == 2) # maximum
      fnobj <- function(par, fix.arg, obs, qdistnam, probs, weights) {
        theoq <- do.call(qdistnam, c(list(p = probs), as.list(par), as.list(fix.arg)))
        max(abs(obs - theoq))
      }
    else if (int.dist == 3) # manhattan
      fnobj <- function(par, fix.arg, obs, qdistnam, probs, weights) {
        theoq <- do.call(qdistnam, c(list(p = probs), as.list(par), as.list(fix.arg)))
        sum(abs(obs - theoq))
      }
  }
  else{
    if (int.dist == 1) # euclidean
      fnobj <- function(par, fix.arg, obs, qdistnam, probs, weights) {
        theoq <- do.call(qdistnam, c(list(p = probs), as.list(par), as.list(fix.arg)))
        sum(weights*(obs - theoq)^2)
      }
    else if (int.dist == 2) # maximum
      fnobj <- function(par, fix.arg, obs, qdistnam, probs, weights) {
        theoq <- do.call(qdistnam, c(list(p = probs), as.list(par), as.list(fix.arg)))
        max(weights*abs(obs - theoq))
      }
    else if (int.dist == 3) # manhattan
      fnobj <- function(par, fix.arg, obs, qdistnam, probs, weights) {
        theoq <- do.call(qdistnam, c(list(p = probs), as.list(par), as.list(fix.arg)))
        sum(weights*abs(obs - theoq))
      }
  }
#   }
#   else stop("Minimum quantile distance estimation is not yet available for censored data.")
  owarn <- getOption("warn")
  if (is.null(custom.optim)) {
    hasbound <- any(is.finite(lower) | is.finite(upper))
    if (optim.method == "default") {
      meth <- ifelse(length(vstart) > 1, "Nelder-Mead", "BFGS")
    }
    else meth <- optim.method
    if (meth == "BFGS" && hasbound && is.null(gradient)) {
      meth <- "L-BFGS-B"
      txt1 <- "The BFGS method cannot be used with bounds without provided the gradient."
      txt2 <- "The method is changed to L-BFGS-B."
      warning(paste(txt1, txt2))
    }
    options(warn = ifelse(silent, -1, 0))
    if (hasbound) {
      if (!is.null(gradient)) {
        opt.fun <- "constrOptim"
      }
      else {
        if (meth == "Nelder-Mead") 
          opt.fun <- "constrOptim"
        else if (meth %in% c("L-BFGS-B", "Brent")) 
          opt.fun <- "optim"
        else {
          txt1 <- paste("The method", meth, "cannot be used by constrOptim() nor optim() without gradient and bounds.")
          txt2 <- "Only optimization methods L-BFGS-B, Brent and Nelder-Mead can be used in such case."
          stop(paste(txt1, txt2))
        }
      }
      if (opt.fun == "constrOptim") {
        npar <- length(vstart)
        lower <- as.double(rep_len(lower, npar))
        upper <- as.double(rep_len(upper, npar))
        haslow <- is.finite(lower)
        Mat <- diag(npar)[haslow, ]
        hasupp <- is.finite(upper)
        Mat <- rbind(Mat, -diag(npar)[hasupp, ])
        colnames(Mat) <- names(vstart)
        rownames(Mat) <- paste0("constr", 1:NROW(Mat))
        Bnd <- c(lower[is.finite(lower)], -upper[is.finite(upper)])
        names(Bnd) <- paste0("constr", 1:length(Bnd))
        initconstr <- Mat %*% vstart - Bnd
        if (any(initconstr < 0)) 
          stop("Starting values must be in the feasible region.")
        opttryerror <- try(opt <- constrOptim(theta = vstart, 
                                              f = fnobj, ui = Mat, ci = Bnd, grad = gradient, 
                                              fix.arg = fix.arg, obs = s, qdistnam = qdistname, 
                                              probs = probs, weights = weights,
                                              hessian = !is.null(gradient), method = meth, 
                                              ...), silent = TRUE)
        if (!inherits(opttryerror, "try-error")) 
          if (length(opt$counts) == 1) 
            opt$counts <- c(opt$counts, NA)
      }
      else {
        opttryerror <- try(opt <- optim(par = vstart, 
                                        fn = fnobj, fix.arg = fix.arg, obs = s, 
                                        gr = gradient, qdistnam = qdistname, 
                                        probs = probs, weights = weights, hessian = TRUE, 
                                        method = meth, lower = lower, upper = upper, 
                                        ...), silent = TRUE)
      }
    }
    else {
      opt.fun <- "optim"
      opttryerror <- try(opt <- optim(par = vstart, fn = fnobj, 
                                      fix.arg = fix.arg, obs = s, gr = gradient, 
                                      qdistnam = qdistname, probs = probs, weights = weights, 
                                      hessian = TRUE, method = meth, 
                                      lower = lower, upper = upper, ...), silent = TRUE)
    }
    options(warn = owarn)
    if (inherits(opttryerror, "try-error")) {
      warnings("The function optim encountered an error and stopped.")
      if (getOption("show.error.messages")) 
        print(attr(opttryerror, "condition"))
      return(list(estimate = rep(NA, length(vstart)), 
                  convergence = 100, loglik = NA, hessian = NA))
    }
    if (opt$convergence > 0) {
      warnings("The function optim failed to converge, with the error code ", 
               opt$convergence)
    }
    if (is.null(names(opt$par))) 
      names(opt$par) <- names(vstart)
    res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                hessian = opt$hessian, probs = probs, dist = dist, 
                optim.function = opt.fun, fix.arg = fix.arg, 
                loglik = .loglik(opt$par, fix.arg, data, ddistname), 
                optim.method = meth, fix.arg.fun = fix.arg.fun, weights = weights, 
                counts = opt$counts, optim.message = opt$message)
  }
  else {
    options(warn = ifelse(silent, -1, 0))
#     if (!cens) 
      opttryerror <- try(opt <- custom.optim(fn = fnobj, 
                                             fix.arg = fix.arg, obs = s, qdistnam = qdistname, 
                                             probs = probs, weights = weights, 
                                             par = vstart, ...), silent = TRUE)
#     else stop("Maximum goodness-of-fit estimation is not yet available for censored data.")
    options(warn = owarn)
    if (inherits(opttryerror, "try-error")) {
      warnings("The customized optimization function encountered an error and stopped.")
      if (getOption("show.error.messages")) 
        print(attr(opttryerror, "condition"))
      return(list(estimate = rep(NA, length(vstart)), 
                  convergence = 100, value = NA, hessian = NA))
    }
    if (opt$convergence > 0) {
      warnings("The customized optimization function failed to converge, with the error code ", 
               opt$convergence)
    }
    if (is.null(names(opt$par))) 
      names(opt$par) <- names(vstart)
    res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                probs = probs, dist = dist, hessian = opt$hessian, 
                optim.function = custom.optim, fix.arg = fix.arg, 
                loglik = .loglik(opt$par, fix.arg, data, ddistname), 
                optim.method = NULL, fix.arg.fun = fix.arg.fun, weights = weights, 
                counts = opt$counts, optim.message = opt$message)
  }
  return(res)
}
