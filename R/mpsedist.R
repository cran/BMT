#' @title Maximum Product of Spacing Fit of Univariate Distributions.
#' @description Fit of univariate distributions for non-censored data using 
#'   maximum product of spacing estimation (mpse), also called maximum spacing 
#'   estimation.
#' 
#' @name mpsedist
#'   
#' @details The \code{mpsedist} function carries out the maximum product of 
#'   spacing estimation numerically, by maximization of the arithmetic mean of
#'   \eqn{\log(F(k) - F(k-1))}.
#'   
#'   The optimization process is the same as 
#'   \code{\link[fitdistrplus]{mledist}}, see the 'details' section of that 
#'   function.
#'   
#'   Optionally, a vector of \code{weights} can be used in the fitting process. 
#'   By default (when \code{weights=NULL}), ordinary mpse is carried out, 
#'   otherwise the specified weights are used to compute a weighted arithmetic
#'   mean.
#'   
#'   We believe this function should be added to the package 
#'   \pkg{fitdistrplus}. Until it is accepted and incorporated into that
#'   package, it will remain in the package \pkg{BMT}. This function is 
#'   internally called in \code{\link{BMTfit.mpse}}.
#'   
#' @param data A numeric vector with the observed values for non-censored data.
#' @param distr A character string \code{"name"} naming a distribution for which
#'   the corresponding density function \code{dname} and the corresponding 
#'   distribution function \code{pname} must be classically defined.
#' @param start A named list giving the initial values of parameters of the 
#'   named distribution or a function of data computing initial values and 
#'   returning a named list. This argument may be omitted (default) for some 
#'   distributions for which reasonable starting values are computed (see the 
#'   'details' section of \code{\link[fitdistrplus]{mledist}}).
#' @param fix.arg An optional named list giving the values of fixed parameters 
#'   of the named distribution or a function of data computing (fixed) parameter
#'   values and returning a named list. Parameters with fixed value are thus NOT
#'   estimated.
#' @param optim.method \code{"default"} (see details) or an optimization method 
#'   to pass to \code{\link{optim}}.
#' @param lower Left bounds on the parameters for the \code{"L-BFGS-B"} method 
#'   (see \code{\link{optim}}) or the \code{\link{constrOptim}} function (as an 
#'   equivalent linear constraint).
#' @param upper Right bounds on the parameters for the \code{"L-BFGS-B"} method 
#'   (see \code{\link{optim}}) or the \code{\link{constrOptim}} function (as an 
#'   equivalent linear constraint).
#' @param custom.optim A function carrying the optimization (see details).
#' @param weights An optional vector of weights to be used in the fitting 
#'   process. Should be \code{NULL} or a numeric vector with strictly positive 
#'   numbers. If non-\code{NULL}, weighted mpse is used, otherwise ordinary 
#'   mpse.
#' @param silent A logical to remove or show warnings when bootstrapping.
#' @param gradient A function to return the gradient of the optimization
#'   objective function for the \code{"BFGS"}, \code{"CG"} and \code{"L-BFGS-B"}
#'   methods. If it is \code{NULL}, a finite-difference approximation will be
#'   used, see \code{\link{optim}}.
#' @param \dots Further arguments passed to the \code{\link{optim}}, 
#'   \code{\link{constrOptim}} or \code{custom.optim} function.
#'   
#' @return \code{mpsedist} returns a list with following components,
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
#'   \item{value}{ the value of the optimization objective function at the solution found. }
#'   
#'   \item{loglik}{ the log-likelihood. }
#'   
#'   \item{hessian}{ a symmetric matrix computed by \code{\link{optim}} as an
#'   estimate of the Hessian at the solution found or computed in the
#'   user-supplied optimization function. }
#'   
#'   \item{optim.function }{ the name of the optimization function used.  }
#'   
#'   \item{fix.arg}{ the named list giving the values of parameters of the named
#'   distribution that must kept fixed rather than estimated by maximum
#'   likelihood or \code{NULL} if there are no such parameters. }
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
#' @references Cheng, R. and N. Amin (1983). \emph{Estimating parameters in 
#'   continuous univariate distributions with a shifted origin}. Journal of the 
#'   Royal Statistical Society. Series B (Methodological), 394-403.
#'   
#'   Ranneby, B. (1984). \emph{The maximum spacing method. an estimation method 
#'   related to the maximum likelihood method}. Scandinavian Journal of 
#'   Statistics, 93-112.
#'   
#' @seealso \code{\link{mqdedist}}, \code{\link[fitdistrplus]{mledist}}, 
#'   \code{\link[fitdistrplus]{mmedist}}, \code{\link[fitdistrplus]{qmedist}}, 
#'   \code{\link[fitdistrplus]{mgedist}}, and \code{\link{optim}}.
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
#' mpse1 <- mpsedist(x1, "norm")
#' mpse1$estimate
#' 
#' # (2) defining your own distribution functions, here for the Gumbel 
#' # distribution for other distributions, see the CRAN task view dedicated 
#' # to probability distributions
#' dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
#' pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
#' qgumbel <- function(p, a, b) a-b*log(-log(p))
#' mpse1 <- mpsedist(x1, "gumbel", start = list(a = 10, b = 5))
#' mpse1$estimate
#' 
#' # (3) fit a discrete distribution (Poisson)
#' set.seed(1234)
#' x2 <- rpois(n = 30, lambda = 2)
#' mpse2 <- mpsedist(x2, "pois")
#' mpse2$estimate
#' 
#' # (4) fit a finite-support distribution (beta)
#' set.seed(1234)
#' x3 <- rbeta(n = 100, shape1 = 5, shape2 = 10)
#' mpse3 <- mpsedist(x3, "beta")
#' mpse3$estimate
#' 
#' # (5) fit frequency distributions on USArrests dataset.
#' x4 <- USArrests$Assault
#' mpse4pois <- mpsedist(x4, "pois")
#' mpse4pois$estimate
#' mpse4nbinom <- mpsedist(x4, "nbinom")
#' mpse4nbinom$estimate
#' 
#' # (6) weighted fit of a normal distribution 
#' set.seed(1234)
#' w1 <- runif(101)
#' mpse1 <- mpsedist(x1, "norm", weights = w1)
#' mpse1$estimate
#' 
#' @keywords distribution


###################
#' @rdname mpsedist
#' @export
mpsedist <- function (data, distr, start = NULL, fix.arg = NULL, optim.method = "default", 
                      lower = -Inf, upper = Inf, custom.optim = NULL, weights = NULL, 
                      silent = TRUE, gradient = NULL, ...) 
{
  if (!is.character(distr)) 
    stop("distr must be a character string naming a distribution")
  else distname <- distr
  ddistname <- paste("d", distname, sep = "")
  if (!exists(ddistname, mode = "function")) 
    stop(paste("The ", ddistname, " function must be defined"))
  pdistname <- paste("p", distname, sep = "")
  if (!exists(pdistname, mode = "function")) 
    stop(paste("The ", pdistname, " function must be defined"))
  if (is.null(custom.optim)) 
    optim.method <- match.arg(optim.method, c("default", 
                                              "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                                              "Brent"))
  start.arg <- start
  if (is.vector(start.arg)) 
    start.arg <- as.list(start.arg)
  txt1 <- "data must be a numeric vector of length greater than 1 for non censored data"
#  txt2 <- "or a dataframe with two columns named left and right and more than one line for censored data"
  if (!is.null(weights)) {
    if (any(weights <= 0)) 
      stop("weights should be a vector of numbers greater than 0")
    if (length(weights) != NROW(data) + 1) 
      stop("weights should be a vector with a length equal to the observation number")
    warning("weights are not taken into account in the default initial values")
  }
  if (is.vector(data)) {
    cens <- FALSE
    if (!(is.numeric(data) & length(data) > 1)) 
      stop(txt1)
  }
  else {
    stop("Maximum product of spacing estimation is not yet available for censored data.")
    #     cens <- TRUE
    #     censdata <- data
    #     if (!(is.vector(censdata$left) & is.vector(censdata$right) & 
    #           length(censdata[, 1]) > 1)) 
    #       stop(paste(txt1, txt2))
  }
  #   if (cens) {
  #     irow.lcens <- is.na(censdata$left)
  #     lcens <- censdata[irow.lcens, ]$right
  #     if (any(is.na(lcens))) 
  #       stop("An observation cannot be both right and left censored, coded with two NA values")
  #     irow.rcens <- is.na(censdata$right)
  #     rcens <- censdata[irow.rcens, ]$left
  #     irow.ncens <- censdata$left == censdata$right & !is.na(censdata$left) & 
  #       !is.na(censdata$right)
  #     ncens <- censdata[irow.ncens, ]$left
  #     irow.icens <- censdata$left != censdata$right & !is.na(censdata$left) & 
  #       !is.na(censdata$right)
  #     icens <- censdata[irow.icens, ]
  #     data <- c(rcens, lcens, ncens, (icens$left + icens$right)/2)
  #   }
  argpdistname <- names(formals(pdistname))
  chfixstt <- checkparam(start.arg = start.arg, fix.arg = fix.arg, 
                         argdistname = argpdistname, errtxt = NULL, 
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
  if (distname == "unif") {
    n <- length(data)
    data <- sort(data)
    par <- c(min = (n*data[1]-data[n])/(n-1), max = (n*data[n] - data[1])/(n-1))
    par <- c(par[!names(par) %in% names(fix.arg)], unlist(fix.arg))
    value <- unname(sum(log(diff(c(par["min"],data,par["max"])))) - (n+1)*log(par["max"]-par["min"]))
    res <- list(estimate = par[!names(par) %in% names(fix.arg)], convergence = 0, 
                value = value,
                loglik = .loglik(par[!names(par) %in% names(fix.arg)], fix.arg, data, ddistname), 
                hessian = NA, optim.function = NA, fix.arg = fix.arg)
    return(res)
  }
  if (!cens && is.null(weights)) {
    fnobj <- function(par, fix.arg, obs, pdistnam, ddistnam) {
      obs <- sort(obs)
      spacing <- diff(c(0, do.call(pdistnam, c(list(obs), as.list(par), as.list(fix.arg))), 1))
      if(any(is.nan(spacing)))
        return(NaN)
      ind <- abs(spacing) < .epsilon
      if(any(ind)){
        aux <- c(obs[1],obs)[ind]
        spacing[ind] <- do.call(ddistnam, c(list(aux), as.list(par), as.list(fix.arg)))
      }
      -sum(log(spacing))
    }
  }
  else if (!cens && !is.null(weights)) {
    fnobj <- function(par, fix.arg, obs, pdistnam, ddistnam) {
      obs <- sort(obs)
      spacing <- diff(c(0, do.call(pdistnam, c(list(obs), as.list(par), as.list(fix.arg))), 1))
      if(any(is.nan(spacing)))
        return(NaN)
      ind <- abs(spacing) < .epsilon
      if(any(ind)){
        aux <- c(obs[1],obs)[ind]
        spacing[ind] <- do.call(ddistnam, c(list(aux), as.list(par), as.list(fix.arg)))
      }
      -sum(weights * log(spacing))
    }
  }
  #   else if (cens && is.null(weights)) {
  #     argpdistname <- names(formals(pdistname))
  #     if (("log" %in% argddistname) & ("log.p" %in% argpdistname)) {
  #       fnobjcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam) {
  #         - sum(do.call(ddistnam, c(list(x = ncens), as.list(par), as.list(fix.arg), list(log = TRUE)))) - 
  #           sum(do.call(pdistnam, c(list(q = lcens), as.list(par), as.list(fix.arg), list(log = TRUE)))) - 
  #           sum(do.call(pdistnam, c(list(q = rcens), as.list(par), as.list(fix.arg), list(lower.tail = FALSE), list(log = TRUE)))) - 
  #           sum(log(do.call(pdistnam, c(list(q = icens$right), as.list(par), as.list(fix.arg))) - 
  #                     do.call(pdistnam, c(list(q = icens$left), as.list(par), as.list(fix.arg)))))
  #       }
  #     }
  #     else {
  #       fnobjcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam) {
  #         -sum(log(do.call(ddistnam, c(list(x = ncens), as.list(par), as.list(fix.arg))))) - 
  #           sum(log(do.call(pdistnam, c(list(q = lcens), as.list(par), as.list(fix.arg))))) - 
  #           sum(log(1 - do.call(pdistnam, c(list(q = rcens), as.list(par), as.list(fix.arg))))) - 
  #           sum(log(do.call(pdistnam, c(list(q = icens$right), as.list(par), as.list(fix.arg))) - 
  #                     do.call(pdistnam, c(list(q = icens$left), as.list(par), as.list(fix.arg)))))
  #       }
  #     }
  #   }
  #   else if (cens && !is.null(weights)) {
  #     fnobjcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam) {
  #       p1 <- log(do.call(ddistnam, c(list(x = ncens), as.list(par), as.list(fix.arg))))
  #       p2 <- log(do.call(pdistnam, c(list(q = lcens), as.list(par), as.list(fix.arg))))
  #       p3 <- log(1 - do.call(pdistnam, c(list(q = rcens), as.list(par), as.list(fix.arg))))
  #       p4 <- log(do.call(pdistnam, c(list(q = icens$right), as.list(par), as.list(fix.arg))) - 
  #                 do.call(pdistnam, c(list(q = icens$left), as.list(par), as.list(fix.arg)))) - 
  #             sum(weights[irow.ncens] * p1) - sum(weights[irow.lcens] * p2) - 
  #             sum(weights[irow.rcens] * p3) - sum(weights[irow.icens] * p4)
  #     }
  #   }
  owarn <- getOption("warn")
  if (is.null(custom.optim)) {
    hasbound <- any(is.finite(lower) | is.finite(upper))
    if (optim.method == "default") {
      meth <- ifelse(length(vstart) > 1, "Nelder-Mead", 
                     "BFGS")
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
        #         if (!cens) {
        opttryerror <- try(opt <- constrOptim(theta = vstart, 
                                              f = fnobj, ui = Mat, ci = Bnd, grad = gradient, 
                                              fix.arg = fix.arg, obs = data, 
                                              pdistnam = pdistname, ddistnam = ddistname, 
                                              hessian = !is.null(gradient), method = meth, 
                                              ...), silent = TRUE)
        #         }
        #         else opttryerror <- try(opt <- constrOptim(theta = vstart, 
        #                                                    f = fnobjcens, ui = Mat, ci = Bnd, grad = gradient, 
        #                                                    ddistnam = ddistname, rcens = rcens, lcens = lcens, 
        #                                                    icens = icens, ncens = ncens, pdistnam = pdistname, 
        #                                                    fix.arg = fix.arg, obs = data, hessian = !is.null(gradient), 
        #                                                    method = meth, ...), silent = TRUE)
        if (!inherits(opttryerror, "try-error")) 
          if (length(opt$counts) == 1) 
            opt$counts <- c(opt$counts, NA)
      }
      else {
        #         if (!cens) 
        opttryerror <- try(opt <- optim(par = vstart, 
                                        fn = fnobj, fix.arg = fix.arg, obs = data, 
                                        pdistnam = pdistname, ddistnam = ddistname, 
                                        gr = gradient, hessian = TRUE, 
                                        method = meth, lower = lower, upper = upper, 
                                        ...), silent = TRUE)
        #         else opttryerror <- try(opt <- optim(par = vstart, 
        #                                              fn = fnobjcens, fix.arg = fix.arg, gr = gradient, 
        #                                              rcens = rcens, lcens = lcens, icens = icens, 
        #                                              ncens = ncens, ddistnam = ddistname, pdistnam = pdistname, 
        #                                              hessian = TRUE, method = meth, lower = lower, 
        #                                              upper = upper, ...), silent = TRUE)
      }
    }
    else {
      opt.fun <- "optim"
      #       if (!cens) 
      opttryerror <- try(opt <- optim(par = vstart, 
                                      fn = fnobj, fix.arg = fix.arg, obs = data, 
                                      pdistnam = pdistname, ddistnam = ddistname,
                                      gr = gradient, hessian = TRUE, 
                                      method = meth, lower = lower, upper = upper, 
                                      ...), silent = TRUE)
      #       else opttryerror <- try(opt <- optim(par = vstart, 
      #                                            fn = fnobjcens, fix.arg = fix.arg, gr = gradient, 
      #                                            rcens = rcens, lcens = lcens, icens = icens, 
      #                                            ncens = ncens, ddistnam = ddistname, pdistnam = pdistname, 
      #                                            hessian = TRUE, method = meth, lower = lower, 
      #                                            upper = upper, ...), silent = TRUE)
    }
    options(warn = owarn)
    if (inherits(opttryerror, "try-error")) {
      warnings("The function optim encountered an error and stopped.")
      if (getOption("show.error.messages")) 
        print(attr(opttryerror, "condition"))
      return(list(estimate = rep(NA, length(vstart)), 
                  convergence = 100, value=NA, loglik = NA, hessian = NA, 
                  optim.function = opt.fun, fix.arg = fix.arg, 
                  optim.method = meth, fix.arg.fun = fix.arg.fun, 
                  counts = c(NA, NA)))
    }
    if (opt$convergence > 0) {
      warnings("The function optim failed to converge, with the error code ", 
               opt$convergence)
    }
    if (is.null(names(opt$par))) 
      names(opt$par) <- names(vstart)
    res <- list(estimate = opt$par, convergence = opt$convergence, value = -opt$value, 
                loglik = .loglik(opt$par, fix.arg, data, ddistname),
                hessian = opt$hessian, optim.function = opt.fun, 
                fix.arg = fix.arg, optim.method = meth, fix.arg.fun = fix.arg.fun, 
                weights = weights, counts = opt$counts, optim.message = opt$message)
  }
  else {
    options(warn = ifelse(silent, -1, 0))
    #     if (!cens) 
    opttryerror <- try(opt <- custom.optim(fn = fnobj, 
                                           fix.arg = fix.arg, obs = data, 
                                           pdistnam = pdistname, ddistnam = ddistname, 
                                           par = vstart, ...), silent = TRUE)
    #     else opttryerror <- try(opt <- custom.optim(fn = fnobjcens, 
    #                                                 fix.arg = fix.arg, rcens = rcens, lcens = lcens, 
    #                                                 icens = icens, ncens = ncens, ddistnam = ddistname, 
    #                                                 pdistnam = pdistname, par = vstart, ...), silent = TRUE)
    options(warn = owarn)
    if (inherits(opttryerror, "try-error")) {
      warnings("The customized optimization function encountered an error and stopped.")
      if (getOption("show.error.messages")) 
        print(attr(opttryerror, "condition"))
      return(list(estimate = rep(NA, length(vstart)), 
                  convergence = 100, value = NA, loglik = NA, hessian = NA, 
                  optim.function = custom.optim, fix.arg = fix.arg, 
                  fix.arg.fun = fix.arg.fun, counts = c(NA, NA)))
    }
    if (opt$convergence > 0) {
      warnings("The customized optimization function failed to converge, with the error code ", 
               opt$convergence)
    }
    argdot <- list(...)
    method.cust <- argdot[argdot == "method"]
    if (length(method.cust) == 0) {
      method.cust <- NULL
    }
    if (is.null(names(opt$par))) 
      names(opt$par) <- names(vstart)
    res <- list(estimate = opt$par, convergence = opt$convergence, value = -opt$value,
                loglik = .loglik(opt$par, fix.arg, data, ddistname),
                hessian = opt$hessian, optim.function = custom.optim, 
                fix.arg = fix.arg, method = method.cust, fix.arg.fun = fix.arg.fun, 
                weights = weights, counts = opt$counts, optim.message = opt$message)
  }
  return(res)
}
