#' @title The BMT-Phi Distribution.
#' @description Density, distribution function, quantile function, random number
#'   generation for the BMT-Phi distribution with mean equal to \code{mean} and 
#'   standard deviation equal to \code{sd}.
#' 
#' @name BMT.Phi
#' @aliases dBMT.Phi
#' @aliases pBMT.Phi
#' @aliases qBMT.Phi
#' @aliases rBMT.Phi
#'   
#' @details If \code{mean} or \code{sd} are not specified they assume the 
#'   default values of 0 and 1, respectively.
#'   
#'   The BMT-Phi distribution is the BMT distribution with \eqn{\kappa_l = 
#'   \kappa_r = 0.58029164978583758}. The BMT-Phi cumulative distribution
#'   function (cdf) is the closest BMT cdf to the normal cdf with the same mean and standard deviation.
#'   
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'   to be the number required
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le 
#'   x]}, otherwise, \eqn{P[X > x]}.
#'   
#' @return \code{dBMT.Phi} gives the density, \code{pBMT.Phi} the distribution 
#'   function, \code{qBMT.Phi} the quantile function, and \code{rBMT.Phi} 
#'   generates random deviates.
#'   
#'   The length of the result is determined by \code{n} for \code{rBMT.Phi}, and
#'   is the maximum of the lengths of the numerical arguments for the other 
#'   functions.
#'   
#'   The numerical arguments other than \code{n} are recycled to the length of 
#'   the result. Only the first elements of the logical arguments are used.
#'   
#'   \code{sd <= 0} is an error and returns \code{NaN}.
#'   
#' @references Torres-Jimenez, C. J. (2018), \emph{The BMT Item Response Theory model: A 
#'   new skewed distribution family with bounded domain and an IRT model based 
#'   on it}, PhD thesis, Doctorado en ciencias - Estadistica, Universidad 
#'   Nacional de Colombia, Sede Bogota.
#'   
#' @seealso \link{Distributions} for other standard distributions.
#'   \code{\link{pBMT}} for the BMT distribution and \code{\link{pBMT.Psi}} for 
#'   the BMT-Psi distribution.
#'   
#' @author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co}
#'   
#' @examples
#' 
#' layout(matrix(1:4,2,2))
#' 
#' curve(pnorm(x), -4, 4, col = "red", lty = 2, ylab = "cdf")
#' curve(pBMT.Phi(x), add = TRUE, col = "blue", lty = 3)
#' legend("topleft", legend = c("norm(0,1)","BMT-Phi(0,1)"), 
#'        bty = "n", col = c("red","blue"), lty = 2:3)
#' 
#' curve(pnorm(x)-pBMT.Phi(x), -4, 4)
#' 
#' curve(qnorm(x), col = "red", lty = 2, xlab = "p", ylab = "qf")
#' curve(qBMT.Phi(x), add = TRUE, col = "blue", lty = 3)
#' 
#' hist(rBMT.Phi(10000), freq = FALSE, breaks = seq(-4,4,0.25), border = "blue")
#' curve(dnorm(x), add = TRUE, col = "red", lty = 2)
#' curve(dBMT.Phi(x), add = TRUE, col = "blue", lty = 3)

#' @rdname BMT.Phi
#' @export
dBMT.Phi <- function(x, mean = 0, sd = 1, log = FALSE) {
  # The length of the result is determined by the maximum of the lengths of the
  # numerical arguments. The numerical arguments are recycled to the length of
  # the result.
  len <- max(length(x), length(mean), length(sd))
  x <- rep(x, len = len)
  mean <- rep(mean, len = len)
  sd <- rep(sd, len = len)
  # Control location-scale parameters
  sd <- replace(sd, sd <= 0, NaN)
  # Transform x to 0,1 given location-scale parameters
  range <- sd / .sd.BMT.DD01
  x <- (x - mean) / range + .mean.BMT.DD01
  # NaNs
  y <- x + mean + sd
  if (any(is.nan(y))) {
    warning("NaNs founded or produced")
  }
  ind <- !is.na(y)
  # For transformed x outside 0,1
  y[ind] <- 0
  # For transformed x inside 0,1
  ind[ind] <- x[ind] > 0 & x[ind] < 1
  if (any(ind)) {
    # inv.x.t.Phi
    t <- .inv.x.t.Phi(x[ind])
    # yf.t.Phi
    y[ind] <- .yf.t.Phi(t) / range[ind]
  }
  # density values y are given as log(y)
  if (log)
    y <- log(y)
  return(y)
}

#' @rdname BMT.Phi
#' @export 
pBMT.Phi <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE){
  # The length of the result is determined by the maximum of the lengths of the 
  # numerical arguments. The numerical arguments are recycled to the length of 
  # the result.
  len <- max(length(q), length(mean), length(sd))
  q <- rep(q, len=len)
  mean <- rep(mean, len=len)
  sd <- rep(sd, len=len)
  # Control location-scale parameters
  sd <- replace(sd, sd <= 0, NaN)
  # Transform q to 0,1 given location-scale parameters
  range <- sd / .sd.BMT.DD01
  q <- (q - mean) / range + .mean.BMT.DD01
  # NaNs
  p <- q + mean + sd
  if(any(is.nan(p)))
    warning("NaNs founded or produced")
  ind <- !is.na(p)
  # For transformed q outside 0,1
  p[ind] <- 0
  ind2 <- ind
  ind2[ind] <- q[ind] >= 1
  p[ind2] <- 1
  # For transformed q inside 0,1
  ind[ind] <- q[ind] > 0 & q[ind] < 1
  if(any(ind)){
    # inv.x.t.Phi
    t <- .inv.x.t.Phi(q[ind])
    # yF.t.Phi
    p[ind] <- .yF.t.Phi(t)
  }
  # probabilities are \eqn{P[X > x]}
  if(!lower.tail)
    p <- 1 - p
  # probabilities p are given as log(p)
  if(log.p)
    p <- log(p)
  return(p)
}

#' @rdname BMT.Phi
#' @export 
qBMT.Phi <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE){
  # probabilities p are given as log(p)
  if(log.p)
    p <- exp(p)
  # probabilities are \eqn{P[X > x]}
  if(!lower.tail)
    p <- 1 - p
  # The length of the result is determined by the maximum of the lengths of the 
  # numerical arguments. The numerical arguments are recycled to the length of 
  # the result.
  len <- max(length(p),length(mean),length(sd))
  p <- rep(p, len=len)
  mean <- rep(mean, len=len)
  sd <- rep(sd, len=len)
  # Control location-scale parameters
  sd <- replace(sd, sd <= 0, NaN)
  # NaNs
  q <- p + mean + sd
  if(any(is.nan(q)))
    warning("NaNs founded or produced")
  ind <- !is.na(q)
  # For p outside (0,1)
  q <- rep(NaN,len)
  ind2 <- ind
  ind2[ind] <- p[ind] == 0
  q[ind2] <- 0
  ind2 <- ind
  ind2[ind] <- p[ind] == 1
  q[ind2] <- 1
  # For p inside (0,1)
  ind[ind] <- p[ind] > 0 & p[ind] < 1
  if(any(ind)){
    # inv.yF.t
    t <- .inv.yF.t.Phi(p[ind])
    # x.t
    q[ind] <- .x.t.Phi(t)
  }
  # Transform q
  range <- sd / .sd.BMT.DD01
  q <- (q - .mean.BMT.DD01) * range + mean
  return(q)
}

#' @rdname BMT.Phi
#' @export 
rBMT.Phi <- function(n, mean = 0, sd = 1){
  # 
  len <- length(n)
  if(len > 1)
    n <- len
  else
    n <- trunc(n)
    if(n < 1)
      stop("invalid arguments")
  # Method of inversion
  p <- runif(n)
  x <- qBMT.Phi(p, mean, sd)
  return(x)
}

## Global constants
.epsilon <- 1e-10
.zero <- 0-.epsilon
.one <- 1+.epsilon
# X \sim BMT.Phi(mean, sd) is the same X \sim BMT(.D, .D, "t w", mean, sd, "l-s")
# .D <- BMTfit.mge(qnorm(1:1e6/(1e6+1)), "KS", start=list(p4=0.5), fix.arg=list(p1=0, p2=1, p3=0), 
#                  type.p.3.4="a-s", type.p.1.2="l-s", custom.optim="nlminb")$estimate
.D <- 0.58029164978583758
# Mean of X \sim BMT(.D, .D, "t w", 0, 1, "c-d")
# .mean.BMT.DD01 <- BMTmean(.D, .D)
.mean.BMT.DD01 <- 0.5
# Standard deviation of X \sim BMT(.D, .D, "t w", 0, 1, "c-d")
# .sd.BMT.DD01 <- BMTsd(.D, .D)
.sd.BMT.DD01 <- 0.17733001242558785
# Maximum of X \sim BMT.Phi(0, 1)
# .d.BMT.Phi <- BMTchangepars(.D, .D, "t w", 0, 1, "l-s")$p2
.d.BMT.Phi <- 2.8196016746449653
# Coefficients of polynomial x.t = a_3 * t^3 + a_2 * t^2 + a_1 * t
# .a_3 <- 6*.D - 2
.a_3 <- 1.4817498987150257
# .a_2 <- -9*.D + 3
.a_2 <- -2.2226248480725381
# .a_1 <- 3*.D
.a_1 <- 1.7408749493575129
# For the real root of x.t = x
# .a <- .a_2/.a_3
.a <- -1.5
# .b <- .a_1/.a_3
.b <- 1.1748777245520317
# .Q <- (.a*.a - 3*.b)/9
.Q <- -0.14162590818401061
# .auxR <- (2*.a*.a - 9*.b)*.a / 54
.auxR <- 0.16871943113800789

# x.t.Phi
.x.t.Phi <- function(t){
  x <- ((.a_3*t + .a_2)*t + .a_1)*t
  return(x)
}

# Inverse function for x.t.Phi
.inv.x.t.Phi <- function(x){
  # Press W.H., Teukolsky S.A., Vetterling W.T. & Flannery B.P. 2007. 
  # Numerical recipes: The art of scientific computing 
  # Section 5.6: Quadratic and Cubic Equations. Page 228.
  len <- length(x)
  c <- -x/.a_3
  R <- .auxR + 0.5*c
  # One real root
  A <- - sign(R) * (abs(R) + sqrt(R*R-.Q*.Q*.Q))^(1/3)
  B <- rep(0, len)
  B[A!=0] <- .Q/A[A!=0]
  r <- (A + B) - .a/3
  r[r < 0] <- 0
  r[r > 1] <- 1
  return(r)
}

# yF.t.Phi (same yF.t)
.yF.t.Phi <- function(t){
  yF <- (-2*t + 3)*t*t
  return(yF)
}

# inv.yF.t.Phi (same inv.yF.t)
.inv.yF.t.Phi <- function(yF){
  t <- 0.5-cos((acos(2*yF-1)-2*pi)/3)
  return(t)
}

# yf.t.Phi
.yf.t.Phi <- function(t){
  yf <- ((-6*t + 6)*t) / ((3*.a_3*t + 2*.a_2)*t + .a_1)
  return(yf)
}