#' @title The BMT-Psi Distribution.
#' @description Density, distribution function, quantile function, random number
#'   generation for the BMT-Psi distribution with mean equal to \code{mean} and 
#'   standard deviation equal to \code{sd}.
#' 
#' @name BMT.Psi
#' @aliases dBMT.Psi
#' @aliases pBMT.Psi
#' @aliases qBMT.Psi
#' @aliases rBMT.Psi
#'   
#' @details If \code{mean} or \code{sd} are not specified they assume the 
#'   default values of 0 and 1, respectively.
#'   
#'   The BMT-Psi distribution is the BMT distribution with \eqn{\kappa_l = 
#'   \kappa_r = 0.63355781127887611515}. The BMT-Psi cumulative distribution 
#'   function (cdf) is the closest BMT cdf to the logistic cdf with scale =
#'   1 / d and d = 1.70174439 (Camilli, 1994, p. 295).
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
#' @return \code{dBMT.Psi} gives the density, \code{pBMT.Psi} the distribution 
#'   function, \code{qBMT.Psi} the quantile function, and \code{rBMT.Psi} 
#'   generates random deviates.
#'   
#'   The length of the result is determined by \code{n} for \code{rBMT.Psi}, and
#'   is the maximum of the lengths of the numerical arguments for the other 
#'   functions.
#'   
#'   The numerical arguments other than \code{n} are recycled to the length of 
#'   the result. Only the first elements of the logical arguments are used.
#'   
#'   \code{sd <= 0} is an error and returns \code{NaN}.
#'   
#' @references Torres-Jimenez, C. J. (2018), \emph{The BMT Item Response Theory 
#'   model: A new skewed distribution family with bounded domain and an IRT 
#'   model based on it}, PhD thesis, Doctorado en ciencias - Estadistica, 
#'   Universidad Nacional de Colombia, Sede Bogota.
#'   
#'   Camilli, G. (1994). Teacher's corner: origin of the scaling constant d= 1.7
#'   in item response theory. Journal of Educational Statistics, 19(3), 293-295.
#'   
#' @seealso \link{Distributions} for other standard distributions. 
#'   \code{\link{pBMT}} for the BMT distribution and \code{\link{pBMT.Phi}} for 
#'   the BMT-Phi distribution.
#'   
#' @author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co}
#'   
#' @examples
#' 
#' layout(matrix(1:4, 2, 2))
#' 
#' curve(plogis(x, scale = 1 / 1.70174439), -4, 4, col = "red", lty = 2, ylab = "cdf")
#' curve(pBMT.Psi(x), add = TRUE, col = "blue", lty = 3)
#' legend("topleft", legend = c("logis(0, 1 / 1.70174439)","BMT-Psi(0,1)"), 
#'        bty = "n", col = c("red","blue"), lty = 2:3)
#' 
#' curve(plogis(x, scale = 1 / 1.70174439)-pBMT.Psi(x), -4, 4)
#' 
#' curve(qlogis(x, scale = 1 / 1.70174439), col = "red", lty = 2, xlab = "p", ylab = "qf")
#' curve(qBMT.Psi(x), add = TRUE, col = "blue", lty = 3)
#' 
#' hist(rBMT.Psi(10000), freq = FALSE, breaks = seq(-4, 4, 0.25), border = "blue")
#' curve(dlogis(x, scale = 1 / 1.70174439), add = TRUE, col = "red", lty = 2)
#' curve(dBMT.Psi(x), add = TRUE, col = "blue", lty = 3)

#' @rdname BMT.Psi
#' @export
dBMT.Psi <- function(x, mean = 0, sd = 1, log = FALSE) {
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
  range <- sd / .sd.BMT.Psi.DD01
  x <- (x - mean) / range + .mean.BMT.Psi.DD01
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
    # inv.x.t.Psi
    t <- .inv.x.t.Psi(x[ind])
    # yf.t.Psi
    y[ind] <- .yf.t.Psi(t) / range[ind]
  }
  # density values y are given as log(y)
  if (log)
    y <- log(y)
  return(y)
}

#' @rdname BMT.Psi
#' @export
pBMT.Psi <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE){
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
  range <- sd / .sd.BMT.Psi.DD01
  q <- (q - mean) / range + .mean.BMT.Psi.DD01
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
    # inv.x.t.Psi
    t <- .inv.x.t.Psi(q[ind])
    # yF.t.Psi
    p[ind] <- .yF.t.Psi(t)
  }
  # probabilities are \eqn{P[X > x]}
  if(!lower.tail)
    p <- 1 - p
  # probabilities p are given as log(p)
  if(log.p)
    p <- log(p)
  return(p)
}

#' @rdname BMT.Psi
#' @export
qBMT.Psi <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE){
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
    t <- .inv.yF.t.Psi(p[ind])
    # x.t
    q[ind] <- .x.t.Psi(t)
  }
  # Transform q
  range <- sd / .sd.BMT.Psi.DD01
  q <- (q - .mean.BMT.Psi.DD01) * range + mean
  return(q)
}

#' @rdname BMT.Psi
#' @export 
rBMT.Psi <- function(n, mean = 0, sd = 1){
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
  x <- qBMT.Psi(p, mean, sd)
  return(x)
}

## Global constants
.epsilon <- 1e-10
.zero <- 0 - .epsilon
.one <- 1 + .epsilon
# X \sim BMT.Psi(mean, sd) is the same X \sim BMT(.D.Psi, .D.Psi, "t w", mean, sd, "l-s")
# .D.Psi <- BMTfit.mge(qlogis(1:1e6/(1e6+1), scale=1/1.70174439),
#                  "KS", start=list(p4=0.5), fix.arg=list(p1=0, p2=1, p3=0), 
#                  type.p.3.4="a-s", type.p.1.2="l-s", custom.optim="nlminb")$estimate
.D.Psi <-             0.63355781127887611515
# Mean of X \sim BMT(.D.Psi, .D.Psi, "t w", 0, 1, "c-d")
# .mean.BMT.Psi.DD01 <- BMTmean(.D.Psi, .D.Psi)
.mean.BMT.Psi.DD01 <- 0.5
# Standard deviation of X \sim BMT(.D.Psi, .D.Psi, "t w", 0, 1, "c-d")
# .sd.BMT.Psi.DD01 <- BMTsd(.D.Psi, .D.Psi)
.sd.BMT.Psi.DD01 <-   0.16771818811837588270
# Maximum of X \sim BMT.Psi(0, 1)
# .d.BMT.Psi <- BMTchangepars(.D.Psi, .D.Psi, "t w", 0, 1, "l-s")$p2
.d.BMT.Psi <-         2.98119128050142556674
# Coefficients of polynomial x.t = a_3.Psi * t^3 + a_2.Psi * t^2 + a_1.Psi * t
# .a_3.Psi <- 6*.D.Psi - 2
.a_3.Psi <-           1.80134686767325646883
# .a_2.Psi <- -9*.D.Psi + 3
.a_2.Psi <-           -2.70202030150988470325
# .a_1.Psi <- 3*.D.Psi
.a_1.Psi <-           1.90067343383662823442
# For the real root of x.t = x
# .a.Psi <- .a_2.Psi/.a_3.Psi = -1.5
.a.Psi <-             -1.5
# .b.Psi <- .a_1.Psi/.a_3.Psi
.b.Psi <-             1.05514016647536013060
# .Q.Psi <- (.a.Psi*.a.Psi - 3*.b.Psi)/9
.Q.Psi <-             -0.10171338882512001578
# .auxR.Psi <- (2*.a.Psi*.a.Psi - 9*.b.Psi)*.a.Psi / 54
.auxR.Psi <-          0.13878504161884000490

# x.t.Psi
.x.t.Psi <- function(t){
  x <- ((.a_3.Psi*t + .a_2.Psi)*t + .a_1.Psi)*t
  return(x)
}

# Inverse function for x.t.Psi
.inv.x.t.Psi <- function(x){
  # Press W.H., Teukolsky S.A., Vetterling W.T. & Flannery B.P. 2007. 
  # Numerical recipes: The art of scientific computing 
  # Section 5.6: Quadratic and Cubic Equations. Page 228.
  len <- length(x)
  c <- -x/.a_3.Psi
  R <- .auxR.Psi + 0.5*c
  # One real root
  A <- - sign(R) * (abs(R) + sqrt(R*R-.Q.Psi*.Q.Psi*.Q.Psi))^(1/3)
  B <- rep(0, len)
  B[A!=0] <- .Q.Psi/A[A!=0]
  r <- (A + B) - .a.Psi/3
  r[r < 0] <- 0
  r[r > 1] <- 1
  return(r)
}

# yF.t.Psi (same yF.t)
.yF.t.Psi <- function(t){
  yF <- (-2*t + 3)*t*t
  return(yF)
}

# inv.yF.t.Psi (same inv.yF.t)
.inv.yF.t.Psi <- function(yF){
  t <- 0.5-cos((acos(2*yF-1)-2*pi)/3)
  return(t)
}

# yf.t.Psi
.yf.t.Psi <- function(t){
  yf <- ((-6*t + 6)*t) / ((3*.a_3.Psi*t + 2*.a_2.Psi)*t + .a_1.Psi)
  return(yf)
}