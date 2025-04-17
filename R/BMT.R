#' @title The BMT Distribution.
#' 
#' @description Density, distribution, quantile function, random number 
#'   generation for the BMT distribution, with \code{p3} and \code{p4} tails 
#'   weights (\eqn{\kappa_l} and \eqn{\kappa_r}) or asymmetry-steepness 
#'   parameters (\eqn{\zeta} and \eqn{\xi}) and \code{p1} and \code{p2} domain 
#'   (minimum and maximum) or location-scale (mean and standard deviation) 
#'   parameters.
#'   
#' @name BMT
#' @aliases dBMT
#' @aliases pBMT
#' @aliases qBMT
#' @aliases rBMT
#'   
#' @details The BMT distribution with tails weights and domain parametrization
#'   (\code{type.p.3.4 = "t w"} and \code{type.p.1.2 = "c-d"}) has quantile
#'   function \deqn{(d - c) [3 t_p ( 1 - t_p )^2 \kappa_l - 3 t_p^2 ( 1 - t_p )
#'   \kappa_r + t_p^2 ( 3 - 2 t_p ) ] + c} where \eqn{0 \le p \le 1}, \eqn{t_p =
#'   1/2 - \cos ( [\arccos ( 2 p - 1 ) - 2 \pi] / 3 )}, and \eqn{0 < \kappa_l <
#'   1} and \eqn{0 < \kappa_r < 1} are, respectively, related to left and right
#'   tail weights or curvatures.
#'   
#'   The BMT coefficient of asymmetry \eqn{-1 < \zeta < 1} is \deqn{\kappa_r - 
#'   \kappa_l}
#'   
#'   The BMT coefficient of steepness \eqn{0 < \xi < 1} is \deqn{(\kappa_r + 
#'   \kappa_l - |\kappa_r - \kappa_l|) / (2 (1 -  |\kappa_r - \kappa_l|))} for 
#'   \eqn{|\kappa_r - \kappa_l| < 1}.
#'   
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'   to be the number required
#' @param p3,p4 tails weights (\eqn{\kappa_l} and \eqn{\kappa_r}) or 
#'   asymmetry-steepness (\eqn{\zeta} and \eqn{\xi}) parameters of the BMT 
#'   distribution.
#' @param type.p.3.4 type of parametrization associated to p3 and p4. "t w" means
#'   tails weights parametrization (default) and "a-s" means asymmetry-steepness
#'   parametrization.
#' @param p1,p2 domain (minimum and maximum) or location-scale (mean and 
#'   standard deviation) parameters of the BMT distribution.
#' @param type.p.1.2 type of parametrization associated to p1 and p2. "c-d" means
#'   domain parametrization (default) and "l-s" means location-scale 
#'   parametrization.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X \le 
#'   x]}, otherwise, \eqn{P[X > x]}.
#'   
#' @return \code{dBMT} gives the density, \code{pBMT} the distribution function,
#'   \code{qBMT} the quantile function, and \code{rBMT} generates random 
#'   deviates.
#'   
#'   The length of the result is determined by \code{n} for \code{rBMT}, and is 
#'   the maximum of the lengths of the numerical arguments for the other 
#'   functions.
#'   
#'   The numerical arguments other than \code{n} are recycled to the length of 
#'   the result. Only the first elements of the logical arguments are used.
#'   
#'   If \code{type.p.3.4 == "t w"}, \code{p3 < 0} and \code{p3 > 1} are errors 
#'   and return \code{NaN}.
#'   
#'   If \code{type.p.3.4 == "a-s"}, \code{p3 < -1} and \code{p3 > 1} are errors 
#'   and return \code{NaN}.
#'   
#'   \code{p4 < 0} and \code{p4 > 1} are errors and return \code{NaN}.
#'   
#'   If \code{type.p.1.2 == "c-d"}, \code{p1 >= p2} is an error and returns 
#'   \code{NaN}.
#'   
#'   If \code{type.p.1.2 == "l-s"}, \code{p2 <= 0} is an error and returns 
#'   \code{NaN}.
#'   
#' @references Torres-Jimenez, C. J. and Montenegro-Diaz, A. M. (2017, September), 
#'   \emph{An alternative to continuous univariate distributions supported on a 
#'   bounded interval: The BMT distribution}. ArXiv e-prints. \url{https://arxiv.org/abs/1709.05534}.
#'   
#'   Torres-Jimenez, C. J. (2017, September), \emph{Comparison of estimation methods 
#'   for the BMT distribution}. ArXiv e-prints.
#'   
#'   Torres-Jimenez, C. J. (2018), \emph{The BMT Item Response Theory model: A 
#'   new skewed distribution family with bounded domain and an IRT model based 
#'   on it}, PhD thesis, Doctorado en ciencias - Estadistica, Universidad 
#'   Nacional de Colombia, Sede Bogota.
#'   
#' @seealso \code{\link{BMTcentral}}, \code{\link{BMTdispersion}}, 
#'   \code{\link{BMTskewness}}, \code{\link{BMTkurtosis}}, 
#'   \code{\link{BMTmoments}} for descriptive measures or moments. 
#'   \code{\link{BMTchangepars}} for parameter conversion between different 
#'   parametrizations.
#'   
#' @author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co} 
#'   and Alvaro Mauricio Montenegro Diaz [ths]
#'   
#' @examples
#' # BMT on [0,1] with left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' z <- seq(0, 1, length.out = 100)
#' F1 <- pBMT(z, 0.25, 0.75, "t w")
#' Q1 <- qBMT(F1, 0.25, 0.75, "t w")
#' max(abs(z - Q1))
#' f1 <- dBMT(z, 0.25, 0.75, "t w")
#' r1 <- rBMT(100, 0.25, 0.75, "t w")
#' layout(matrix(c(1,2,1,3), 2, 2))
#' hist(r1, freq = FALSE, xlim = c(0,1))
#' lines(z, f1)
#' plot(z, F1, type="l")
#' plot(F1, Q1, type="l")
#' 
#' # BMT on [0,1] with asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.5
#' F2 <- pBMT(z, 0.5, 0.5, "a-s")
#' Q2 <- qBMT(F2, 0.5, 0.5, "a-s")
#' f2 <- dBMT(z, 0.5, 0.5, "a-s")
#' r2 <- rBMT(100, 0.5, 0.5, "a-s")
#' max(abs(f1 - f2))
#' max(abs(F1 - F2))
#' max(abs(Q1 - Q2))
#' 
#' # BMT on [-1.783489, 3.312195] with 
#' # left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' x <- seq(-1.783489, 3.312195, length.out = 100)
#' F3 <- pBMT(x, 0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' Q3 <- qBMT(F3, 0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' max(abs(x - Q3))
#' f3 <- dBMT(x, 0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' r3 <- rBMT(100, 0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' layout(matrix(c(1,2,1,3), 2, 2))
#' hist(r3, freq = FALSE, xlim = c(-1.783489,3.312195))
#' lines(x, f3)
#' plot(x, F3, type="l")
#' plot(F3, Q3, type="l")
#' 
#' # BMT with mean equal to 0, standard deviation equal to 1, 
#' # asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.5
#' f4 <- dBMT(x, 0.5, 0.5, "a-s", 0, 1, "l-s")
#' F4 <- pBMT(x, 0.5, 0.5, "a-s", 0, 1, "l-s")
#' Q4 <- qBMT(F4, 0.5, 0.5, "a-s", 0, 1, "l-s")
#' r4 <- rBMT(100, 0.5, 0.5, "a-s", 0, 1, "l-s")
#' max(abs(f3 - f4))
#' max(abs(F3 - F4))
#' max(abs(Q3 - Q4))
#' 
#' @keywords distribution

#' @rdname BMT
#' @export 
dBMT <- function(x, p3, p4, type.p.3.4 = "t w",
                 p1 = 0, p2 = 1, type.p.1.2 = "c-d",
                 log = FALSE) {
  # Control type.p.3.4
  TYPE.P.3.4 <- c("t w", "a-s") # tail weights or asymmetry-steepness
  int.type.p.3.4 <- pmatch(type.p.3.4, TYPE.P.3.4)
  if (is.na(int.type.p.3.4))
    stop("invalid type of parametrization for parameters 3 and 4")
  if (int.type.p.3.4 == -1)
    stop("ambiguous type of parametrization for parameters 3 and 4")
  # Control type.p.1.2
  TYPE.P.1.2 <- c("c-d", "l-s") # domain or location-scale
  int.type.p.1.2 <- pmatch(type.p.1.2, TYPE.P.1.2)
  if (is.na(int.type.p.1.2))
    stop("invalid type of parametrization for parameters 1 and 2")
  if (int.type.p.1.2 == -1)
    stop("ambiguous type of parametrization for parameters 1 and 2")
  # The length of the result is determined by the maximum of the lengths of the
  # numerical arguments. The numerical arguments are recycled to the length of
  # the result.
  len <- max(length(x),length(p1),length(p2),length(p3),length(p4))
  x <- rep(x, len = len)
  p1 <- rep(p1, len = len)
  p2 <- rep(p2, len = len)
  p3 <- rep(p3, len = len)
  p4 <- rep(p4, len = len)
  # Transform x to 0,1 given domain or location-scale parameters
  if (int.type.p.1.2 == 1) {
    # domain parametrization
    # Control domain parameters
    min <- replace(p1, p1 >= p2, NaN)
    max <- replace(p2, p1 >= p2, NaN)
    # Transform x
    range <- max - min
    x <- (x - min) / range
  }
  else{
    # location-scale parametrization
    # Control location-scale parameters
    mu <- p1
    sigma <- replace(p2, p2 <= 0, NaN)
    # Transform x
    range <- sigma / BMTsd(p3, p4, type.p.3.4)
    x <- (x - mu) / range + BMTmean(p3, p4, type.p.3.4)
  }
  # Obtain coefficients of polynomials x.t and yf.t given tail weights or
  # asymmetry-steepness parameters
  if (int.type.p.3.4 == 1) {
    # tail weights parametrization
    # Control tail weights parameters
    kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
    kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # Coefficients a_3*t^3+a_2*t^2+a_1*t+a_0
    a_3 <- 3 * kappa_l + 3 * kappa_r - 2
    a_2 <- (-6 * kappa_l - 3 * kappa_r + 3)
    a_1 <- (3 * kappa_l)
  }
  else{
    # asymmetry-steepness parametrization
    # Control asymmetry-steepness parameters
    zeta <- replace(p3, p3 < -1 | p3 > 1, NaN)
    xi <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # Coefficients a_3*t^3+a_2*t^2+a_1*t+a_0
    abs.zeta <- abs(zeta)
    aux1 <- 0.5 - xi
    a_3 <- 6 * (xi + abs.zeta * aux1) - 2
    a_2 <- -9 * (xi + abs.zeta * aux1) + 1.5 * zeta + 3
    a_1 <- 3 * (xi + abs.zeta * aux1) - 1.5 * zeta
  }
  # NaNs
  y <- x + a_3 + a_2 + a_1
  if (any(is.nan(y))) {
    warning("NaNs founded or produced")
  }
  ind <- !is.na(y)
  # Transformed x outside 0,1
  y[ind] <- 0
  # Transformed x inside 0,1
  ind[ind] <- x[ind] > 0 & x[ind] < 1
  if (any(ind)) {
    # inv.x.t
    t <- .inv.x.t(x[ind], a_3[ind], a_2[ind], a_1[ind])
    # yf.t
    y[ind] <- .yf.t(t, a_3[ind], a_2[ind], a_1[ind]) / range[ind]
  }
  # density values y are given as log(y)
  if (log)
    y <- log(y)
  return(y)
}

#' @rdname BMT
#' @export 
pBMT <- function(q, p3, p4, type.p.3.4 = "t w", 
                 p1 = 0, p2 = 1, type.p.1.2 = "c-d", 
                 lower.tail = TRUE, log.p = FALSE){
  # Control type.p.3.4
  TYPE.P.3.4 <- c("t w", "a-s") # tail weights or asymmetry-steepness
  int.type.p.3.4 <- pmatch(type.p.3.4, TYPE.P.3.4)
  if (is.na(int.type.p.3.4)) 
    stop("invalid type of parametrization for parameters 3 and 4")
  if (int.type.p.3.4 == -1) 
    stop("ambiguous type of parametrization for parameters 3 and 4")
  # Control type.p.1.2
  TYPE.P.1.2 <- c("c-d", "l-s") # domain or location-scale
  int.type.p.1.2 <- pmatch(type.p.1.2, TYPE.P.1.2)
  if (is.na(int.type.p.1.2)) 
    stop("invalid type of parametrization for parameters 1 and 2")
  if (int.type.p.1.2 == -1) 
    stop("ambiguous type of parametrization for parameters 1 and 2")
  # The length of the result is determined by the maximum of the lengths of the 
  # numerical arguments. The numerical arguments are recycled to the length of 
  # the result.
  len <- max(length(q),length(p1),length(p2),length(p3),length(p4))
  q <- rep(q, len=len)
  p1 <- rep(p1, len=len)
  p2 <- rep(p2, len=len)
  p3 <- rep(p3, len=len)
  p4 <- rep(p4, len=len)
  # Transform q to 0,1 given domain or location-scale parameters
  if(int.type.p.1.2 == 1){ # domain parametrization
    # Control domain parameters
    min <- replace(p1, p1 >= p2, NaN)
    max <- replace(p2, p1 >= p2, NaN)
    # Transform q
    range <- max - min
    q <- (q - min)/range
  }
  else{ # location-scale parametrization
    # Control location-scale parameters
    mu <- p1
    sigma <- replace(p2, p2 <= 0, NaN)
    # Transform q
    range <- sigma/BMTsd(p3, p4, type.p.3.4)
    q <- (q - mu)/range + BMTmean(p3, p4, type.p.3.4)
  }
  # Obtain coefficients of polynomials x.t and yf.t given tail weights or
  # asymmetry-steepness parameters
  if(int.type.p.3.4 == 1){ # tail weights parametrization
    # Control tail weights parameters
    kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
    kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # Coefficients a_3*t^3+a_2*t^2+a_1*t+a_0
    a_3 <- 3*kappa_l+3*kappa_r-2
    a_2 <- (-6*kappa_l-3*kappa_r+3)
    a_1 <- (3*kappa_l)
  }
  else{ # asymmetry-steepness parametrization
    # Control asymmetry-steepness parameters
    zeta <- replace(p3, p3 < -1 | p3 > 1, NaN)
    xi <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # Coefficients a_3*t^3+a_2*t^2+a_1*t+a_0
    abs.zeta <- abs(zeta)
    aux1 <- 0.5-xi
    a_3 <- 6*(xi+abs.zeta*aux1)-2
    a_2 <- -9*(xi+abs.zeta*aux1)+1.5*zeta+3
    a_1 <- 3*(xi+abs.zeta*aux1)-1.5*zeta
  }
  # NaNs
  p <- q+a_3+a_2+a_1
  if(any(is.nan(p)))
    warning("NaNs founded or produced")
  ind <- !is.na(p)
  # Transformed q outside 0,1
  p[ind] <- 0
  ind2 <- ind
  ind2[ind] <- q[ind] >= 1
  p[ind2] <- 1
  # Transformed q inside 0,1
  ind[ind] <- q[ind] > 0 & q[ind] < 1
  if(any(ind)){
    # inv.x.t
    t <- .inv.x.t(q[ind], a_3[ind], a_2[ind], a_1[ind])
    # yF.t
    p[ind] <- .yF.t(t)
  }
  # probabilities are \eqn{P[X > x]}
  if(!lower.tail)
    p <- 1 - p
  # probabilities p are given as log(p)
  if(log.p)
    p <- log(p)
  return(p)
}

#' @rdname BMT
#' @export 
qBMT <- function(p, p3, p4, type.p.3.4 = "t w", 
                 p1 = 0, p2 = 1, type.p.1.2 = "c-d", 
                 lower.tail = TRUE, log.p = FALSE){
  # probabilities p are given as log(p)
  if(log.p)
    p <- exp(p)
  # probabilities are \eqn{P[X > x]}
  if(!lower.tail)
    p <- 1 - p
  # Control type.p.3.4
  TYPE.P.3.4 <- c("t w", "a-s") # tail weights or asymmetry-steepness
  int.type.p.3.4 <- pmatch(type.p.3.4, TYPE.P.3.4)
  if (is.na(int.type.p.3.4)) 
    stop("invalid type of parametrization for parameters 3 and 4")
  if (int.type.p.3.4 == -1) 
    stop("ambiguous type of parametrization for parameters 3 and 4")
  # Control type.p.1.2
  TYPE.P.1.2 <- c("c-d", "l-s") # domain or location-scale
  int.type.p.1.2 <- pmatch(type.p.1.2, TYPE.P.1.2)
  if (is.na(int.type.p.1.2)) 
    stop("invalid type of parametrization for parameters 1 and 2")
  if (int.type.p.1.2 == -1) 
    stop("ambiguous type of parametrization for parameters 1 and 2")
  # The length of the result is determined by the maximum of the lengths of the 
  # numerical arguments. The numerical arguments are recycled to the length of 
  # the result.
  len <- max(length(p),length(p1),length(p2),length(p3),length(p4))
  p <- rep(p, len=len)
  p1 <- rep(p1, len=len)
  p2 <- rep(p2, len=len)
  p3 <- rep(p3, len=len)
  p4 <- rep(p4, len=len)
  # Obtain coefficients of polynomials x.t and yf.t given tail weights or
  # asymmetry-steepness parameters
  if(int.type.p.3.4 == 1){ # tail weights parametrization
    # Control tail weights parameters
    kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
    kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # Coefficients a_3*t^3+a_2*t^2+a_1*t+a_0
    a_3 <- 3*kappa_l+3*kappa_r-2
    a_2 <- (-6*kappa_l-3*kappa_r+3)
    a_1 <- (3*kappa_l)
  }
  else{ # asymmetry-steepness parametrization
    # Control asymmetry-steepness parameters
    zeta <- replace(p3, p3 < -1 | p3 > 1, NaN)
    xi <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # Coefficients a_3*t^3+a_2*t^2+a_1*t+a_0
    abs.zeta <- abs(zeta)
    aux1 <- 0.5-xi
    a_3 <- 6*(xi+abs.zeta*aux1)-2
    a_2 <- -9*(xi+abs.zeta*aux1)+1.5*zeta+3
    a_1 <- 3*(xi+abs.zeta*aux1)-1.5*zeta
  }
  # NaNs
  q <- p+a_3+a_2+a_1
  if(any(is.nan(q)))
    warning("NaNs founded or produced")
  ind <- !is.na(q)
  # q outside (0,1)
  q <- rep(NaN,len)
  ind2 <- ind
  ind2[ind] <- p[ind] == 0
  q[ind2] <- 0
  ind2 <- ind
  ind2[ind] <- p[ind] == 1
  q[ind2] <- 1
  # q inside (0,1)
  ind[ind] <- p[ind] > 0 & p[ind] < 1
  if(any(ind)){
    # inv.yF.t
    t <- .inv.yF.t(p[ind])
    # x.t
    q[ind] <- .x.t(t, a_3[ind], a_2[ind], a_1[ind])
  }
  # Transform q to [c,d] given domain or location-scale parameters
  if(int.type.p.1.2 == 1){ # domain parametrization
    # Control domain parameters
    min <- replace(p1, p1 >= p2, NaN)
    max <- replace(p2, p1 >= p2, NaN)
    # Transform q
    range <- max - min
    q <- q * range + min
  }
  else{ # location-scale parametrization
    # Control location-scale parameters
    mu <- p1
    sigma <- replace(p2, p2 <= 0, NaN)
    # Transform q
    range <- sigma/BMTsd(p3, p4, type.p.3.4)
    q <- (q - BMTmean(p3, p4, type.p.3.4)) * range + mu
  }
  return(q)
}

#' @rdname BMT
#' @export 
rBMT <- function(n, p3, p4, type.p.3.4 = "t w", 
                 p1 = 0, p2 = 1, type.p.1.2 = "c-d"){
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
  x <- qBMT(p, p3, p4, type.p.3.4, p1, p2, type.p.1.2)
  return(x)
}

# Global constants
.epsilon <- 1e-10
.zero <- 0-.epsilon
.one <- 1+.epsilon

# x.t
.x.t <- function(t, a_3, a_2, a_1){
  x <- ((a_3*t + a_2)*t + a_1)*t
  return(x)
}

# Inverse function for x.t
.inv.x.t <- function(x, a_3, a_2, a_1){
  # Press W.H., Teukolsky S.A., Vetterling W.T. & Flannery B.P. 2007. 
  # Numerical recipes: The art of scientific computing 
  # Section 5.6: Quadratic and Cubic Equations. Page 228.
  len <- length(x)
  a <- a_2/a_3
  b <- a_1/a_3
  c <- -x/a_3
  Q <- (a*a - 3*b)/9
  R <- ((2*a*a - 9*b)*a + 27*c)/54
  r <- rep(0,len)
  # All real roots
  ind.r <- Q*Q*Q-R*R > 0
  a.v <- a[ind.r]
  Q.v <- Q[ind.r]
  R.v <- R[ind.r]
  theta <- acos(R.v/sqrt(Q.v*Q.v*Q.v))
  aux1 <- -2*sqrt(Q.v)
  aux2 <- a.v/3
  r.v <- aux1*cos(theta/3)-aux2
  ind.no01 <- r.v < .zero | r.v > .one
  r.v[ind.no01] <- aux1[ind.no01]*cos((theta[ind.no01]+2*pi)/3)-aux2[ind.no01]
  ind.no01 <- r.v < .zero | r.v > .one
  r.v[ind.no01] <- aux1[ind.no01]*cos((theta[ind.no01]-2*pi)/3)-aux2[ind.no01]
  r[ind.r] <- r.v
  # Two complex roots
  a.v <- a[!ind.r]
  Q.v <- Q[!ind.r]
  R.v <- R[!ind.r]
  aux2 <- a.v/3
  A <- -sign(R.v)*(abs(R.v)+sqrt(R.v*R.v-Q.v*Q.v*Q.v))^(1/3)
  B <- rep(0,length(A))
  B[A!=0] <- Q.v[A!=0]/A[A!=0]
  r.v <- (A+B)-aux2
  ind.no01 <- r.v < .zero | r.v > .one
  r.v[ind.no01] <- -0.5*(r.v[ind.no01])-1.5*aux2[ind.no01]
  r[!ind.r] <- r.v
  # Considering an epsilon, all roots in [0,1]
  r[r >= .zero & r < 0] <- 0
  r[r > 1 & r <= .one] <- 1
  return(r)
}

# yF.t
.yF.t <- function(t){
  yF <- (-2*t + 3)*t*t
  return(yF)
}

# inv.yF.t
.inv.yF.t <- function(yF){
  t <- 0.5-cos((acos(2*yF-1)-2*pi)/3)
  return(t)
}

# yf.t
.yf.t <- function(t, a_3, a_2, a_1){
  yf <- ((-6*t + 6)*t) / ((3*a_3*t + 2*a_2)*t + a_1)
  return(yf)
}
