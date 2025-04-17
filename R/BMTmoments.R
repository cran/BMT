#' @title The BMT Distribution Moments, Moment-Generating Function and 
#'   Characteristic Function.
#' @description Any raw, central or standardized moment, the moment-generating 
#'   function and the characteristic function for the BMT distribution, with 
#'   \code{p3} and \code{p4} tails weights (\eqn{\kappa_l} and \eqn{\kappa_r}) 
#'   or asymmetry-steepness parameters (\eqn{\zeta} and \eqn{\xi}) and \code{p1}
#'   and \code{p2} domain (minimum and maximum) or location-scale (mean and 
#'   standard deviation) parameters.
#' 
#' @name BMTmoments
#' @aliases BMTmoment
#' @aliases BMTmgf
#' @aliases BMTchf
#'   
#' @details See References.
#'   
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
#' @param order order of the moment.
#' @param type type of the moment: raw, central or standardized (default).
#' @param method method to obtain the moment: exact formula or Chebyshev-Gauss 
#'   quadrature (default).
#' @param s variable for the moment-generating and characteristic functions.
#'   
#' @return \code{BMTmoment} gives any raw, central or standardized moment, 
#'   \code{BMTmgf} the moment-generating function and \code{BMTchf} the 
#'   characteristic function
#'   
#'   The arguments are recycled to the length of the result. Only the first 
#'   elements of \code{type.p.3.4}, \code{type.p.1.2}, \code{type} and
#'   \code{method} are used.
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
#'   Torres-Jimenez, C. J. (2018), \emph{The BMT Item Response Theory model: A 
#'   new skewed distribution family with bounded domain and an IRT model based 
#'   on it}, PhD thesis, Doctorado en ciencias - Estadistica, Universidad 
#'   Nacional de Colombia, Sede Bogota.
#'   
#' @seealso \code{\link{BMTcentral}}, \code{\link{BMTdispersion}}, 
#'   \code{\link{BMTskewness}}, \code{\link{BMTkurtosis}} for specific 
#'   descriptive measures or moments.
#'   
#' @author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co}
#'   
#' @examples
#' layout(matrix(1:4, 2, 2, TRUE))
#' s <- seq(-1, 1, length.out = 100)
#' 
#' # BMT on [0,1] with left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' BMTmoment(0.25, 0.75, order = 5) # hyperskewness by Gauss-Legendre quadrature
#' BMTmoment(0.25, 0.75, order = 5, method = "exact") # hyperskewness by exact formula
#' mgf <- BMTmgf(s, 0.25, 0.75) # moment-generation function
#' plot(s, mgf, type="l")
#' chf <- BMTchf(s, 0.25, 0.75) # characteristic function
#' 
#' # BMT on [0,1] with asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.5
#' BMTmoment(0.5, 0.5, "a-s", order = 5)
#' BMTmoment(0.5, 0.5, "a-s", order = 5, method = "exact")
#' mgf <- BMTmgf(s, 0.5, 0.5, "a-s")
#' plot(s, mgf, type="l")
#' chf <- BMTchf(s, 0.5, 0.5, "a-s")
#' 
#' # BMT on [-1.783489, 3.312195] with 
#' # left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' BMTmoment(0.25, 0.75, "t w", -1.783489, 3.312195, "c-d", order = 5)
#' BMTmoment(0.25, 0.75, "t w", -1.783489, 3.312195, "c-d", order = 5, method = "exact")
#' mgf <- BMTmgf(s, 0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' plot(s, mgf, type="l")
#' chf <- BMTchf(s, 0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' 
#' # BMT with mean equal to 0, standard deviation equal to 1, 
#' # asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.5
#' BMTmoment(0.5, 0.5, "a-s", 0, 1, "l-s", order = 5)
#' BMTmoment(0.5, 0.5, "a-s", 0, 1, "l-s", order = 5, method = "exact")
#' mgf <- BMTmgf(s, 0.5, 0.5, "a-s", 0, 1, "l-s")
#' plot(s, mgf, type="l")
#' chf <- BMTchf(s, 0.5, 0.5, "a-s", 0, 1, "l-s")

#' @rdname BMTmoments
#' @export
BMTmoment <- function(p3, p4, type.p.3.4 = "t w",
                 p1 = 0, p2 = 1, type.p.1.2 = "c-d",
                 order, type = "standardised", method = "quadrature"){
  # Control order
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  if (any(!is.wholenumber(order)) || any(order < 1)) 
    stop("order should be a vector of integers greater or equal than 1")
  # Control type
  type <- match.arg(type, c("raw","central","standardised"))
  # Control method
  method <- match.arg(method, c("quadrature","exact"))
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
  len1 <- max(length(p3),length(p4))
  p3 <- rep(p3, len=len1)
  p4 <- rep(p4, len=len1)
  len2 <- max(length(p1),length(p2))
  p1 <- rep(p1, len=len2)
  p2 <- rep(p2, len=len2)
  # domain or location-scale parametrization
  if(int.type.p.1.2 == 1){ # domain parametrization
    # Control domain parameters
    min <- replace(p1, p1 >= p2, NaN)
    max <- replace(p2, p1 >= p2, NaN)
    # scale
    a <- max - min
    # shift
    b <- min
  }
  else{ # location-scale parametrization
    # Control location-scale parameters
    mu <- p1
    sigma <- replace(p2, p2 <= 0, NaN)
    # scale
    a <- sigma/BMTsd(p3, p4, type.p.3.4)
    # shift
    b <- mu - a * BMTmean(p3, p4, type.p.3.4)
  }
  # Obtain moments
  if(method=="quadrature"){ # by quadrature
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
    # function 
    funct1 <- function(order,a_3,a_2,a_1,a,b){
      # 10 points for the Gauss-Legendre quadrature over [0,1] (22 digits)
      t <- 0.5*.GL.10.points + 0.5
      # x.t
      x.t <- .x.t(t, a_3, a_2, a_1)
      # scaling and shifting (or minus mean for central or standardised)
      if(a!=1)
        x.t <- a * x.t
      if(b!=0)
        x.t <- x.t + b
      # Derivative of yF.t
      yFp.t <- 6*t*(1-t)
      # Gauss-Legendre quadrature over [0,1]
      return(0.5*sum(.GL.10.weights*(x.t^order)*yFp.t))
    }
    # by type of moment
    if(type=="raw"){ # raw
      # moments (vectorised form)
      m <- mapply(funct1,order,a_3,a_2,a_1,a,b) 
    }
    else{
      # mean
      mean <- BMTmean(p3, p4, type.p.3.4)
      # moments (vectorised form)
      m <- mapply(funct1,order,a_3,a_2,a_1,rep(1,len=len2),-mean)
      if(type=="central"){ # central
        # scaled moments
        m <- a^order * m
      }
      else{ # standardised
        # standard deviation
        sigma <- BMTsd(p3, p4, type.p.3.4)
        # standardised moments
        m <- m / ( sigma^order )
      }
    }
  }
  else{ # by exact formula
    # tail weights or asymmetry-steepness parametrization
    if(int.type.p.3.4 == 1){ # tail weights parametrization
      # Control tail weights parameters
      kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
      kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
    }
    else{ # asymmetry-steepness parametrization
      # change parametrization
      p <- BMTchangepars(p3, p4, type.p.3.4)
      kappa_l <- p$p3
      kappa_r <- p$p4
    }
    # function
    funct2 <- function(kappa_l,kappa_r,order,a,b){
      # Order 4 composition of order including zero
      K <- partitions::compositions(order, 4, include.zero=TRUE)
      # function for each term of the sum
      term4 <- function(v,kappa_l,kappa_r,order,a,b){
        term4 <- factorial(order) * 3^(v[2]+v[3]) * 
          ifelse(v[1]==0,1,(b)^v[1]) * 
          ifelse(v[2]==0,1,(a*kappa_l+b)^v[2]) * 
          ifelse(v[3]==0,1,(a*(1-kappa_r)+b)^v[3]) * 
          ifelse(v[4]==0,1,(a+b)^v[4]) /
          factorial(v[1]) / 
          factorial(v[2]) / 
          factorial(v[3]) / 
          factorial(v[4]) / 
          choose((3*order+2),(1+v[2]+2*v[3]+3*v[4]))
        return(term4)
      }
      return(2/(order+1) * sum(apply(K,2,term4,kappa_l,kappa_r,order,a,b)))
    }
    # by type of moment    
    if(type=="raw"){
      # moment 
      m <- mapply(funct2,kappa_l,kappa_r,order,a,b)
    }
    else{
      # mean
      mean <- BMTmean(kappa_l, kappa_r)
      # moment
      m <- mapply(funct2,kappa_l,kappa_r,order,rep(1,len=len2),-mean)
      if(type=="central"){ # central
        # scaled moments
        m <- a^order * m
      }
      else{ # standardised
        # standard deviation
        sigma <- BMTsd(kappa_l, kappa_r)
        # standardised moments
        m <- m / ( sigma^order )
      }
    }
  }
  return(m)
}

#' @rdname BMTmoments
#' @export
BMTmgf <- function(s, p3, p4, type.p.3.4 = "t w",
                   p1 = 0, p2 = 1, type.p.1.2 = "c-d"){
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
  len1 <- max(length(p3),length(p4))
  p3 <- rep(p3, len=len1)
  p4 <- rep(p4, len=len1)
  len2 <- max(length(s),length(p1),length(p2))
  s <- rep(s, len=len2)
  p1 <- rep(p1, len=len2)
  p2 <- rep(p2, len=len2)
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
  #
  funct3 <- function(s,a_3,a_2,a_1){
    # 10 points for the Gauss-Legendre quadrature over [0,1] (22 digits)
    t <- 0.5*.GL.10.points + 0.5
    # x.t
    x.t <- .x.t(t, a_3, a_2, a_1)
    # Derivative of yF.t
    yFp.t <- 6*t*(1-t)
    # Gauss-Legendre quadrature over [0,1]
    return(0.5*sum(.GL.10.weights*exp(s*x.t)*yFp.t))
  }
  # domain or location-scale parametrization
  if(int.type.p.1.2 == 1){ # domain parametrization
    # Control domain parameters
    min <- replace(p1, p1 >= p2, NaN)
    max <- replace(p2, p1 >= p2, NaN)
    # range
    range <- max - min
    # scaled and shifted
    y <- mapply(funct3,range*s,a_3,a_2,a_1)*exp(min*s)
  }
  else{ # location-scale parametrization
    # Control location-scale parameters
    mu <- p1
    sigma <- replace(p2, p2 <= 0, NaN)
    # range
    range <- sigma/BMTsd(p3, p4, type.p.3.4)
    # scaled and shifted
    y <- mapply(funct3,range*s,a_3,a_2,a_1)*exp((mu-range*BMTmean(p3, p4, type.p.3.4))*s)
  }
  return(y)
}

#' @rdname BMTmoments
#' @export
BMTchf <- function(s, p3, p4, type.p.3.4 = "t w",
                   p1 = 0, p2 = 1, type.p.1.2 = "c-d"){
  y <- BMTmgf(1i*s, p3, p4, type.p.3.4, p1, p2, type.p.1.2)
  return(y)
}

#' @rdname BMTmoments
#' @export
mBMT <- function(order, p3, p4, type.p.3.4, p1, p2, type.p.1.2){
  fun <- switch(order,BMTmean,BMTsd,BMTskew,BMTkurt)
  return(fun(p3, p4, type.p.3.4, p1, p2, type.p.1.2))
}

# Global constants
# 10 points for the Gauss-Legendre quadrature over [-1,1] (22 digits)
.GL.10.points <- c(-0.973906528517171720078,
                   -0.8650633666889845107321,
                   -0.6794095682990244062343,
                   -0.4333953941292471907993,
                   -0.148874338981631210885,
                   0.1488743389816312108848,
                   0.433395394129247190799,
                   0.6794095682990244062343,
                   0.8650633666889845107321,
                   0.973906528517171720078)
# Weights for 10 points of the Gauss-Legendre quadrature
.GL.10.weights <- c(0.0666713443086881375936,
                    0.149451349150580593146,
                    0.2190863625159820439955,
                    0.2692667193099963550912,
                    0.295524224714752870174,
                    0.295524224714752870174,
                    0.2692667193099963550913,
                    0.219086362515982043995,
                    0.1494513491505805931458,
                    0.0666713443086881375936)
