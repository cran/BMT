#' @title The BMT Distribution Descriptive Measures - Central Tendency.
#' @description Mean, median and mode for the BMT distribution, with \code{p3} 
#'   and \code{p4} tails weights (\eqn{\kappa_l} and \eqn{\kappa_r}) or 
#'   asymmetry-steepness parameters (\eqn{\zeta} and \eqn{\xi}) and \code{p1} 
#'   and \code{p2} domain (minimum and maximum) or location-scale (mean and 
#'   standard deviation) parameters.
#' 
#' @name BMTcentral
#' @aliases BMTmean
#' @aliases BMTmedian
#' @aliases BMTmode
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
#'   
#' @return \code{BMTmean} gives the mean, \code{BMTmedian} the median and 
#'   \code{BMTmode} the mode for the BMT distribution.
#'   
#'   The arguments are recycled to the length of the result. Only the first 
#'   elements of \code{type.p.3.4} and \code{type.p.1.2} are used.
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
#' @seealso \code{\link{BMTdispersion}}, \code{\link{BMTskewness}}, 
#'   \code{\link{BMTkurtosis}}, \code{\link{BMTmoments}} for other descriptive 
#'   measures or moments.
#'   
#' @author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co}
#'   
#' @examples
#' # BMT on [0,1] with left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' BMTmean(0.25, 0.75, "t w")
#' BMTmedian(0.25, 0.75, "t w")
#' BMTmode(0.25, 0.75, "t w")
#' 
#' # BMT on [0,1] with asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.75
#' BMTmean(0.5, 0.5, "a-s")
#' BMTmedian(0.5, 0.5, "a-s")
#' BMTmode(0.5, 0.5, "a-s")
#' 
#' # BMT on [-1.783489,3.312195] with 
#' # left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' BMTmean(0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' BMTmedian(0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' BMTmode(0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' 
#' # BMT with mean equal to 0, standard deviation equal to 1, 
#' # asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.75
#' BMTmean(0.5, 0.5, "a-s", 0, 1, "l-s")
#' BMTmedian(0.5, 0.5, "a-s", 0, 1, "l-s")
#' BMTmode(0.5, 0.5, "a-s", 0, 1, "l-s")

#' @rdname BMTcentral
#' @export
BMTmean <- function(p3, p4, type.p.3.4 = "t w", 
                    p1 = 0, p2 = 1, type.p.1.2 = "c-d"){
  # The length of the result is determined by the maximum of the lengths of the
  # numerical arguments. The numerical arguments are recycled to the length of
  # the result.
  len <- max(length(p1),length(p2),length(p3),length(p4))
  p1 <- rep(p1, len=len)
  p2 <- rep(p2, len=len)
  p3 <- rep(p3, len=len)
  p4 <- rep(p4, len=len)
  # Control type.p.3.4
  TYPE.P.3.4 <- c("t w", "a-s")
  int.type.p.3.4 <- pmatch(type.p.3.4, TYPE.P.3.4)
  if (is.na(int.type.p.3.4)) 
    stop("invalid type of parametrization for parameters 3 and 4")
  if (int.type.p.3.4 == -1) 
    stop("ambiguous type of parametrization for parameters 3 and 4")
  # Control type.p.1.2
  TYPE.P.1.2 <- c("c-d", "l-s")
  int.type.p.1.2 <- pmatch(type.p.1.2, TYPE.P.1.2)
  if (is.na(int.type.p.1.2)) 
    stop("invalid type of parametrization for parameters 1 and 2")
  if (int.type.p.1.2 == -1) 
    stop("ambiguous type of parametrization for parameters 1 and 2")
  # domain or location-scale parametrization
  if(int.type.p.1.2 == 1){ # domain parametrization
    # Control domain parameters
    min <- replace(p1, p1 >= p2, NaN)
    max <- replace(p2, p1 >= p2, NaN)
    # range
    range <- max - min
    # tail weights or asymmetry-steepness parametrization
    if(int.type.p.3.4 == 1){ # tail weights parametrization
      # Control tail weights parameters
      kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
      kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
      # mean
      m <- 0.3*(kappa_l - kappa_r) + 0.5
    }
    else{ # asymmetry-steepness parametrization
      # Control asymmetry-steepness parameters
      zeta <- replace(p3, p3 < -1 | p3 > 1, NaN)
      xi <- replace(p4, p4 < 0 | p4 > 1, NaN)
      # mean
      m <- -0.3*zeta + 0.5
    }
    # scaled and shifted mean
    m <- range*m + min
  }
  else{ # location-scale parametrization
    # mean
    m <- p1
  }
  return(m)
}

#' @rdname BMTcentral
#' @export 
BMTmedian <- function(p3, p4, type.p.3.4 = "t w", 
                      p1 = 0, p2 = 1, type.p.1.2 = "c-d"){
  # The length of the result is determined by the maximum of the lengths of the
  # numerical arguments. The numerical arguments are recycled to the length of
  # the result.
  len <- max(length(p1),length(p2),length(p3),length(p4))
  p1 <- rep(p1, len=len)
  p2 <- rep(p2, len=len)
  p3 <- rep(p3, len=len)
  p4 <- rep(p4, len=len)
  # Control type.p.3.4
  TYPE.P.3.4 <- c("t w", "a-s")
  int.type.p.3.4 <- pmatch(type.p.3.4, TYPE.P.3.4)
  if (is.na(int.type.p.3.4)) 
    stop("invalid type of parametrization for parameters 3 and 4")
  if (int.type.p.3.4 == -1) 
    stop("ambiguous type of parametrization for parameters 3 and 4")
  # Control type.p.1.2
  TYPE.P.1.2 <- c("c-d", "l-s")
  int.type.p.1.2 <- pmatch(type.p.1.2, TYPE.P.1.2)
  if (is.na(int.type.p.1.2)) 
    stop("invalid type of parametrization for parameters 1 and 2")
  if (int.type.p.1.2 == -1) 
    stop("ambiguous type of parametrization for parameters 1 and 2")
  # tail weights or asymmetry-steepness parametrization
  if(int.type.p.3.4 == 1){ # tail weights parametrization
    # Control tail weights parameters
    kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
    kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # median
    m <- 0.375*(kappa_l - kappa_r) + 0.5
  }
  else{ # asymmetry-steepness parametrization
    # Control asymmetry-steepness parameters
    zeta <- replace(p3, p3 < -1 | p3 > 1, NaN)
    xi <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # median
    m <- -0.375*zeta + 0.5
  }
  # domain or location-scale parametrization
  if(int.type.p.1.2 == 1){ # domain parametrization
    # Control domain parameters
    min <- replace(p1, p1 >= p2, NaN)
    max <- replace(p2, p1 >= p2, NaN)
    # range
    range <- max - min
    # scaled and shifted median
    m <- range*m + min
  }
  else{ # location-scales parametrization
    # Control location-scale parameters
    mu <- p1
    sigma <- replace(p2, p2 <= 0, NaN)
    # range
    range <- sigma/BMTsd(p3, p4, type.p.3.4)
    # scaled and shifted median
    m <- range*(m - BMTmean(p3, p4, type.p.3.4)) + mu
  }
  return(m)
}

#' @rdname BMTcentral
#' @export
BMTmode <- function(p3, p4, type.p.3.4 = "t w", 
                    p1 = 0, p2 = 1, type.p.1.2 = "c-d"){
  # The length of the result is determined by the maximum of the lengths of the
  # numerical arguments. The numerical arguments are recycled to the length of
  # the result.
  len <- max(length(p1),length(p2),length(p3),length(p4))
  p1 <- rep(p1, len=len)
  p2 <- rep(p2, len=len)
  p3 <- rep(p3, len=len)
  p4 <- rep(p4, len=len)
  # Control type.p.3.4
  TYPE.P.3.4 <- c("t w", "a-s")
  int.type.p.3.4 <- pmatch(type.p.3.4, TYPE.P.3.4)
  if (is.na(int.type.p.3.4)) 
    stop("invalid type of parametrization for parameters 3 and 4")
  if (int.type.p.3.4 == -1) 
    stop("ambiguous type of parametrization for parameters 3 and 4")
  # Control type.p.1.2
  TYPE.P.1.2 <- c("c-d", "l-s")
  int.type.p.1.2 <- pmatch(type.p.1.2, TYPE.P.1.2)
  if (is.na(int.type.p.1.2)) 
    stop("invalid type of parametrization for parameters 1 and 2")
  if (int.type.p.1.2 == -1) 
    stop("ambiguous type of parametrization for parameters 1 and 2")
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
  # mode
  aux1 <- sqrt(kappa_l*kappa_r)
  m <- ifelse(kappa_l==kappa_r,0.5,
              (kappa_l^2 - 5*kappa_l*kappa_r + kappa_l*aux1 + 3*kappa_r*aux1 + 
               9*kappa_l*kappa_r^2 + 3*kappa_l^2*kappa_r - 3*kappa_r^2*aux1 - 
               9*kappa_l*kappa_r*aux1)*(kappa_l-aux1)/(kappa_l-kappa_r)^3)
  # domain or location-scale parametrization
  if(int.type.p.1.2 == 1){ # domain parametrization
    # Control domain parameters
    min <- replace(p1, p1 >= p2, NaN)
    max <- replace(p2, p1 >= p2, NaN)
    # range
    range <- max - min
    # scaled and shifted mode
    m <- range*m + min
  }
  else{ # location-scale parametrization
    # Control location-scale parameters
    mu <- p1
    sigma <- replace(p2, p2 <= 0, NaN)
    # range
    range <- sigma/BMTsd(p3, p4, type.p.3.4)
    # scaled and shifted mode
    m <- range*(m - BMTmean(p3, p4, type.p.3.4)) + mu
  }
  return(m)
}
