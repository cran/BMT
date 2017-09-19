#' @title The BMT Distribution Descriptive Measures - Skewness.
#' @description Skewness and an asymmetry coefficient for the BMT distribution,
#'   with \code{p3} and \code{p4} tails weights (\eqn{\kappa_l} and
#'   \eqn{\kappa_r}) or asymmetry-steepness parameters (\eqn{\zeta} and
#'   \eqn{\xi}) and \code{p1} and \code{p2} domain (minimum and maximum) or
#'   location-scale (mean and standard deviation) parameters.
#' @rdname BMTskewness
#' @name BMTskewness
#' @aliases BMTskew
#' @aliases BMTasymm
#'   
#' @details See References.
#'   
#' @param p3,p4 tails weights (\eqn{\kappa_l} and \eqn{\kappa_r}) or 
#'   asymmetry-steepness (\eqn{\zeta} and \eqn{\xi}) parameters of the BMT 
#'   distribution.
#' @param type.p.3.4 type of parametrization asociated to p3 and p4. "t w" means
#'   tails weights parametrization (default) and "a-s" means asymmetry-steepness
#'   parametrization.
#' @param p1,p2 domain (minimum and maximum) or location-scale (mean and 
#'   standard deviation) parameters of the BMT ditribution.
#' @param type.p.1.2 type of parametrization asociated to p1 and p2. "c-d" means
#'   domain parametrization (default) and "l-s" means location-scale 
#'   parametrization.
#'   
#' @return \code{BMTskew} gives the Pearson's skewness and \code{BMTasymm} the 
#'   proposed asymmetry coefficient for the BMT distribution.
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
#'   bounded interval: The BMT distribution}. ArXiv e-prints.
#'   
#'   Torres-Jimenez, C. J. (2018), \emph{The BMT Item Response Theory model: A 
#'   new skewed distribution family with bounded domain and an IRT model based 
#'   on it}, PhD thesis, Doctorado en ciencias - Estadistica, Universidad 
#'   Nacional de Colombia, Sede Bogota.
#'   
#' @seealso \code{\link{BMTcentral}}, \code{\link{BMTdispersion}}, 
#'   \code{\link{BMTkurtosis}}, \code{\link{BMTmoments}} for other descriptive 
#'   measures or moments.
#'   
#' @author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co}
#'   
#' @examples
#' # BMT on [0,1] with left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' BMTskew(0.25, 0.75, "t w")
#' BMTasymm(0.25, 0.75, "t w")
#' 
#' # BMT on [0,1] with asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.75
#' BMTskew(0.5, 0.5, "a-s")
#' BMTasymm(0.5, 0.5, "a-s")
#' 
#' # domain or location-scale parameters do not affect 
#' # the skewness and the asymmetry coefficient
#' 
#' # BMT on [-1.783489,3.312195] with 
#' # left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' BMTskew(0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' BMTasymm(0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' 
#' # BMT with mean equal to 0, standard deviation equal to 1, 
#' # asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.75
#' BMTskew(0.5, 0.5, "a-s", 0, 1, "l-s")
#' BMTasymm(0.5, 0.5, "a-s", 0, 1, "l-s")

#' @rdname BMTskewness
#' @export BMTskew
BMTskew <- function(p3, p4, type.p.3.4 = "t w", 
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
  # tail weigths or asymmetry-steepness parametrization
  if(int.type.p.3.4 == 1){ # tail weights parametrization
    # Control tail weights parameters
    kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
    kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # skewness
    m <- 9*(kappa_r-kappa_l)*((13*kappa_l-65+4*kappa_r)*kappa_l+(13*kappa_r-65)*kappa_r+150) /
         (77000*BMTvar(kappa_l,kappa_r,type.p.3.4)^1.5)
  }
  else{ # asymmetry-steepness parametrization
    # Control asymmetry-steepness parameters
    zeta <- replace(p3, p3 < -1 | p3 > 1, NaN)
    xi <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # skewness
    abs.zeta <- abs(zeta)
    m <- 9*zeta*
         ((((30*xi-30)*xi+13)*abs.zeta+
           ((-60*xi+160)*xi-65))*abs.zeta+
           ((30*xi-130)*xi+150)) /
         (77000*BMTvar(zeta,xi,type.p.3.4)^1.5)
  }
  # domain or location-scale parameters do not affect the skewness
  return(m)
}

#' @rdname BMTskewness
#' @export BMTasymm
BMTasymm <- function(p3, p4, type.p.3.4 = "t w", 
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
  # tail weigths or asymmetry-steepness parametrization
  if(int.type.p.3.4 == 1){ # tail weights parametrization
    # Control tail weights parameters
    kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
    kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # asymmetry coefficient
    m <- kappa_r - kappa_l
  }
  else{ # asymmetry-steepness parametrization
    # Control asymmetry-steepness parameters
    zeta <- replace(p3, p3 < -1 | p3 > 1, NaN)
    # asymmetry coefficient
    m <- zeta
  }
  # domain or location-scale parameters do not affect the asymmetry coefficient
  return(m)
}