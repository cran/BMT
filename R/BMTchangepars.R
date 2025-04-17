#' @title The BMT Distribution Parameter Conversion.
#' @description Parameter conversion for different parameterizations for the BMT 
#'   distribution, with \code{p3} and \code{p4} tails weights (\eqn{\kappa_l} 
#'   and \eqn{\kappa_r}) or asymmetry-steepness parameters (\eqn{\zeta} and 
#'   \eqn{\xi}) and \code{p1} and \code{p2} domain (minimum and maximum) or 
#'   location-scale (mean and standard deviation) parameters.
#' 
#' @name BMTchangepars
#' @aliases BMTchangepars
#'   
#' @details The BMT coefficient of asymmetry \eqn{-1 < \zeta < 1} is 
#'   \deqn{\kappa_r - \kappa_l}
#'   
#'   The BMT coefficient of steepness \eqn{0 < \xi < 1} is \deqn{(\kappa_r + 
#'   \kappa_l - |\kappa_r - \kappa_l|) / (2 (1 -  |\kappa_r - \kappa_l|))} for 
#'   \eqn{|\kappa_r - \kappa_l| < 1}.
#'   
#'   The BMT distribution has mean \eqn{( d - c ) BMTmean(\kappa_l, \kappa_r) + 
#'   c} and standard deviation \eqn{( d - c ) BMTsd(\kappa_l, \kappa_r)}
#'   
#'   From these equations, we can go back and forth with each parameterization.
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
#' @return \code{BMTchangepars} reparametrize \code{p3}, \code{p4}, \code{p1}, 
#'   \code{p2} according to the alternative parameterizations from the given 
#'   \code{type.p.3.4} and \code{type.p.1.2}. \code{BMTchangepars} returns a 
#'   list with the alternative arguments to those received.
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
#' @references Torres-Jimenez, C. J. (2018), \emph{The BMT Item Response Theory 
#'   model: A new skewed distribution family with bounded domain and an IRT 
#'   model based on it}, PhD thesis, Doctorado en ciencias - Estadistica, 
#'   Universidad Nacional de Colombia, Sede Bogota.
#'   
#' @seealso \code{\link{BMT}} for the BMT density, distribution, quantile 
#'   function and random deviates.
#'   
#' @author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co} 
#'   and Alvaro Mauricio Montenegro Diaz [ths]
#'   
#' @examples
#' # BMT on [0,1] with left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' parameters <- BMTchangepars(0.25, 0.75, "t w")
#' parameters # Parameters of the BMT in the asymmetry-steepness parametrization
#' 
#' # BMT with mean equal to 0, standard deviation equal to 1, 
#' # asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.75
#' parameters <- BMTchangepars(0.5, 0.5, "a-s", 0, 1, "l-s")
#' parameters # Parameters of the BMT in the tail weight and domain parametrization

#' @rdname BMTchangepars
#' @export 
BMTchangepars <- function(p3, p4, type.p.3.4 = "t w", 
                          p1 = NULL, p2 = NULL, type.p.1.2 = NULL){
  # The length of the result is determined by the maximum of the lengths of the
  # numerical arguments. The numerical arguments are recycled to the length of
  # the result.
  if(is.null(p1) || is.null(p2)){
    len <- max(length(p3),length(p4))
    p3 <- rep(p3, len=len)
    p4 <- rep(p4, len=len)
  }
  else{
    len <- max(length(p1),length(p2),length(p3),length(p4))
    p1 <- rep(p1, len=len)
    p2 <- rep(p2, len=len)
    p3 <- rep(p3, len=len)
    p4 <- rep(p4, len=len)
  }
  # Control type.p.3.4
  TYPE.P.3.4 <- c("t w", "a-s")
  int.type.p.3.4 <- pmatch(type.p.3.4, TYPE.P.3.4)
  if (is.na(int.type.p.3.4)) 
    stop("invalid type of parametrization for parameters 3 and 4")
  if (int.type.p.3.4 == -1) 
    stop("ambiguous type of parametrization for parameters 3 and 4")
  # tail weights or asymmetry-steepness parametrization
  if(int.type.p.3.4 == 1){ # tail weights parametrization
    # Control tail weights parameters
    kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
    kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # Parameter conversion
    zeta <- BMTasymm(kappa_l, kappa_r, type.p.3.4)
    xi <- BMTsteep(kappa_l, kappa_r, type.p.3.4)
    #
    zeta <- replace(zeta, zeta > 1 & zeta < .one, 1)
    zeta <- replace(zeta, zeta < -1 & zeta > -.one, -1)
    xi <- replace(xi, xi > 1 & xi < .one, 1)
    xi <- replace(xi, xi < 0 & xi > .zero, 0)
    # 
    p <- list(p3=zeta, p4=xi, type.p.3.4="a-s")
    if(is.null(type.p.1.2)){
      return(p)
    }
  }
  else{ # asymmetry-steepness parametrization
    # Control asymmetry-steepness parameters
    zeta <- replace(p3, p3 < -1 | p3 > 1, NaN)
    xi <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # Parameter conversion
    kappa_r <- xi + abs(zeta)*(0.5 - xi) + 0.5*zeta
    kappa_l <- kappa_r - zeta
    #
    kappa_r <- replace(kappa_r, kappa_r > 1 & kappa_r < .one, 1)
    kappa_r <- replace(kappa_r, kappa_r < 0 & kappa_r > .zero, 0)
    kappa_l <- replace(kappa_l, kappa_l > 1 & kappa_l < .one, 1)
    kappa_l <- replace(kappa_l, kappa_l < 0 & kappa_l > .zero, 0)
    # 
    p <- list(p3=kappa_l, p4=kappa_r, type.p.3.4="t w")
    if(is.null(type.p.1.2)){
      return(p)
    }
  }
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
    # 
    p$p1 <- range*BMTmean(p3, p4, type.p.3.4) + min
    p$p2 <- range*BMTsd(p3, p4, type.p.3.4)
    p$type.p.1.2 <- "l-s"
  }
  else{ # location-scale parametrization
    # Control location-scale parameters
    mu <- p1
    sigma <- replace(p2, p2 <= 0, NaN)
    # range
    range <- sigma/BMTsd(p3, p4, type.p.3.4)
    # 
    p$p1 <- mu - range*BMTmean(p3, p4, type.p.3.4)
    p$p2 <- range + p$p1
    p$type.p.1.2 <- "c-d"
  }
  return(p)
}
