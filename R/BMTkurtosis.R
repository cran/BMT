#' @title The BMT Distribution Descriptive Measures - Kurtosis.
#' @description Kurtosis and steepness coefficient for the BMT distribution with
#'   \code{p3} and \code{p4} tails weights (\eqn{\kappa_l} and \eqn{\kappa_r})
#'   or asymmetry-steepness parameters (\eqn{\zeta} and \eqn{\xi}) and \code{p1}
#'   and \code{p2} domain (minimum and maximum) or location-scale (mean and 
#'   standard deviation) parameters.
#' 
#' @name BMTkurtosis
#' @aliases BMTkurt
#' @aliases BMTsteep
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
#' @return \code{BMTkurt} gives the Pearson's kurtosis and \code{BMTsteep} the 
#'   proposed steepness coefficient for the BMT distribution.
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
#' @seealso \code{\link{BMTcentral}}, \code{\link{BMTdispersion}}, 
#'   \code{\link{BMTskewness}}, \code{\link{BMTmoments}} for other descriptive 
#'   measures or moments.
#'   
#' @author Camilo Jose Torres-Jimenez [aut,cre] \email{cjtorresj@unal.edu.co}
#'   
#' @examples
#' # BMT on [0,1] with left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' BMTkurt(0.25, 0.75, "t w")
#' BMTsteep(0.25, 0.75, "t w")
#' 
#' # BMT on [0,1] with asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.75
#' BMTkurt(0.5, 0.5, "a-s")
#' BMTsteep(0.5, 0.5, "a-s")
#' 
#' # domain or location-scale parameters do not affect 
#' # the skewness and the asymmetry coefficient
#' 
#' # BMT on [-1.783489,3.312195] with 
#' # left tail weight equal to 0.25 and 
#' # right tail weight equal to 0.75
#' BMTkurt(0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' BMTsteep(0.25, 0.75, "t w", -1.783489, 3.312195, "c-d")
#' 
#' # BMT with mean equal to 0, standard deviation equal to 1, 
#' # asymmetry coefficient equal to 0.5 and 
#' # steepness coefficient equal to 0.75
#' BMTkurt(0.5, 0.5, "a-s", 0, 1, "l-s")
#' BMTsteep(0.5, 0.5, "a-s", 0, 1, "l-s")

#' @rdname BMTkurtosis
#' @export 
BMTkurt <- function(p3, p4, type.p.3.4 = "t w", 
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
  # tail weights or asymmetry-steepness parametrization
  if(int.type.p.3.4 == 1){ # tail weights parametrization
    # Control tail weights parameters
    kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
    kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # kurtosis
    m <- ((((6507*kappa_l - 43380)*kappa_l + 135900)*kappa_l - 150000)*kappa_l + 
          (((6507*kappa_r - 43380)*kappa_r + 135900)*kappa_r - 150000)*kappa_r + 
          ((432*kappa_l - 28620 + 13122*kappa_r)*kappa_l + 
          (432*kappa_r - 28620)*kappa_r + 29700)*kappa_l*kappa_r + 125125) / 
         (10010000*BMTvar(kappa_l,kappa_r,type.p.3.4)^2)
  }
  else{ # asymmetry-steepness parametrization
    # Control asymmetry-steepness parameters
    zeta <- replace(p3, p3 < -1 | p3 > 1, NaN)
    xi <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # kurtosis
    abs.zeta <- abs(zeta)
    m <- ((((((((27000*xi-54000)*xi+53460)*xi-26460)*xi+6507)*abs.zeta+ 
            ((((-108000*xi+306000)*xi-322920)*xi+185220)*xi-43380))*abs.zeta+
            ((((162000*xi-594000)*xi+786960)*xi-460260)*xi+135900))*abs.zeta+
            ((((-108000*xi+486000)*xi-819000)*xi+601500)*xi-150000))*abs.zeta+
            ((((27000*xi-144000)*xi+301500)*xi-300000)*xi+125125)) /
         (10010000*BMTvar(zeta,xi,type.p.3.4)^2)
  }
  # domain or location-scale parameters do not affect the kurtosis
  return(m)
}

#' @rdname BMTkurtosis
#' @export
BMTsteep <- function(p3, p4, type.p.3.4 = "t w", 
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
  # tail weights or asymmetry-steepness parametrization
  if(int.type.p.3.4 == 1){ # tail weights parametrization
    # Control tail weights parameters
    kappa_l <- replace(p3, p3 < 0 | p3 > 1, NaN)
    kappa_r <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # steepness coefficient
    abs.dif <- abs(kappa_r - kappa_l)
    m <- ifelse(abs.dif == 1, 1, (kappa_r + kappa_l - abs.dif) / (2*(1 - abs.dif)))
  }
  else{
    # Control asymmetry-steepness parameters
    xi <- replace(p4, p4 < 0 | p4 > 1, NaN)
    # steepness coefficient
    m <- xi
  }
  # domain or location-scale parameters do not affect the steepness coefficient
  return(m)
}
