#' @title f21hyper
#' @aliases   f21hyper
#' @description Computes the value of a Gaussian hypergeometric function \eqn{ F(a,b,c,z) } for \eqn{-1 \leq z \leq 1} and \eqn{a,b,c \geq 0}
#' @author Martin Feldkircher and Stefan Zeugner
#' @param a The parameter \code{a} of the Gaussian hypergeometric function, must be a positive scalar here
#' @param b  The parameter \code{b} of the Gaussian hypergeometric function, must be a positive scalar here
#' @param c  The parameter \code{c} of the Gaussian hypergeometric function, must be a positive scalar here
#' @param z  The parameter \code{z} of the Gaussian hypergeometric function, must be between -1 and 1 here
#'
#'
#' @references
#'
#'   Liang F., Paulo R., Molina G., Clyde M., Berger J.(2008): Mixtures of g-priors for Bayesian variable selection. J. Am. Statist. Assoc. 103, p. 410-423
#'
#' Yuancheng Si and Saralees Nadarajah and Xiaodong Song, (2020). On the distribution of quotient of random variables conditioned to the positive quadrant. \emph{Communications in Statistics - Theory and Methods}, \bold{49}, pp2514-2528.
#'
#'
#' Saralees Nadarajah and Y.Si (2020) A note on the “L-logistic regression models: Prior sensitivity analysis, robustness to outliers and applications”. \emph{Brazilian Journal of Probability and Statistics},\bold{34},p. 183-187.
#'
#'
#' @details
#' The function \code{f21hyper} complements the analysis of the 'hyper-g prior' introduced by Liang et al. (2008)..
#'
#' @returns the value of the Gaussian hypergeometric function \eqn{ F(a,b,c,z)
#' @returns Invalid arguments will return an error message.
#'
#' @examples
#' f21hyper(30,1,20,.8) #returns about 165.8197
#' f21hyper(30,10,20,0) #returns one
#' f21hyper(10,15,20,-0.1) # returns about 0.4872972
#'
#'
#'
#'
#'
#' @rdname f21hyper
#'
#' @export




f21hyper <-
  function (a, b, c, z)
  {
    if ((length(a) != 1) | (length(b) != 1) | (length(c) != 1) |
        (length(z) != 1))
      stop("All function arguments need to be scalars")
    if ((a < 0) | (b < 0) | (c < 0))
      stop("Arguments a, b, and c need to be non-negative")
    if ((z > 1) | (z <= (-1)))
      stop("Argument z needs to be between -1 and 1")
    nmx = max(100, 3 * floor(((a + b) * z - c - 1)/(1 - z)))
    if (nmx > 10000)
      warning("Power series probably does not converge")
    serie = 0:nmx
    return(1 + sum(cumprod((a + serie)/(c + serie) * (b + serie)/(1 + serie) * z)))
  }
