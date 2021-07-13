#' @title Bibs_expPR
#' @aliases   dBibs_expPR
#' @description probability density function of quotient of Balakrishna and Shiji's bivariate exponential random variables  conditioned to the positive quadrant.For more detailed information please read the first reference paper.
#' @param x vector of positive quantiles.
#' @param a  parameter for  Balakrishna and Shiji's bivariate exponential distribution
#' @param r  parameter for  Balakrishna and Shiji's bivariate exponential distribution
#' @references
#'
#' Yuancheng Si and Saralees Nadarajah and Xiaodong Song, (2020). On the distribution of quotient of random variables conditioned to the positive quadrant. \emph{Communications in Statistics - Theory and Methods}, \bold{49}, pp2514-2528.
#'
#'
#' Balakrishnan, N. and Lai, C. -D. (2009).\emph{Continuous Bivariate Distributions}.Springer Verlag, New York.
#'
#' Balakrishna, N. and Shiji, K. (2014).On a class of bivariate exponential distributions.\emph{Statistics and Probability Letters}, \bold{85}, pp153-160.
#'
#'@details
#'
#' Probability density function
#' \deqn{f_R (r \mid X > 0, Y > 0) = \frac {a}{2 \sqrt{r}} \left( r + \frac {a^2}{4 r} \right)^{-3 / 2}}
#'
#' For \eqn{r > 0},\eqn{a > 0}
#'
#' @return \code{dBibs_expPR} gives the probability density function for  quotient of Balakrishna and Shiji's bivariate exponential random variables conditioned to the positive quadrant.
#' @return Invalid arguments will return an error message.
#'
#' @examples
#'
#'x <- seq(0.1,5,0.1)
#'y <- dBibs_expPR(x, 2, 2)
#'plot(x,y,type = 'l')


#' @rdname Bibs_expPR
#'
#' @export
dBibs_expPR <- function(x, a, r)
{
  stopifnot(x > 0, r >0, a > 0)
  pdf = a / (2 * sqrt(x)) * (x + a^2 / (4 * x))^(-3 / 2)
  return(pdf)
}

