#' @title BiMG_expPR
#' @aliases   dBiMG_expPR
#' @description probability density function of quotient of Morgenstern type bivariate exponential random variables  conditioned to the positive quadrant.For more detailed information please read the first reference paper.
#' @author Saralees Nadarajah & Yuancheng Si \email{siyuanchengman@gmail.com}
#' @param x vector of positive quantiles.
#' @param a  parameter for Morgenstern type bivariate exponential distribution
#' @param b  parameter for Morgenstern type bivariate exponential distribution
#' @param alpha parameter for Morgenstern type bivariate exponential distribution
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
#' \deqn{f_R (r \mid X > 0, Y > 0) = \frac  {(1 + \alpha) \exp (a + b)}{\Pr (X > 0, Y > 0) (1 + r)^2} - \frac  {2 \alpha \exp (a + 2 b)}{\Pr (X > 0, Y > 0) (2 + r)^2} - \frac  {2 \alpha \exp (2 a + b)}{\Pr (X > 0, Y > 0) (1 + 2 r)^2} + \frac  {\alpha \exp (2 a + 2 b)}{\Pr (X > 0, Y > 0) (1 + r)^2}}
#'
#' For \eqn{r > 0},\eqn{-1 \leq \alpha \leq 1,  a > -\infty,  b > -\infty } These correlated exponential random variables can also be used to model the stress and strength components of a system, hence the quotient distribution can be used to estimate the probability of failure of the system
#'
#' @return \code{dBiMG_expPR} gives the probability density function for  quotient of Morgenstern type bivariate exponential random variables  conditioned to the positive quadrant
#' @return Invalid arguments will return an error message.
#'
#'
#' @examples
#' x <- seq(0.1,5,0.1)
#' y <- dBiMG_expPR(x, 3, 2, 0.5)
#' plot(x,y,type = 'l')
#'
#'


#' @rdname BiMG_expPR
#'
#' @export
dBiMG_expPR <- function(x, a, b, alpha)
{
  stopifnot(x > 0, alpha >= -1 && alpha <= 1)
  s = exp(a + b) * (1 + alpha * (1 - exp(a))) * (1 - exp(b))
  pdf = (1 + alpha) * exp(a + b) / (s * (1 + x)^2) - 2 * alpha * exp(2 * a + b) / (s * (1 + 2 * x)^2) - 2 * alpha * exp(2 * b + a) / (s * (2 + x)^2) + alpha * exp(2 * a + 2 * b) / (s * (1 + x)^2)
  return(pdf)
}

