#' @title BiparetoPR
#' @aliases   dBiparetoPR
#' @description probability density function of quotient of Bivariate Pareto random variables  conditioned to the positive quadrant.For more detailed information please read the first reference paper.
#' @author Saralees Nadarajah & Yuancheng Si \email{siyuanchengman@gmail.com}
#' @param x vector of positive quantiles.
#' @references
#'
#' Yuancheng Si and Saralees Nadarajah and Xiaodong Song, (2020). On the distribution of quotient of random variables conditioned to the positive quadrant. \emph{Communications in Statistics - Theory and Methods}, \bold{49}, pp2514-2528.
#'
#'
#' Mardia, K. V. (1962).Multivariate Pareto distributions.\emph{Annals of Mathematical Statistics}, \bold{33}, 1008-1015.
#'
#' Nadarajah, S. (2009) A bivariate Pareto model for drought.\emph{Stochastic Environmental Research and Risk Assessment}, \bold{23}, pp811-822.
#'
#'
#'
#'@details
#'
#' Probability density function
#' \deqn{f_R (r \mid X > 0, Y > 0) = (r + 1)^{-2}}
#'
#' For \eqn{r > 0},Nadarajah (2009) used this distribution to model the proportion of droughts defined as a quotient of drought durations and non-drought durations.
#'
#' @return \code{dBiparetoPR} gives the probability density function for  quotient of Bivariate Pareto random variables conditioned to the positive quadrant.
#' @return Invalid arguments will return an error message.
#'
#' @examples
#'
#' x <- seq(0.1,5,0.1)
#' y <- dBiparetoPR(x)
#' plot(x,y,type = 'l')
#'


#' @rdname BiparetoPR
#'
#' @export
dBiparetoPR <- function(x)
{
  stopifnot(x > 0)
  pdf = (1 + x)^(-2)
  return(pdf)
}

