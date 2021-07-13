#' @title BiexpweightedPR
#' @aliases   dBiexpweightedPR
#' @description probability density function of quotient of Bivariate exponential random variables resulting from weighted linear combinations conditioned to the positive quadrant.For more detailed information please read the first reference paper.
#' @author Saralees Nadarajah & Yuancheng Si \email{siyuanchengman@gmail.com}
#' @param x vector of positive quantiles.
#' @param a  parameter for Bivariate exponential random variables resulting from weighted linear combinations
#' @param b  parameter for Bivariate exponential random variables resulting from weighted linear combinations
#' @param c  parameter for Bivariate exponential random variables resulting from weighted linear combinations
#' @references
#'
#' Yuancheng Si and Saralees Nadarajah and Xiaodong Song, (2020). On the distribution of quotient of random variables conditioned to the positive quadrant. \emph{Communications in Statistics - Theory and Methods}, \bold{49}, pp2514-2528.
#'
#'
#' Balakrishnan, N. and Lai, C. -D. (2009).\emph{Continuous Bivariate Distributions}.Springer Verlag, New York.
#'
#' Caginalp, C. and Caginalp, G. (2018).The quotient of normal random variables and application to asset price fat tails.\emph{Physica A---Statistical Mechanics and Its Applications}, \bold{499}, pp457-471.
#'
#'Louzada, F., Ara, A. and Fernandes, G. (2017).The bivariate alpha-skew-normal distribution.\emph{Communications in Statistics - Theory and Methods}, \bold{46}, pp7147-7156.
#'
#'Nadarajah, S. (2009).A bivariate Pareto model for drought.\emph{Stochastic Environmental Research and Risk Assessment}, \bold{23}, pp811-822.
#'
#'Nadarajah, S. and Kotz, S. (2006).Reliability models based on bivariate exponential distributions.\emph{Probabilistic Engineering Mechanics}, \bold{21},  pp338-351.
#'
#'Nadarajah, S. and Kotz, S. (2007).Financial Pareto ratios.\emph{Quantitative Finance}, \bold{7}, pp257-260.
#'
#'
#'@details
#'
#' Probability density function
#' \deqn{f_R (r \mid X > 0, Y > 0) = \frac {(1 - 2 c) \exp \left[ (1 - 2 c) a + b \right]} {\Pr (X > 0, Y > 0) \left[ 1 + (1 - 2 c) r \right]^2}}
#'
#' For \eqn{x > a > -\infty},\eqn{y > b > -\infty,r > 0,0 < c < 1},These correlated exponential random variables can be used to model the stress and strength components of a system, hence the quotient distribution can be used to estimate the probability of failure of the system
#'
#' @return \code{dBiexpweightedPR} gives the probability density function for  quotient of Bivariate exponential random variables resulting from weighted linear combinations conditioned to the positive quadrant.
#' @return Invalid arguments will return an error message.
#'
#' @examples
#' x <- seq(0.1,5,0.1)
#' y <- dBiexpweightedPR(x, 4, 2, 0.2)
#' plot(x,y,type = 'l')
#'
#'


#' @rdname BiexpweightedPR
#'
#' @export
dBiexpweightedPR <- function(x, a, b, c)
{
  stopifnot(x > 0, 0 < c && c < 1)
  s = exp((1 - 2 * c) * a + b)
  #Pr(x>0,y>0)
  pdf = (1 - 2 * c) / (1 + (1 - 2 * c) * x)^2
  return(pdf)
}

