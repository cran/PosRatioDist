#' @title BicauchyPR
#' @aliases   dBicauchyPR
#' @description probability density function of quotient of Bivariate cauchy random variables conditioned to the positive quadrant.For more detailed information please read the first reference paper.
#' @author Saralees Nadarajah & Yuancheng Si \email{siyuanchengman@gmail.com}
#' @param x single real positive scalar
#' @param a  parameter for bivaraite cauchy distribution
#' @param b  parameter for bivaraite cauchy distribution
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
#' \deqn{f_R (r \mid X > 0, Y > 0) =\frac {1}{2 \pi \Pr (X > 0, Y > 0)}J_1 \left( r^2 + 1, A r + B, C, \frac {3}{2} \right)}
#'
#' For \eqn{-\infty < x < \infty},\eqn{-\infty < y < \infty,r > 0,-\infty < a < \infty,-\infty < b < \infty},where \eqn{A = -2 a, B = -2 b,C = 1 + a^2 + b^2} and \eqn{J_1} is given by first reference paper section (2.5).
#'
#' @return \code{dBicauchyPR} gives the probability density function for  quotient of Bivariate cauchy random variables conditioned to the positive quadrant.
#' @return Invalid arguments will return an error message.
#'
#' @examples
#' x <- seq(0.1,5,0.1)
#' y <- c()
#' for (i in x){y=c(y,dBicauchyPR(i,1,2))}
#' plot(x,y,type = 'l')
#'
#'
#'
#'

#'
#'
#' @rdname BicauchyPR
#'
#' @export
#'
dBicauchyPR <- function(x, a, b)
{
  stopifnot(x > 0)
  s = 1 / 4 + 1 / (2 * pi) * (-atan(-a) - atan(-b) + atan(a * b / (sqrt(1 + a^2 + b^2))))
  #Pr(x>0,y>0)
  J_1 <- function(a, b, c, alpha)
  {
    stopifnot(a > 0,b^2 < 4 * a * c,-1 < 1 && 1 < 2 * alpha - 1)
    a^(-1) * c^(1 - alpha) * beta(2, 2 * alpha -2) * f21hyper(1, alpha - 1, alpha + 1 / 2, 1 - b^2 / (4 * a * c))
  }
  # J_1(a,b,c,alpha)
  A = -2 * a
  B = -2 * b
  C = 1 + a^2 + b^2
  x1 = x^2 + 1
  x2 = A * x + B
  x3 = C
  x4 = 3 / 2
  #PARAMETERS FOR J_1
  pdf = 1 / (2 * pi * s) * J_1(x1, x2, x3, x4)
  return(pdf)
}

