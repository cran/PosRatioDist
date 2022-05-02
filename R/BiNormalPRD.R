#' @title BinormalPR
#' @aliases   dBinormalPR
#' @description probability density function of quotient of Bivariate normal random variables conditioned to the positive quadrant.For more detailed information please read the first reference paper.
#' @author Saralees Nadarajah & Yuancheng Si \email{siyuanchengman@gmail.com}
#' @param x vector of positive quantiles.
#' @param a  parameter
#' @param b  parameter
#' @param rho correlation coefficient,\eqn{-1<\rho<1}
#' @references
#'
#' Yuancheng Si and Saralees Nadarajah and Xiaodong Song, (2020). On the distribution of quotient of random variables conditioned to the positive quadrant. \emph{Communications in Statistics - Theory and Methods}, \bold{49}, pp2514-2528.
#'
#' Balakrishna, N. and Shiji, K. (2014). On a class of bivariate exponential distributions.\emph{Statistics and Probability Letters}, \bold{85}, pp153-160.
#'
#' Arnold, B. C. and Strauss, D. (1988).Pseudolikelihood estimation.\emph{Sankhya B} , \bold{53}, pp233-243.
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
#' \deqn{f_R (r \mid X > 0, Y > 0) =\frac {1}{2 \pi \sqrt{1 - \rho^2} \Pr (X > 0, Y > 0)}\exp \left[ -\frac {a^2 + b^2 - 2 \rho a b}{2 \left( 1 - \rho^2 \right)} \right]I_1 \left( \frac {1 + C r + r^2}{2 \left( 1 - \rho^2 \right)},\frac {A r + B}{2 \left( 1 - \rho^2 \right)} \right)}
#'
#' For \eqn{-\infty < x < \infty},\eqn{-\infty < y < \infty,r > 0,-\infty < a < \infty,-\infty < b < \infty,-1 < \rho < 1},where \eqn{A = -2 a + 2 \rho b,B = -2 b + 2 \rho a,C = -2 \rho}
#'
#' @return \code{dBinormalPR} gives the probability density function for  quotient of Bivariate normal random variables conditioned to the positive quadrant.
#' @return Invalid arguments will return an error message.
#'
#' @examples
#' x <- seq(0.1,5,0.1)
#' y <- dBinormalPR(x, 2, 1, 0.5)
#' plot(x,y,type = 'l')
#'
#'


#' @rdname BinormalPR
#'
#' @export
dBinormalPR <- function(x, a, b, rho)
{
  stopifnot(x > 0, rho < 1 && rho > -1)
  n = 2
  mean = rep(0 , 2)
  lower = rep(-Inf , 2)
  upper = c(-a , -b)
  corr = diag(n)
  corr[lower.tri(corr)] = rho
  corr[upper.tri(corr)] = rho
  Survival00 = 1 - stats::pnorm(-a) - stats::pnorm(-b) + mvtnorm::pmvnorm(lower, upper, mean, corr)
  s = Survival00[1]
  #Pr(x>0,y>0)
  I_1 <- function(a , b){- sqrt(pi) * b / (4 * a^(3 / 2)) * exp(b^2 / (4 * a)) * 2 * stats::pnorm(-sqrt(2) * (b / (2 * sqrt(a))) + 1  / (2 * a))}
  # I_1(a,b)
  A = -2 * a + 2 * rho * b
  B = -2 * b + 2 * rho * a
  C = -2 * rho
  x1 = (1 + C * x + x^2) / (2 * (1 - rho^2))
  x2 = (A * x + B) / (2 * (1 - rho^2))
  #PARAMETERS FOR I_1
  pdf = 1 / (2 * pi * sqrt(1 - rho^2 ) * s) * exp(- (a^2 + b^2 - 2 * rho * a * b) / (2 * (1 - rho^2))) * I_1(x1, x2)
  return(pdf)
  }

