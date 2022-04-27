#' @title BitPR
#' @aliases   dBitPR
#' @description probability density function of quotient of Bivariate t random variables conditioned to the positive quadrant.For more detailed information please read the first reference paper.
#' @author Saralees Nadarajah & Yuancheng Si \email{siyuanchengman@gmail.com}
#' @param x single positive scalar,for quotient of Bivariate t random variables conditioned to the positive quadrant
#' @param a  parameter for Bivariate t distribution
#' @param b  parameter for Bivariate t distribution
#' @param v  parameter, degree of freedom of Bivariate t distribution
#' @param rho correlation coefficient,\eqn{-1<\rho<1}
#' @references
#'
#' Yuancheng Si and Saralees Nadarajah and Xiaodong Song, (2020). On the distribution of quotient of random variables conditioned to the positive quadrant. \emph{Communications in Statistics - Theory and Methods}, \bold{49}, pp2514-2528.
#'
#' Balakrishnan, N. and Lai, C. -D. (2009).\emph{Continuous Bivariate Distributions}.Springer Verlag, New York.
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
#' \deqn{f_R (r \mid X > 0, Y > 0) =\frac {\Gamma \left( \frac {\nu + 2}{2} \right) \nu^{\frac {\nu}{2}}\left( 1 - \rho^2 \right)^{\frac {\nu + 1}{2}}}{\Gamma \left( \frac {\nu}{2} \right) \pi \Pr (X > 0, Y > 0)}J_1 \left( r^2 - 2 \rho r + 1, A r + B, C + \nu \left( 1 - \rho^2 \right),\frac {\nu}{2} + 1 \right)}
#'
#' For \eqn{-\infty < x < \infty},\eqn{-\infty < y < \infty,r > 0,-\infty < a < \infty,-\infty < b < \infty,-1 < \rho < 1},where \eqn{A = -2 a + 2 \rho b,B = -2 b + 2 \rho a,C = a^2 + b^2 - 2 \rho a b} and \eqn{J_1} is given by first reference paper section (2.5).
#'
#' @return \code{dBitPR} gives the probability density function for  quotient of Bivariate t random variables conditioned to the positive quadrant.
#' @return Invalid arguments will return an error message.
#'
#' @examples
#' x <- seq(0.1,5,0.1)
#' y <- c()
#' for (i in x){y=c(y,dBitPR(i,1,2,0.5,2))}
#' plot(x,y,type = 'l')
#'
#'


#' @rdname BitPR
#'
#' @export
#'
#'
#''
#'
#'
#'
dBitPR <- function(x, a, b, rho, v)
{
  stopifnot(x > 0, rho < 1 && rho > -1, v %% 1 == 0 && v >= 0)
  n = 2
  mean = rep(0 , 2)
  lower = rep(-Inf , 2)
  upper = c(-a , -b)
  corr = diag(n)
  corr[lower.tri(corr)] = rho
  corr[upper.tri(corr)] = rho
  delta = rep(0, 2)
  Survival00 = 1 - stats::pt(-a, v) - stats::pt(-b, v) + mvtnorm::pmvt(lower, upper, delta, v, corr)
  s = Survival00[1]
  #Pr(x>0,y>0)
  J_1 <- function(a, b, c, alpha){a^(-1) * c^(1 - alpha) * beta(2, 2 * alpha -2) *
      f21hyper(1, alpha - 1, alpha + 1 / 2, 1 - b^2 / (4 * a * c))}
  # J_1(a,b,c,alpha)
  A = -2 * a + 2 * rho * b
  B = -2 * b + 2 * rho * a
  C = -2 * rho * a * b + a^2 + b^2
  x1 = x^2 + 1 - 2 * x * rho
  x2 = A * x + B
  x3 = C + v * (1 - rho^2)
  x4 = v / 2 + 1
  #PARAMETERS FOR J_1
  pdf = gamma((v + 2) / 2) * v^(v / 2) * (1 - rho^2)^((v + 1) / 2)/ (gamma(v / 2) * v * pi * s) * J_1(x1, x2, x3, x4)
  return(pdf)
}

