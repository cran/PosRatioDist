#' @title BilomaxPR
#' @aliases   dBilomaxPR
#' @description probability density function of quotient of Bivariate Lomax random variables  conditioned to the positive quadrant.For more detailed information please read the first reference paper.
#' @author Saralees Nadarajah & Yuancheng Si \email{siyuanchengman@gmail.com}
#' @param x single positive scalar for quotient
#' @param alpha  parameter for Bivariate lomax distribution
#' @param beta  parameter for Bivariate lomax distribution
#' @param theta  parameter for Bivariate lomax distribution
#' @param a parameter for Bivariate lomax distribution
#' @param b parameter for Bivariate lomax distribution
#' @param c parameter for Bivariate lomax distribution
#'
#' @references
#'
#' Yuancheng Si and Saralees Nadarajah and Xiaodong Song, (2020). On the distribution of quotient of random variables conditioned to the positive quadrant. \emph{Communications in Statistics - Theory and Methods}, \bold{49}, pp2514-2528.
#'
#'
#' Balakrishnan, N. and Lai, C. -D. (2009).\emph{Continuous Bivariate Distributions}.Springer Verlag, New York.
#'
#'
#'@details
#'
#' Probability density function
#' \deqn{f_R (r \mid X > 0, Y > 0) = \frac {c^2 \theta^2 r}{\Pr (X > 0, Y > 0)} J_3 \left( \theta r, \beta - \theta a + \left( \alpha - \theta b \right) r, 1 - \alpha a - \beta b + \theta a b, c + 2 \right) +\frac {c^2 \theta \left[ (\alpha - \theta b) r + \beta - \theta a \right]} {\Pr (X > 0, Y > 0)} J_2 \left( \theta r, \beta - \theta a + \left( \alpha - \theta b \right) r, 1 - \alpha a - \beta b + \theta a b, c + 2 \right) +\frac {c \left[ c (\alpha - \theta b) (\beta - \theta a) + \alpha \beta - \theta \right]}{\Pr (X > 0, Y > 0)}J_1 \left( \theta r, \beta - \theta a + \left( \alpha - \theta b \right) r,1 - \alpha a - \beta b + \theta a b, c + 2 \right)}
#'
#' For \eqn{r > 0},\eqn{\alpha > 0}, \eqn{\beta > 0}, \eqn{\theta > 0}, \eqn{0 \leq \theta \leq (c + 1) \alpha \beta} where \eqn{J_1,J_2,J_3} are given by  first reference paper section (2.5)
#'
#' @return \code{dBilomaxPR} gives the probability density function for bivariate lomax random variables conditioned to the positive quadrant.
#' @return Invalid arguments will return an error message.
#'
#'
#'
#'
#'


#' @rdname dBilomaxPR
#'
#' @export
#'
#'
#'
#'
#'
#'
dBilomaxPR <- function(x, a, b, c, alpha, beta, theta)
{
  stopifnot(x > 0, alpha >0, beta > 0, theta >= 0 && theta <= (c + 1) * alpha * beta)
  s = (1 - alpha * a - beta * b + theta * a * b)^(-c)
  x1 = theta * x
  x2 = beta - theta * a + (alpha - theta * b) * x
  x3 = 1 - alpha * a - beta * b + theta * a * b
  x4 = c + 2

  J_1 <- function(a, b, c, alpha)
  {
    stopifnot(a > 0,b^2 < 4 * a * c,-1 < 1 && 1 < 2 * alpha - 1)
    a^(-1) * c^(1 - alpha) * beta(2, 2 * alpha -2) * f21hyper(1, alpha - 1, alpha + 1 / 2, 1 - b^2 / (4 * a * c))
  }
  # J_1(a,b,c,alpha)
  J_2 <- function(a, b, c, alpha)
  {
    stopifnot(a > 0,b^2 < 4 * a * c,-1 < 2 && 2 < 2 * alpha - 1)
    a^(-3 / 2) * c^(3 / 2 - alpha) * beta(3, 2 * alpha - 3) * f21hyper(3 / 2, alpha - 3 / 2, alpha + 1 / 2, 1 - b^2 / (4 * a * c))
  }
  # J_2(a,b,c,alpha)
  J_3 <- function(a, b, c, alpha)
  {
    stopifnot(a > 0,b^2 < 4 * a * c,-1 < 3 && 3 < 2 * alpha - 1)
    a^(-2) * c^(2 - alpha) * beta(4, 2 * alpha - 4) * f21hyper(2, alpha - 2, alpha + 1 / 2, 1 - b^2 / (4 * a * c))
  }
  # J_3(a,b,c,alpha)

  pdf = c^2 * theta^2 * x / s * J_3(x1,x2,x3,x4) + c^2 * theta * ((alpha - theta * b) * x + beta - theta * a) / s * J_2(x1,x2,x3,x4) + c * (c * (alpha - theta * b) * (beta - theta * a) + alpha * beta - theta) / s * J_1(x1,x2,x3,x4)
  return(pdf)
}

