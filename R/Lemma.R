#' @title Lemma
#' @aliases  I_1 I_2 I_3 J_1 J_2 J_3
#' @description Technical Lemmas for calculating quotient of  random variables conditioned to the positive quadrant.For more detailed information please read the first reference paper section 2.2.
#' @author Saralees Nadarajah & Yuancheng Si \email{siyuanchengman@gmail.com}
#' @param a  parameter
#' @param b  parameter
#' @param c parameter
#' @param alpha parameter
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
#' \eqn{I_n} Type I Integration
#' \deqn{I_n (a, b) = \int_0^\infty y^n \exp \left( -a y^2 - b y \right) dy}
#'
#' For \eqn{-\infty < a < \infty,-\infty < b < \infty},where n is positive integer.
#'
#' In particular,for \eqn{a > 0},we have expressions below
#' \deqn{I_1  (a, b) = -\frac {\sqrt{\pi} b}{4 a^{3 / 2}}\exp \left( \frac {b^2}{4 a} \right)\ {\rm erfc} \left( \frac {b}{2 \sqrt{a}} \right)  + \frac {1}{2 a}}
#' \deqn{I_2  (a, b) = \frac {\sqrt{\pi}}{4 a^{3 / 2}}\exp \left( \frac {b^2}{4 a} \right)\ {\rm erfc} \left( \frac {b}{2 \sqrt{a}} \right) +\frac {\sqrt{\pi} b^2}{8 a^{5 / 2}} \exp \left( \frac {b^2}{4 a} \right)\ {\rm erfc} \left( \frac {b}{2 \sqrt{a}} \right) - \frac {b}{4 a^2}}
#' \deqn{I_3  (a, b) = -\frac {3 \sqrt{\pi} b}{8 a^{5 / 2}}\exp \left( \frac {b^2}{4 a} \right)\ {\rm erfc} \left( \frac {b}{2 \sqrt{a}} \right) -\frac {\sqrt{\pi} b^3}{16 a^{7 / 2}}\exp \left( \frac {b^2}{4 a} \right)\ {\rm erfc} \left( \frac {b}{2 \sqrt{a}} \right) + \frac {1}{2 a^2} + \frac {b^2}{8 a^3}}
#'
#'
#' \eqn{J_n} Type J Integration
#' \deqn{J_n  (a, b, c, \alpha) = \int_0^\infty y^n \left( a y^2 + b y + c \right)^{-\alpha} dy}
#' In particular,for \eqn{a > 0,b^2 < 4ac, -1 < n < 2\alpha - 1},we have expressions below
#'
#' \deqn{J_1  (a, b, c, \alpha) = a^{-1} c^{1 - \alpha} B \left( 2, 2 \alpha - 2 \right) \ {}_2F_1 \left( 1, \alpha - 1; \alpha + \frac {1}{2}; 1 - \frac {b^2}{4 a c} \right)}
#' \deqn{J_2  (a, b, c, \alpha) = a^{-\frac {3}{2}} c^{\frac {3}{2} - \alpha} B \left( 3, 2 \alpha - 3 \right) \ {}_2F_1 \left( \frac {3}{2}, \alpha - \frac {3}{2}; \alpha + \frac {1}{2}; 1 - \frac {b^2}{4 a c} \right)}
#' \deqn{J_3  (a, b, c, \alpha) =  a^{-2} c^{2 - \alpha} B \left( 4, 2 \alpha - 4 \right) \ {}_2F_1 \left( 2, \alpha - 2; \alpha + \frac {1}{2}; 1 - \frac {b^2}{4 a c} \right)}
#'
#'
#'
#'
#' @return \code{I_1} gives value of Type I integration with \eqn{n = 1}
#' @return \code{I_2} gives value of Type I integration with \eqn{n = 2}
#' @return \code{I_3} gives value of Type I integration with \eqn{n = 3}
#' @return \code{J_1} gives value of Type J integration with \eqn{n = 1}
#' @return \code{J_2} gives value of Type J integration with \eqn{n = 2}
#' @return \code{J_3} gives value of Type J integration with \eqn{n = 3}
#'
#' @return Invalid arguments will return an error message.
#' @examples
#' I_1(1,2)
#' I_2(1,2)
#' I_3(1,2)
#' J_1(1,2,3,3)
#' J_2(1,2,3,3)
#' J_3(1,2,3,3)
#'
#'
#'
#'
#'
#'
#'
#'
#'


#' @rdname Lemma
#' @export
I_1 <- function(a , b)
{
  stopifnot(a > 0)
  - sqrt(pi) * b / (4 * a^(3 / 2)) * exp(b^2 / (4 * a)) * 2 * stats::pnorm(-sqrt(2) * (b / (2 * sqrt(a)) + 1  / (2 * a)))
}
# I_1(a,b)


#' @rdname Lemma
#' @export
I_2 <- function(a , b)
{
  stopifnot(a > 0)
   sqrt(pi) * b / (4 * a^(3 / 2)) * exp(b^2 / (4 * a)) * 2 * stats::pnorm(-sqrt(2) * (b / (2 * sqrt(a)))) + sqrt(pi) * b^2 / (8 * a^(5 / 2)) * exp(b^2 / (4 * a)) * 2 * stats::pnorm(-sqrt(2) *(b / (2 * sqrt(a)))) - b / (4 * a^2)
}
# I_2(a,b)


#' @rdname Lemma
#' @export
I_3 <- function(a , b)
{
  stopifnot(a > 0)
  -3 * sqrt(pi) * b / (8 * a^(5 / 2)) * exp(b^2 / (4 * a)) * 2 * stats::pnorm(-sqrt(2) *(b / (2 * sqrt(a)))) - sqrt(pi) * b^3 / (16 * a^(7 / 2)) * exp(b^2 / (4 * a)) * 2 * stats::pnorm(-sqrt(2) *(b / (2 * sqrt(a)))) + 1 / (2 * a^2) + b^2 / (8 * a^3)
}
# I_3(a,b)


#' @rdname Lemma
#' @export
J_1 <- function(a, b, c, alpha)
{
  stopifnot(a > 0,b^2 < 4 * a * c,-1 < 1 && 1 < 2 * alpha - 1)
  a^(-1) * c^(1 - alpha) * beta(2, 2 * alpha -2) * f21hyper(1, alpha - 1, alpha + 1 / 2, 1 - b^2 / (4 * a * c))
}
# J_1(a,b,c,alpha)


#' @rdname Lemma
#' @export
J_2 <- function(a, b, c, alpha)
{
  stopifnot(a > 0,b^2 < 4 * a * c,-1 < 2 && 2 < 2 * alpha - 1)
  a^(-3 / 2) * c^(3 / 2 - alpha) * beta(3, 2 * alpha - 3) * f21hyper(3 / 2, alpha - 3 / 2, alpha + 1 / 2, 1 - b^2 / (4 * a * c))
}
# J_2(a,b,c,alpha)


#' @rdname Lemma
#' @export
J_3 <- function(a, b, c, alpha)
{
  stopifnot(a > 0,b^2 < 4 * a * c,-1 < 3 && 3 < 2 * alpha - 1)
  a^(-2) * c^(2 - alpha) * beta(4, 2 * alpha - 4) * f21hyper(2, alpha - 2, alpha + 1 / 2, 1 - b^2 / (4 * a * c))
}
# J_3(a,b,c,alpha)

