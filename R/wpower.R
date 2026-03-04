
#' Estimate p1, p2, and p3 from a location shift model by simulation
#'
#'Simulates independent samples from a chosen distribution under the model
#'\eqn{X} and \eqn{X'} are i.i.d., and \eqn{Y = X' + \theta}.
#'
#' Estimates:
#' \itemsize{
#'   \item \eqn{p_1 = P(Y >= X)}
#'   \item \eqn{p_2 = P(Y >= X_1 AND Y >= X_2)}
#'   \item \eqn{p_3 = P(Y_1 >= X AND Y_2 >= X)}
#'}
#'
#' @param theta Numeric Treatment Effect. Location Shift
#' @param rdist Function. Random generator for the error distribution. The default is \code{rnorm}
#' @param nsim Number of simulation
#'
#' @return A list with p1, p2, p3
#'
#' @examples
#' estimate_p123(theta = 0.7, rdist = rnorm)
#' estimate_p123(theta = 0.7, rdist= rlogis, nsim=5000)
#'


estimate_p123 <- function(theta, rdist, nsim = 30000) {

  #P1
  x <- rdist(nsim)
  y <- rdist(nsim) + theta
  p1 <- mean(y  > x)

  #P2 = P(Y > X1 and Y > X2)
  #
  x1 <- rdist(nsim)
  x2 <- rdist(nsim)
  Y <- rdist(nsim) + theta
  p2 <- mean((Y  > x1) & (Y  > x2))

  #P3
  y1 <- rdist(nsim) + theta
  y2 <- rdist(nsim) + theta
  X <-  rdist(nsim)
  p3 <- mean((y1 > X)  & (y2 > X))

  list(p1 = p1, p2 = p2, p3 = p3)
}






#' Estimate p1, p2, p3 from pilot data
#'
#' Given observed data from two groups, estimates p1, p2, p3 using
#' pairwise comparisons. This follows estimators described on page 363 of
#' in the Sample Size Calculations in Clinical Research Book.
#'
#' The implementation is based on the indicator form estimators.
#'
#' @param x Numeric Vector. Control group observations
#' @param y Numeric vector. Treatment group observations
#'
#' @return A list with p1, p2, p3
#'
#'
#'


estimate_p123_data <- function(x, y) {

  n1 <- length(x)
  n2 <- length(y)

  A <- rowSums(outer(y,x, ">="))
  B <- colSums(outer(y,x, ">="))

  p1 <- mean(outer(y, x, ">="))
  p2 <- sum(A *(A-1))/ (n1 * n2 *(n1-1))
  p3 <- sum(B *(B-1))/ (n1 * n2 *(n2-1))


  list(p1 = p1, p2 = p2, p3 = p3)
}





#' Noether (1987) sample size approximation for the Wilcoxon rank-sum test
#'
#' Computes the total sample size, using the Noether(1987) approximation
#' for the two sample Wilcoxon test in the continuous case. This method depends
#' only on the effect size.
#'
#'
#' @param alpha Numeric. significance level
#' @param beta Numeric. Type II error (1-Power)
#' @param p1 Numeric. Effect size probability
#' @param k Numeric. Allocation ratio. Default is 1
#'
#' @return Integer. Total Sample Size N
#'
#' @examples
#' noether_N(alpha = 0.05, beta=0.2, p1=0.7, k=1)

noether_N <- function(alpha, beta, p1, k = 1) {

  rho <- k / (k + 1)

  z_a <- qnorm(1 - alpha/2)
  z_b <- qnorm(1 - beta)

  N <- (z_a + z_b)^2 / (12 * rho * (1 - rho) * (p1 - 0.5)^2)
  ceiling(N)
}



#' Simulate power of the Wilcoxon rank-sum test
#'
#' Estimates the empirical power of \code{\link[stats]{wilcox.test}}
#' under a location shift model by Monte Carlo simulation.
#'
#'
#'
#' @param n1 Integer. Sample size for group 1 (control)
#' @param n2 Integer. Sample size for group 2 (treatment)
#' @param theta Numeric. Location Shift
#' @param alpha Numeric Significance level Default: 0.05
#' @param rdist Function. Random Generator. Default: rnorm
#' @param nsim Integer. Number of simulation replicates Default is 2000
#' @param alternative Character Type of alternative hypothesis
#' one of "two.sided", "less", "greater". Default: two.sided
#' @param exact logical Computes exact p-value. Default: FALSE
#' @param correct Logical. Whether to apply continuity correction.
#' Default: FALSE
#'
#'
#' @return Numeric. Estimated Power
#'
#'@examples
#' wilcox_power(n1 = 20, n2 = 20, theta = 0.7, rdist = rnorm, nsim = 1000)
#'





wilcox_power <- function(n1,n2,theta,alpha=0.05, rdist=rnorm, nsim=2000,
                         alternative = "two.sided",exact = FALSE,
                         correct = FALSE){


  # Power Simulation
  rejections <- replicate(nsim, {
    x <- rdist(n1)
    y <- theta + rdist(n2)


    p_val <- stats::wilcox.test(x, y, alternative = alternative,
                      exact = exact, correct = correct)$p.value

     p_val < alpha })

  mean(rejections)}



#' Compares Noether and TrialSize sample size calculations for Wilcoxon
#' rank sum test.
#'
#' Computes sample sizes using two approaches for the two sample Wilcoxon rank
#' sum test
#'
#'
#'
#' The function then estimates achieved power for both designs by simulation
#' under the specified location shift \code{theta} and distribution \code{rdist}.
#'
#'
#' @param theta Numeric. Location shift (treatment effect)
#' @param alpha Numeric. Significance Level. Default is 0.05
#' @param power Numeric. Target power. Default is 0.8
#' @param k Numeric. Allocation ration \eqn{k=n_1/n_2}. Default is 1
#'
#' @param p1, p2, p3 Optional numeric. If supplied used as effect size inputs
#' for the TrialSize calculations.
#'
#' @param rdist Function. Random generator. Default rnorm
#' @param pilot_x, pilot_y Optional numeric vectors. If supplied used to
#' estimate \eqn{p_1,p_2,p_3} from pilot data
#'
#' @param alternative Character. Passed to \code{\link[stats]{wilcox.test}}.
#' @param exact Logical Passed to \code{\link[stats]{wilcox.test}}.
#' @param correct Logical. Passed to \code{\link[stats]{wilcox.test}}.
#'
#' @param nsim_est Integer Number of simulations for estimating \eqn{p_1,p_2,p_3}
#' @param nsim_pow Integer Number of simulations for estimating power
#' @param ... Additional arguments passed to \code{\link[stats]{wilcox.test}}.
#'
#' @return An object of class \code{"wpower}: a list with components
#' \describe{
#'   \item{trialsize}{list with \code{n1}, \code{n2}, \code{N_total}, and simulated \code{power}.}
#'   \item{noether}{list with \code{n1}, \code{n2}, \code{N_total}, and simulated \code{power}.}
#'   \item{theta, alpha, target_power, k}{Input Settings}
#'   \item{rdist}{Generator function used for simulation}
#'   \item{alternative, exact, correct}{Wilxcoxon test options used in simulation }
#' }
#'
#' @examples
#' # Example 1: Simulated estimated of p1,p2,p3 under a normal shift model
#' Example1 <- wpower(theta = 0.7, rdist = rnorm, nsim_est = 50000, nsim_pow = 1000)
#' Example1$trialsize
#' Example1$noether
#'
#' @export
wpower <- function(theta, alpha = 0.05, power = 0.80, k = 1,
                 p1 = NULL, p2 = NULL, p3 = NULL,
                 rdist = rnorm,
                 pilot_x = NULL, pilot_y = NULL,
                 alternative = "two.sided",
                 exact = FALSE, correct = FALSE,
                 nsim_est = 100000, nsim_pow = 5000, ...) {

  beta <- 1-power


  # p1 p2 p3
  if (!is.null(pilot_x) && !is.null(pilot_y)) {
    # Estimate from pilot data

    probs <- estimate_p123_data(pilot_x, pilot_y)
    if (is.null(p1)) p1 <- probs$p1
    if (is.null(p2)) p2 <- probs$p2
    if (is.null(p3)) p3 <- probs$p3

  } else if (is.null(p1) || is.null(p2) || is.null(p3)) {
    # Estimate missing values from simulation

    probs <- estimate_p123(theta, rdist, nsim = nsim_est)
    if (is.null(p1)) p1 <- probs$p1
    if (is.null(p2)) p2 <- probs$p2
    if (is.null(p3)) p3 <- probs$p3
  }

  # Sample Sizes
  # Trial Size returns n2
  n2_ts <- ceiling(Nonpara.Two.Sample(alpha, beta, k,p1,p2,p3))
  n1_ts <- ceiling(k * n2_ts)

  # Noether: returns N
  N_no <- noether_N(alpha, beta, p1, k)
  rho <- k/(k+1)
  n1_no <- round(N_no * rho)
  n2_no <- N_no - n1_no

  # Power Simulation
  power_ts <- wilcox_power(n1_ts, n2_ts, theta, alpha, rdist, nsim_pow,
                           alternative = alternative,
                           exact = exact, correct = correct)

  power_no <- wilcox_power(n1_no, n2_no, theta, alpha, rdist, nsim_pow,
                           alternative = alternative,
                           exact = exact, correct = correct)

  ## Result Object
  result <- list(
    trialsize = list(n1 = n1_ts, n2=n2_ts,
                     N_total = n1_ts+n2_ts, power = power_ts),

    noether = list(n1 = n1_no, n2=n2_no,
                   N_total = N_no, power = power_no),

    theta = theta, alpha = alpha, target_power = power, k=k,
    p1=p1, p2=p2, p3=p3,
    rdist = rdist,
    alternative = alternative, exact = exact, correct = correct
  )


  class(result) <- "wpower"
  return(result)

}
