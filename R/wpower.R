##############################################################
# ESTIMATION OF p1, p2, p3
##############################################################

#' Estimate p1, p2, p3 from a location shift model via simulation
#'
#'
#' p1 = P(Y >= X)
#' p2 = P(Y >= X1 AND Y >= X2)
#' p3 = P(Y1 >= X AND Y2 >= X)
#'
#'
#' @param theta Numeric Treatment Effect. Location Shift
#' @param rdist Function. Default: rnorm
#' @param nsim Number of simulation
#' @return A list with p1, p2, p3
#'
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
###########################################################

#' Estimate p1, p2, p3 from pilot data
#'
#' Given observed data from two groups, estimates p1, p2, p3 using
#' pairwise comparisons. This follows estimators described on page 363 of
#' T=Book
#'
#'
#'
#' @param x Numeric Vector. Control group observations
#' @param y Numeric vector. Treatment group observations
#' @return A list with p1, p2, p3


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



##################################################
# SAMPLE SIZE NOETHER
##################################################

#' Noether (1987) sample size formula
#'
#' simpler only needs p1
#' @param alpha significance level
#' @param beta Type II error (1-Power)
#' @param p1
#' @return Total Sample Size N
#'
#'

noether_N <- function(alpha, beta, p1, k = 1) {

  rho <- k / (k + 1)

  z_a <- qnorm(1 - alpha/2)
  z_b <- qnorm(1 - beta)

  N <- (z_a + z_b)^2 / (12 * rho * (1 - rho) * (p1 - 0.5)^2)
  ceiling(N)
}


##############################################################
# POWER SIMULATION
##############################################################

#' Simulate power of the Wilcoxon rank-sum test
#'
#'
#'
#' @param n1 Integer. Sample size for group 1
#' @param n2 Integer. Sample size for group 2
#' @param thera Numeric. Location Shift
#' @param alpha Numeric Significance level Default: 0.05
#' @param rdist Default: rnorm
#' @param nsim Integer. Number of simulation replicates Default
#' @param alternative Character Type of alternative hypothesis
#' one if "two.sided", "less", "greater". Default: two.sided
#' @param exact logical Compuyer exact p-value. Default: FALSE
#' @param correct Logical. Whether to apply continuity correction.
#' Default: FALSE
#'
#'
#' @return Numeric, Estimated Power
#'
#'


wilcox_power <- function(n1,n2,theta,alpha=0.05, rdist=rnorm, nsim=2000,
                         alternative = "two.sided",exact = FALSE,
                         correct = FALSE){


  # Power Simulation
  rejections <- replicate(nsim, {
    x <- rdist(n1)
    y <- theta + rdist(n2)

    p_val <- stats::wilcox.test(x, y,
                         alternative = alternative,
                         exact = exact,
                         correct = correct)$p.value


    p_val < alpha
  })

  mean(rejections)}



#####################################################################
# MAIN FUNCTION
#####################################################################


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

















