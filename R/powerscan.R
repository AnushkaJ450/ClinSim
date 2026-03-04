#' Power and sample size curves across a range of effect sizes
#'
#' Runs \code{\link{wpower}} over a grid of \code{theta} values and returns a
#' data frame comparing Noether and TrialSize sample sizes and simulated power.
#' Produces two base plots: power vs theta and total sample size vs theta.
#'
#' @param theta_range Numeric vector of effect sizes (\eqn{\theta}) to evaluate
#' @param alpha Significance Level. Default is 0.05
#' @param power Target power. Default is 0.8
#' @param k Allocation ratio.
#' @param rdist Random generator.
#'
#' @param alternative Passed to \code{\link[stats]{wilcox.test}}.
#' @param exact Passed to \code{\link[stats]{wilcox.test}}.
#' @param correct Passed to \code{\link[stats]{wilcox.test}}.
#'
#' @param nsim_est Integer Number of simulations for estimating \eqn{p_1,p_2,p_3}
#' @param nsim_pow Integer Number of simulations for estimating power.
#'
#' @examples
#' results <- powerscan(theta_range = seq(0.3, 1.2, by = 0.1))
#' results2 <- powerscan(theta_range = c(0.4, 0.6, 0.8, 1.0, 1.2),nsim_pow=3000)
#'@export
#'
#'

powerscan <- function(theta_range, alpha = 0.05, power = 0.80,
                             k = 1, rdist = rnorm,
                             alternative = "two.sided",
                             exact = FALSE, correct = FALSE,
                             nsim_est = 100000, nsim_pow = 5000) {


  results <- data.frame(
    theta = numeric(),
    p1 = numeric(), p2 = numeric(), p3 = numeric(),
    N_noether = numeric(), power_noether = numeric(),
    N_trialsize = numeric(), power_trialsize = numeric()
  )

  # Run for each theta
  for (th in theta_range) {

    res <- wpower(theta = th, alpha = alpha, power = power,
                              k = k, rdist = rdist,
                              alternative = alternative,
                              exact = exact, correct = correct,
                              nsim_est = nsim_est, nsim_pow = nsim_pow)

    results <- rbind(results, data.frame(
      theta = th,
      p1 = res$p1, p2 = res$p2, p3 = res$p3,
      N_noether = res$noether$N_total,
      power_noether = res$noether$power,
      N_trialsize = res$trialsize$N_total,
      power_trialsize = res$trialsize$power
    ))
  }

  # Plot
  oldpar <- par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  on.exit(par(oldpar))

  # Power curves
  plot(results$theta, results$power_noether, type = "b", pch = 16,
       col = "navy", ylim = c(0, 1),
       xlab = expression(theta), ylab = "Simulated Power",
       main = "Power vs Effect Size")
  lines(results$theta, results$power_trialsize, type = "b", pch = 17,
        col = "magenta")
  abline(h = power, lty = 2, col = "red")
  legend("bottomright", legend = c("Noether", "TrialSize"),
         col = c("navy", "magenta"), pch = c(16, 17), lty = 1)

  # Sample size curves
  plot(results$theta, results$N_noether, type = "b", pch = 16,
       col = "navy",
       xlab = expression(theta), ylab = "Total N",
       main = "Sample Size vs Effect Size")
  lines(results$theta, results$N_trialsize, type = "b", pch = 17,
        col = "magenta")
  legend("topright", legend = c("Noether", "TrialSize"),
         col = c("navy", "magenta"), pch = c(16, 17), lty = 1)


  invisible(results)
}

