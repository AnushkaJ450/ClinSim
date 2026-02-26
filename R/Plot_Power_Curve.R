################################################################################
# power_curve.R
#
#  generating power curves across a range of
# theta values, comparing Noether and TrialSize methods.
################################################################################


#' Plot power curves across a range of effect sizes
#'

plot_power_curve <- function(theta_range, alpha = 0.05, power = 0.80,
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

results <- plot_power_curve(theta_range = seq(0.3, 1.2, by = 0.1))
results2 <- plot_power_curve(theta_range = c(0.4, 0.6, 0.8, 1.0, 1.2),nsim_pow=3000)
