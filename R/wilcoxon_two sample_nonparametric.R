#' @param n1 Sample size in group 1
#' @param n2 Sample size in group 2
#' @param theta Treatment effect  (0 = no effect, >0 = treatment effect)
#' @param alpha Significance level
#' @param n_sim Number of simulation replications


simulate_wilcoxon <- function(n1, n2, theta, alpha = 0.05, n_sim = 2000){

  rejections <- replicate(n_sim, {
    x <- rnorm(n1, 0, 1)  # control group
    y <- rnorm(n2, mean = theta, sd = 1) # Treatment group
    p <- wilcox.test(x, y, exact = FALSE)$p.value
    p < alpha
  })

  mean(rejections)
}
