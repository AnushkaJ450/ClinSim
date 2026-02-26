#' Simulates Two Sample t-test for Power Estimation
#'
#'
#' @param alpha Significance level
#' @param margin (µ2 - µ1)
#' @param sigma SD
#' @param n1 sample size in group 1
#' @param n2 Sample size in group 2
#' @param n_sum Number of simulations
#'
#' @return Estimated Power
#'
#'



simulate_two_sample <- function(n1, n2, margin, sigma, alpha = 0.05, n_sim=10000 ){

  rejections <- replicate(n_sim, {

    g1 <- rnorm(n1, mean =0, sd = sigma)
    g2 <- rnorm(n2, mean = margin, sd= sigma)

    p_val <- t.test(g1, g2, var.equal = TRUE)$p.value
    p_val < alpha

  })


  mean(rejections)


}


















