#' Validate TrialSize two-sample mean sample size using simulation
#'
#' @description
#' Computes the required sample size for a two-sample mean equality design using
#' TrialSize: TwoSampleMean.Equality() and then checks the achieved power
#' empirically via Monte Carlo simulation using a two-sample t-test with equal
#' variances.
#'
#' This function is intended as a simulation-based check of the analytical
#' calculation, returning both the target (theoretical) power and the simulated
#' (empirical) power.
#'
#' @details
#' Data are generated under a normal model with common standard deviation
#' sigma. Group 1 outcomes are generated from \eqn{N(0, \sigma^2)} and
#' group 2 outcomes from \eqn{N(\mathrm{margin}, \sigma^2)}, where
#' \code{margin = mu2 - mu1}.
#'
#' The simulation applies \code{\link[stats]{t.test}} with \code{var.equal = TRUE}
#' and rejects \eqn{H_0: \mu_1 = \mu_2} when the p-value is below \code{alpha_test}.
#'
#' @param alpha Numeric. Significance level used in the TrialSize sample size
#'   calculation (typically 0.05).
#' @param beta Numeric. Type II error rate used in the TrialSize calculation,
#'   where target power is \code{1 - beta}.
#' @param sigma Numeric. Common standard deviation (must be positive).
#' @param k Numeric. Allocation ratio \eqn{k = n1/n2}. For equal allocation use
#'   \code{k = 1}.
#' @param margin Numeric. True mean difference \eqn{(\mu_2 - \mu_1)} under the
#'   alternative used for the power calculation.
#' @param alpha_test Numeric. Significance level used for the t-test in the
#'   simulation. Default is 0.05. (Often set equal to \code{alpha}.)
#' @param n_sim Integer. Number of simulations. Default is 10000.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{n1}: Required sample size for group 1 from TrialSize.
#'   \item \code{n2}: Corresponding sample size for group 2 implied by \code{k}.
#'   \item \code{theoretical_power}: Target power \code{1 - beta}.
#'   \item \code{empirical_power}: Simulated power estimate (proportion of rejections).
#' }
#'
#' @examples
#' # Example: check if TrialSize sample size achieves ~80% power
#' # (requires TrialSize package)
#' # res <- two_sample_mean_B(alpha = 0.05, beta = 0.2, sigma = 1, k = 1, margin = 0.3,
#' #                          alpha_test = 0.05, n_sim = 2000)
#' # res$empirical_power
#'
#' @export


two_sample_mean_B <-  function(alpha, beta, sigma, k, margin,
                               alpha_test = 0.05, n_sim = 10000){


  n1 <- TwoSampleMean.Equality(alpha, beta, sigma, k, margin)
  n2 <- n1/k


  rejections <- replicate(n_sim, {

    g1 <- rnorm(n1, 0, sigma)
    g2 <- rnorm(n2, margin, sigma)

    p_val <- t.test(g1, g2, var.equal = TRUE)$p.value
    p_val < alpha_test

  })

  empirical_power <- mean(rejections)


  theoretical_power <- 1 - beta
  list(
    n1 = n1,
    n2 = n2,
    theoretical_power = theoretical_power,
    empirical_power = empirical_power
  )

}







