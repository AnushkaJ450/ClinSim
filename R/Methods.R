################################################################
# METHODS
#   - print.wpower()   Short output od results
#   - summary.wpower() Detailed Summary
################################################################


#' Print method for
#'
#' Short Summary of key results
#' @parmam x A clinsim object
#' @param ... Additional arguments


print.wpower <- function(x,...){
  cat("Wilcoxon-Mann-Whitney Sample Size/Power Calculation\n")
  cat("Call: ")
  print(x$call)
  cat("\n")
  cat(sprintf("  theta = %.3f | alpha = %.3f | target power = %.2f | k = %.1f\n",
              x$theta, x$alpha, x$target_power, x$k))


  cat(sprintf("  Noether:    n1 = %d, n2 = %d, total = %d, power = %.3f\n",
              x$noether$n1, x$noether$n2, x$noether$N_total, x$noether$power))

  cat(sprintf("  TrialSize:  n1 = %d, n2 = %d, total = %d, power = %.3f\n",
              x$trialsize$n1, x$trialsize$n2, x$trialsize$N_total, x$trialsize$power))

  invisible(x)
}



#' Summary method for
#'
#' Detailed summary Method with explanation of results
#'
#'
#' @param object A wpower object
#' @param ... additional arguments
#'
#'
#'

summary.wpower<- function(object, ...) {

  cat("========================================================\n")
  cat("  Wilcoxon-Mann-Whitney Sample Size & Power Calculation\n")
  cat("========================================================\n\n")

  cat("Call:\n  ")
  print(object$call)
  cat("\n")

  # Settings
  cat("Settings:\n")
  cat(sprintf("  Significance level (alpha):  %.3f\n", object$alpha))
  cat(sprintf("  Target power (1 - beta):     %.2f\n", object$target_power))
  cat(sprintf("  Location shift (theta):      %.4f\n", object$theta))
  cat(sprintf("  Allocation ratio (k=n1/n2):  %.1f\n", object$k))
  cat(sprintf("  wilcox.test settings:        alternative = %s, exact = %s, correct = %s\n",
              object$wt_args$alternative,
              as.character(object$wt_args$exact),
              as.character(object$wt_args$correct)))
  cat(sprintf("  Estimation replicates:       %d\n", object$nsim_est))
  cat(sprintf("  Power simulation replicates: %d\n\n", object$nsim_pow))

  cat("\n")
  # Probability parameters
  cat("Probability parameters:\n")
  cat(sprintf("  p1 = P(Y > X)             = %.4f\n", object$p1))
  cat(sprintf("  p2 = P(Y > X1, Y > X2)    = %.4f\n", object$p2))
  cat(sprintf("  p3 = P(Y1 > X, Y2 > X)    = %.4f\n", object$p3))

  cat("\n")
  # Results table
  cat("Results:\n")
  cat("  Method       n1    n2    Total   Simulated Power\n")
  cat("  -------------------------------------------------\n")

  print_row <- function(method,x){
    cat(sprintf(
      "  %-12s %-5d %-5d %-7d %.3f\n",
      method, x$n1, x$n2, x$N_total, x$power))

  }

  print_row("Noether",   object$noether)
  print_row("TrialSize", object$trialsize)

  cat("\n")

  # Interpretation
  diff_n <- object$noether$N_total - object$trialsize$N_total
  cat("Interpretation:\n")
  cat(sprintf("  Noether recommends %d %s subjects than TrialSize.\n",
              abs(diff_n), if (diff_n > 0) "more" else "fewer"))

  no_ok <- object$noether$power >= (object$target_power - 0.02)
  ts_ok <- object$trialsize$power >= (object$target_power - 0.02)

  if (no_ok) {
    cat("  Noether achieves approximately the target power.\n")
  } else {
    cat(sprintf("  Noether power (%.3f) is below the target (%.2f).\n",
                object$noether$power, object$target_power))
  }
  if (ts_ok) {
    cat("  TrialSize achieves approximately the target power.\n")
  } else {
    cat(sprintf("  TrialSize power (%.3f) is below the target (%.2f).\n",
                object$trialsize$power, object$target_power))
  }


  invisible(object)
}




