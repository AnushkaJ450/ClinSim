#' ClinSim: Sample size checking by simulation
#'
#' @description clinsim provides helper functions for planning and checking
#' sample size calculations using simulation. It includes a non parametric two
#' sample Wilcoxon rank-sum workflow and a simpler parametric example, designed
#' for teaching and demonstration.
#'
#' @details The main functions are:
#' Main functions:
#'  - `mean2power()` — parametric two-sample mean sample size + simulation power check
#'  - `wpower()` — Wilcoxon rank-sum sample size (Noether vs TrialSize) + simulation power check
#'  - `powerscan()` — run `wpower()` across a grid of effect sizes and plot power / sample size curves
#'
#'
#' The following example datasets are provided.
#'
#' @docType package
#' @aliases ClinSim-package
#' @importFrom graphics plot
"_PACKAGE"
