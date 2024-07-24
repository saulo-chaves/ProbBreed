#' The 'ProbBreed' package.
#'
#' @description
#' ProbBreed uses probability theory under the Bayesian framework for calculating
#' the risk of selecting candidates in a multi-environment context.
#' Contained are functions used to fit a Bayesian multi-environment model
#' (based on the available presets), extract posterior values and maximum posterior values,
#' compute the variance components, check the model’s convergence, and calculate the probabilities.
#' For both across and within-environments scopes, the package computes the probability of superior performance and the pairwise probability of superior performance.
#' Furthermore, the probability of superior stability and the pairwise probability of superior stability across environments is estimated.
#' @docType package
#' @name ProbBreed-package
#' @aliases ProbBreed
#' @useDynLib ProbBreed, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import ggplot2
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rlang .data
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.32.6. https://mc-stan.org
#'
#' Dias, K. O. G, Santos J. P. R., Krause, M. D., Piepho H. -P., Guimarães, L. J. M.,
#' Pastina, M. M., and Garcia, A. A. F. (2022). Leveraging probability concepts
#'  for cultivar recommendation in multi-environment trials. \emph{Theoretical and
#' Applied Genetics}, 133(2):443-455. \doi{10.1007/s00122-022-04041-y}
#'
#' @keywords internal
"_PACKAGE"
## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL
