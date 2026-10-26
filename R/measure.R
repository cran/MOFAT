#' @title
#' Screening Measures
#'
#' @description
#' This function can be used for computing screening measures.
#'
#' @details
#' The \code{measure} function computes the screening measures such as the total Sobol' indices (Sobol' 1993)
#' and \eqn{\mu^*} measure of Campolongo et al. (2007). The design matrix should have the Sobol' design structure.
#' Please see Xiao et al. (2022) and Yu and Joseph (2026) for further details.
#'
#' @author Wei-Yang Yu, Difan Song, Qian Xiao, and V. Roshan Joseph
#'
#' @importFrom stats var
#'
#' @param design design matrix, which should have the Sobol' design structure
#' @param y response vector
#'
#' @return
#' \item{mustar}{\eqn{\mu^*} measure}
#' \item{t}{Total Sobol' index}
#'
#' @references
#'
#' Sobol’, I. M. (1993). On Sensitivity Estimation for Nonlinear Mathematical Models. Mathematical Modeling and Computational Experiments, 1, 407–414.
#'
#' Campolongo, F., Cariboni, J., and Saltelli, A. (2007). An Effective Screening Design for Sensitivity Analysis of Large Models. Environmental Modelling and Software, 22, 1509–1518.
#'
#' Xiao, Q., Joseph, V. R., and Ray, D. M. (2022). Maximum One-Factor-At-A-Time Designs for Screening in Computer Experiments. Technometrics, 65, 261--270.
#'
#' Yu, W.-Y. and Joseph, V. R. (2026). Space-Filling One-Factor-At-A-Time Designs. \doi{doi:10.48550/arXiv.2606.02533}
#'
#' @seealso \code{\link{mofat}}, \code{\link{mofatQQ}}, \code{\link{soft}}
#'
#' @export
#'
#' @examples
#' #Friedman function
#' fun <- function (X)
#' {
#'  Y <- 10*sin(pi*X[1]*X[2]) + 20*(X[3] - 0.5)^2 + 10*X[4] + 5*X[5]
#'  return(Y)
#' }
#' design = mofat(p=10, l=3)
#' y = apply(design, 1, fun)
#'
#' #Screening measures
#' measure(design, y)

# sensitivity measure for SOFT designs
measure <- function(design, y) {
  p <- ncol(design)
  l <- nrow(design) / (p + 1)
  design_ls <- lapply(1:(p + 1), function(i) {
    design[((i - 1) * l + 1):(i * l), , drop = FALSE]
  })
  y_ls <- split(y, rep(1:(p + 1), each = l))
  y_diff <- rep(list(NA), p)
  delta_ls <- rep(list(NA), p)
  for (i in 1:p) {
    A_info <- find_A(design_ls, i + 1, p)
    y_diff[[A_info$which_factor]] <- abs(
      y_ls[[i + 1]] - y_ls[[A_info$which_matrix]]
    )
    delta_ls[[A_info$which_factor]] <- abs(
      design_ls[[i + 1]][, A_info$which_factor] -
        design_ls[[A_info$which_matrix]][, A_info$which_factor]
    )
  }

  # mustar
  mustar <- lapply(1:p, function(i) {
    sum(y_diff[[i]] / delta_ls[[i]]) / l
  })
  mustar <- unlist(mustar)

  # t
  var_y <- var(y)
  vtot <- unlist(lapply(y_diff, function(i) {
    sum(i^2)
  })) /
    (2 * l)
  t <- vtot / var_y

  result <- list(mustar = mustar, t = t)
  return(result)
}
