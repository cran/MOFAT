############################################ simple functions for soft.R
# Transform the levels of "design" to "levels". Here, "design" means the A0 matrix (matrix of the l base runs).
#' @noRd
level_trans <- function(design, l, levels) {
  b <- design
  b[design < ((l+1)/2)] <- 2*design[design < ((l+1)/2)]-1
  b[design > ((l+1)/2)] <- 2*design[design > ((l+1)/2)]
  design <- matrix(levels[b], nrow = l)
  return(design)
}

# Generate representative vertex points. Here, "design" means the Monte Carlo samples (qmc).
#' @noRd
approx_corner <- function(design) {
  a <- design
  a[design < 0.5] <- 1
  a[design > 0.5] <- 0
  return(a)
}

#' Evaluation Set for the Q Criterion
#'
#' Generate Sobol points and augment them with corner points to evaluate
#' the \eqn{Q} criterion. Duplicate points are removed from the resulting
#' evaluation set.
#' See Yu and Joseph (2026) for further details.
#' 
#' @author Wei-Yang Yu and V. Roshan Joseph
#'
#' @param n number of Sobol points.
#' @param p number of factors.
#' @param threshold dimension cutoff for corner augmentation. When
#'   \code{p <= threshold}, all \eqn{2^p} binary corners are added. Otherwise,
#'   representative corners derived from the Sobol points are added.
#'
#' @return A matrix with \code{p} columns containing the unique Sobol
#'   and corner evaluation points. The number of rows may differ from \code{n}
#'   because corner points are added and duplicate rows are removed.
#'
#' @references
#' Yu, W.-Y. and Joseph, V. R. (2026). Space-Filling One-Factor-At-A-Time Designs. \doi{10.48550/arXiv.2606.02533}
#'
#' @seealso \code{\link{Q}}
#'
#' @export
#'
#' @examples
#' qmc <- qmc_generate(n = 1000, p = 2)

qmc_generate <- function(n,p,threshold = 10) {
  qmc <- spacefillr::generate_sobol_set(n = n, dim = p, seed = 0)
  qmc <- matrix(qmc, ncol = p)
  if(p <= threshold) {
    qmc <- rbind(qmc, SFDesign::full.factorial(p,2))
  } else {
    qmc <- rbind(qmc, approx_corner(qmc))
  }
  qmc <- unique(qmc)
  return(qmc)
}

# Generate A0 matrix. 
#' @noRd
A_fun <- function(p,l) {
  L <- factorial(l)
  if(p < L) {
    A <- matrix(SFDesign::maximinLHD(l,p)$design, ncol = p)
    A <- A*l + 0.5
    a <- NULL
    b <- NULL
  } else {
    a <- p%/%(L)
    b <- p%%L
    P_l <- t(gtools::permutations(l,l,1:l))
    A1 <- do.call(cbind, replicate(a, P_l, simplify = FALSE))
    if(b > 0) {
      A2 <- P_l[,1:b]
    } else {
      A2 = NULL
    }
    A <- cbind(A1, A2)
  }
  ans <- list(A = A, L = L, quotient = a, mod = b)
  return(ans)
}

############################################ simple functions for measure.R
# Given A(i), find A(oi). 
#' @noRd
find_A <- function(X_ls, idx, p) {
  B <- X_ls[[idx]]
  for(j in seq_len(idx-1)) {
    A <- X_ls[[j]]
    col_same <- vapply(seq_len(p), function(k) { # check each column one by one and return a vector of TRUE and FALSE.
      isTRUE(all.equal(A[, k], B[, k], check.attributes = FALSE))
    }, logical(1))
    if(sum(!col_same) == 1) { # if exactly one column differs, return.
      return(list(which_matrix = j, which_factor = which(!col_same))) # Compared to X_ls[[idx]], which factor is changed from which matrix.
    }
  }
}
