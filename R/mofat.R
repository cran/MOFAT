#' @title
#' MOFAT
#'
#' @description
#' This function can be used for generating MOFAT designs.
#'
#' @details
#' The \code{mofat} function generates the MOFAT design
#' for a given number of factors (\eqn{p\ge 2}) and
#' number of base runs (\eqn{l \ge 3}). The total number of runs in the MOFAT design will be \eqn{l(p+1)}.
#' A MOFAT design can be viewed as an optimized version of Morris screening design (Morris 1991) by exploiting
#' its connections with the Monte Carlo-based design of Sobol' (1993).
#' Please see Xiao et al. (2022) for details.
#'
#' Three choices for the \code{method} are given: "uniform", "projection", and "best". Option "uniform" gives \code{l} equally-spaced levels
#' for the entire design, which are also balanced. "projection" option adjusts the levels of the two base matrices A and B such that
#' there are \eqn{2l} or \eqn{2l-1} levels in the design depending on \code{l} is even or odd. Option "best" (default) chooses the best
#' among the first two options using maximin distance criterion.
#'
#' @author Qian Xiao and V. Roshan Joseph
#'
#' @importFrom SLHD maximinSLHD
#' @importFrom stats dist
#' @importFrom utils combn
#'
#' @param p number of factors
#' @param l number of base runs
#' @param method choose among "uniform", "projection", and "best"
#' @param nstart number of random starts used when \code{l > 3}
#' @param init initialization method for the base matrix when \code{l > 3};
#' choose among \code{"rLHD"}, \code{"MaxPro"}, and \code{"MmLHD"}
#'
#' @return
#' \item{design}{MOFAT design}
#'
#' @references
#'
#' Morris, M. D. (1991), “Factorial sampling plans for preliminary computational experiments,”
#' Technometrics, 33, 161–174.
#'
#' Sobol’, I. M. (1993), “On sensitivity estimation for nonlinear mathematical models,” Mathematical
#' Modeling and Computational Experiments, 1, 407–414.
#'
#' Xiao, Q., Joseph, V. R., and Ray, D. M. (2022). “Maximum One-Factor-At-A-Time  Designs for Screening in Computer Experiments”.
#' Technometrics, 65, 261--270.
#'
#' @export
#'
#' @examples
#' #MOFAT with three base runs
#' mofat(p=10, l=3, method="uniform")
#' mofat(p=10, l=3, method="projection")
#'
#' #MOFAT with five base runs
#' mofat(p=10,l=5)
#' dim(mofat(p=125,l=5))

#library(SLHD)

###########################################
####MOFAT design with l levels

mofat <- function(p, l, method = "best", nstart = 1, init = "MaxPro") {
  #p is the factor size; l is the level size.
  if (
    length(p) != 1 || !is.numeric(p) || is.na(p) || p < 2 || p != as.integer(p)
  ) {
    stop("p must be a single integer greater than or equal to 2.")
  }
  if (
    length(l) != 1 || !is.numeric(l) || is.na(l) || l < 3 || l != as.integer(l)
  ) {
    stop("l must be a single integer greater than or equal to 3.")
  }
  if (
    length(nstart) != 1 ||
      !is.numeric(nstart) ||
      is.na(nstart) ||
      nstart < 1 ||
      nstart != as.integer(nstart)
  ) {
    stop("nstart must be a single integer greater than or equal to 1.")
  }
  if (!method %in% c("uniform", "projection", "best")) {
    stop("Invalid method. Choose from 'uniform', 'projection', or 'best'.")
  }
  if (!init %in% c("rLHD", "MaxPro", "MmLHD")) {
    stop(
      "Invalid initialization method. Choose from 'rLHD', 'MaxPro', or 'MmLHD'."
    )
  }

  #number of runs needed
  n <- l * (p + 1)

  ##################function needed
  criterion <- function(D) {
    p <- ncol(D)
    l <- nrow(D) / (p + 1)
    pairs <- combn(p, 2)
    out <- mapply(
      function(i, j) {
        mean(
          1 /
            (dist(
              D[
                c(1:l, (i * l + 1):((i + 1) * l), (j * l + 1):((j + 1) * l)),
                c(i, j)
              ]
            )^2 +
              1 / l^2)
        )
      },
      pairs[1, ],
      pairs[2, ]
    )
    mean(out)
  }

  #function for level permutation
  permut <- function(n) {
    if (n == 1) {
      return(matrix(1))
    } else {
      sp <- permut(n - 1)
      p <- nrow(sp)
      A <- matrix(nrow = n * p, ncol = n)
      for (i in 1:n) {
        A[(i - 1) * p + 1:p, ] <- cbind(i, sp + (sp >= i))
      }
      return(A)
    }
  }

  #generate B matrix based on A matrix
  B.gen <- function(A, l) {
    B <- A
    for (j in 1:dim(A)[2]) {
      B[, j] <- (rank(A[, j] - .5 + as.numeric(A[, j] < .5)) - 0.5) / l
    }
    return(B)
  }

  #generate design based on A and B
  sobol_design <- function(A, B) {
    d <- dim(A)[2]
    res <- A
    for (i in 1:d) {
      interm <- A
      interm[, i] <- B[, i]
      res <- rbind(res, interm)
    }
    return(res)
  }

  adjust <- function(A, B) {
    A0 <- A
    B0 <- B
    A[A < .5] <- A[A < .5] - .25 / l
    A[A >= .5] <- A[A >= .5] + .25 / l
    B[B < .5] <- B[B < .5] + .25 / l
    B[B >= .5] <- B[B >= .5] - .25 / l
    D1 <- sobol_design(A, B)
    d1 <- min(dist(D1, method = "manhattan"))

    A <- A0
    B <- B0
    A[A < .5] <- A[A < .5] + .25 / l
    A[A >= .5] <- A[A >= .5] - .25 / l
    B[B < .5] <- B[B < .5] - .25 / l
    B[B >= .5] <- B[B >= .5] + .25 / l
    D2 <- sobol_design(A, B)
    d2 <- min(dist(D2, method = "manhattan"))
    if (d2 < d1) {
      design <- D1
    } else {
      design <- D2
    }
    return(design)
  }

  choose_design <- function(A, B) {
    design <- sobol_design(A, B)
    d0 <- min(dist(design, method = "manhattan"))
    if (method == "projection") {
      design <- adjust(A, B)
    }
    if (method == "best") {
      design_adj <- adjust(A, B)
      d_adj <- min(dist(design_adj, method = "manhattan"))
      if (d_adj > d0) design <- design_adj
    }
    return(design)
  }

  randomize_design <- function(design) {
    randomness <- sample(p)
    design[
      c(1:l, sapply(randomness, function(x) x * l + 1:l)),
      randomness
    ]
  }

  generate_base_matrix <- function(k) {
    if (k == 0) {
      return(NULL)
    }
    if (k == 1) {
      return(matrix(((1:l) - 0.5) / l, ncol = 1))
    }
    if (init == "rLHD") {
      return(
        (vapply(seq_len(k), function(i) sample.int(l), integer(l)) - 0.5) / l
      )
    }
    if (init == "MaxPro") {
      return(as.matrix(MaxPro::MaxProLHD(l, k)$Design))
    }
    if (init == "MmLHD") {
      return(maximinSLHD(t = 1, m = l, k = k)$StandDesign)
    }
  }

  choose_base_matrix <- function(generator) {
    if (nstart > 1) {
      crit <- numeric(nstart)
      seed <- stats::runif(1, 0, 1e4)
      for (i in seq_len(nstart)) {
        A <- withr::with_seed(seed + i, generator())
        B <- B.gen(A, l)
        design <- choose_design(A, B)
        crit[i] <- criterion(design)
      }
      return(withr::with_seed(seed + which.min(crit), generator()))
    }
    generator()
  }

  ###Mofat with l=3 case:
  if (l == 3) {
    ###when the number of factors k is no larger than l!
    if (p <= 6) {
      #creation of matrix P_l
      allperm <- t(permut(l))
      A <- allperm[, 1:p]
      ###standardize to (0,1) range.
      A <- (A - 0.5) / l
      B <- B.gen(A, l)
      design <- choose_design(A, B)
      return(randomize_design(design))
    } else {
      #creation of matrix P_l
      allperm <- t(permut(l))
      #find the best A from Algebraic construction
      if (p %% 6 == 0) {
        A <- matrix(
          rep(allperm, floor(p / 6)),
          l,
          floor(p / 6) * dim(allperm)[2]
        )
      } else {
        A <- matrix(
          rep(allperm, floor(p / 6)),
          l,
          floor(p / 6) * dim(allperm)[2]
        )
        A <- cbind(A, allperm[, 1:(p %% 6)])
      }
      ###standardize to (0,1) range.
      A <- (A - 0.5) / l
      B <- B.gen(A, l)
      design <- choose_design(A, B)
      return(randomize_design(design))
    }
  } else {
    ###when the number of factors k is no larger than l!
    if (p <= factorial(l)) {
      A <- choose_base_matrix(function() generate_base_matrix(p))

      #matrix B
      B <- B.gen(A, l)
      design <- choose_design(A, B)
      return(design)
    } else {
      #creation of matrix P_l
      allperm <- t(permut(l))
      #column number in allperm
      cn <- dim(allperm)[2]
      #rescale allperm
      allperm <- (allperm - 0.5) / l

      A <- choose_base_matrix(function() {
        A_p <- generate_base_matrix(p %% cn)
        A <- matrix(rep(allperm, floor(p / cn)), l, floor(p / cn) * cn)
        if (!is.null(A_p)) {
          A <- cbind(A, A_p)
        }
        A
      })
      #matrix B
      B <- B.gen(A, l)
      design <- choose_design(A, B)
      return(design)
    }
  }
}
