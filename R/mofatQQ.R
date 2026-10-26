#' @title
#' MOFAT Designs for Quantitative and Qualitative Factors
#'
#' @description
#' This function generates MOFAT-based designs for continuous, discrete
#' numeric, and nominal factors.
#'
#' @details
#' The \code{mofatQQ} function first generates a base MOFAT design with
#' \code{\link{mofat}}, then transforms the columns according to factor type.
#' The total number of factors must be at least two.
#' Continuous factors are kept on the unit scale or rescaled to user-supplied
#' bounds. Discrete numeric factors are mapped either to integer levels or to
#' user-supplied numeric levels. Nominal factors are encoded by integer level
#' labels, and the nominal columns are placed to the right of the numeric
#' columns in the returned design.
#' Each discrete numeric and nominal factor must have at least two levels.
#'
#' When \code{nstart > 1}, the function generates multiple candidate mixed
#' designs and selects the one with the smallest projected space-filling
#' criterion. In that criterion, continuous and discrete numeric distances are
#' measured on a rescaled numeric scale, while nominal differences are treated
#' as binary mismatches.
#'
#' @author Difan Song and V. Roshan Joseph
#'
#' @param l number of base runs
#' @param p_continuous number of continuous factors
#' @param b_continuous optional list of length \code{p_continuous}; each
#' element should be a length-two numeric vector giving the lower and upper
#' bounds used to rescale a continuous factor
#' @param p_disnum number of discrete numeric factors
#' @param l_disnum specification of the discrete numeric levels. Supply either
#' a numeric vector giving the number of levels for each factor, or a list whose
#' \eqn{j}th element contains the numeric levels for factor \eqn{j}
#' @param p_nom number of nominal factors
#' @param l_nom specification of the nominal levels. Supply either a numeric
#' vector giving the number of levels for each nominal factor, or a list whose
#' element lengths define the number of nominal levels. The returned design
#' always encodes nominal factors as integer level labels
#' @param nstart number of candidate mixed-factor designs to generate
#' @param ... additional arguments passed to \code{\link{mofat}}, such as
#' \code{method} or \code{init}
#'
#' @return
#' A design matrix with columns ordered as continuous, discrete numeric, and
#' nominal factors. Nominal factors are returned as integer level labels.
#'
#' @references
#'
#' Xiao, Q., Joseph, V. R., and Ray, D. M. (2022). Maximum One-Factor-At-A-Time Designs for Screening in Computer Experiments. Technometrics, 65, 261--270.
#'
#' Song, D., and Joseph, V. R. (2026). Efficient Screening Designs for Expensive Black-Box Models with Qualitative and Quantitative Factors. \doi{10.48550/arXiv.2609.02087}
#'
#' @seealso \code{\link{mofat}}, \code{\link{measure}}
#'
#' @export
#'
#' @examples
#' # Mixed continuous and nominal factors
#' D <- mofatQQ(
#'   l = 5,
#'   p_continuous = 2,
#'   b_continuous = list(c(-1, 1), c(0, 10)),
#'   p_nom = 1,
#'   l_nom = 3,
#'   nstart = 2,
#'   init = "rLHD"
#' )
#' dim(D)
mofatQQ <- function(
  l,
  p_continuous = 0,
  b_continuous = NULL,
  p_disnum = 0,
  l_disnum = NULL,
  p_nom = 0,
  l_nom = NULL,
  nstart = 1,
  ...
) {
  if (
    length(l) != 1 || !is.numeric(l) || is.na(l) || l < 3 || l != as.integer(l)
  ) {
    stop("l must be a single integer greater than or equal to 3.")
  }
  for (arg in c("p_continuous", "p_disnum", "p_nom")) {
    value <- get(arg)
    if (
      length(value) != 1 ||
        !is.numeric(value) ||
        is.na(value) ||
        value < 0 ||
        value != as.integer(value)
    ) {
      stop(paste0(arg, " must be a single nonnegative integer."))
    }
  }
  if (p_continuous + p_disnum + p_nom < 2) {
    stop("At least two factors must be specified.")
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
  if (!is.null(b_continuous)) {
    if (!is.list(b_continuous) || length(b_continuous) != p_continuous) {
      stop("Input error: Length of b_continuous does not match p_continuous.")
    }
    for (j in seq_len(p_continuous)) {
      if (
        !is.numeric(b_continuous[[j]]) ||
          length(b_continuous[[j]]) != 2 ||
          anyNA(b_continuous[[j]])
      ) {
        stop(
          "Each element of b_continuous must be a numeric vector of length 2."
        )
      }
    }
  }
  if (p_disnum > 0) {
    if (is.numeric(l_disnum)) {
      if (length(l_disnum) != p_disnum) {
        stop("Input error: Length of l_disnum does not match p_disnum.")
      }
      if (
        anyNA(l_disnum) ||
          any(l_disnum < 2) ||
          any(l_disnum != as.integer(l_disnum))
      ) {
        stop(
          "Each discrete numeric factor must have at least two integer levels."
        )
      }
    } else if (is.list(l_disnum)) {
      if (length(l_disnum) != p_disnum) {
        stop("Input error: Length of l_disnum does not match p_disnum.")
      }
      if (any(vapply(l_disnum, function(x) length(x) < 2, logical(1)))) {
        stop("Each discrete numeric factor must have at least two levels.")
      }
    } else {
      stop("Input error: l_disnum must be either numeric or a list.")
    }
  }
  if (p_nom > 0) {
    if (is.numeric(l_nom)) {
      if (length(l_nom) != p_nom) {
        stop("Input error: Length of l_nom does not match p_nom.")
      }
      if (anyNA(l_nom) || any(l_nom < 2) || any(l_nom != as.integer(l_nom))) {
        stop("Each nominal factor must have at least two integer levels.")
      }
    } else if (is.list(l_nom)) {
      if (length(l_nom) != p_nom) {
        stop("Input error: Length of l_nom does not match p_nom.")
      }
      if (any(vapply(l_nom, function(x) length(x) < 2, logical(1)))) {
        stop("Each nominal factor must have at least two levels.")
      }
    } else {
      stop("Input error: l_nom must be either numeric or a list.")
    }
  }

  U_l <- function(x, m) {
    (x - 1) %/% (length(unique(c(x))) / m) + 1
  }

  V_l <- function(x, m) {
    (U_l(x, m) +
      ((stats::ave(U_l(x, m), U_l(x, m), FUN = seq_along) - 1) %% (m - 1))) %%
      m +
      1
  }

  find_max_derangement <- function(values, k) {
    m <- length(values)
    h <- ifelse(m == 3, 2, ceiling((m - 1) / 2))
    mat <- matrix(NA, min(k, h), m)

    cost <- max(values) - sapply(values, function(x) abs(x - values))
    diag(cost) <- sum(cost) + 1
    for (i in 1:min(k, h)) {
      pairs <- RcppHungarian::HungarianSolver(cost)$pairs
      mat[i, ] <- values[pairs[, 2]]
      cost[pairs] <- sum(cost) + 1
    }

    if (i == k) {
      return(mat)
    }
    do.call(
      rbind,
      replicate(k %/% ceiling((m - 1) / 2) + 1, mat, simplify = FALSE)
    )[1:k, ]
  }

  criterion <- function(D, p_nom = 0, scale_continuous = TRUE) {
    n <- nrow(D)
    p <- ncol(D)
    l <- n / (p + 1)
    eps <- 1 / l^2

    p_numeric <- p - p_nom
    numeric_idx <- if (p_numeric > 0) seq_len(p_numeric) else integer(0)
    X <- D

    if (p_numeric > 0 && isTRUE(scale_continuous)) {
      for (j in numeric_idx) {
        r <- range(X[, j])
        X[, j] <- (X[, j] - r[1]) / (r[2] - r[1])
      }
    }

    pairs <- combn(p, 2)
    out <- numeric(ncol(pairs))

    for (k in seq_len(ncol(pairs))) {
      i <- pairs[1, k]
      j <- pairs[2, k]
      ridx <- c(1:l, (i * l + 1):((i + 1) * l), (j * l + 1):((j + 1) * l))

      xi <- X[ridx, i]
      xj <- X[ridx, j]

      if (i <= p_numeric) {
        di <- as.numeric(stats::dist(xi))
        di2 <- di * di
      } else {
        di2 <- as.numeric(stats::dist(xi) > 0)
      }
      if (j <= p_numeric) {
        dj <- as.numeric(stats::dist(xj))
        dj2 <- dj * dj
      } else {
        dj2 <- as.numeric(stats::dist(xj) > 0)
      }

      d2 <- di2 + dj2
      out[k] <- mean(1 / (d2 + eps))
    }

    mean(out)
  }

  gen_design <- function() {
    raw <- mofat(p = p_continuous + p_disnum + p_nom, l = l, nstart = 1, ...)

    if (p_continuous > 0) {
      continuous <- raw[, 1:p_continuous, drop = FALSE]
      if (!is.null(b_continuous)) {
        for (j in 1:p_continuous) {
          continuous[, j] <- b_continuous[[j]][1] +
            (b_continuous[[j]][2] - b_continuous[[j]][1]) * continuous[, j]
        }
      }
    } else {
      continuous <- NULL
    }

    if (p_disnum > 0) {
      disnum <- round(
        raw[, (p_continuous + 1):(p_continuous + p_disnum), drop = FALSE] *
          l +
          0.5
      )
      for (j in 1:p_disnum) {
        j_ <- j + p_continuous
        if (is.numeric(l_disnum)) {
          A_j <- disnum[1:l, j]
          disnum[setdiff(1:nrow(raw), (j_ * l + 1):(j_ * l + l)), j] <- U_l(
            A_j,
            l_disnum[j]
          )

          derangement <- find_max_derangement(
            1:l_disnum[j],
            ceiling(l / l_disnum[j])
          )

          disnum[(j_ * l + 1):(j_ * l + l), j] <- derangement[cbind(
            U_l(A_j, l_disnum[j]),
            (stats::ave(
              U_l(A_j, l_disnum[j]),
              U_l(A_j, l_disnum[j]),
              FUN = seq_along
            ) -
              1)
          )]
        } else if (is.list(l_disnum)) {
          A_j <- disnum[1:l, j]
          disnum[
            setdiff(1:nrow(raw), (j_ * l + 1):(j_ * l + l)),
            j
          ] <- l_disnum[[j]][U_l(A_j, length(l_disnum[[j]]))]

          derangement <- find_max_derangement(
            l_disnum[[j]],
            ceiling(l / length(l_disnum[[j]]))
          )

          disnum[(j_ * l + 1):(j_ * l + l), j] <- derangement[cbind(
            stats::ave(
              U_l(A_j, length(l_disnum[[j]])),
              U_l(A_j, length(l_disnum[[j]])),
              FUN = seq_along
            ),
            U_l(A_j, length(l_disnum[[j]]))
          )]
        }
      }
    } else {
      disnum <- NULL
    }

    if (p_nom > 0) {
      nom <- round(
        raw[,
          (p_continuous + p_disnum + 1):(p_continuous + p_disnum + p_nom),
          drop = FALSE
        ] *
          l +
          0.5
      )
      if (is.numeric(l_nom)) {
        for (j in 1:p_nom) {
          j_ <- j + p_continuous + p_disnum
          A_j <- nom[1:l, j]
          nom[setdiff(1:nrow(raw), (j_ * l + 1):(j_ * l + l)), j] <- U_l(
            A_j,
            l_nom[j]
          )
          nom[(j_ * l + 1):(j_ * l + l), j] <- V_l(A_j, l_nom[j])
        }
      } else if (is.list(l_nom)) {
        for (j in 1:p_nom) {
          j_ <- j + p_continuous + p_disnum
          A_j <- nom[1:l, j]
          nom[setdiff(1:nrow(raw), (j_ * l + 1):(j_ * l + l)), j] <- U_l(
            A_j,
            length(l_nom[[j]])
          )
          nom[(j_ * l + 1):(j_ * l + l), j] <- V_l(A_j, length(l_nom[[j]]))
        }
      }
    } else {
      nom <- NULL
    }

    cbind(continuous, disnum, nom)
  }

  crit <- numeric(nstart)
  seed <- stats::runif(1, 0, 1e4)
  for (i in 1:nstart) {
    D <- withr::with_seed(seed + i, gen_design())
    crit[i] <- criterion(D, p_nom = p_nom)
  }

  withr::with_seed(seed + which.min(crit), gen_design())
}
