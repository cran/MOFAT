`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.soft_positive_integer <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != floor(x) || x < minimum) {
    stop(sprintf("%s must be an integer of at least %d.", name, minimum),
         call. = FALSE)
  }
  as.integer(x)
}

.soft_positive_scalar <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x <= 0) {
    stop(sprintf("%s must be a positive finite scalar.", name),
         call. = FALSE)
  }
  as.numeric(x)
}

.soft_numeric_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || !length(x) ||
      nrow(x) < 1L || ncol(x) < 1L || any(!is.finite(x))) {
    stop(sprintf("%s must be a nonempty finite numeric matrix.", name),
         call. = FALSE)
  }
  matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x),
         dimnames = dimnames(x))
}

.soft_structure <- function(structure) {
  if (length(structure) == 1L && !is.na(structure) &&
      is.character(structure) && structure %in% c("standard", "strict")) {
    structure
  } else {
    "standard"
  }
}

.soft_backend_capabilities <- function(refresh = FALSE) {
  soft_backend_capabilities_rcpp(refresh)
}

#' Inspect MOFAT Computational Backends
#'
#' Report whether the CPU, Metal, and CUDA backends were compiled and whether
#' they are usable in the current R process. CUDA probing includes helper ABI
#' validation and a small runtime self-test. Results are cached for the R
#' session unless \code{refresh = TRUE} is requested.
#'
#' @param refresh logical scalar. If \code{TRUE}, discard the cached CUDA
#'   probe and repeat helper loading, device discovery, and the self-test.
#'
#' @return A named list with \code{cpu}, \code{metal}, and \code{cuda}
#'   components. Every component contains \code{compiled}, \code{available},
#'   \code{device}, and \code{reason}. CUDA additionally reports helper,
#'   toolkit, runtime, driver, and compute-capability details when available.
#'
#' @examples
#' info <- mofat_backend_info()
#' vapply(info, function(x) isTRUE(x$available), logical(1))
#' info$cuda$reason
#'
#' @export
mofat_backend_info <- function(refresh = FALSE) {
  if (length(refresh) != 1L || !is.logical(refresh) || is.na(refresh)) {
    stop("refresh must be TRUE or FALSE.", call. = FALSE)
  }
  .soft_backend_capabilities(refresh)
}

.soft_resolve_backend <- function(backend) {
  requested <- if (length(backend) == 1L && is.character(backend) &&
                   !is.na(backend)) tolower(backend) else ""
  valid <- c("auto", "cpu", "cuda", "metal")

  if (!requested %in% valid) {
    return(list(backend = "cpu", device = "CPU",
                fallback_reason = "The requested backend is invalid."))
  }
  if (identical(requested, "cpu")) {
    return(list(backend = "cpu", device = "CPU",
                fallback_reason = NULL))
  }
  capabilities <- .soft_backend_capabilities()
  if (!identical(requested, "auto")) {
    info <- capabilities[[requested]]
    if (!isTRUE(info$available)) {
      return(list(
        backend = "cpu",
        device = capabilities$cpu$device %||% "CPU",
        fallback_reason = info$reason %||%
          sprintf("The requested %s backend is unavailable.", requested)
      ))
    }
    return(list(backend = requested, device = info$device,
                fallback_reason = NULL))
  }

  platform <- Sys.info()[["sysname"]]
  if (identical(platform, "Darwin") && isTRUE(capabilities$metal$available)) {
    return(list(backend = "metal", device = capabilities$metal$device,
                fallback_reason = NULL))
  }
  if (identical(platform, "Linux") && isTRUE(capabilities$cuda$available)) {
    return(list(backend = "cuda", device = capabilities$cuda$device,
                fallback_reason = NULL))
  }
  list(backend = "cpu", device = capabilities$cpu$device %||% "CPU",
       fallback_reason = NULL)
}

.soft_metal_kernel_path <- function() {
  system.file("metal", "mmse_kernels.metal", package = "MOFAT")
}

.soft_prepare <- function(p, l, m, total_levels = 2L * l) {
  qmc <- qmc_generate(m, p)
  A_info <- A_fun(p, l)
  x <- 0.25 / total_levels
  levels <- seq(x, 1 - x, length.out = 2L * l)
  A <- level_trans(A_info$A, l, levels)
  list(
    p = p,
    l = l,
    m = m,
    qmc = qmc,
    A = A,
    levels = levels,
    delta = levels[l + 1L] - levels[1L],
    factorial_l = A_info$L,
    remainder = A_info$mod,
    quotient = A_info$quotient,
    alpha = 1,
    beta = 1 / (2L * l)^2,
    nugget = 1e-6
  )
}

.soft_shift_block <- function(block, factor, delta) {
  block[, factor] <- ifelse(
    block[, factor] < 0.5,
    block[, factor] + delta,
    block[, factor] - delta
  )
  block
}

.soft_standard_design <- function(A, p, delta) {
  blocks <- vector("list", p + 1L)
  blocks[[1L]] <- A
  for (factor in seq_len(p)) {
    blocks[[factor + 1L]] <- .soft_shift_block(A, factor, delta)
  }
  do.call(rbind, blocks)
}

.soft_q_backend <- function(backend, qmc, design, alpha, beta, nugget = 1e-6) {
  switch(
    backend,
    cpu = Q_score_cpu_rcpp(qmc, design, alpha, beta, nugget),
    metal = soft_q_metal_rcpp(qmc, design, alpha, beta, nugget,
                              .soft_metal_kernel_path()),
    cuda = soft_q_cuda_rcpp(qmc, design, alpha, beta, nugget),
    stop("Internal backend dispatch error.", call. = FALSE)
  )
}

.soft_log_correlations <- function(X, evaluation, alpha, beta) {
  answer <- matrix(0, nrow(X), nrow(evaluation))
  for (row in seq_len(nrow(X))) {
    differences <- sweep(evaluation, 2L, X[row, ], FUN = "-")
    answer[row, ] <- rowSums(-alpha * log1p(differences^2 / beta))
  }
  answer
}

.soft_incremental_state <- function(design, evaluation, alpha, beta, nugget) {
  log_R <- .soft_log_correlations(design, design, alpha, beta)
  R <- exp(log_R)
  diag(R) <- diag(R) + nugget
  L <- t(chol(R))
  log_rhs <- .soft_log_correlations(design, evaluation, alpha, beta)
  scales <- apply(log_rhs, 2L, max)
  rhs <- exp(sweep(log_rhs, 2L, scales, FUN = "-"))
  solutions <- forwardsolve(L, rhs)
  norms <- colSums(solutions^2)
  list(
    design = design,
    evaluation = evaluation,
    alpha = alpha,
    beta = beta,
    nugget = nugget,
    L = L,
    solutions = solutions,
    scales = scales,
    norms = norms,
    log_values = 2 * scales + log(pmax(norms, .Machine$double.xmin))
  )
}

.soft_append_state <- function(state, block) {
  C <- exp(.soft_log_correlations(
    state$design, block, state$alpha, state$beta
  ))
  D <- exp(.soft_log_correlations(block, block, state$alpha, state$beta))
  diag(D) <- diag(D) + state$nugget
  W <- forwardsolve(state$L, C)
  G <- t(chol(D - crossprod(W)))

  log_rhs <- .soft_log_correlations(
    block, state$evaluation, state$alpha, state$beta
  )
  block_scales <- apply(log_rhs, 2L, max)
  new_scales <- pmax(state$scales, block_scales)
  old_factors <- exp(state$scales - new_scales)
  scaled_solutions <- sweep(
    state$solutions, 2L, old_factors, FUN = "*"
  )
  rhs <- exp(sweep(log_rhs, 2L, new_scales, FUN = "-")) -
    crossprod(W, scaled_solutions)
  Z <- forwardsolve(G, rhs)
  new_norms <- state$norms * old_factors^2 + colSums(Z^2)

  old_n <- nrow(state$L)
  block_n <- nrow(block)
  new_L <- matrix(0, old_n + block_n, old_n + block_n)
  new_L[seq_len(old_n), seq_len(old_n)] <- state$L
  new_L[old_n + seq_len(block_n), seq_len(old_n)] <- t(W)
  new_L[old_n + seq_len(block_n), old_n + seq_len(block_n)] <- G
  state$design <- rbind(state$design, block)
  state$L <- new_L
  state$solutions <- rbind(scaled_solutions, Z)
  state$scales <- new_scales
  state$norms <- new_norms
  state$log_values <- 2 * new_scales +
    log(pmax(new_norms, .Machine$double.xmin))
  state
}

.soft_strict_cpu <- function(A, p, delta, qmc, alpha, beta, nugget) {
  blocks <- list(A, .soft_shift_block(A, 1L, delta))
  state <- .soft_incremental_state(
    do.call(rbind, blocks), qmc, alpha, beta, nugget
  )
  if (p >= 2L) {
    for (factor in 2L:p) {
      left_block <- .soft_shift_block(blocks[[1L]], factor, delta)
      right_block <- .soft_shift_block(blocks[[length(blocks)]], factor, delta)
      left_state <- .soft_append_state(state, left_block)
      right_state <- .soft_append_state(state, right_block)
      left_score <- min(left_state$log_values)
      right_score <- min(right_state$log_values)
      if (left_score >= right_score) {
        blocks <- c(list(left_block), blocks)
        state <- left_state
      } else {
        blocks <- c(blocks, list(right_block))
        state <- right_state
      }
    }
  }
  design <- do.call(rbind, blocks)
  list(design = design, log_Q = min(state$log_values), device = "CPU")
}

.soft_initial_cpu <- function(A, prepared, structure) {
  if (identical(structure, "standard")) {
    design <- .soft_standard_design(A, prepared$p, prepared$delta)
    score <- Q_score_cpu_rcpp(prepared$qmc, design, prepared$alpha,
                              prepared$beta, prepared$nugget)
    list(design = design, log_Q = score$log_Q, device = "CPU")
  } else {
    .soft_strict_cpu(A, prepared$p, prepared$delta, prepared$qmc,
                     prepared$alpha, prepared$beta, prepared$nugget)
  }
}

.soft_optimize_cpu_reference <- function(prepared, structure) {
  A <- prepared$A
  current <- .soft_initial_cpu(A, prepared, structure)
  p0 <- if (is.null(prepared$remainder)) prepared$p else prepared$remainder
  p1 <- if (is.null(prepared$quotient)) 0L else
    prepared$factorial_l * prepared$quotient

  if (p0 > 1L) {
    for (i in seq_len(prepared$l - 1L)) {
      for (j in (i + 1L):prepared$l) {
        for (k in 2L:p0) {
          column <- p1 + k
          A_new <- A
          A_new[c(i, j), column] <- A_new[c(j, i), column]
          proposed <- .soft_initial_cpu(A_new, prepared, structure)
          if (proposed$log_Q > current$log_Q) {
            A <- A_new
            current <- proposed
          }
        }
      }
    }
  }
  current
}

.soft_optimize_cpu <- function(prepared, structure) {
  remainder <- if (is.null(prepared$remainder)) -1L else prepared$remainder
  quotient <- if (is.null(prepared$quotient)) -1L else prepared$quotient
  soft_opt_cpu_rcpp(
    prepared$qmc, prepared$A, prepared$p, prepared$l, prepared$delta,
    remainder, quotient, structure, prepared$nugget
  )
}

.soft_optimize_backend <- function(prepared, structure, backend) {
  if (identical(backend, "cpu")) {
    return(.soft_optimize_cpu(prepared, structure))
  }
  remainder <- if (is.null(prepared$remainder)) -1L else prepared$remainder
  quotient <- if (is.null(prepared$quotient)) -1L else prepared$quotient
  if (identical(backend, "metal")) {
    return(soft_opt_metal_rcpp(
      prepared$qmc, prepared$A, prepared$p, prepared$l, prepared$delta,
      remainder, quotient, structure, prepared$nugget,
      .soft_metal_kernel_path()
    ))
  }
  soft_opt_cuda_rcpp(
    prepared$qmc, prepared$A, prepared$p, prepared$l, prepared$delta,
    remainder, quotient, structure, prepared$nugget
  )
}

.soft_augment_cpu_reference <- function(design, candidates, n_aug, beta,
                                        nugget = 1e-6) {
  eligible <- rep(TRUE, nrow(candidates))
  state <- .soft_incremental_state(design, candidates, 1, beta, nugget)
  for (step in seq_len(n_aug)) {
    repeat {
      if (!any(eligible)) {
        stop("No numerically eligible augmentation candidate remains.",
             call. = FALSE)
      }
      values <- state$log_values
      values[!eligible] <- Inf
      selected <- which.min(values)
      proposed <- tryCatch(
        .soft_append_state(state, candidates[selected, , drop = FALSE]),
        error = function(error) NULL
      )
      if (!is.null(proposed)) {
        state <- proposed
        eligible[selected] <- FALSE
        break
      }
      eligible[selected] <- FALSE
    }
  }
  list(design = state$design, device = "CPU")
}

.soft_augment_cpu <- function(design, candidates, n_aug, beta,
                              nugget = 1e-6) {
  soft_augment_cpu_rcpp(design, candidates, n_aug, beta, nugget)
}

.soft_augment_backend <- function(design, candidates, n_aug, beta, nugget,
                                  backend) {
  if (identical(backend, "cpu")) {
    return(.soft_augment_cpu(design, candidates, n_aug, beta, nugget))
  }
  if (identical(backend, "metal")) {
    return(soft_augment_metal_backend_rcpp(
      design, candidates, n_aug, beta, nugget, .soft_metal_kernel_path()
    ))
  }
  soft_augment_cuda_backend_rcpp(design, candidates, n_aug, beta, nugget)
}
