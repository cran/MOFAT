#' SOFT Designs
#'
#' Generate a SOFT design under either the standard or strict OFAT structure.
#' For \eqn{p} factors and \eqn{l} base runs, the resulting design contains
#' \eqn{l(p+1)} runs. The design is constructed by maximizing the
#' \eqn{Q} criterion under either the standard or strict OFAT
#' structure. Computation can use an available GPU backend, with CPU execution
#' used when GPU acceleration is unavailable. See Yu and Joseph (2026) for
#' further details.
#' 
#' @author Wei-Yang Yu and V. Roshan Joseph
#'
#' @param p number of factors; an integer of at least 2.
#' @param l even number of base runs; an integer of at least 2.
#' @param structure indicating the OFAT structure. Choose from
#'   \code{"standard"} or \code{"strict"}. The default is \code{"standard"}.
#' @param m number of Sobol evaluation points generated before vertex
#'   augmentation. The default is 500.
#' @param backend one of \code{"auto"}, \code{"cpu"}, \code{"cuda"}, or \code{"metal"}. Invalid,
#'   unavailable, or failed GPU backends silently use CPU.
#'
#'
#' @details
#' When optimizing \eqn{Q}, we set \eqn{\alpha = 1} and \eqn{\beta = (2*l)^{-2}} internally.
#'
#'
#' @return
#' A list with the following components:
#'   \item{design}{The SOFT design matrix.}
#'   \item{Q}{The \eqn{Q} criterion of the design.}
#'   \item{normalized_Q}{The normalized criterion defined as \eqn{Q^{1/p}}.}
#'   \item{qmc_size}{The number of unique evaluation points after adding corner points and removing duplicates.
#'     This value may differ from the input \code{m}.}
#'   \item{backend_used}{The computational backend actually used.}
#'   \item{device}{The CPU or GPU device used for computation.}
#'   \item{fallback_reason}{The reason for falling back to the CPU backend,
#'     or \code{NULL} if no fallback occurred.}
#' 
#' @references 
#' Yu, W.-Y. and Joseph, V. R. (2026). Space-Filling One-Factor-At-A-Time Designs. \doi{10.48550/arXiv.2606.02533}
#' 
#' @seealso \code{\link{Q}}, \code{\link{qmc_generate}}
#' 
#' @export
#' @examples
#' result <- soft(p = 2, l = 4)
#' result$design


soft <- function(p, l, structure = "standard", m = 500L,
                 backend = "auto") {
  p <- .soft_positive_integer(p, "p", 2L)
  l <- .soft_positive_integer(l, "l", 2L)
  if (l %% 2L != 0L) {
    stop("l must be even.", call. = FALSE)
  }
  m <- .soft_positive_integer(m, "m")
  structure <- .soft_structure(structure)
  prepared <- .soft_prepare(p, l, m, total_levels = 2L * l)
  selected <- .soft_resolve_backend(backend)
  fallback_reason <- selected$fallback_reason

  if (identical(selected$backend, "cpu")) {
    result <- .soft_optimize_cpu(prepared, structure)
  } else {
    result <- tryCatch(
      .soft_optimize_backend(prepared, structure, selected$backend),
      error = function(error) {
        fallback_reason <<- sprintf(
          "%s backend execution failed during design generation: %s",
          selected$backend, conditionMessage(error)
        )
        .soft_optimize_cpu(prepared, structure)
      }
    )
  }
  backend_used <- if (identical(result$device, "CPU")) "cpu" else
    selected$backend
  log_Q <- as.numeric(result$log_Q)
  list(
    design = result$design,
    Q = exp(log_Q),
    normalized_Q = exp(log_Q / p),
    qmc_size = nrow(prepared$qmc),
    backend_used = backend_used,
    device = result$device %||% selected$device,
    fallback_reason = fallback_reason
  )
}
