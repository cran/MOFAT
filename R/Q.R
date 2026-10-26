#' Q Criterion
#'
#' Evaluate the \eqn{Q} criterion for a given design. This criterion measures the design's space-filling performance
#' and serves as the objective function used to construct SOFT designs.
#' See Yu and Joseph (2026) for further details.
#' 
#' @author Wei-Yang Yu and V. Roshan Joseph
#'
#' @param alpha positive MIM kernel exponent.
#' @param beta positive MIM kernel scale.
#' @param design the design matrix.
#' @param qmc the evaluation set.
#' @param backend one of \code{"auto"}, \code{"cpu"}, \code{"cuda"}, or \code{"metal"}. Invalid,
#'   unavailable, or failed GPU backends silently use CPU.
#'
#' @details The \eqn{Q} criterion for a design \eqn{\mathbf{D}=[\mathbf{x}_1,\ldots,\mathbf{x}_n]^T} is defined as
#' \deqn{
#' \min_{\mathbf{x}\in\mathbf{E}}
#' \mathbf{r}(\mathbf{x})^T
#' \mathbf{R}^{-1}\mathbf{r}(\mathbf{x}),
#' }
#' where \eqn{\mathbf{E}} is the evaluation set provided by \code{qmc}. 
#' Here, \eqn{\mathbf{r}(\mathbf{x}) = [R(\mathbf{x}_1,\mathbf{x}), \ldots, R(\mathbf{x}_n,\mathbf{x})]^T} and
#' \eqn{\mathbf{R}_{i,j} = R(\mathbf{x}_i,\mathbf{x}_j)}, where
#' \deqn{
#' R(\mathbf{x},\mathbf{x}') =
#' \prod_{j=1}^{p}
#' \left\{1+\frac{(x_j-x'_j)^2}{\beta}\right\}^{-\alpha}.
#' }
#' Larger values of \eqn{Q} indicate better space-filling performance.
#'
#' @return 
#' A list with the following components:
#'   \item{Q}{the \eqn{Q} criterion of the design.}
#'   \item{normalized_Q}{the normalized criterion defined as \eqn{Q^{1/p}}.}
#'   \item{backend_used}{The computational backend actually used.}
#' 
#' @references 
#' 
#' Yu, W.-Y. and Joseph, V. R. (2026). Space-Filling One-Factor-At-A-Time Designs. \doi{10.48550/arXiv.2606.02533}
#' 
#' @seealso \code{\link{qmc_generate}}
#' 
#' @export
#' @examples
#' design <- soft(p = 2, l = 4)$design
#' qmc <- qmc_generate(n = 1000, p = 2)
#' Q(1, 1/(2*4)^2, design, qmc)


Q <- function(alpha, beta, design, qmc, backend = "auto") {
  alpha <- .soft_positive_scalar(alpha, "alpha")
  beta <- .soft_positive_scalar(beta, "beta")
  design <- .soft_numeric_matrix(design, "design")
  qmc <- .soft_numeric_matrix(qmc, "qmc")
  if (ncol(design) != ncol(qmc)) {
    stop("design and qmc must have equal column counts.", call. = FALSE)
  }
  selected <- .soft_resolve_backend(backend)
  backend_used <- selected$backend
  if (identical(backend_used, "cpu")) {
    result <- .soft_q_backend("cpu", qmc, design, alpha, beta)
  } else {
    result <- tryCatch(
      .soft_q_backend(backend_used, qmc, design, alpha, beta),
      error = function(error) {
        backend_used <<- "cpu"
        .soft_q_backend("cpu", qmc, design, alpha, beta)
      }
    )
  }
  list(
    Q = as.numeric(result$Q),
    normalized_Q = as.numeric(result$normalized_Q),
    backend_used = backend_used
  )
}
