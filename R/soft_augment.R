#' Augmented SOFT (A-SOFT) Designs
#' 
#' Construct an augmented SOFT (A-SOFT) design by combining a SOFT component
#' with sequentially selected unrestricted space-filling points. The resulting
#' design has the form
#' \deqn{
#'   D =
#'   \begin{pmatrix}
#'     D_{\mathrm{soft}} \\
#'     D_{\mathrm{aug}}
#'   \end{pmatrix},
#' }
#' where \eqn{D_{\mathrm{soft}}} is the initial SOFT design and
#' \eqn{D_{\mathrm{aug}}} contains the augmented points.
#' 
#' When \code{ofat_design = NULL}, the SOFT component contains
#' \eqn{n_{\mathrm{soft}} = l(p+1)} runs, and the final design contains
#' \eqn{l(p+1) + n_{\mathrm{aug}}} runs. Alternatively, a user-supplied
#' initial design can be provided through \code{ofat_design}, in which case
#' the supplied design is held fixed and augmented directly.
#' 
#' 
#' The augmented points are selected sequentially from a Sobol candidate pool
#' by maximizing the \eqn{Q} criterion. The augmented points themselves are not constrained to
#' follow the OFAT structure. CPU and optional GPU computation are supported.
#' See Yu and Joseph (2026) for further details.
#' 
#' @author Wei-Yang Yu and V. Roshan Joseph
#' 
#' @param p number of factors; an integer of at least 2. Not required when
#'   \code{ofat_design} is supplied; inferred as \code{ncol(ofat_design)},
#'   overriding any supplied value.
#' @param l even number of base runs; an integer of at least 2. Not required
#'   when \code{ofat_design} is supplied; inferred as
#'   \code{nrow(ofat_design) / (p + 1)}, overriding any supplied value.
#'   Also determines the augmentation kernel scale.
#' @param structure indicating the OFAT structure. Choose from
#'   \code{"standard"} or \code{"strict"}. The default is \code{"standard"}.
#'   Ignored when \code{ofat_design} is supplied.
#' @param n_aug number of points to augment.
#' @param m number of Sobol evaluation points for initial SOFT
#'   optimization. The default is 500. Ignored when \code{ofat_design} is
#'   supplied.
#' @param N size of the Sobol candidate pool for the sequential augmentation. The default is 10000.
#' @param backend one of \code{"auto"}, \code{"cpu"}, \code{"cuda"}, or \code{"metal"}. Invalid,
#'   unavailable, or failed GPU backends silently use CPU.
#' @param ofat_design \code{NULL}, or a nonempty numeric matrix with
#'   \code{p} columns, \code{l * (p + 1)} runs for an even integer \code{l},
#'   and values in \code{[0, 1]}.
#'   When supplied, \code{p} and \code{l} are inferred from its dimensions, and
#'   this matrix is augmented directly instead of generating an initial SOFT
#'   design. Its range does not affect the Sobol candidate bounds.
#'
#' @details
#' The initial SOFT design and the sequential augmentation use different kernel
#' scales for optimizing the \eqn{Q} criterion. For the initial SOFT optimization, the
#' kernel parameters are fixed internally at \eqn{\alpha=1} and
#' \eqn{\beta=(2l)^{-2}}. During sequential augmentation, they are fixed at
#' \eqn{\alpha=1} and
#' \eqn{\beta=(2l+n_{\mathrm{aug}})^{-2}}.
#'
#' Regardless of whether \code{ofat_design} is supplied, candidate bounds for the augmentation
#' are fixed at
#' \eqn{[0.25/(2l+n_{\mathrm{aug}}), 1-0.25/(2l+n_{\mathrm{aug}})]}. 
#'
#' @return 
#' A list with the following components:
#'   \item{design}{The design matrix containing the generated or supplied
#'     initial design followed by the sequentially augmented points.}
#'   \item{backend_used}{The computational backend actually used.}
#'   \item{device}{The CPU or GPU device used for computation.}
#'   \item{fallback_reason}{The reason for falling back to the CPU backend,
#'     or \code{NULL} if no fallback occurred.}
#' 
#' @references 
#' 
#' Yu, W.-Y. and Joseph, V. R. (2026). Space-Filling One-Factor-At-A-Time Designs. \doi{10.48550/arXiv.2606.02533}
#' 
#' @seealso \code{\link{soft}}, \code{\link{Q}}, \code{\link{qmc_generate}}
#' 
#' @export
#' 
#' @examples
#' result <- a_soft(p = 2, l = 4, n_aug = 8)
#' result$design
#'
#' existing <- mofat(p = 2, l = 4, method = "projection")
#' extended <- a_soft(n_aug = 8, ofat_design = existing)
#' extended$design


a_soft <- function(p, l, structure = "standard", n_aug,
                         m = 500L, N = 10000L, backend = "auto",
                         ofat_design = NULL) {
  has_ofat_design <- !is.null(ofat_design)
  if (has_ofat_design) {
    ofat_design <- .soft_numeric_matrix(ofat_design, "ofat_design")
    p <- ncol(ofat_design)
    l <- nrow(ofat_design) / (p + 1L)
    if (l %% 2L != 0L) {
      stop("The OFAT design should contain l \u00d7 (p + 1) runs.",
           call. = FALSE)
    }
    if (any(ofat_design < 0 | ofat_design > 1)) {
      stop("ofat_design values must lie in [0, 1].", call. = FALSE)
    }
  }
  p <- .soft_positive_integer(p, "p", 2L)
  l <- .soft_positive_integer(l, "l", 2L)
  if (l %% 2L != 0L) {
    stop("l must be even.", call. = FALSE)
  }
  n_aug <- .soft_positive_integer(n_aug, "n_aug")
  N <- .soft_positive_integer(N, "N")
  if (n_aug > N) {
    stop("n_aug must not exceed N.", call. = FALSE)
  }

  total_levels <- 2L * l + n_aug
  boundary <- 0.25 / total_levels
  candidate_bounds <- c(boundary, 1 - boundary)
  if (!has_ofat_design) {
    m <- .soft_positive_integer(m, "m")
    structure <- .soft_structure(structure)
    prepared <- .soft_prepare(p, l, m, total_levels = total_levels)
  }

  candidates <- spacefillr::generate_sobol_set(n = N, dim = p, seed = 0)
  candidates <- matrix(candidates, ncol = p)
  candidates <- scales::rescale(
    candidates,
    to = candidate_bounds,
    from = c(0, 1)
  )
  candidates <- unique(candidates)
  if (nrow(candidates) < n_aug) {
    stop("The deduplicated candidate pool contains fewer than n_aug rows.",
         call. = FALSE)
  }

  selected <- .soft_resolve_backend(backend)
  fallback_reason <- selected$fallback_reason
  beta_aug <- 1 / total_levels^2
  nugget <- 1e-6

  run_workflow <- function(chosen_backend) {
    if (has_ofat_design) {
      base_design <- ofat_design
      base_device <- NULL
    } else {
      stage <- "initial design generation"
      base <- tryCatch(
        .soft_optimize_backend(prepared, structure, chosen_backend),
        error = function(error) {
          attr(error, "soft_stage") <- stage
          stop(error)
        }
      )
      base_design <- base$design
      base_device <- base$device
    }
    stage <- "augmentation"
    augmented <- tryCatch(
      .soft_augment_backend(base_design, candidates, n_aug, beta_aug,
                            nugget, chosen_backend),
      error = function(error) {
        attr(error, "soft_stage") <- stage
        stop(error)
      }
    )
    list(design = augmented$design,
         device = augmented$device %||% base_device)
  }

  if (identical(selected$backend, "cpu")) {
    result <- run_workflow("cpu")
    backend_used <- "cpu"
  } else {
    result <- tryCatch(
      run_workflow(selected$backend),
      error = function(error) {
        stage <- attr(error, "soft_stage") %||% "workflow"
        fallback_reason <<- sprintf(
          "%s backend execution failed during %s: %s",
          selected$backend, stage, conditionMessage(error)
        )
        run_workflow("cpu")
      }
    )
    backend_used <- if (identical(result$device, "CPU")) "cpu" else
      selected$backend
  }

  if (has_ofat_design && !is.null(colnames(ofat_design))) {
    colnames(result$design) <- colnames(ofat_design)
  }

  list(
    design = result$design,
    backend_used = backend_used,
    device = result$device %||% selected$device,
    fallback_reason = fallback_reason
  )
}
