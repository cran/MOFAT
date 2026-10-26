direct_q <- function(alpha, beta, design, qmc, eta = 1e-6) {
  kernel <- function(a, b) {
    prod((1 + (a - b)^2 / beta)^(-alpha))
  }
  n <- nrow(design)
  R <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) R[i, j] <- kernel(design[i, ], design[j, ])
  }
  R <- R + diag(eta, n)
  values <- apply(qmc, 1, function(point) {
    r <- apply(design, 1, kernel, b = point)
    drop(crossprod(r, solve(R, r)))
  })
  min(values)
}

test_that("Q agrees with a direct small-matrix calculation", {
  design <- matrix(c(
    0.1, 0.2,
    0.7, 0.3,
    0.4, 0.9
  ), ncol = 2, byrow = TRUE)
  qmc <- matrix(c(
    0.0, 0.0,
    0.3, 0.6,
    1.0, 1.0
  ), ncol = 2, byrow = TRUE)
  expected <- direct_q(1.25, 0.07, design, qmc)
  result <- Q(1.25, 0.07, design, qmc, backend = "cpu")
  expect_named(result, c("Q", "normalized_Q", "backend_used"),
               ignore.order = FALSE)
  expect_equal(result$Q, expected, tolerance = 1e-11)
  expect_equal(result$normalized_Q, expected^(1 / ncol(design)),
               tolerance = 1e-12)
  expect_identical(result$backend_used, "cpu")
})

test_that("Q uses the supplied evaluation points without adding corners", {
  design <- matrix(c(0.2, 0.2, 0.8, 0.8), ncol = 2, byrow = TRUE)
  qmc <- matrix(c(0.5, 0.5), nrow = 1)
  expected <- direct_q(1, 0.05, design, qmc)
  result <- Q(1, 0.05, design, qmc, backend = "cpu")
  expect_equal(result$Q, expected, tolerance = 1e-12)
})

test_that("Q validates scalars and matrices", {
  design <- matrix(1:4, ncol = 2)
  qmc <- matrix(1:4, ncol = 2)
  expect_error(Q(0, 1, design, qmc), "alpha")
  expect_error(Q(1, Inf, design, qmc), "beta")
  expect_error(Q(1, 1, as.data.frame(design), qmc), "design")
  expect_error(Q(1, 1, design, matrix(1:3, ncol = 1)), "column")
})

test_that("the alpha-one fast path retains log-space underflow protection", {
  p <- 200L
  design <- rbind(rep(0, p), rep(1, p))
  qmc <- matrix(0.5, nrow = 1L, ncol = p)
  result <- Q(1, 1e-4, design, qmc, backend = "cpu")

  expect_identical(result$Q, 0)
  expect_true(is.finite(result$normalized_Q))
  expect_gt(result$normalized_Q, 0)
})
