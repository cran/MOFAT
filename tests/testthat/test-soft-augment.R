test_that("a_soft has the exact public contract and row count", {
  result <- a_soft(3, 2, n_aug = 2, m = 16, N = 64,
                         backend = "cpu")
  expect_named(result, c("design", "backend_used", "device",
                         "fallback_reason"), ignore.order = FALSE)
  expect_equal(dim(result$design), c(2 * (3 + 1) + 2, 3))
  expect_identical(result$backend_used, "cpu")
  expect_identical(result$device, "CPU")
  expect_null(result$fallback_reason)
})

test_that("augmentation uses the widened range for base and candidate rows", {
  p <- 3
  l <- 2
  n_aug <- 2
  boundary <- 0.25 / (2 * l + n_aug)
  result <- a_soft(p, l, n_aug = n_aug, m = 16, N = 64,
                         backend = "cpu")
  expect_true(all(result$design >= boundary - 1e-15))
  expect_true(all(result$design <= 1 - boundary + 1e-15))
  expect_true(any(abs(result$design[seq_len(l * (p + 1)), ] - boundary) < 1e-15))
  expect_true(any(abs(result$design[seq_len(l * (p + 1)), ] -
                      (1 - boundary)) < 1e-15))
})

test_that("augmentation defaults and validation match the public API", {
  expect_identical(formals(a_soft)$m, 500L)
  expect_identical(formals(a_soft)$N, 10000L)
  expect_null(formals(a_soft)$ofat_design)
  expect_identical(tail(names(formals(a_soft)), 1L), "ofat_design")
  expect_error(a_soft(3, 2, n_aug = 3, m = 8, N = 2),
               "must not exceed")
  expect_error(a_soft(3, 3, n_aug = 1, m = 8, N = 8), "even")
})

test_that("ofat_design NULL preserves legacy calls", {
  omitted <- a_soft(3, 2, "standard", 1, 8, 16, "cpu")
  explicit <- a_soft(3, 2, "standard", 1, 8, 16, "cpu", NULL)
  expect_identical(explicit, omitted)
})

test_that("a supplied design is augmented directly", {
  existing <- matrix(c(
    0.10, 0.20,
    0.35, 0.90,
    0.80, 0.40,
    0.70, 0.60,
    0.25, 0.30,
    0.50, 0.75
  ), ncol = 2, byrow = TRUE,
  dimnames = list(NULL, c("factor_a", "factor_b")))

  result <- a_soft(
    structure = stop("structure must not be evaluated"), n_aug = 2,
    m = stop("m must not be evaluated"), N = 32, backend = "cpu",
    ofat_design = existing
  )

  expect_named(result, c("design", "backend_used", "device",
                         "fallback_reason"), ignore.order = FALSE)
  expect_identical(dim(result$design), c(8L, 2L))
  expect_identical(result$design[seq_len(nrow(existing)), , drop = FALSE],
                   existing)
  expect_identical(colnames(result$design), colnames(existing))
  expect_identical(result$backend_used, "cpu")
  expect_identical(result$device, "CPU")
  expect_null(result$fallback_reason)
})

test_that("ofat_design dimensions override supplied p and l", {
  existing <- matrix(0.5, nrow = 6, ncol = 2)
  args <- list(n_aug = 1, N = 8, backend = "cpu", ofat_design = existing)
  expected <- do.call(a_soft, args)
  for (supplied in list(list(p = 2, l = 2), list(p = 5), list(l = 7),
                       list(p = 5, l = 7))) {
    expect_identical(do.call(a_soft, c(supplied, args)), expected)
  }
  expect_identical(
    a_soft(p = stop("p must not be evaluated"),
           l = stop("l must not be evaluated"), n_aug = 1, N = 8,
           backend = "cpu", ofat_design = existing),
    expected
  )
})

test_that("candidate bounds are fixed across initial-design modes", {
  p <- 2
  l <- 2
  n_aug <- 2
  N <- 32
  narrow <- matrix(c(
    0.40, 0.45,
    0.55, 0.60
  ), nrow = l * (p + 1), ncol = p, byrow = TRUE)
  wide <- matrix(c(
    0.00, 0.10,
    0.90, 1.00
  ), nrow = l * (p + 1), ncol = p, byrow = TRUE)
  captured <- list()

  local_mocked_bindings(
    .soft_optimize_backend = function(prepared, structure, backend) {
      list(design = matrix(0.5, nrow = l * (p + 1L), ncol = p),
           device = "CPU")
    },
    .soft_augment_backend = function(design, candidates, n_aug, beta,
                                     nugget, backend) {
      captured[[length(captured) + 1L]] <<- list(
        design = design, candidates = candidates, n_aug = n_aug,
        beta = beta, nugget = nugget, backend = backend
      )
      list(design = rbind(design, candidates[seq_len(n_aug), , drop = FALSE]),
           device = "CPU")
    },
    .package = "MOFAT"
  )

  a_soft(p, l, n_aug = n_aug, m = 8, N = N, backend = "cpu")
  a_soft(n_aug = n_aug, N = N, backend = "cpu",
               ofat_design = narrow)
  a_soft(n_aug = n_aug, N = N, backend = "cpu",
               ofat_design = wide)

  raw_candidates <- spacefillr::generate_sobol_set(n = N, dim = p, seed = 0)
  raw_candidates <- matrix(raw_candidates, ncol = p)
  boundary <- 0.25 / (2 * l + n_aug)
  expected_candidates <- scales::rescale(
    raw_candidates, to = c(boundary, 1 - boundary), from = c(0, 1)
  )
  expected_candidates <- unique(expected_candidates)

  expect_length(captured, 3L)
  for (call in captured) {
    expect_identical(call$candidates, expected_candidates)
    expect_identical(call$n_aug, as.integer(n_aug))
    expect_equal(call$beta, 1 / (2 * l + n_aug)^2)
    expect_identical(call$nugget, 1e-6)
    expect_identical(call$backend, "cpu")
  }
  expect_identical(captured[[2L]]$design, narrow)
  expect_identical(captured[[3L]]$design, wide)
})

test_that("ofat_design validation is strict", {
  augment_existing <- function(ofat_design) {
    a_soft(n_aug = 1, N = 8, backend = "cpu",
                 ofat_design = ofat_design)
  }

  expect_error(augment_existing(c(0.1, 0.9)),
               "nonempty finite numeric matrix")
  expect_error(augment_existing(data.frame(x = c(0.1, 0.9),
                                           y = c(0.2, 0.8))),
               "nonempty finite numeric matrix")
  expect_error(augment_existing(matrix(numeric(), nrow = 0, ncol = 2)),
               "nonempty finite numeric matrix")
  expect_error(augment_existing(matrix(letters[1:4], ncol = 2)),
               "nonempty finite numeric matrix")
  for (bad_value in c(NA_real_, NaN, Inf, -Inf)) {
    invalid <- matrix(c(0.1, 0.2, 0.8, bad_value), ncol = 2)
    expect_error(augment_existing(invalid),
                 "nonempty finite numeric matrix")
  }
  expect_error(augment_existing(matrix(0.5, nrow = 4, ncol = 1)),
               "p must be an integer of at least 2")
  for (runs in c(1, 2, 3, 4, 5, 7, 9, 15)) {
    expect_error(augment_existing(matrix(0.5, nrow = runs, ncol = 2)),
                 "The OFAT design should contain l \u00d7 (p + 1) runs.",
                 fixed = TRUE)
  }
  expect_error(augment_existing(matrix(c(-0.1, 0.2, 0.8, 0.9),
                                        nrow = 6, ncol = 2)),
               "values must lie")
  expect_error(augment_existing(matrix(c(0.1, 0.2, 0.8, 1.1),
                                        nrow = 6, ncol = 2)),
               "values must lie")
  constant <- matrix(0.5, nrow = 6, ncol = 2)
  result <- expect_silent(augment_existing(constant))
  expect_identical(dim(result$design), c(7L, 2L))
  expect_identical(result$design[seq_len(nrow(constant)), , drop = FALSE],
                   constant)
})

test_that("invalid augmentation structure silently uses standard", {
  standard <- a_soft(3, 2, structure = "standard", n_aug = 1,
                           m = 16, N = 32, backend = "cpu")
  invalid <- expect_silent(
    a_soft(3, 2, structure = "invalid", n_aug = 1,
                 m = 16, N = 32, backend = "cpu")
  )
  expect_equal(invalid$design, standard$design)
})

test_that("candidate deduplication is row-based and independent of design", {
  candidates <- rbind(c(0.2, 0.3), c(0.2, 0.3), c(0.7, 0.8))
  expect_equal(nrow(unique(candidates)), 2)
  design <- matrix(c(0.2, 0.3, 0.8, 0.7), ncol = 2, byrow = TRUE)
  augmented <- MOFAT:::.soft_augment_cpu(
    design, unique(candidates), 1, beta = 1 / 25
  )$design
  expect_equal(nrow(augmented), nrow(design) + 1)
})

test_that("compiled CPU augmentation matches the retained R reference", {
  prepared <- MOFAT:::.soft_prepare(4L, 4L, 80L, total_levels = 11L)
  initial <- MOFAT:::.soft_optimize_cpu(prepared, "standard")$design
  set.seed(11)
  candidates <- matrix(runif(300L * 4L), ncol = 4L)
  compiled <- MOFAT:::.soft_augment_cpu(
    initial, candidates, 5L, beta = 1 / 11^2
  )
  reference <- MOFAT:::.soft_augment_cpu_reference(
    initial, candidates, 5L, beta = 1 / 11^2
  )
  expect_identical(compiled$design, reference$design)
})

test_that("the initial augmented design keeps the base beta", {
  p <- 3
  l <- 2
  n_aug <- 2
  prepared <- MOFAT:::.soft_prepare(
    p, l, m = 16, total_levels = 2 * l + n_aug
  )
  expect_equal(prepared$beta, 1 / (2 * l)^2)
  initial <- MOFAT:::.soft_optimize_cpu(prepared, "standard")$design
  result <- a_soft(p, l, n_aug = n_aug, m = 16, N = 64,
                         backend = "cpu")$design
  expect_equal(result[seq_len(nrow(initial)), , drop = FALSE], initial)
})

test_that("a GPU augmentation failure restarts the complete workflow on CPU", {
  expected <- a_soft(3, 2, n_aug = 1, m = 8, N = 16,
                           backend = "cpu")$design
  calls <- character()
  original_augment <- MOFAT:::.soft_augment_backend
  local_mocked_bindings(
    .soft_resolve_backend = function(backend) {
      list(backend = "metal", device = "Test GPU", fallback_reason = NULL)
    },
    .soft_optimize_backend = function(prepared, structure, backend) {
      calls <<- c(calls, paste0("opt-", backend))
      result <- MOFAT:::.soft_optimize_cpu(prepared, structure)
      if (backend == "metal") result$device <- "Test GPU"
      result
    },
    .soft_augment_backend = function(design, candidates, n_aug, beta,
                                     nugget, backend) {
      calls <<- c(calls, paste0("augment-", backend))
      if (backend == "metal") stop("simulated GPU error")
      original_augment(design, candidates, n_aug, beta, nugget, backend)
    },
    .package = "MOFAT"
  )
  result <- expect_silent(
    a_soft(3, 2, n_aug = 1, m = 8, N = 16, backend = "metal")
  )
  expect_identical(calls, c("opt-metal", "augment-metal", "opt-cpu",
                            "augment-cpu"))
  expect_equal(result$design, expected)
  expect_identical(result$backend_used, "cpu")
  expect_match(result$fallback_reason, "augmentation")
})

test_that("a supplied-design GPU failure retries only augmentation on CPU", {
  existing <- matrix(c(
    0.10, 0.20,
    0.35, 0.90,
    0.80, 0.40,
    0.70, 0.60,
    0.25, 0.30,
    0.50, 0.75
  ), ncol = 2, byrow = TRUE)
  calls <- character()
  original_augment <- MOFAT:::.soft_augment_backend

  local_mocked_bindings(
    .soft_resolve_backend = function(backend) {
      list(backend = "metal", device = "Test GPU", fallback_reason = NULL)
    },
    .soft_optimize_backend = function(...) {
      stop("initial optimization must not run")
    },
    .soft_augment_backend = function(design, candidates, n_aug, beta,
                                     nugget, backend) {
      calls <<- c(calls, paste0("augment-", backend))
      if (backend == "metal") stop("simulated GPU error")
      original_augment(design, candidates, n_aug, beta, nugget, backend)
    },
    .package = "MOFAT"
  )

  result <- expect_silent(
    a_soft(n_aug = 2, N = 32, backend = "metal",
                 ofat_design = existing)
  )

  expect_identical(calls, c("augment-metal", "augment-cpu"))
  expect_identical(result$design[seq_len(nrow(existing)), , drop = FALSE],
                   existing)
  expect_identical(result$backend_used, "cpu")
  expect_identical(result$device, "CPU")
  expect_match(result$fallback_reason, "augmentation")
})
