test_that("soft validates inputs and has the exact public contract", {
  result <- soft(3, 2, m = 16, backend = "cpu")
  expect_named(result, c(
    "design", "Q", "normalized_Q", "qmc_size", "backend_used",
    "device", "fallback_reason"
  ), ignore.order = FALSE)
  expect_equal(dim(result$design), c(2 * (3 + 1), 3))
  expect_equal(result$normalized_Q, result$Q^(1 / 3), tolerance = 1e-12)
  expect_equal(result$qmc_size, 16 + 2^3)
  expect_identical(result$backend_used, "cpu")
  expect_identical(result$device, "CPU")
  expect_null(result$fallback_reason)

  expect_error(soft(1, 2, m = 8), "p")
  expect_error(soft(2, 3, m = 8), "even")
  expect_error(soft(2, 2, m = 0), "m")
})

test_that("compiled CPU optimization matches the retained R reference", {
  prepared <- MOFAT:::.soft_prepare(4L, 4L, 80L)
  for (structure in c("standard", "strict")) {
    compiled <- MOFAT:::.soft_optimize_cpu(prepared, structure)
    reference <- MOFAT:::.soft_optimize_cpu_reference(prepared, structure)
    expect_identical(compiled$design, reference$design)
    expect_equal(compiled$log_Q, reference$log_Q, tolerance = 1e-12)
  }
})

test_that("invalid structures silently use standard", {
  standard <- soft(3, 2, structure = "standard", m = 16, backend = "cpu")
  invalid <- expect_silent(
    soft(3, 2, structure = "not-a-structure", m = 16, backend = "cpu")
  )
  expect_equal(invalid$design, standard$design)
  expect_equal(invalid$Q, standard$Q)
})

test_that("invalid and unavailable backend requests silently use CPU", {
  invalid <- expect_silent(soft(3, 2, m = 16, backend = "invalid"))
  expect_identical(invalid$backend_used, "cpu")
  expect_identical(invalid$device, "CPU")
  expect_match(invalid$fallback_reason, "invalid")

  capabilities <- MOFAT:::.soft_backend_capabilities()
  unavailable <- if (!isTRUE(capabilities$metal$available)) "metal" else "cuda"
  forced <- expect_silent(soft(3, 2, m = 16, backend = unavailable))
  expect_identical(forced$backend_used, "cpu")
  expect_true(is.character(forced$fallback_reason))
})

test_that("standard blocks change one specified factor from the base", {
  p <- 4
  l <- 2
  design <- soft(p, l, m = 16, backend = "cpu")$design
  base <- design[seq_len(l), , drop = FALSE]
  for (factor in seq_len(p)) {
    rows <- (factor * l + 1L):((factor + 1L) * l)
    block <- design[rows, , drop = FALSE]
    changed <- which(vapply(seq_len(p), function(column) {
      !isTRUE(all.equal(base[, column], block[, column], tolerance = 0))
    }, logical(1)))
    expect_identical(changed, factor)
  }
})

test_that("neighboring strict blocks differ in exactly one factor", {
  p <- 4
  l <- 2
  design <- soft(p, l, structure = "strict", m = 16,
                 backend = "cpu")$design
  blocks <- lapply(seq_len(p + 1L), function(index) {
    rows <- ((index - 1L) * l + 1L):(index * l)
    design[rows, , drop = FALSE]
  })
  for (index in seq_len(p)) {
    differences <- vapply(seq_len(p), function(column) {
      !isTRUE(all.equal(blocks[[index]][, column],
                        blocks[[index + 1L]][, column], tolerance = 0))
    }, logical(1))
    expect_equal(sum(differences), 1)
  }
})
