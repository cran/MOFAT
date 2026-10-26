test_that("mofat_backend_info exposes the stable backend contract", {
  info <- mofat_backend_info()
  expect_named(info, c("cpu", "metal", "cuda"), ignore.order = FALSE)
  for (backend in info) {
    expect_true(all(c("compiled", "available", "device", "reason") %in%
                    names(backend)))
    expect_type(backend$compiled, "logical")
    expect_length(backend$compiled, 1L)
    expect_type(backend$available, "logical")
    expect_length(backend$available, 1L)
  }
  expect_true(info$cpu$compiled)
  expect_true(info$cpu$available)
  expect_identical(info$cpu$device, "CPU")
  expect_null(info$cpu$reason)
  expect_true(all(c("helper_abi", "toolkit_version", "runtime_version",
                    "driver_version", "compute_capability", "helper_path") %in%
                  names(info$cuda)))
})

test_that("backend refresh validates its public argument", {
  expect_error(mofat_backend_info(NA), "TRUE or FALSE")
  expect_error(mofat_backend_info(1), "TRUE or FALSE")
  expect_error(mofat_backend_info(c(TRUE, FALSE)), "TRUE or FALSE")
  expect_type(mofat_backend_info(TRUE), "list")
})

test_that("a broken CUDA helper is diagnosed without crashing R", {
  skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
  fake_helper <- tempfile(fileext = ".so")
  writeBin(charToRaw("not a shared library"), fake_helper)
  old_helper <- Sys.getenv("MOFAT_CUDA_HELPER", unset = NA_character_)
  on.exit({
    if (is.na(old_helper)) {
      Sys.unsetenv("MOFAT_CUDA_HELPER")
    } else {
      Sys.setenv(MOFAT_CUDA_HELPER = old_helper)
    }
    mofat_backend_info(refresh = TRUE)
    unlink(fake_helper)
  }, add = TRUE)

  Sys.setenv(MOFAT_CUDA_HELPER = fake_helper)
  info <- mofat_backend_info(refresh = TRUE)$cuda
  expect_true(info$compiled)
  expect_false(info$available)
  expect_null(info$device)
  expect_match(info$reason, "load CUDA helper|Cannot load CUDA helper")

  fallback <- expect_silent(soft(3, 2, m = 16, backend = "cuda"))
  expect_identical(fallback$backend_used, "cpu")
  expect_match(fallback$fallback_reason, "CUDA helper|cuda backend")
})
