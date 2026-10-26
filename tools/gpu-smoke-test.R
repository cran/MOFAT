args <- commandArgs(trailingOnly = TRUE)
backend <- if (length(args) >= 1L) tolower(args[[1L]]) else "cuda"
mode <- if (length(args) >= 2L) tolower(args[[2L]]) else "all"

if (!backend %in% c("cuda", "metal")) {
  stop("backend must be 'cuda' or 'metal'.", call. = FALSE)
}

suppressPackageStartupMessages(library(MOFAT))

relative_error <- function(actual, expected) {
  abs(actual - expected) / max(1, abs(actual), abs(expected))
}

assert_gpu <- function(result) {
  if (!identical(result$backend_used, backend)) {
    reason <- result$fallback_reason
    if (is.null(reason)) reason <- "no fallback reason was returned"
    stop(sprintf("Expected backend '%s' but used '%s': %s",
                 backend, result$backend_used, reason), call. = FALSE)
  }
}

info <- mofat_backend_info(refresh = TRUE)
selected_info <- info[[backend]]
if (!isTRUE(selected_info$available)) {
  stop(sprintf("%s backend is unavailable: %s", backend,
               selected_info$reason), call. = FALSE)
}

report <- list(backend = backend, mode = mode, info = info,
               started = Sys.time())

run_info <- function() {
  print(info)
  invisible(TRUE)
}

run_q <- function() {
  design <- matrix(c(0.1, 0.2, 0.7, 0.3, 0.4, 0.9),
                   ncol = 2, byrow = TRUE)
  qmc <- matrix(c(0.0, 0.0, 0.3, 0.6, 1.0, 1.0),
                ncol = 2, byrow = TRUE)
  cpu <- Q(1.25, 0.07, design, qmc, backend = "cpu")
  gpu <- Q(1.25, 0.07, design, qmc, backend = backend)
  stopifnot(identical(gpu$backend_used, backend))
  if (relative_error(gpu$normalized_Q, cpu$normalized_Q) > 4e-4) {
    stop("GPU Q result exceeds the 4e-4 relative tolerance.", call. = FALSE)
  }
  list(cpu = cpu, gpu = gpu)
}

run_soft <- function(structure) {
  result <- soft(4, 4, structure = structure, m = 80, backend = backend)
  assert_gpu(result)
  stopifnot(identical(dim(result$design), c(20L, 4L)))
  stopifnot(all(is.finite(result$design)))
  stopifnot(all(result$design >= 0), all(result$design <= 1))
  result
}

run_augment <- function() {
  generated <- a_soft(4, 4, structure = "strict", n_aug = 3,
                            m = 80, N = 256, backend = backend)
  assert_gpu(generated)
  stopifnot(identical(dim(generated$design), c(23L, 4L)))
  stopifnot(all(is.finite(generated$design)))

  existing <- generated$design[seq_len(20L), , drop = FALSE]
  colnames(existing) <- paste0("x", seq_len(4L))
  supplied <- a_soft(n_aug = 2, N = 256, backend = backend,
                           ofat_design = existing)
  assert_gpu(supplied)
  stopifnot(identical(dim(supplied$design), c(22L, 4L)))
  stopifnot(identical(
    supplied$design[seq_len(nrow(existing)), , drop = FALSE], existing
  ))
  stopifnot(identical(colnames(supplied$design), colnames(existing)))
  stopifnot(all(is.finite(supplied$design)))

  list(generated = generated, supplied = supplied)
}

run_lifecycle <- function() {
  design <- soft(4, 4, m = 80, backend = backend)
  assert_gpu(design)
  qmc <- qmc_generate(128, 4)
  for (iteration in seq_len(25L)) {
    score <- Q(1, 1 / 64, design$design, qmc, backend = backend)
    stopifnot(identical(score$backend_used, backend))
    stopifnot(is.finite(score$normalized_Q), score$normalized_Q > 0)
  }
  TRUE
}

run_benchmark <- function() {
  benchmark_strict <- function(chosen_backend) {
    result <- NULL
    elapsed <- unname(system.time(
      result <- soft(20, 6, structure = "strict", m = 500,
                     backend = chosen_backend)
    )[["elapsed"]])
    if (identical(chosen_backend, backend)) assert_gpu(result)
    elapsed
  }
  benchmark_augmentation <- function(chosen_backend) {
    result <- NULL
    elapsed <- unname(system.time(
      result <- a_soft(8, 4, structure = "strict", n_aug = 8,
                             m = 500, N = 10000,
                             backend = chosen_backend)
    )[["elapsed"]])
    if (identical(chosen_backend, backend)) assert_gpu(result)
    elapsed
  }
  strict_cpu <- median(replicate(3L, benchmark_strict("cpu")))
  strict_gpu <- median(replicate(3L, benchmark_strict(backend)))
  augment_cpu <- median(replicate(3L, benchmark_augmentation("cpu")))
  augment_gpu <- median(replicate(3L, benchmark_augmentation(backend)))
  require_speedup <- identical(Sys.getenv("MOFAT_REQUIRE_GPU_SPEEDUP"), "1")
  if (require_speedup &&
      (!(strict_gpu < strict_cpu) || !(augment_gpu < augment_cpu))) {
    stop(sprintf(paste(
      "GPU did not beat CPU on both required L4 workloads:",
      "strict %.3fs/%.3fs, augmentation %.3fs/%.3fs (GPU/CPU)."
    ), strict_gpu, strict_cpu, augment_gpu, augment_cpu), call. = FALSE)
  }
  list(
    strict = list(cpu_seconds = strict_cpu, gpu_seconds = strict_gpu,
                  speedup = strict_cpu / strict_gpu),
    augmentation_n10000 = list(cpu_seconds = augment_cpu,
                               gpu_seconds = augment_gpu,
                               speedup = augment_cpu / augment_gpu)
  )
}

if (mode %in% c("all", "info")) report$diagnostics <- run_info()
if (mode %in% c("all", "q")) report$q <- run_q()
if (mode %in% c("all", "standard")) report$standard <- run_soft("standard")
if (mode %in% c("all", "strict")) report$strict <- run_soft("strict")
if (mode %in% c("all", "augment")) report$augment <- run_augment()
if (mode %in% c("all", "lifecycle")) report$lifecycle <- run_lifecycle()
if (mode == "benchmark") report$benchmark <- run_benchmark()

known_modes <- c("all", "info", "q", "standard", "strict", "augment",
                 "lifecycle", "benchmark")
if (!mode %in% known_modes) stop("Unknown GPU smoke-test mode.", call. = FALSE)

report$finished <- Sys.time()
report_path <- Sys.getenv("MOFAT_GPU_REPORT")
if (nzchar(report_path)) saveRDS(report, report_path)
print(report[c("backend", "mode", "started", "finished")])
