# GPU installation and verification

MOFAT uses optional GPU backends. CPU support is always compiled. GPU build,
load, device, or runtime failures fall back to CPU and are reported by
`mofat_backend_info()` and, where present, `fallback_reason`.

## Check the current runtime

```r
library(MOFAT)
mofat_backend_info(refresh = TRUE)
```

`compiled = TRUE` means the platform helper was installed. `available = TRUE`
means it loaded, passed ABI validation, found a device, and completed a small
self-test. A CUDA result also reports its helper ABI, runtime and driver
versions, compute capability, and helper path.

## Linux and Google Colab CUDA

Install an NVIDIA driver, CUDA Toolkit 12.x or 13.x, and an NVCC-compatible
GCC/Clang host compiler. `nvcc` must be on `PATH`, or set `CUDA_HOME` or
`CUDA_PATH` to the toolkit root.

```sh
R CMD INSTALL MOFAT_2.0.tar.gz
```

For a specific NVCC host compiler or GPU target:

```sh
NVCC_CCBIN=/usr/bin/g++ MOFAT_CUDA_ARCH=native \
  R CMD INSTALL MOFAT_2.0.tar.gz
```

In an R-enabled Colab runtime after uploading the archive:

```r
Sys.setenv(MOFAT_CUDA_ARCH = "native")
install.packages("/content/MOFAT_2.0.tar.gz",
                 repos = NULL, type = "source")
library(MOFAT)
mofat_backend_info(refresh = TRUE)
soft(4, 4, structure = "strict", backend = "cuda")$backend_used
```

The distributed CUDA helper contains no Rcpp or Armadillo code. It exchanges
primitive buffers through a versioned C ABI, avoiding the NVCC/R compiler ABI
mismatch that caused the earlier Colab segmentation fault.

## Windows CPU

Windows installations use the CPU backend. CUDA acceleration is supported
only on Linux, so Windows does not require CUDA, Visual Studio Build Tools, or
an NVIDIA driver beyond the normal requirements of R and Rtools source
packages.

## macOS Metal

Metal is enabled automatically when the Metal and Metal Performance Shaders
frameworks are available. No CUDA Toolkit is required.

```r
mofat_backend_info(refresh = TRUE)$metal
soft(4, 4, backend = "metal")$backend_used
```

Some headless macOS processes cannot obtain a Metal device. Those processes
fall back to CPU even though the Metal backend was compiled.

## Force CPU-only installation

```sh
MOFAT_DISABLE_GPU=1 R CMD INSTALL MOFAT_2.0.tar.gz
```

## Troubleshooting

- Inspect `config.log` in the extracted source tree after an optional CUDA
  helper build fails.
- Confirm `mofat_backend_info(refresh = TRUE)$cuda$helper_path` exists.
- Confirm the CUDA driver is new enough for the reported runtime.
- Reinstall after changing `CUDA_HOME`, `CUDA_PATH`, `NVCC_CCBIN`, or
  `MOFAT_CUDA_ARCH`.
- Use `backend_used` to confirm that a requested operation actually used the
  GPU; silent fallback never labels a CPU computation as CUDA or Metal.
