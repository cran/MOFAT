#ifndef MOFAT_CUDA_HELPER_API_H
#define MOFAT_CUDA_HELPER_API_H

#include <stddef.h>

#define MOFAT_CUDA_HELPER_ABI_VERSION 1
#define MOFAT_CUDA_ERROR_CAPACITY 1024
#define MOFAT_CUDA_DEVICE_CAPACITY 256

enum MofatCudaStatus {
  MOFAT_CUDA_SUCCESS = 0,
  MOFAT_CUDA_INVALID_ARGUMENT = 1,
  MOFAT_CUDA_UNAVAILABLE = 2,
  MOFAT_CUDA_RUNTIME_ERROR = 3,
  MOFAT_CUDA_NUMERICAL_ERROR = 4,
  MOFAT_CUDA_ABI_ERROR = 5
};

#if defined(__GNUC__)
#define MOFAT_CUDA_EXPORT __attribute__((visibility("default")))
#else
#define MOFAT_CUDA_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

MOFAT_CUDA_EXPORT int mofat_cuda_helper_abi_version(void);

MOFAT_CUDA_EXPORT int mofat_cuda_helper_info(
    int* toolkit_version, int* runtime_version, int* driver_version,
    int* compute_major, int* compute_minor,
    char* device, size_t device_capacity,
    char* error, size_t error_capacity);

MOFAT_CUDA_EXPORT int mofat_cuda_helper_self_test(
    char* error, size_t error_capacity);

MOFAT_CUDA_EXPORT int mofat_cuda_helper_q(
    const double* qmc, int qmc_rows, int p,
    const double* design, int design_rows,
    double alpha, double beta, double nugget,
    double* log_q, char* error, size_t error_capacity);

MOFAT_CUDA_EXPORT int mofat_cuda_helper_optimize(
    const double* qmc, int qmc_rows,
    const double* initial_a, int p, int l,
    double delta, int remainder, int quotient,
    int strict_structure, double nugget,
    double* output_design, double* log_q,
    char* error, size_t error_capacity);

MOFAT_CUDA_EXPORT int mofat_cuda_helper_augment(
    const double* design, int design_rows,
    const double* candidates, int candidate_rows,
    int p, int n_aug, double beta, double nugget,
    double* output_design,
    char* error, size_t error_capacity);

#ifdef __cplusplus
}
#endif

#endif
