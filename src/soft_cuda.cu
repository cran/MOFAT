#include "cuda_helper_api.h"

#include <cuda_runtime.h>
#include <math_constants.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

class Matrix {
 public:
  Matrix() : rows_(0), cols_(0) {}
  Matrix(int rows, int cols)
      : rows_(rows), cols_(cols), values_(checked_size(rows, cols), 0.0) {}
  Matrix(const double* values, int rows, int cols)
      : rows_(rows), cols_(cols),
        values_(values, values + checked_size(rows, cols)) {}

  int rows() const { return rows_; }
  int cols() const { return cols_; }
  std::size_t size() const { return values_.size(); }
  double* data() { return values_.data(); }
  const double* data() const { return values_.data(); }

  double& operator()(int row, int col) {
    return values_[static_cast<std::size_t>(row) * cols_ + col];
  }
  double operator()(int row, int col) const {
    return values_[static_cast<std::size_t>(row) * cols_ + col];
  }

 private:
  int rows_, cols_;
  std::vector<double> values_;

  static std::size_t checked_size(int rows, int cols) {
    if (rows < 0 || cols < 0)
      throw std::invalid_argument("Matrix dimensions must be nonnegative.");
    return static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
  }
};

void copy_message(char* destination, std::size_t capacity,
                  const std::string& message) {
  if (!destination || capacity == 0) return;
  const std::size_t count = std::min(capacity - 1, message.size());
  std::memcpy(destination, message.data(), count);
  destination[count] = '\0';
}

void clear_message(char* destination, std::size_t capacity) {
  if (destination && capacity > 0) destination[0] = '\0';
}

void cuda_check(cudaError_t status, const char* operation) {
  if (status != cudaSuccess)
    throw std::runtime_error(std::string(operation) + ": " +
                             cudaGetErrorString(status));
}

void cublas_check(cublasStatus_t status, const char* operation) {
  if (status != CUBLAS_STATUS_SUCCESS)
    throw std::runtime_error(std::string(operation) +
                             " failed with status " +
                             std::to_string(static_cast<int>(status)));
}

void cusolver_check(cusolverStatus_t status, const char* operation) {
  if (status != CUSOLVER_STATUS_SUCCESS)
    throw std::runtime_error(std::string(operation) +
                             " failed with status " +
                             std::to_string(static_cast<int>(status)));
}

void require(bool condition, const char* message) {
  if (!condition) throw std::invalid_argument(message);
}

bool finite_matrix(const double* values, std::size_t count) {
  if (!values) return false;
  for (std::size_t i = 0; i < count; ++i)
    if (!std::isfinite(values[i])) return false;
  return true;
}

std::vector<float> float_values(const Matrix& matrix) {
  std::vector<float> result(matrix.size());
  for (std::size_t i = 0; i < matrix.size(); ++i)
    result[i] = static_cast<float>(matrix.data()[i]);
  return result;
}

Matrix bind_rows(const Matrix& top, const Matrix& bottom) {
  require(top.cols() == bottom.cols(), "Matrix column counts must match.");
  Matrix result(top.rows() + bottom.rows(), top.cols());
  std::copy(top.data(), top.data() + top.size(), result.data());
  std::copy(bottom.data(), bottom.data() + bottom.size(),
            result.data() + top.size());
  return result;
}

Matrix bind_blocks(const std::vector<Matrix>& blocks) {
  require(!blocks.empty(), "At least one matrix block is required.");
  int rows = 0;
  const int cols = blocks.front().cols();
  for (const Matrix& block : blocks) {
    require(block.cols() == cols, "Matrix block column counts must match.");
    rows += block.rows();
  }
  Matrix result(rows, cols);
  std::size_t offset = 0;
  for (const Matrix& block : blocks) {
    std::copy(block.data(), block.data() + block.size(),
              result.data() + offset);
    offset += block.size();
  }
  return result;
}

double correlation(const Matrix& left, int left_row,
                   const Matrix& right, int right_row, double beta) {
  double value = 1.0;
  const double inverse_beta = 1.0 / beta;
  for (int column = 0; column < left.cols(); ++column) {
    const double difference = left(left_row, column) - right(right_row, column);
    value *= 1.0 / (1.0 + difference * difference * inverse_beta);
  }
  return value;
}

bool cholesky_lower(const Matrix& input, Matrix& output) {
  require(input.rows() == input.cols(), "Cholesky input must be square.");
  const int n = input.rows();
  output = Matrix(n, n);
  for (int row = 0; row < n; ++row) {
    for (int col = 0; col <= row; ++col) {
      double value = input(row, col);
      for (int k = 0; k < col; ++k)
        value -= output(row, k) * output(col, k);
      if (row == col) {
        if (!(value > 0.0) || !std::isfinite(value)) return false;
        output(row, col) = std::sqrt(value);
      } else {
        output(row, col) = value / output(col, col);
      }
    }
  }
  return true;
}

Matrix solve_lower(const Matrix& lower, const Matrix& rhs) {
  require(lower.rows() == lower.cols(), "Triangular matrix must be square.");
  require(lower.rows() == rhs.rows(), "Triangular solve dimensions mismatch.");
  Matrix result(rhs.rows(), rhs.cols());
  for (int col = 0; col < rhs.cols(); ++col) {
    for (int row = 0; row < rhs.rows(); ++row) {
      double value = rhs(row, col);
      for (int k = 0; k < row; ++k)
        value -= lower(row, k) * result(k, col);
      result(row, col) = value / lower(row, row);
    }
  }
  return result;
}

Matrix shifted_block(const Matrix& block, int factor, double delta) {
  Matrix result = block;
  for (int row = 0; row < result.rows(); ++row)
    result(row, factor) += result(row, factor) < 0.5 ? delta : -delta;
  return result;
}

Matrix standard_design(const Matrix& initial, int p, double delta) {
  Matrix result(initial.rows() * (p + 1), initial.cols());
  for (int row = 0; row < initial.rows(); ++row)
    for (int col = 0; col < initial.cols(); ++col)
      result(row, col) = initial(row, col);
  for (int factor = 0; factor < p; ++factor) {
    const Matrix shifted = shifted_block(initial, factor, delta);
    for (int row = 0; row < shifted.rows(); ++row)
      for (int col = 0; col < shifted.cols(); ++col)
        result((factor + 1) * initial.rows() + row, col) = shifted(row, col);
  }
  return result;
}

bool gpu_close(double left, double right) {
  const double tolerance = std::sqrt(std::numeric_limits<float>::epsilon()) *
    std::max({1.0, std::abs(left), std::abs(right)});
  return std::abs(left - right) <= tolerance;
}

__global__ void correlation_matrix_kernel(const float* x, float* r, int n,
                                           int p, float alpha, float beta,
                                           float nugget) {
  const int row = blockIdx.y * blockDim.y + threadIdx.y;
  const int col = blockIdx.x * blockDim.x + threadIdx.x;
  if (row >= n || col >= n) return;
  float log_value = 0.0f;
  for (int k = 0; k < p; ++k) {
    const float difference = x[row * p + k] - x[col * p + k];
    log_value -= alpha * log1pf(difference * difference / beta);
  }
  float value = expf(log_value);
  if (row == col) value += nugget;
  r[row + col * n] = value;
}

__global__ void log_rhs_kernel(const float* x, const float* qmc,
                               float* rhs, int n, int m, int p,
                               float alpha, float beta) {
  const int q = blockIdx.x * blockDim.x + threadIdx.x;
  const int row = blockIdx.y * blockDim.y + threadIdx.y;
  if (row >= n || q >= m) return;
  float value = 0.0f;
  for (int k = 0; k < p; ++k) {
    const float difference = x[row * p + k] - qmc[q * p + k];
    value -= alpha * log1pf(difference * difference / beta);
  }
  rhs[row + q * n] = value;
}

__global__ void scale_rhs_kernel(float* rhs, float* scales, int n, int m) {
  const int q = blockIdx.x * blockDim.x + threadIdx.x;
  if (q >= m) return;
  float maximum = -CUDART_INF_F;
  for (int row = 0; row < n; ++row)
    maximum = fmaxf(maximum, rhs[row + q * n]);
  scales[q] = maximum;
  for (int row = 0; row < n; ++row)
    rhs[row + q * n] = expf(rhs[row + q * n] - maximum);
}

__global__ void log_norm_kernel(const float* solved, const float* scales,
                                float* logs, int n, int m) {
  const int q = blockIdx.x * blockDim.x + threadIdx.x;
  if (q >= m) return;
  float norm = 0.0f;
  for (int row = 0; row < n; ++row) {
    const float value = solved[row + q * n];
    norm = fmaf(value, value, norm);
  }
  logs[q] = 2.0f * scales[q] + logf(fmaxf(norm, FLT_MIN));
}

__global__ void reduce_min_kernel(const float* values, float* minimum, int n) {
  if (blockIdx.x != 0 || threadIdx.x != 0) return;
  float answer = CUDART_INF_F;
  for (int i = 0; i < n; ++i) answer = fminf(answer, values[i]);
  minimum[0] = answer;
}

__global__ void incremental_initialize_kernel(
    const float* x, const float* qmc, const float* lower, float* solutions,
    float* scales, float* norms, float* logs, int n, int m, int p,
    float beta) {
  const int q = blockIdx.x * blockDim.x + threadIdx.x;
  if (q >= m) return;
  float scale = -CUDART_INF_F;
  for (int row = 0; row < n; ++row) {
    float log_r = 0.0f;
    for (int dim = 0; dim < p; ++dim) {
      const float difference = x[row * p + dim] - qmc[q * p + dim];
      log_r -= log1pf(difference * difference / beta);
    }
    solutions[row * m + q] = log_r;
    scale = fmaxf(scale, log_r);
  }
  float norm = 0.0f;
  for (int row = 0; row < n; ++row) {
    float value = expf(solutions[row * m + q] - scale);
    for (int col = 0; col < row; ++col)
      value -= lower[row * n + col] * solutions[col * m + q];
    value /= lower[row * n + row];
    solutions[row * m + q] = value;
    norm = fmaf(value, value, norm);
  }
  scales[q] = scale;
  norms[q] = norm;
  logs[q] = 2.0f * scale + logf(fmaxf(norm, FLT_MIN));
}

__global__ void incremental_candidate_kernel(
    const float* block, const float* qmc, const float* w, const float* g,
    const float* old_solutions, const float* old_scales,
    const float* old_norms, float* new_scales, float* new_norms,
    float* new_z, float* logs, int n, int k, int m, int p, float beta) {
  const int q = blockIdx.x * blockDim.x + threadIdx.x;
  if (q >= m) return;
  float scale = old_scales[q];
  for (int row = 0; row < k; ++row) {
    float log_r = 0.0f;
    for (int dim = 0; dim < p; ++dim) {
      const float difference = block[row * p + dim] - qmc[q * p + dim];
      log_r -= log1pf(difference * difference / beta);
    }
    new_z[row * m + q] = log_r;
    scale = fmaxf(scale, log_r);
  }
  const float old_factor = expf(old_scales[q] - scale);
  for (int row = 0; row < k; ++row) {
    float value = expf(new_z[row * m + q] - scale);
    for (int col = 0; col < n; ++col)
      value -= w[col * k + row] * old_solutions[col * m + q] * old_factor;
    for (int col = 0; col < row; ++col)
      value -= g[row * k + col] * new_z[col * m + q];
    value /= g[row * k + row];
    new_z[row * m + q] = value;
  }
  float norm = old_norms[q] * old_factor * old_factor;
  for (int row = 0; row < k; ++row) {
    const float value = new_z[row * m + q];
    norm = fmaf(value, value, norm);
  }
  new_scales[q] = scale;
  new_norms[q] = norm;
  logs[q] = 2.0f * scale + logf(fmaxf(norm, FLT_MIN));
}

__global__ void incremental_accept_kernel(
    float* solutions, float* scales, float* norms, float* logs,
    const float* candidate_scales, const float* candidate_norms,
    const float* candidate_z, int n, int k, int m) {
  const int q = blockIdx.x * blockDim.x + threadIdx.x;
  if (q >= m) return;
  const float factor = expf(scales[q] - candidate_scales[q]);
  for (int row = 0; row < n; ++row)
    solutions[row * m + q] *= factor;
  for (int row = 0; row < k; ++row)
    solutions[(n + row) * m + q] = candidate_z[row * m + q];
  scales[q] = candidate_scales[q];
  norms[q] = candidate_norms[q];
  logs[q] = 2.0f * scales[q] + logf(fmaxf(norms[q], FLT_MIN));
}

struct ScoreResult { double minimum; };

class CudaEvaluator {
 public:
  CudaEvaluator(const Matrix& qmc, double alpha, double beta, double nugget)
      : qmc_(nullptr), x_(nullptr), factor_(nullptr), rhs_(nullptr),
        scales_(nullptr), logs_(nullptr), minimum_(nullptr), workspace_(nullptr),
        info_(nullptr), stream_(nullptr), cublas_(nullptr), solver_(nullptr),
        m_(qmc.rows()), p_(qmc.cols()), n_(0), workspace_size_(0),
        alpha_(static_cast<float>(alpha)), beta_(static_cast<float>(beta)),
        nugget_(static_cast<float>(nugget)) {
    try {
      int count = 0;
      cuda_check(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
      if (count < 1) throw std::runtime_error("No CUDA device is available.");
      cuda_check(cudaGetDevice(&device_), "cudaGetDevice");
      cudaDeviceProp properties;
      cuda_check(cudaGetDeviceProperties(&properties, device_),
                 "cudaGetDeviceProperties");
      device_name_ = properties.name;
      cuda_check(cudaStreamCreate(&stream_), "cudaStreamCreate");
      cublas_check(cublasCreate(&cublas_), "cublasCreate");
      cusolver_check(cusolverDnCreate(&solver_), "cusolverDnCreate");
      cublas_check(cublasSetStream(cublas_, stream_), "cublasSetStream");
      cusolver_check(cusolverDnSetStream(solver_, stream_),
                     "cusolverDnSetStream");
      const std::vector<float> host = float_values(qmc);
      allocate(qmc_, host.size(), "cudaMalloc QMC");
      cuda_check(cudaMemcpyAsync(qmc_, host.data(), host.size() * sizeof(float),
                                 cudaMemcpyHostToDevice, stream_),
                 "cudaMemcpy QMC");
      cuda_check(cudaStreamSynchronize(stream_), "cudaStreamSynchronize");
    } catch (...) {
      cleanup();
      throw;
    }
  }

  ~CudaEvaluator() { cleanup(); }

  void cleanup() noexcept {
    if (info_) cudaFree(info_);
    if (workspace_) cudaFree(workspace_);
    if (minimum_) cudaFree(minimum_);
    if (logs_) cudaFree(logs_);
    if (scales_) cudaFree(scales_);
    if (rhs_) cudaFree(rhs_);
    if (factor_) cudaFree(factor_);
    if (x_) cudaFree(x_);
    if (qmc_) cudaFree(qmc_);
    if (solver_) cusolverDnDestroy(solver_);
    if (cublas_) cublasDestroy(cublas_);
    if (stream_) cudaStreamDestroy(stream_);
    info_ = nullptr;
    workspace_ = minimum_ = logs_ = scales_ = rhs_ = factor_ = x_ = qmc_ =
      nullptr;
    solver_ = nullptr;
    cublas_ = nullptr;
    stream_ = nullptr;
  }

  const std::string& device_name() const { return device_name_; }

  ScoreResult score(const Matrix& design) {
    require(design.rows() > 0 && design.cols() == p_,
            "Invalid design dimensions for CUDA Q evaluation.");
    const int n = design.rows();
    ensure_resources(n);
    const std::vector<float> host_x = float_values(design);
    cuda_check(cudaMemcpyAsync(x_, host_x.data(), host_x.size() * sizeof(float),
                               cudaMemcpyHostToDevice, stream_),
               "cudaMemcpy design");
    const dim3 threads2(16, 16);
    correlation_matrix_kernel<<<dim3((n + 15) / 16, (n + 15) / 16),
                                threads2, 0, stream_>>>(
      x_, factor_, n, p_, alpha_, beta_, nugget_);
    cuda_check(cudaGetLastError(), "CUDA correlation-matrix kernel launch");
    log_rhs_kernel<<<dim3((m_ + 15) / 16, (n + 15) / 16), threads2,
                     0, stream_>>>(x_, qmc_, rhs_, n, m_, p_, alpha_, beta_);
    cuda_check(cudaGetLastError(), "CUDA log-RHS kernel launch");
    const int threads = 256;
    scale_rhs_kernel<<<(m_ + threads - 1) / threads, threads, 0, stream_>>>(
      rhs_, scales_, n, m_);
    cuda_check(cudaGetLastError(), "CUDA RHS-scaling kernel launch");
    cusolver_check(cusolverDnSpotrf(
      solver_, CUBLAS_FILL_MODE_LOWER, n, factor_, n, workspace_,
      workspace_size_, info_), "cusolverDnSpotrf");
    const float one = 1.0f;
    cublas_check(cublasStrsm(
      cublas_, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER,
      CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, m_, &one, factor_, n, rhs_, n),
      "cublasStrsm");
    log_norm_kernel<<<(m_ + threads - 1) / threads, threads, 0, stream_>>>(
      rhs_, scales_, logs_, n, m_);
    cuda_check(cudaGetLastError(), "CUDA log-norm kernel launch");
    reduce_min_kernel<<<1, 1, 0, stream_>>>(logs_, minimum_, m_);
    cuda_check(cudaGetLastError(), "CUDA minimum-reduction kernel launch");
    int host_info = 0;
    float host_minimum = std::numeric_limits<float>::infinity();
    cuda_check(cudaMemcpyAsync(&host_info, info_, sizeof(int),
                               cudaMemcpyDeviceToHost, stream_),
               "cudaMemcpy factor status");
    cuda_check(cudaMemcpyAsync(&host_minimum, minimum_, sizeof(float),
                               cudaMemcpyDeviceToHost, stream_),
               "cudaMemcpy minimum Q");
    cuda_check(cudaStreamSynchronize(stream_), "cudaStreamSynchronize");
    if (host_info != 0)
      throw std::runtime_error("CUDA Cholesky factorization failed with info " +
                               std::to_string(host_info));
    if (!std::isfinite(host_minimum))
      throw std::runtime_error("CUDA Q evaluation returned a non-finite value.");
    return ScoreResult{static_cast<double>(host_minimum)};
  }

 private:
  float *qmc_, *x_, *factor_, *rhs_, *scales_, *logs_, *minimum_, *workspace_;
  int* info_;
  cudaStream_t stream_;
  cublasHandle_t cublas_;
  cusolverDnHandle_t solver_;
  int device_, m_, p_, n_, workspace_size_;
  float alpha_, beta_, nugget_;
  std::string device_name_;

  static void allocate(float*& pointer, std::size_t count,
                       const char* operation) {
    cuda_check(cudaMalloc(reinterpret_cast<void**>(&pointer),
                          count * sizeof(float)), operation);
  }

  void ensure_resources(int n) {
    if (n_ != 0) {
      if (n_ != n)
        throw std::runtime_error("CUDA evaluator design size changed.");
      return;
    }
    n_ = n;
    allocate(x_, static_cast<std::size_t>(n) * p_, "cudaMalloc design");
    allocate(factor_, static_cast<std::size_t>(n) * n, "cudaMalloc factor");
    allocate(rhs_, static_cast<std::size_t>(n) * m_, "cudaMalloc RHS");
    allocate(scales_, m_, "cudaMalloc scales");
    allocate(logs_, m_, "cudaMalloc logs");
    allocate(minimum_, 1, "cudaMalloc minimum");
    cuda_check(cudaMalloc(reinterpret_cast<void**>(&info_), sizeof(int)),
               "cudaMalloc info");
    cusolver_check(cusolverDnSpotrf_bufferSize(
      solver_, CUBLAS_FILL_MODE_LOWER, n, factor_, n, &workspace_size_),
      "cusolverDnSpotrf_bufferSize");
    allocate(workspace_, workspace_size_, "cudaMalloc cuSOLVER workspace");
  }
};

class CudaIncrementalEvaluator {
 public:
  struct Slot {
    float *block_device = nullptr, *w_device = nullptr, *g_device = nullptr;
    float *scales = nullptr, *norms = nullptr, *z = nullptr;
    float *logs = nullptr, *minimum = nullptr;
    Matrix w, g, block;
    double min_log = -std::numeric_limits<double>::infinity();
  };

  CudaIncrementalEvaluator(const Matrix& evaluation, double beta,
                           double nugget, int maximum_n, int block_size)
      : qmc_(nullptr), solutions_(nullptr), scales_(nullptr), norms_(nullptr),
        logs_(nullptr), minimum_(nullptr), init_x_(nullptr), init_l_(nullptr),
        stream_(nullptr), m_(evaluation.rows()), p_(evaluation.cols()),
        max_n_(maximum_n), k_(block_size), n_(0), beta_(beta), nugget_(nugget),
        current_min_(-std::numeric_limits<double>::infinity()),
        pending_slot_(-1), pending_n_(0) {
    try {
      int count = 0;
      cuda_check(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
      if (count < 1) throw std::runtime_error("No CUDA device is available.");
      int device = 0;
      cuda_check(cudaGetDevice(&device), "cudaGetDevice");
      cudaDeviceProp properties;
      cuda_check(cudaGetDeviceProperties(&properties, device),
                 "cudaGetDeviceProperties");
      device_name_ = properties.name;
      cuda_check(cudaStreamCreate(&stream_), "cudaStreamCreate");
      const std::vector<float> qmc_host = float_values(evaluation);
      allocate(qmc_, qmc_host.size(), "CUDA incremental QMC");
      allocate(solutions_, static_cast<std::size_t>(max_n_) * m_,
               "CUDA incremental solutions");
      allocate(scales_, m_, "CUDA incremental scales");
      allocate(norms_, m_, "CUDA incremental norms");
      allocate(logs_, m_, "CUDA incremental logs");
      allocate(minimum_, 1, "CUDA incremental minimum");
      allocate(init_x_, static_cast<std::size_t>(max_n_) * p_,
               "CUDA initial design workspace");
      allocate(init_l_, static_cast<std::size_t>(max_n_) * max_n_,
               "CUDA initial factor workspace");
      for (Slot& slot : slots_) {
        allocate(slot.block_device, static_cast<std::size_t>(k_) * p_,
                 "CUDA candidate block workspace");
        allocate(slot.w_device, static_cast<std::size_t>(max_n_) * k_,
                 "CUDA candidate W workspace");
        allocate(slot.g_device, static_cast<std::size_t>(k_) * k_,
                 "CUDA candidate G workspace");
        allocate(slot.scales, m_, "CUDA candidate scales");
        allocate(slot.norms, m_, "CUDA candidate norms");
        allocate(slot.z, static_cast<std::size_t>(k_) * m_,
                 "CUDA candidate solutions");
        allocate(slot.logs, m_, "CUDA candidate logs");
        allocate(slot.minimum, 1, "CUDA candidate minimum");
      }
      cuda_check(cudaMemcpyAsync(qmc_, qmc_host.data(),
                                 qmc_host.size() * sizeof(float),
                                 cudaMemcpyHostToDevice, stream_),
                 "cudaMemcpy incremental QMC");
      cuda_check(cudaStreamSynchronize(stream_), "cudaStreamSynchronize");
    } catch (...) {
      cleanup();
      throw;
    }
  }

  ~CudaIncrementalEvaluator() { cleanup(); }

  void cleanup() noexcept {
    for (Slot& slot : slots_) {
      if (slot.minimum) cudaFree(slot.minimum);
      if (slot.logs) cudaFree(slot.logs);
      if (slot.z) cudaFree(slot.z);
      if (slot.norms) cudaFree(slot.norms);
      if (slot.scales) cudaFree(slot.scales);
      if (slot.g_device) cudaFree(slot.g_device);
      if (slot.w_device) cudaFree(slot.w_device);
      if (slot.block_device) cudaFree(slot.block_device);
    }
    if (init_l_) cudaFree(init_l_);
    if (init_x_) cudaFree(init_x_);
    if (minimum_) cudaFree(minimum_);
    if (logs_) cudaFree(logs_);
    if (norms_) cudaFree(norms_);
    if (scales_) cudaFree(scales_);
    if (solutions_) cudaFree(solutions_);
    if (qmc_) cudaFree(qmc_);
    if (stream_) cudaStreamDestroy(stream_);
    for (Slot& slot : slots_) {
      slot.minimum = slot.logs = slot.z = slot.norms = slot.scales = nullptr;
      slot.g_device = slot.w_device = slot.block_device = nullptr;
    }
    init_l_ = init_x_ = minimum_ = logs_ = norms_ = scales_ = solutions_ =
      qmc_ = nullptr;
    stream_ = nullptr;
  }

  double current_min() const { return current_min_; }
  double candidate_min(int slot_index) const { return slots_[slot_index].min_log; }

  void initialize(const Matrix& design) {
    require(design.rows() > 0 && design.cols() == p_,
            "Invalid initial CUDA design dimensions.");
    require(design.rows() <= max_n_, "CUDA incremental design exceeds capacity.");
    pending_slot_ = -1;
    x_ = design;
    n_ = design.rows();
    build_factor(x_, lower_);
    const std::vector<float> x_host = float_values(x_);
    const std::vector<float> l_host = float_values(lower_);
    cuda_check(cudaMemcpyAsync(init_x_, x_host.data(),
                               x_host.size() * sizeof(float),
                               cudaMemcpyHostToDevice, stream_),
               "cudaMemcpy initial design");
    cuda_check(cudaMemcpyAsync(init_l_, l_host.data(),
                               l_host.size() * sizeof(float),
                               cudaMemcpyHostToDevice, stream_),
               "cudaMemcpy initial factor");
    const int threads = 256;
    incremental_initialize_kernel<<<(m_ + threads - 1) / threads, threads,
                                    0, stream_>>>(
      init_x_, qmc_, init_l_, solutions_, scales_, norms_, logs_, n_, m_, p_,
      static_cast<float>(beta_));
    cuda_check(cudaGetLastError(),
               "CUDA incremental-initialization kernel launch");
    reduce_min_kernel<<<1, 1, 0, stream_>>>(logs_, minimum_, m_);
    cuda_check(cudaGetLastError(),
               "CUDA incremental minimum-reduction kernel launch");
    float minimum = 0.0f;
    cuda_check(cudaMemcpyAsync(&minimum, minimum_, sizeof(float),
                               cudaMemcpyDeviceToHost, stream_),
               "cudaMemcpy initial minimum");
    cuda_check(cudaStreamSynchronize(stream_), "cudaStreamSynchronize");
    if (!std::isfinite(minimum))
      throw std::runtime_error("CUDA initial incremental score is non-finite.");
    current_min_ = minimum;
  }

  bool evaluate(const Matrix& block, int slot_index) {
    if (!prepare_slot(block, slot_index)) return false;
    flush_pending_accept();
    launch_candidate(slot_index);
    Slot& slot = slots_[slot_index];
    float minimum = 0.0f;
    cuda_check(cudaMemcpyAsync(&minimum, slot.minimum, sizeof(float),
                               cudaMemcpyDeviceToHost, stream_),
               "cudaMemcpy candidate minimum");
    cuda_check(cudaStreamSynchronize(stream_), "cudaStreamSynchronize");
    if (!std::isfinite(minimum)) return false;
    slot.min_log = minimum;
    return true;
  }

  bool evaluate_pair(const Matrix& left, const Matrix& right) {
    if (!prepare_slot(left, 0) || !prepare_slot(right, 1)) return false;
    flush_pending_accept();
    launch_candidate(0);
    launch_candidate(1);
    float minima[2] = {0.0f, 0.0f};
    cuda_check(cudaMemcpyAsync(&minima[0], slots_[0].minimum, sizeof(float),
                               cudaMemcpyDeviceToHost, stream_),
               "cudaMemcpy left candidate minimum");
    cuda_check(cudaMemcpyAsync(&minima[1], slots_[1].minimum, sizeof(float),
                               cudaMemcpyDeviceToHost, stream_),
               "cudaMemcpy right candidate minimum");
    cuda_check(cudaStreamSynchronize(stream_), "cudaStreamSynchronize");
    if (!std::isfinite(minima[0]) || !std::isfinite(minima[1])) return false;
    slots_[0].min_log = minima[0];
    slots_[1].min_log = minima[1];
    return true;
  }

  void accept(int slot_index, bool needs_gpu_state) {
    if (needs_gpu_state) {
      pending_slot_ = slot_index;
      pending_n_ = n_;
    } else {
      pending_slot_ = -1;
    }
    update_host_state(slot_index);
  }

  void accept_now(int slot_index) {
    launch_accept(slot_index, n_);
    cuda_check(cudaStreamSynchronize(stream_), "cudaStreamSynchronize");
    pending_slot_ = -1;
    update_host_state(slot_index);
  }

  std::vector<float> log_values() const {
    std::vector<float> answer(m_);
    cuda_check(cudaMemcpy(answer.data(), logs_, m_ * sizeof(float),
                          cudaMemcpyDeviceToHost),
               "cudaMemcpy incremental log values");
    return answer;
  }

 private:
  float *qmc_, *solutions_, *scales_, *norms_, *logs_, *minimum_;
  float *init_x_, *init_l_;
  cudaStream_t stream_;
  Slot slots_[2];
  int m_, p_, max_n_, k_, n_;
  double beta_, nugget_, current_min_;
  int pending_slot_, pending_n_;
  std::string device_name_;
  Matrix x_, lower_;

  static void allocate(float*& pointer, std::size_t count,
                       const char* operation) {
    cuda_check(cudaMalloc(reinterpret_cast<void**>(&pointer),
                          count * sizeof(float)), operation);
  }

  void build_factor(const Matrix& x, Matrix& lower) const {
    Matrix correlation_matrix(x.rows(), x.rows());
    for (int row = 0; row < x.rows(); ++row) {
      correlation_matrix(row, row) = 1.0 + nugget_;
      for (int col = row + 1; col < x.rows(); ++col) {
        const double value = correlation(x, row, x, col, beta_);
        correlation_matrix(row, col) = value;
        correlation_matrix(col, row) = value;
      }
    }
    if (!cholesky_lower(correlation_matrix, lower))
      throw std::runtime_error("Initial CUDA incremental Cholesky failed.");
  }

  bool build_block_factor(const Matrix& block, Matrix& w, Matrix& g) const {
    const int block_n = block.rows();
    Matrix c(n_, block_n);
    Matrix d(block_n, block_n);
    for (int row = 0; row < n_; ++row)
      for (int col = 0; col < block_n; ++col)
        c(row, col) = correlation(x_, row, block, col, beta_);
    for (int row = 0; row < block_n; ++row) {
      d(row, row) = 1.0 + nugget_;
      for (int col = row + 1; col < block_n; ++col) {
        const double value = correlation(block, row, block, col, beta_);
        d(row, col) = value;
        d(col, row) = value;
      }
    }
    w = solve_lower(lower_, c);
    Matrix schur(block_n, block_n);
    for (int row = 0; row < block_n; ++row) {
      for (int col = 0; col < block_n; ++col) {
        double value = d(row, col);
        for (int k = 0; k < n_; ++k) value -= w(k, row) * w(k, col);
        schur(row, col) = value;
      }
    }
    return cholesky_lower(schur, g);
  }

  bool prepare_slot(const Matrix& block, int slot_index) {
    require(slot_index >= 0 && slot_index < 2, "Invalid CUDA slot index.");
    require(block.rows() == k_ && block.cols() == p_,
            "Invalid CUDA candidate block dimensions.");
    require(n_ + block.rows() <= max_n_, "CUDA candidate exceeds capacity.");
    Slot& slot = slots_[slot_index];
    slot.block = block;
    if (!build_block_factor(block, slot.w, slot.g)) return false;
    const std::vector<float> block_host = float_values(block);
    const std::vector<float> w_host = float_values(slot.w);
    const std::vector<float> g_host = float_values(slot.g);
    cuda_check(cudaMemcpyAsync(slot.block_device, block_host.data(),
                               block_host.size() * sizeof(float),
                               cudaMemcpyHostToDevice, stream_),
               "cudaMemcpy candidate block");
    cuda_check(cudaMemcpyAsync(slot.w_device, w_host.data(),
                               w_host.size() * sizeof(float),
                               cudaMemcpyHostToDevice, stream_),
               "cudaMemcpy candidate W");
    cuda_check(cudaMemcpyAsync(slot.g_device, g_host.data(),
                               g_host.size() * sizeof(float),
                               cudaMemcpyHostToDevice, stream_),
               "cudaMemcpy candidate G");
    return true;
  }

  void launch_candidate(int slot_index) {
    Slot& slot = slots_[slot_index];
    const int threads = 256;
    incremental_candidate_kernel<<<(m_ + threads - 1) / threads, threads,
                                   0, stream_>>>(
      slot.block_device, qmc_, slot.w_device, slot.g_device, solutions_,
      scales_, norms_, slot.scales, slot.norms, slot.z, slot.logs, n_,
      slot.block.rows(), m_, p_, static_cast<float>(beta_));
    cuda_check(cudaGetLastError(), "CUDA candidate kernel launch");
    reduce_min_kernel<<<1, 1, 0, stream_>>>(slot.logs, slot.minimum, m_);
    cuda_check(cudaGetLastError(),
               "CUDA candidate minimum-reduction kernel launch");
  }

  void launch_accept(int slot_index, int accepted_n) {
    Slot& slot = slots_[slot_index];
    const int threads = 256;
    incremental_accept_kernel<<<(m_ + threads - 1) / threads, threads,
                                0, stream_>>>(
      solutions_, scales_, norms_, logs_, slot.scales, slot.norms, slot.z,
      accepted_n, slot.block.rows(), m_);
    cuda_check(cudaGetLastError(), "CUDA incremental accept kernel");
  }

  void flush_pending_accept() {
    if (pending_slot_ < 0) return;
    launch_accept(pending_slot_, pending_n_);
    pending_slot_ = -1;
  }

  void update_host_state(int slot_index) {
    Slot& slot = slots_[slot_index];
    const int block_n = slot.block.rows();
    Matrix updated(n_ + block_n, n_ + block_n);
    for (int row = 0; row < n_; ++row)
      for (int col = 0; col < n_; ++col)
        updated(row, col) = lower_(row, col);
    for (int row = 0; row < block_n; ++row) {
      for (int col = 0; col < n_; ++col)
        updated(n_ + row, col) = slot.w(col, row);
      for (int col = 0; col < block_n; ++col)
        updated(n_ + row, n_ + col) = slot.g(row, col);
    }
    lower_ = std::move(updated);
    x_ = bind_rows(x_, slot.block);
    n_ += block_n;
    current_min_ = slot.min_log;
  }
};

struct StrictResult { Matrix design; double score; };

StrictResult strict_design(const Matrix& initial, int p, double delta,
                           CudaIncrementalEvaluator& evaluator) {
  std::vector<Matrix> blocks{initial, shifted_block(initial, 0, delta)};
  evaluator.initialize(bind_rows(blocks[0], blocks[1]));
  for (int factor = 1; factor < p; ++factor) {
    const Matrix left_block = shifted_block(blocks.front(), factor, delta);
    const Matrix right_block = shifted_block(blocks.back(), factor, delta);
    if (!evaluator.evaluate_pair(left_block, right_block))
      throw std::runtime_error("CUDA strict incremental factorization failed.");
    const double left_score = evaluator.candidate_min(0);
    const double right_score = evaluator.candidate_min(1);
    if (left_score > right_score || gpu_close(left_score, right_score)) {
      blocks.insert(blocks.begin(), left_block);
      evaluator.accept(0, factor + 1 < p);
    } else {
      blocks.push_back(right_block);
      evaluator.accept(1, factor + 1 < p);
    }
  }
  return StrictResult{bind_blocks(blocks), evaluator.current_min()};
}

int exception_status(const std::exception& error) {
  if (dynamic_cast<const std::invalid_argument*>(&error))
    return MOFAT_CUDA_INVALID_ARGUMENT;
  const std::string message = error.what();
  if (message.find("Cholesky") != std::string::npos ||
      message.find("factorization") != std::string::npos ||
      message.find("numerically") != std::string::npos)
    return MOFAT_CUDA_NUMERICAL_ERROR;
  if (message.find("No CUDA device") != std::string::npos ||
      message.find("cudaGetDeviceCount") != std::string::npos)
    return MOFAT_CUDA_UNAVAILABLE;
  return MOFAT_CUDA_RUNTIME_ERROR;
}

template <typename Function>
int guarded_call(Function function, char* error, std::size_t error_capacity) {
  clear_message(error, error_capacity);
  try {
    function();
    return MOFAT_CUDA_SUCCESS;
  } catch (const std::exception& caught) {
    copy_message(error, error_capacity, caught.what());
    return exception_status(caught);
  } catch (...) {
    copy_message(error, error_capacity, "Unknown CUDA helper failure.");
    return MOFAT_CUDA_RUNTIME_ERROR;
  }
}

}  // namespace

extern "C" MOFAT_CUDA_EXPORT int mofat_cuda_helper_abi_version(void) {
  return MOFAT_CUDA_HELPER_ABI_VERSION;
}

extern "C" MOFAT_CUDA_EXPORT int mofat_cuda_helper_info(
    int* toolkit_version, int* runtime_version, int* driver_version,
    int* compute_major, int* compute_minor,
    char* device, size_t device_capacity,
    char* error, size_t error_capacity) {
  return guarded_call([&]() {
    require(toolkit_version && runtime_version && driver_version &&
              compute_major && compute_minor,
            "CUDA info output pointers must not be null.");
    require(device && device_capacity > 0,
            "CUDA device output buffer must not be null or empty.");
    int count = 0;
    cuda_check(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
    if (count < 1) throw std::runtime_error("No CUDA device is available.");
    int current = 0;
    cuda_check(cudaGetDevice(&current), "cudaGetDevice");
    cudaDeviceProp properties;
    cuda_check(cudaGetDeviceProperties(&properties, current),
               "cudaGetDeviceProperties");
    *toolkit_version = CUDART_VERSION;
    cuda_check(cudaRuntimeGetVersion(runtime_version), "cudaRuntimeGetVersion");
    cuda_check(cudaDriverGetVersion(driver_version), "cudaDriverGetVersion");
    *compute_major = properties.major;
    *compute_minor = properties.minor;
    copy_message(device, device_capacity, properties.name);
  }, error, error_capacity);
}

extern "C" MOFAT_CUDA_EXPORT int mofat_cuda_helper_self_test(
    char* error, size_t error_capacity) {
  return guarded_call([&]() {
    const double qmc_values[] = {0.25, 0.25, 0.75, 0.75};
    const double design_values[] = {0.20, 0.30, 0.80, 0.70};
    Matrix qmc(qmc_values, 2, 2);
    Matrix design(design_values, 2, 2);
    CudaEvaluator evaluator(qmc, 1.0, 0.05, 1e-6);
    const ScoreResult result = evaluator.score(design);
    if (!std::isfinite(result.minimum))
      throw std::runtime_error("CUDA self-test produced a non-finite result.");
  }, error, error_capacity);
}

extern "C" MOFAT_CUDA_EXPORT int mofat_cuda_helper_q(
    const double* qmc, int qmc_rows, int p,
    const double* design, int design_rows,
    double alpha, double beta, double nugget,
    double* log_q, char* error, size_t error_capacity) {
  return guarded_call([&]() {
    require(qmc_rows > 0 && design_rows > 0 && p > 0,
            "CUDA Q dimensions must be positive.");
    require(log_q, "CUDA Q output pointer must not be null.");
    require(std::isfinite(alpha) && alpha > 0.0,
            "alpha must be positive and finite.");
    require(std::isfinite(beta) && beta > 0.0,
            "beta must be positive and finite.");
    require(std::isfinite(nugget) && nugget >= 0.0,
            "nugget must be finite and nonnegative.");
    require(finite_matrix(qmc, static_cast<std::size_t>(qmc_rows) * p),
            "qmc must contain finite values.");
    require(finite_matrix(design, static_cast<std::size_t>(design_rows) * p),
            "design must contain finite values.");
    Matrix qmc_matrix(qmc, qmc_rows, p);
    Matrix design_matrix(design, design_rows, p);
    CudaEvaluator evaluator(qmc_matrix, alpha, beta, nugget);
    *log_q = evaluator.score(design_matrix).minimum;
  }, error, error_capacity);
}

extern "C" MOFAT_CUDA_EXPORT int mofat_cuda_helper_optimize(
    const double* qmc, int qmc_rows,
    const double* initial_a, int p, int l,
    double delta, int remainder, int quotient,
    int strict_structure, double nugget,
    double* output_design, double* log_q,
    char* error, size_t error_capacity) {
  return guarded_call([&]() {
    require(qmc_rows > 0 && p >= 2 && l >= 2,
            "CUDA optimization dimensions are invalid.");
    require(output_design && log_q,
            "CUDA optimization output pointers must not be null.");
    require(std::isfinite(delta) && delta > 0.0,
            "delta must be positive and finite.");
    require(std::isfinite(nugget) && nugget >= 0.0,
            "nugget must be finite and nonnegative.");
    require(finite_matrix(qmc, static_cast<std::size_t>(qmc_rows) * p),
            "qmc must contain finite values.");
    require(finite_matrix(initial_a, static_cast<std::size_t>(l) * p),
            "initial_A must contain finite values.");
    const Matrix qmc_matrix(qmc, qmc_rows, p);
    Matrix a(initial_a, l, p);
    const int p0 = remainder < 0 ? p : remainder;
    const int p1 = quotient < 0 ? 0 :
      static_cast<int>(std::tgamma(static_cast<double>(l) + 1.0)) * quotient;
    require(p0 >= 0 && p1 >= 0 && p1 + std::max(0, p0 - 1) < p,
            "CUDA optimization swap columns are out of range.");
    const double beta = 1.0 / std::pow(2.0 * l, 2.0);
    StrictResult current;
    if (strict_structure) {
      CudaIncrementalEvaluator evaluator(qmc_matrix, beta, nugget,
                                         l * (p + 1), l);
      current = strict_design(a, p, delta, evaluator);
      if (p0 > 1) {
        for (int i = 0; i < l - 1; ++i) {
          for (int j = i + 1; j < l; ++j) {
            for (int k = 1; k < p0; ++k) {
              Matrix candidate_a = a;
              std::swap(candidate_a(i, p1 + k), candidate_a(j, p1 + k));
              StrictResult proposed = strict_design(candidate_a, p, delta,
                                                    evaluator);
              if (proposed.score > current.score &&
                  !gpu_close(proposed.score, current.score)) {
                a = std::move(candidate_a);
                current = std::move(proposed);
              }
            }
          }
        }
      }
    } else {
      CudaEvaluator evaluator(qmc_matrix, 1.0, beta, nugget);
      current.design = standard_design(a, p, delta);
      current.score = evaluator.score(current.design).minimum;
      if (p0 > 1) {
        for (int i = 0; i < l - 1; ++i) {
          for (int j = i + 1; j < l; ++j) {
            for (int k = 1; k < p0; ++k) {
              Matrix candidate_a = a;
              std::swap(candidate_a(i, p1 + k), candidate_a(j, p1 + k));
              StrictResult proposed;
              proposed.design = standard_design(candidate_a, p, delta);
              proposed.score = evaluator.score(proposed.design).minimum;
              if (proposed.score > current.score &&
                  !gpu_close(proposed.score, current.score)) {
                a = std::move(candidate_a);
                current = std::move(proposed);
              }
            }
          }
        }
      }
    }
    std::copy(current.design.data(), current.design.data() + current.design.size(),
              output_design);
    *log_q = current.score;
  }, error, error_capacity);
}

extern "C" MOFAT_CUDA_EXPORT int mofat_cuda_helper_augment(
    const double* design, int design_rows,
    const double* candidates, int candidate_rows,
    int p, int n_aug, double beta, double nugget,
    double* output_design,
    char* error, size_t error_capacity) {
  return guarded_call([&]() {
    require(design_rows > 0 && candidate_rows > 0 && p > 0,
            "CUDA augmentation dimensions must be positive.");
    require(n_aug > 0 && n_aug <= candidate_rows,
            "n_aug must be between 1 and the number of candidates.");
    require(output_design, "CUDA augmentation output pointer must not be null.");
    require(std::isfinite(beta) && beta > 0.0,
            "beta must be positive and finite.");
    require(std::isfinite(nugget) && nugget >= 0.0,
            "nugget must be finite and nonnegative.");
    require(finite_matrix(design, static_cast<std::size_t>(design_rows) * p),
            "design must contain finite values.");
    require(finite_matrix(candidates,
                          static_cast<std::size_t>(candidate_rows) * p),
            "candidates must contain finite values.");
    Matrix current(design, design_rows, p);
    const Matrix candidate_matrix(candidates, candidate_rows, p);
    CudaIncrementalEvaluator evaluator(candidate_matrix, beta, nugget,
                                       design_rows + n_aug, 1);
    evaluator.initialize(current);
    std::vector<bool> eligible(candidate_rows, true);
    for (int step = 0; step < n_aug; ++step) {
      while (true) {
        const std::vector<float> scores = evaluator.log_values();
        int selected = candidate_rows;
        float best = std::numeric_limits<float>::infinity();
        for (int row = 0; row < candidate_rows; ++row) {
          if (eligible[row] && std::isfinite(scores[row]) && scores[row] < best) {
            selected = row;
            best = scores[row];
          }
        }
        if (selected == candidate_rows)
          throw std::runtime_error(
            "No numerically eligible CUDA candidate remains.");
        Matrix point(1, p);
        for (int col = 0; col < p; ++col)
          point(0, col) = candidate_matrix(selected, col);
        if (evaluator.evaluate(point, 0)) {
          evaluator.accept_now(0);
          current = bind_rows(current, point);
          eligible[selected] = false;
          break;
        }
        eligible[selected] = false;
      }
    }
    std::copy(current.data(), current.data() + current.size(), output_design);
  }, error, error_capacity);
}
