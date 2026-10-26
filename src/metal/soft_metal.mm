// Copied into src and compiled as Objective-C++ on supported macOS hosts.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
#include "soft_gpu_api.h"
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct CorrParams {
  uint32_t n;
  uint32_t m;
  uint32_t p;
  uint32_t stride;
  float beta;
  float nugget;
  float alpha;
};

struct IncrementalParams {
  uint32_t n;
  uint32_t k;
  uint32_t m;
  uint32_t p;
  uint32_t y_stride;
  float beta;
  float alpha;
};

struct AugmentParams {
  uint32_t n;
  uint32_t m;
  uint32_t p;
  uint32_t y_stride;
  float beta;
  float g;
  float alpha;
};

struct ScoreResult {
  double min_log;
  double normalized(int p) const { return std::exp(min_log / p); }
};

bool gpu_close(double a, double b) {
  const double tolerance = std::sqrt(std::numeric_limits<float>::epsilon()) *
    std::max({1.0, std::abs(a), std::abs(b)});
  return std::abs(a - b) <= tolerance;
}

std::string read_text_file(const std::string& path) {
  std::ifstream input(path);
  if (!input) Rcpp::stop("Cannot read Metal kernel file: " + path);
  std::ostringstream text;
  text << input.rdbuf();
  return text.str();
}

std::vector<float> row_major_float(const arma::mat& x) {
  std::vector<float> out(x.n_elem);
  for (arma::uword i = 0; i < x.n_rows; ++i)
    for (arma::uword j = 0; j < x.n_cols; ++j)
      out[i * x.n_cols + j] = static_cast<float>(x(i, j));
  return out;
}

double correlation_double(const arma::rowvec& a, const arma::rowvec& b,
                          double beta, double alpha = 1.0) {
  // beta = theta^2
  const double inv_beta = 1.0 / beta;
  if (alpha == 1.0) {
    double value = 1.0;
    for (arma::uword k = 0; k < a.n_elem; ++k) {
      const double d = a[k] - b[k];
      value *= 1.0 / (1.0 + d * d * inv_beta);
    }
    return value;
  }
  double log_value = 0.0;
  for (arma::uword k = 0; k < a.n_elem; ++k) {
    const double d = a[k] - b[k];
    log_value -= alpha * std::log1p(d * d * inv_beta);
  }
  return std::exp(log_value);
}

class MetalEvaluator {
 public:
  MetalEvaluator(const arma::mat& qmc, double alpha, double beta, double nugget,
                 const std::string& kernel_path)
      : device_(nil), queue_(nil), library_(nil), corr_xx_(nil),
        corr_xqmc_(nil), scale_columns_(nil), mmse_(nil), reduce_(nil),
        qmc_buffer_(nil), x_buffer_(nil), r_buffer_(nil), factor_buffer_(nil),
        rhs_buffer_(nil), solution_buffer_(nil), values_buffer_(nil),
        scale_buffer_(nil), status_buffer_(nil), minimum_buffer_(nil),
        r_matrix_(nil), factor_matrix_(nil), rhs_matrix_(nil),
        solution_matrix_(nil), chol_(nil), solve_(nil), n_(0),
        r_row_bytes_(0), rhs_row_bytes_(0),
        m_(qmc.n_rows), p_(qmc.n_cols), alpha_(alpha), beta_(beta), nugget_(nugget),
        evaluations_(0) {
    device_ = MTLCreateSystemDefaultDevice();
    if (device_ == nil)
      Rcpp::stop("Metal GPU is unavailable in this process.");
    queue_ = [device_ newCommandQueue];
    if (queue_ == nil) Rcpp::stop("Could not create Metal command queue.");

    const std::string source = read_text_file(kernel_path);
    NSString* source_ns = [NSString stringWithUTF8String:source.c_str()];
    NSError* error = nil;
    MTLCompileOptions* options = [[MTLCompileOptions alloc] init];
    library_ = [device_ newLibraryWithSource:source_ns options:options error:&error];
    [options release];
    if (library_ == nil) {
      const char* message = error ? [[error localizedDescription] UTF8String]
                                  : "unknown Metal compilation error";
      Rcpp::stop(std::string("Metal kernel compilation failed: ") + message);
    }

    corr_xx_ = make_pipeline("correlation_xx");
    corr_xqmc_ = make_pipeline("log_correlation_xqmc");
    scale_columns_ = make_pipeline("scale_log_columns");
    mmse_ = make_pipeline("log_mmse_values");
    reduce_ = make_pipeline("reduce_min_log_values");

    const std::vector<float> qmc_data = row_major_float(qmc);
    qmc_buffer_ = [device_ newBufferWithBytes:qmc_data.data()
                                       length:qmc_data.size() * sizeof(float)
                                      options:MTLResourceStorageModeShared];
    if (qmc_buffer_ == nil) Rcpp::stop("Could not allocate QMC Metal buffer.");
  }

  ~MetalEvaluator() {
    [solve_ release];
    [chol_ release];
    [solution_matrix_ release];
    [rhs_matrix_ release];
    [factor_matrix_ release];
    [r_matrix_ release];
    [minimum_buffer_ release];
    [status_buffer_ release];
    [scale_buffer_ release];
    [values_buffer_ release];
    [solution_buffer_ release];
    [rhs_buffer_ release];
    [factor_buffer_ release];
    [r_buffer_ release];
    [x_buffer_ release];
    [qmc_buffer_ release];
    [reduce_ release];
    [mmse_ release];
    [scale_columns_ release];
    [corr_xqmc_ release];
    [corr_xx_ release];
    [library_ release];
    [queue_ release];
    [device_ release];
  }

  std::string device_name() const {
    return std::string([[device_ name] UTF8String]);
  }

  std::size_t evaluations() const { return evaluations_; }

  ScoreResult score(const arma::mat& X) {
    if (X.n_cols != p_) Rcpp::stop("X and QMC column counts differ.");
    if (X.n_rows == 0 || m_ == 0) Rcpp::stop("X and QMC must be nonempty.");
    ++evaluations_;

    const NSUInteger n = X.n_rows;
    const NSUInteger m = m_;
    ensure_resources(n);
    const std::vector<float> x_data = row_major_float(X);
    std::memcpy([x_buffer_ contents], x_data.data(),
                x_data.size() * sizeof(float));
    *static_cast<int32_t*>([status_buffer_ contents]) = 0;

    CorrParams rpar{static_cast<uint32_t>(n), static_cast<uint32_t>(m),
                    static_cast<uint32_t>(p_), static_cast<uint32_t>(r_row_bytes_ / sizeof(float)),
                    static_cast<float>(beta_), static_cast<float>(nugget_),
                    static_cast<float>(alpha_)};
    CorrParams qpar = rpar;
    qpar.stride = static_cast<uint32_t>(rhs_row_bytes_ / sizeof(float));

    id<MTLCommandBuffer> command = [queue_ commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:corr_xx_];
    [encoder setBuffer:x_buffer_ offset:0 atIndex:0];
    [encoder setBuffer:r_buffer_ offset:0 atIndex:1];
    [encoder setBytes:&rpar length:sizeof(rpar) atIndex:2];
    dispatch_2d(encoder, corr_xx_, n, n);

    [encoder setComputePipelineState:corr_xqmc_];
    [encoder setBuffer:x_buffer_ offset:0 atIndex:0];
    [encoder setBuffer:qmc_buffer_ offset:0 atIndex:1];
    [encoder setBuffer:rhs_buffer_ offset:0 atIndex:2];
    [encoder setBytes:&qpar length:sizeof(qpar) atIndex:3];
    dispatch_2d(encoder, corr_xqmc_, m, n);

    [encoder setComputePipelineState:scale_columns_];
    [encoder setBuffer:rhs_buffer_ offset:0 atIndex:0];
    [encoder setBuffer:scale_buffer_ offset:0 atIndex:1];
    [encoder setBytes:&qpar length:sizeof(qpar) atIndex:2];
    dispatch_1d(encoder, scale_columns_, m);
    [encoder endEncoding];

    [chol_ encodeToCommandBuffer:command sourceMatrix:r_matrix_
                    resultMatrix:factor_matrix_ status:status_buffer_];
    [solve_ encodeToCommandBuffer:command sourceMatrix:factor_matrix_
             rightHandSideMatrix:rhs_matrix_ solutionMatrix:solution_matrix_];

    encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:mmse_];
    [encoder setBuffer:rhs_buffer_ offset:0 atIndex:0];
    [encoder setBuffer:solution_buffer_ offset:0 atIndex:1];
    [encoder setBuffer:scale_buffer_ offset:0 atIndex:2];
    [encoder setBuffer:values_buffer_ offset:0 atIndex:3];
    [encoder setBytes:&qpar length:sizeof(qpar) atIndex:4];
    dispatch_1d(encoder, mmse_, m);
    [encoder setComputePipelineState:reduce_];
    [encoder setBuffer:values_buffer_ offset:0 atIndex:0];
    [encoder setBuffer:minimum_buffer_ offset:0 atIndex:1];
    [encoder setBytes:&qpar length:sizeof(qpar) atIndex:2];
    dispatch_1d(encoder, reduce_, 1);
    [encoder endEncoding];
    [command commit];
    [command waitUntilCompleted];

    if ([command status] == MTLCommandBufferStatusError) {
      const char* message = [[[command error] localizedDescription] UTF8String];
      Rcpp::stop(std::string("Metal command failed: ") + (message ? message : "unknown error"));
    }
    const int32_t status = *static_cast<int32_t*>([status_buffer_ contents]);
    if (status != MPSMatrixDecompositionStatusSuccess) {
      Rcpp::stop("Metal Cholesky failed with status " + std::to_string(status));
    }

    const float minimum = *static_cast<const float*>([minimum_buffer_ contents]);
    return ScoreResult{static_cast<double>(minimum)};
  }

 private:
  id<MTLDevice> device_;
  id<MTLCommandQueue> queue_;
  id<MTLLibrary> library_;
  id<MTLComputePipelineState> corr_xx_;
  id<MTLComputePipelineState> corr_xqmc_;
  id<MTLComputePipelineState> scale_columns_;
  id<MTLComputePipelineState> mmse_;
  id<MTLComputePipelineState> reduce_;
  id<MTLBuffer> qmc_buffer_;
  id<MTLBuffer> x_buffer_, r_buffer_, factor_buffer_, rhs_buffer_;
  id<MTLBuffer> solution_buffer_, values_buffer_, scale_buffer_;
  id<MTLBuffer> status_buffer_, minimum_buffer_;
  MPSMatrix* r_matrix_;
  MPSMatrix* factor_matrix_;
  MPSMatrix* rhs_matrix_;
  MPSMatrix* solution_matrix_;
  MPSMatrixDecompositionCholesky* chol_;
  MPSMatrixSolveCholesky* solve_;
  NSUInteger n_, r_row_bytes_, rhs_row_bytes_;
  arma::uword m_;
  arma::uword p_;
  double alpha_;
  double beta_;
  double nugget_;
  std::size_t evaluations_;

  void ensure_resources(NSUInteger n) {
    if (n_ != 0) {
      if (n != n_) Rcpp::stop("Metal evaluator design size changed.");
      return;
    }
    n_ = n;
    r_row_bytes_ =
      [MPSMatrixDescriptor rowBytesForColumns:n dataType:MPSDataTypeFloat32];
    rhs_row_bytes_ =
      [MPSMatrixDescriptor rowBytesForColumns:m_ dataType:MPSDataTypeFloat32];
    const NSUInteger r_bytes = r_row_bytes_ * n;
    const NSUInteger rhs_bytes = rhs_row_bytes_ * n;
    x_buffer_ = [device_ newBufferWithLength:n * p_ * sizeof(float)
                                      options:MTLResourceStorageModeShared];
    r_buffer_ = [device_ newBufferWithLength:r_bytes
                                      options:MTLResourceStorageModeShared];
    factor_buffer_ = [device_ newBufferWithLength:r_bytes
                                           options:MTLResourceStorageModeShared];
    rhs_buffer_ = [device_ newBufferWithLength:rhs_bytes
                                        options:MTLResourceStorageModeShared];
    solution_buffer_ = [device_ newBufferWithLength:rhs_bytes
                                             options:MTLResourceStorageModeShared];
    values_buffer_ = [device_ newBufferWithLength:m_ * sizeof(float)
                                           options:MTLResourceStorageModeShared];
    scale_buffer_ = [device_ newBufferWithLength:m_ * sizeof(float)
                                          options:MTLResourceStorageModeShared];
    status_buffer_ = [device_ newBufferWithLength:sizeof(int32_t)
                                           options:MTLResourceStorageModeShared];
    minimum_buffer_ = [device_ newBufferWithLength:sizeof(float)
                                            options:MTLResourceStorageModeShared];
    if (!x_buffer_ || !r_buffer_ || !factor_buffer_ || !rhs_buffer_ ||
        !solution_buffer_ || !values_buffer_ || !scale_buffer_ ||
        !status_buffer_ || !minimum_buffer_)
      Rcpp::stop("Metal buffer allocation failed.");

    MPSMatrixDescriptor* r_desc =
      [MPSMatrixDescriptor matrixDescriptorWithRows:n columns:n
                                            rowBytes:r_row_bytes_
                                            dataType:MPSDataTypeFloat32];
    MPSMatrixDescriptor* rhs_desc =
      [MPSMatrixDescriptor matrixDescriptorWithRows:n columns:m_
                                            rowBytes:rhs_row_bytes_
                                            dataType:MPSDataTypeFloat32];
    r_matrix_ = [[MPSMatrix alloc] initWithBuffer:r_buffer_ descriptor:r_desc];
    factor_matrix_ = [[MPSMatrix alloc] initWithBuffer:factor_buffer_
                                            descriptor:r_desc];
    rhs_matrix_ = [[MPSMatrix alloc] initWithBuffer:rhs_buffer_
                                          descriptor:rhs_desc];
    solution_matrix_ = [[MPSMatrix alloc] initWithBuffer:solution_buffer_
                                               descriptor:rhs_desc];
    chol_ = [[MPSMatrixDecompositionCholesky alloc]
      initWithDevice:device_ lower:YES order:n];
    solve_ = [[MPSMatrixSolveCholesky alloc]
      initWithDevice:device_ upper:NO order:n numberOfRightHandSides:m_];
    if (!r_matrix_ || !factor_matrix_ || !rhs_matrix_ || !solution_matrix_ ||
        !chol_ || !solve_)
      Rcpp::stop("Metal MPS resource allocation failed.");
  }

  id<MTLComputePipelineState> make_pipeline(const char* name) {
    NSString* ns_name = [NSString stringWithUTF8String:name];
    id<MTLFunction> function = [library_ newFunctionWithName:ns_name];
    if (function == nil) Rcpp::stop(std::string("Missing Metal function: ") + name);
    NSError* error = nil;
    id<MTLComputePipelineState> pipeline =
      [device_ newComputePipelineStateWithFunction:function error:&error];
    [function release];
    if (pipeline == nil) {
      const char* message = error ? [[error localizedDescription] UTF8String] : "unknown error";
      Rcpp::stop(std::string("Could not create Metal pipeline ") + name + ": " + message);
    }
    return pipeline;
  }

  static void dispatch_1d(id<MTLComputeCommandEncoder> encoder,
                          id<MTLComputePipelineState> pipeline, NSUInteger width) {
    const NSUInteger threads = std::min<NSUInteger>(
      [pipeline maxTotalThreadsPerThreadgroup], 256);
    [encoder dispatchThreads:MTLSizeMake(width, 1, 1)
       threadsPerThreadgroup:MTLSizeMake(threads, 1, 1)];
  }

  static void dispatch_2d(id<MTLComputeCommandEncoder> encoder,
                          id<MTLComputePipelineState> pipeline,
                          NSUInteger width, NSUInteger height) {
    const NSUInteger w = std::min<NSUInteger>(16, [pipeline threadExecutionWidth]);
    const NSUInteger max_threads = [pipeline maxTotalThreadsPerThreadgroup];
    const NSUInteger h = std::max<NSUInteger>(1, std::min<NSUInteger>(16, max_threads / w));
    [encoder dispatchThreads:MTLSizeMake(width, height, 1)
       threadsPerThreadgroup:MTLSizeMake(w, h, 1)];
  }

};

class IncrementalStrictEvaluator {
 public:
  struct Slot {
    id<MTLBuffer> block_buffer = nil;
    id<MTLBuffer> w_buffer = nil;
    id<MTLBuffer> g_buffer = nil;
    id<MTLBuffer> scale = nil;
    id<MTLBuffer> norm = nil;
    id<MTLBuffer> z = nil;
    id<MTLBuffer> logs = nil;
    double min_log = -std::numeric_limits<double>::infinity();
    arma::mat W;
    arma::mat G;
    arma::mat block;
  };

  IncrementalStrictEvaluator(const arma::mat& qmc, double beta, double nugget,
                             arma::uword maximum_n, arma::uword block_size,
                             const std::string& kernel_path)
      : device_(nil), queue_(nil), library_(nil), init_(nil), candidate_(nil),
        accept_(nil), qmc_buffer_(nil), y_(nil), scale_(nil), norm_(nil),
        init_x_(nil), init_l_(nil), init_logs_(nil),
        m_(qmc.n_rows), p_(qmc.n_cols), max_n_(maximum_n), k_(block_size),
        beta_(beta), nugget_(nugget), n_(0), evaluations_(0),
        pending_slot_(-1), pending_n_(0) {
    device_ = MTLCreateSystemDefaultDevice();
    if (!device_) Rcpp::stop("Metal GPU is unavailable.");
    queue_ = [device_ newCommandQueue];
    const std::string source = read_text_file(kernel_path);
    NSError* error = nil;
    library_ = [device_ newLibraryWithSource:[NSString stringWithUTF8String:source.c_str()]
                                     options:nil error:&error];
    if (!library_) Rcpp::stop(std::string("Metal compilation failed: ") +
      (error ? [[error localizedDescription] UTF8String] : "unknown error"));
    init_ = make_pipeline("strict_initialize");
    candidate_ = make_pipeline("strict_incremental_candidate");
    accept_ = make_pipeline("strict_accept_candidate");

    const std::vector<float> qdata = row_major_float(qmc);
    qmc_buffer_ = [device_ newBufferWithBytes:qdata.data()
                                       length:qdata.size() * sizeof(float)
                                      options:MTLResourceStorageModeShared];
    y_ = [device_ newBufferWithLength:max_n_ * m_ * sizeof(float)
                               options:MTLResourceStorageModeShared];
    scale_ = [device_ newBufferWithLength:m_ * sizeof(float)
                                   options:MTLResourceStorageModeShared];
    norm_ = [device_ newBufferWithLength:m_ * sizeof(float)
                                  options:MTLResourceStorageModeShared];
    init_x_ = [device_ newBufferWithLength:max_n_ * p_ * sizeof(float)
                                    options:MTLResourceStorageModeShared];
    init_l_ = [device_ newBufferWithLength:max_n_ * max_n_ * sizeof(float)
                                    options:MTLResourceStorageModeShared];
    init_logs_ = [device_ newBufferWithLength:m_ * sizeof(float)
                                       options:MTLResourceStorageModeShared];
    for (auto& slot : slots_) {
      slot.block_buffer = [device_ newBufferWithLength:k_ * p_ * sizeof(float)
                                               options:MTLResourceStorageModeShared];
      slot.w_buffer = [device_ newBufferWithLength:max_n_ * k_ * sizeof(float)
                                           options:MTLResourceStorageModeShared];
      slot.g_buffer = [device_ newBufferWithLength:k_ * k_ * sizeof(float)
                                           options:MTLResourceStorageModeShared];
      slot.scale = [device_ newBufferWithLength:m_ * sizeof(float)
                                         options:MTLResourceStorageModeShared];
      slot.norm = [device_ newBufferWithLength:m_ * sizeof(float)
                                        options:MTLResourceStorageModeShared];
      slot.z = [device_ newBufferWithLength:k_ * m_ * sizeof(float)
                                     options:MTLResourceStorageModeShared];
      slot.logs = [device_ newBufferWithLength:m_ * sizeof(float)
                                        options:MTLResourceStorageModeShared];
    }
  }

  ~IncrementalStrictEvaluator() {
    for (auto& slot : slots_) {
      [slot.logs release]; [slot.z release]; [slot.norm release]; [slot.scale release];
      [slot.g_buffer release]; [slot.w_buffer release]; [slot.block_buffer release];
    }
    [init_logs_ release]; [init_l_ release]; [init_x_ release];
    [norm_ release]; [scale_ release]; [y_ release]; [qmc_buffer_ release];
    [accept_ release]; [candidate_ release]; [init_ release];
    [library_ release]; [queue_ release]; [device_ release];
  }

  void initialize(const arma::mat& X) {
    pending_slot_ = -1;
    X_ = X;
    n_ = X.n_rows;
    build_cpu_factor(X_, L_);
    const std::vector<float> xdata = row_major_float(X_);
    const std::vector<float> ldata = row_major_float(L_);
    std::memcpy([init_x_ contents], xdata.data(), xdata.size() * sizeof(float));
    std::memcpy([init_l_ contents], ldata.data(), ldata.size() * sizeof(float));
    IncrementalParams par{static_cast<uint32_t>(n_), 0,
      static_cast<uint32_t>(m_), static_cast<uint32_t>(p_),
      static_cast<uint32_t>(m_), static_cast<float>(beta_), 1.0f};
    id<MTLCommandBuffer> command = [queue_ commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:init_];
    [encoder setBuffer:init_x_ offset:0 atIndex:0];
    [encoder setBuffer:qmc_buffer_ offset:0 atIndex:1];
    [encoder setBuffer:init_l_ offset:0 atIndex:2];
    [encoder setBuffer:y_ offset:0 atIndex:3];
    [encoder setBuffer:scale_ offset:0 atIndex:4];
    [encoder setBuffer:norm_ offset:0 atIndex:5];
    [encoder setBuffer:init_logs_ offset:0 atIndex:6];
    [encoder setBytes:&par length:sizeof(par) atIndex:7];
    dispatch(encoder, init_, m_);
    [encoder endEncoding]; [command commit]; [command waitUntilCompleted];
    check(command);
    current_log_ = minimum(init_logs_);
    ++evaluations_;
  }

  void evaluate_pair(const arma::mat& block0, const arma::mat& block1,
                     double& score0, double& score1) {
    const arma::mat* blocks[2] = {&block0, &block1};
    for (int slot_index = 0; slot_index < 2; ++slot_index) {
      Slot& slot = slots_[slot_index];
      slot.block = *blocks[slot_index];
      build_block_factor(slot.block, slot.W, slot.G);
      const std::vector<float> bdata = row_major_float(slot.block);
      const std::vector<float> wdata = row_major_float(slot.W);
      const std::vector<float> gdata = row_major_float(slot.G);
      std::memcpy([slot.block_buffer contents], bdata.data(),
                  bdata.size() * sizeof(float));
      std::memcpy([slot.w_buffer contents], wdata.data(),
                  wdata.size() * sizeof(float));
      std::memcpy([slot.g_buffer contents], gdata.data(),
                  gdata.size() * sizeof(float));
    }
    IncrementalParams par{static_cast<uint32_t>(n_), static_cast<uint32_t>(k_),
      static_cast<uint32_t>(m_), static_cast<uint32_t>(p_),
      static_cast<uint32_t>(m_), static_cast<float>(beta_), 1.0f};
    id<MTLCommandBuffer> command = [queue_ commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    if (pending_slot_ >= 0) {
      // Apply the previous selection immediately before the evaluations that
      // depend on it, so one command-buffer wait covers both operations.
      encode_accept(encoder, pending_slot_, pending_n_);
      pending_slot_ = -1;
    }
    [encoder setComputePipelineState:candidate_];
    for (int slot_index = 0; slot_index < 2; ++slot_index) {
      Slot& slot = slots_[slot_index];
      [encoder setBuffer:slot.block_buffer offset:0 atIndex:0];
      [encoder setBuffer:qmc_buffer_ offset:0 atIndex:1];
      [encoder setBuffer:slot.w_buffer offset:0 atIndex:2];
      [encoder setBuffer:slot.g_buffer offset:0 atIndex:3];
      [encoder setBuffer:y_ offset:0 atIndex:4];
      [encoder setBuffer:scale_ offset:0 atIndex:5];
      [encoder setBuffer:norm_ offset:0 atIndex:6];
      [encoder setBuffer:slot.scale offset:0 atIndex:7];
      [encoder setBuffer:slot.norm offset:0 atIndex:8];
      [encoder setBuffer:slot.z offset:0 atIndex:9];
      [encoder setBuffer:slot.logs offset:0 atIndex:10];
      [encoder setBytes:&par length:sizeof(par) atIndex:11];
      dispatch(encoder, candidate_, m_);
    }
    [encoder endEncoding]; [command commit]; [command waitUntilCompleted];
    check(command);
    for (int slot_index = 0; slot_index < 2; ++slot_index) {
      slots_[slot_index].min_log = minimum(slots_[slot_index].logs);
    }
    evaluations_ += 2;
    score0 = slots_[0].min_log;
    score1 = slots_[1].min_log;
  }

  void accept(int slot_index, bool needs_gpu_state) {
    Slot& slot = slots_[slot_index];
    // The CPU factor is needed to build the next pair immediately. The GPU
    // update can wait until that pair is submitted; the final update is unused.
    pending_slot_ = needs_gpu_state ? slot_index : -1;
    pending_n_ = n_;

    arma::mat updated(n_ + k_, n_ + k_, arma::fill::zeros);
    updated.submat(0, 0, n_ - 1, n_ - 1) = L_;
    updated.submat(n_, 0, n_ + k_ - 1, n_ - 1) = slot.W.t();
    updated.submat(n_, n_, n_ + k_ - 1, n_ + k_ - 1) = slot.G;
    L_ = std::move(updated);
    X_ = arma::join_cols(X_, slot.block);
    n_ += k_;
    current_log_ = slot.min_log;
  }

  double current_log() const { return current_log_; }
  std::size_t evaluations() const { return evaluations_; }
  std::string device_name() const { return std::string([[device_ name] UTF8String]); }

 private:
  id<MTLDevice> device_; id<MTLCommandQueue> queue_; id<MTLLibrary> library_;
  id<MTLComputePipelineState> init_, candidate_, accept_;
  id<MTLBuffer> qmc_buffer_, y_, scale_, norm_;
  id<MTLBuffer> init_x_, init_l_, init_logs_;
  Slot slots_[2];
  arma::uword m_, p_, max_n_, k_, n_;
  double beta_, nugget_, current_log_;
  std::size_t evaluations_;
  int pending_slot_;
  arma::uword pending_n_;
  arma::mat X_, L_;

  id<MTLComputePipelineState> make_pipeline(const char* name) {
    id<MTLFunction> fn = [library_ newFunctionWithName:[NSString stringWithUTF8String:name]];
    NSError* error = nil;
    id<MTLComputePipelineState> result =
      [device_ newComputePipelineStateWithFunction:fn error:&error];
    [fn release];
    if (!result) Rcpp::stop(std::string("Pipeline creation failed: ") + name);
    return result;
  }
  static void dispatch(id<MTLComputeCommandEncoder> enc,
                       id<MTLComputePipelineState> pipe, NSUInteger count) {
    const NSUInteger width = std::min<NSUInteger>(256, [pipe maxTotalThreadsPerThreadgroup]);
    [enc dispatchThreads:MTLSizeMake(count, 1, 1)
     threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
  }
  void encode_accept(id<MTLComputeCommandEncoder> encoder, int slot_index,
                     arma::uword accepted_n) {
    Slot& slot = slots_[slot_index];
    IncrementalParams par{static_cast<uint32_t>(accepted_n),
      static_cast<uint32_t>(k_), static_cast<uint32_t>(m_),
      static_cast<uint32_t>(p_), static_cast<uint32_t>(m_),
      static_cast<float>(beta_), 1.0f};
    [encoder setComputePipelineState:accept_];
    [encoder setBuffer:y_ offset:0 atIndex:0];
    [encoder setBuffer:scale_ offset:0 atIndex:1];
    [encoder setBuffer:norm_ offset:0 atIndex:2];
    [encoder setBuffer:slot.scale offset:0 atIndex:3];
    [encoder setBuffer:slot.norm offset:0 atIndex:4];
    [encoder setBuffer:slot.z offset:0 atIndex:5];
    [encoder setBytes:&par length:sizeof(par) atIndex:6];
    dispatch(encoder, accept_, m_);
  }
  static void check(id<MTLCommandBuffer> command) {
    if ([command status] == MTLCommandBufferStatusError)
      Rcpp::stop(std::string("Metal strict command failed: ") +
        [[[[command error] localizedDescription] description] UTF8String]);
  }
  double minimum(id<MTLBuffer> values) const {
    const float* data = static_cast<const float*>([values contents]);
    float answer = std::numeric_limits<float>::infinity();
    for (arma::uword j = 0; j < m_; ++j) answer = std::min(answer, data[j]);
    return answer;
  }
  void build_cpu_factor(const arma::mat& X, arma::mat& L) const {
    arma::mat R(X.n_rows, X.n_rows, arma::fill::eye);
    for (arma::uword i = 0; i < X.n_rows; ++i)
      for (arma::uword j = i + 1; j < X.n_rows; ++j)
        R(i, j) = R(j, i) = correlation_double(X.row(i), X.row(j), beta_);
    R.diag() += nugget_;
    if (!arma::chol(L, R, "lower")) Rcpp::stop("Initial strict Cholesky failed.");
  }
  void build_block_factor(const arma::mat& block, arma::mat& W, arma::mat& G) const {
    arma::mat C(n_, k_);
    arma::mat D(k_, k_, arma::fill::eye);
    for (arma::uword i = 0; i < n_; ++i)
      for (arma::uword j = 0; j < k_; ++j)
        C(i, j) = correlation_double(X_.row(i), block.row(j), beta_);
    for (arma::uword i = 0; i < k_; ++i)
      for (arma::uword j = i + 1; j < k_; ++j)
        D(i, j) = D(j, i) = correlation_double(block.row(i), block.row(j), beta_);
    D.diag() += nugget_;
    W = arma::solve(arma::trimatl(L_), C, arma::solve_opts::fast);
    arma::mat S = D - W.t() * W;
    if (!arma::chol(G, S, "lower")) Rcpp::stop("Incremental strict Cholesky failed.");
  }
};

class SequentialAugmentEvaluator {
 public:
  struct Selection {
    arma::uword index;
    double log_value;
  };

  SequentialAugmentEvaluator(const arma::mat& candidates, double beta,
                             double nugget, arma::uword maximum_n,
                             const std::string& kernel_path)
      : device_(nil), queue_(nil), library_(nil), init_(nil), accept_(nil),
        candidates_buffer_(nil), y_(nil), scale_(nil), norm_(nil), logs_(nil),
        init_x_(nil), init_l_(nil), point_buffer_(nil), w_buffer_(nil),
        candidates_(candidates), eligible_(candidates.n_rows, true),
        m_(candidates.n_rows), p_(candidates.n_cols), max_n_(maximum_n), n_(0),
        beta_(beta), nugget_(nugget) {
    device_ = MTLCreateSystemDefaultDevice();
    if (!device_) Rcpp::stop("Metal GPU is unavailable.");
    queue_ = [device_ newCommandQueue];
    if (!queue_) Rcpp::stop("Could not create Metal command queue.");

    const std::string source = read_text_file(kernel_path);
    NSError* error = nil;
    library_ = [device_ newLibraryWithSource:
      [NSString stringWithUTF8String:source.c_str()] options:nil error:&error];
    if (!library_) Rcpp::stop(std::string("Metal compilation failed: ") +
      (error ? [[error localizedDescription] UTF8String] : "unknown error"));
    init_ = make_pipeline("augment_initialize");
    accept_ = make_pipeline("augment_accept");

    const std::vector<float> candidate_data = row_major_float(candidates_);
    candidates_buffer_ = buffer(candidate_data);
    y_ = [device_ newBufferWithLength:max_n_ * m_ * sizeof(float)
                               options:MTLResourceStorageModeShared];
    scale_ = [device_ newBufferWithLength:m_ * sizeof(float)
                                   options:MTLResourceStorageModeShared];
    norm_ = [device_ newBufferWithLength:m_ * sizeof(float)
                                  options:MTLResourceStorageModeShared];
    logs_ = [device_ newBufferWithLength:m_ * sizeof(float)
                                  options:MTLResourceStorageModeShared];
    init_x_ = [device_ newBufferWithLength:max_n_ * p_ * sizeof(float)
                                    options:MTLResourceStorageModeShared];
    init_l_ = [device_ newBufferWithLength:max_n_ * max_n_ * sizeof(float)
                                    options:MTLResourceStorageModeShared];
    point_buffer_ = [device_ newBufferWithLength:p_ * sizeof(float)
                                          options:MTLResourceStorageModeShared];
    w_buffer_ = [device_ newBufferWithLength:max_n_ * sizeof(float)
                                      options:MTLResourceStorageModeShared];
    if (!candidates_buffer_ || !y_ || !scale_ || !norm_ || !logs_ ||
        !init_x_ || !init_l_ || !point_buffer_ || !w_buffer_)
      Rcpp::stop("Metal augmentation buffer allocation failed.");
  }

  ~SequentialAugmentEvaluator() {
    [w_buffer_ release]; [point_buffer_ release];
    [init_l_ release]; [init_x_ release];
    [logs_ release]; [norm_ release]; [scale_ release]; [y_ release];
    [candidates_buffer_ release]; [accept_ release]; [init_ release];
    [library_ release]; [queue_ release]; [device_ release];
  }

  void initialize(const arma::mat& X) {
    X_ = X;
    n_ = X_.n_rows;
    build_cpu_factor();
    const std::vector<float> x_data = row_major_float(X_);
    const std::vector<float> l_data = row_major_float(L_);
    std::memcpy([init_x_ contents], x_data.data(),
                x_data.size() * sizeof(float));
    std::memcpy([init_l_ contents], l_data.data(),
                l_data.size() * sizeof(float));
    AugmentParams par{static_cast<uint32_t>(n_), static_cast<uint32_t>(m_),
      static_cast<uint32_t>(p_), static_cast<uint32_t>(m_),
      static_cast<float>(beta_), 0.0f, 1.0f};
    id<MTLCommandBuffer> command = [queue_ commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    [encoder setComputePipelineState:init_];
    [encoder setBuffer:init_x_ offset:0 atIndex:0];
    [encoder setBuffer:candidates_buffer_ offset:0 atIndex:1];
    [encoder setBuffer:init_l_ offset:0 atIndex:2];
    [encoder setBuffer:y_ offset:0 atIndex:3];
    [encoder setBuffer:scale_ offset:0 atIndex:4];
    [encoder setBuffer:norm_ offset:0 atIndex:5];
    [encoder setBuffer:logs_ offset:0 atIndex:6];
    [encoder setBytes:&par length:sizeof(par) atIndex:7];
    dispatch(encoder, init_, m_);
    [encoder endEncoding]; [command commit]; [command waitUntilCompleted];
    check(command);
  }

  Selection select_and_accept() {
    while (true) {
      const float* logs = static_cast<const float*>([logs_ contents]);
      arma::uword selected = m_;
      float best = std::numeric_limits<float>::infinity();
      for (arma::uword j = 0; j < m_; ++j) {
        if (eligible_[j] && std::isfinite(logs[j]) && logs[j] < best) {
          selected = j;
          best = logs[j];
        }
      }
      if (selected == m_)
        Rcpp::stop("No numerically augmentable candidates remain.");

      arma::vec c(n_);
      for (arma::uword row = 0; row < n_; ++row)
        c[row] = correlation_double(X_.row(row), candidates_.row(selected), beta_);
      arma::vec w = arma::solve(arma::trimatl(L_), c, arma::solve_opts::fast);
      const double g_squared = 1.0 + nugget_ - arma::dot(w, w);
      if (!(g_squared > 0.0) || !std::isfinite(g_squared)) {
        eligible_[selected] = false;
        continue;
      }

      const double selected_log = static_cast<double>(best);
      const double g = std::sqrt(g_squared);
      const arma::rowvec point = candidates_.row(selected);
      const std::vector<float> point_data = row_major_float(arma::mat(point));
      const std::vector<float> w_data = row_major_float(arma::mat(w));
      std::memcpy([point_buffer_ contents], point_data.data(),
                  point_data.size() * sizeof(float));
      std::memcpy([w_buffer_ contents], w_data.data(),
                  w_data.size() * sizeof(float));
      AugmentParams par{static_cast<uint32_t>(n_), static_cast<uint32_t>(m_),
        static_cast<uint32_t>(p_), static_cast<uint32_t>(m_),
        static_cast<float>(beta_), static_cast<float>(g), 1.0f};
      id<MTLCommandBuffer> command = [queue_ commandBuffer];
      id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
      [encoder setComputePipelineState:accept_];
      [encoder setBuffer:point_buffer_ offset:0 atIndex:0];
      [encoder setBuffer:candidates_buffer_ offset:0 atIndex:1];
      [encoder setBuffer:w_buffer_ offset:0 atIndex:2];
      [encoder setBuffer:y_ offset:0 atIndex:3];
      [encoder setBuffer:scale_ offset:0 atIndex:4];
      [encoder setBuffer:norm_ offset:0 atIndex:5];
      [encoder setBuffer:logs_ offset:0 atIndex:6];
      [encoder setBytes:&par length:sizeof(par) atIndex:7];
      dispatch(encoder, accept_, m_);
      [encoder endEncoding]; [command commit]; [command waitUntilCompleted];
      check(command);

      arma::mat updated(n_ + 1, n_ + 1, arma::fill::zeros);
      updated.submat(0, 0, n_ - 1, n_ - 1) = L_;
      updated.submat(n_, 0, n_, n_ - 1) = w.t();
      updated(n_, n_) = g;
      L_ = std::move(updated);
      X_ = arma::join_cols(X_, point);
      ++n_;
      eligible_[selected] = false;
      return Selection{selected, selected_log};
    }
  }

  const arma::mat& design() const { return X_; }
  std::string device_name() const { return std::string([[device_ name] UTF8String]); }

 private:
  id<MTLDevice> device_; id<MTLCommandQueue> queue_; id<MTLLibrary> library_;
  id<MTLComputePipelineState> init_, accept_;
  id<MTLBuffer> candidates_buffer_, y_, scale_, norm_, logs_;
  id<MTLBuffer> init_x_, init_l_, point_buffer_, w_buffer_;
  arma::mat candidates_;
  std::vector<bool> eligible_;
  arma::uword m_, p_, max_n_, n_;
  double beta_, nugget_;
  arma::mat X_, L_;

  id<MTLComputePipelineState> make_pipeline(const char* name) {
    id<MTLFunction> function = [library_ newFunctionWithName:
      [NSString stringWithUTF8String:name]];
    NSError* error = nil;
    id<MTLComputePipelineState> result =
      [device_ newComputePipelineStateWithFunction:function error:&error];
    [function release];
    if (!result) Rcpp::stop(std::string("Pipeline creation failed: ") + name);
    return result;
  }

  id<MTLBuffer> buffer(const std::vector<float>& data) {
    return [device_ newBufferWithBytes:data.data()
                                length:data.size() * sizeof(float)
                               options:MTLResourceStorageModeShared];
  }

  static void dispatch(id<MTLComputeCommandEncoder> encoder,
                       id<MTLComputePipelineState> pipeline, NSUInteger count) {
    const NSUInteger width = std::min<NSUInteger>(
      256, [pipeline maxTotalThreadsPerThreadgroup]);
    [encoder dispatchThreads:MTLSizeMake(count, 1, 1)
       threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
  }

  static void check(id<MTLCommandBuffer> command) {
    if ([command status] == MTLCommandBufferStatusError)
      Rcpp::stop(std::string("Metal augmentation command failed: ") +
        [[[[command error] localizedDescription] description] UTF8String]);
  }

  void build_cpu_factor() {
    arma::mat R(n_, n_, arma::fill::eye);
    for (arma::uword i = 0; i < n_; ++i)
      for (arma::uword j = i + 1; j < n_; ++j)
        R(i, j) = R(j, i) = correlation_double(X_.row(i), X_.row(j), beta_);
    R.diag() += nugget_;
    if (!arma::chol(L_, R, "lower"))
      Rcpp::stop("Initial augmentation Cholesky failed.");
  }
};

arma::mat standard_design(const arma::mat& A, arma::uword p, double delta) {
  arma::mat design(A.n_rows * (p + 1), A.n_cols);
  design.rows(0, A.n_rows - 1) = A;
  for (arma::uword col = 0; col < p; ++col) {
    arma::mat block = A;
    for (arma::uword row = 0; row < block.n_rows; ++row)
      block(row, col) += block(row, col) < 0.5 ? delta : -delta;
    design.rows((col + 1) * A.n_rows, (col + 2) * A.n_rows - 1) = block;
  }
  return design;
}

arma::mat shifted_block(const arma::mat& block, arma::uword col, double delta) {
  arma::mat out = block;
  for (arma::uword row = 0; row < out.n_rows; ++row)
    out(row, col) += out(row, col) < 0.5 ? delta : -delta;
  return out;
}

arma::mat bind_blocks(const std::vector<arma::mat>& blocks) {
  arma::uword total_rows = 0;
  for (const arma::mat& block : blocks) total_rows += block.n_rows;
  arma::mat out(total_rows, blocks.front().n_cols);
  arma::uword first_row = 0;
  for (const arma::mat& block : blocks) {
    out.rows(first_row, first_row + block.n_rows - 1) = block;
    first_row += block.n_rows;
  }
  return out;
}

struct StrictResult { arma::mat design; double score; };

StrictResult strict_construct(const arma::mat& A, arma::uword p, double delta,
                              IncrementalStrictEvaluator& evaluator) {
  std::vector<arma::mat> blocks;
  blocks.push_back(A);
  blocks.push_back(shifted_block(A, 0, delta));
  arma::mat old_X = bind_blocks(blocks);
  evaluator.initialize(old_X);
  double current_score = evaluator.current_log();

  for (arma::uword col = 1; col < p; ++col) {
    const arma::mat candidate1 = shifted_block(blocks[0], col, delta);
    const arma::mat candidate2 = shifted_block(blocks[col], col, delta);
    double score1, score2;
    evaluator.evaluate_pair(candidate1, candidate2, score1, score2);
    if (score1 > score2 || gpu_close(score1, score2)) {
      blocks.insert(blocks.begin(), candidate1);
      evaluator.accept(0, col + 1 < p);
      current_score = score1;
    } else {
      blocks.push_back(candidate2);
      evaluator.accept(1, col + 1 < p);
      current_score = score2;
    }
  }
  return StrictResult{bind_blocks(blocks), current_score};
}

} // namespace

Rcpp::List metal_info_rcpp() {
  @autoreleasepool {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (device == nil)
      return Rcpp::List::create(Rcpp::Named("available") = false);
    Rcpp::List result = Rcpp::List::create(
      Rcpp::Named("available") = true,
      Rcpp::Named("device") = std::string([[device name] UTF8String]),
      Rcpp::Named("unified_memory") = static_cast<bool>([device hasUnifiedMemory]),
      Rcpp::Named("recommended_max_working_set_bytes") =
        static_cast<double>([device recommendedMaxWorkingSetSize]));
    [device release];
    return result;
  }
}

Rcpp::List mmse_metal_log_rcpp(const arma::mat& qmc, const arma::mat& X,
                               double alpha, double beta, double nugget,
                               const std::string& kernel_path) {
  @autoreleasepool {
    MetalEvaluator evaluator(qmc, alpha, beta, nugget, kernel_path);
    const ScoreResult score = evaluator.score(X);
    return Rcpp::List::create(
      Rcpp::Named("log_objective") = score.min_log,
      Rcpp::Named("normalized_objective") = score.normalized(X.n_cols),
      Rcpp::Named("original_objective") = std::exp(score.min_log));
  }
}

Rcpp::List soft_aug_metal_rcpp(const arma::mat& design,
                               const arma::mat& candidates,
                               int n_aug, double beta, double nugget,
                               const std::string& kernel_path) {
  @autoreleasepool {
    if (design.n_rows == 0 || design.n_cols == 0)
      Rcpp::stop("design must be a nonempty matrix.");
    if (candidates.n_rows == 0)
      Rcpp::stop("candidates must contain at least one point.");
    if (design.n_cols != candidates.n_cols)
      Rcpp::stop("design and candidates must have the same number of columns.");
    if (!design.is_finite() || !candidates.is_finite())
      Rcpp::stop("design and candidates must contain only finite values.");
    if (n_aug < 1 || static_cast<arma::uword>(n_aug) > candidates.n_rows)
      Rcpp::stop("n_aug must be between 1 and the number of candidates.");
    if (!(beta > 0.0) || !std::isfinite(beta))
      Rcpp::stop("beta must be positive and finite.");
    if (!(nugget >= 0.0) || !std::isfinite(nugget))
      Rcpp::stop("nugget must be finite and nonnegative.");

    SequentialAugmentEvaluator evaluator(
      candidates, beta, nugget, design.n_rows + n_aug, kernel_path);
    evaluator.initialize(design);

    arma::uvec selected_indices(n_aug);
    arma::vec selected_logs(n_aug);
    for (int step = 0; step < n_aug; ++step) {
      const auto selected = evaluator.select_and_accept();
      selected_indices[step] = selected.index;
      selected_logs[step] = selected.log_value;
    }
    arma::mat selected_points(n_aug, candidates.n_cols);
    for (int step = 0; step < n_aug; ++step)
      selected_points.row(step) = candidates.row(selected_indices[step]);
    return Rcpp::List::create(
      Rcpp::Named("design") = evaluator.design(),
      Rcpp::Named("selected_points") = selected_points,
      Rcpp::Named("selected_indices") = selected_indices + 1,
      Rcpp::Named("selected_log_values") = selected_logs,
      Rcpp::Named("selected_values") = arma::exp(selected_logs),
      Rcpp::Named("beta") = beta,
      Rcpp::Named("device") = evaluator.device_name());
  }
}

Rcpp::List stand_opt_metal_rcpp(
    const arma::mat& qmc, arma::mat A, int p, int l, double delta,
    int p0, int p1, double nugget,
    const std::string& kernel_path) {
  @autoreleasepool {
    const double beta = 1.0 / std::pow(2.0 * l, 2.0);
    MetalEvaluator evaluator(qmc, 1.0, beta, nugget, kernel_path);
    arma::mat design = standard_design(A, p, delta);
    double old_score = evaluator.score(design).min_log;

    if (p0 > 1) {
      for (int i = 0; i < l - 1; ++i) {
        for (int j = i + 1; j < l; ++j) {
          for (int k = 1; k < p0; ++k) {
            arma::mat candidate_A = A;
            const int col = p1 + k;
            std::swap(candidate_A(i, col), candidate_A(j, col));
            const arma::mat candidate_design = standard_design(candidate_A, p, delta);
            double new_score = evaluator.score(candidate_design).min_log;
            double comparison_old = old_score;
            if (new_score > comparison_old &&
                !gpu_close(new_score, comparison_old)) {
              A = std::move(candidate_A);
              design = candidate_design;
              old_score = new_score;
            }
          }
        }
      }
    }
    return Rcpp::List::create(
      Rcpp::Named("design") = design,
      Rcpp::Named("log_objective") = old_score,
      Rcpp::Named("normalized_objective") = std::exp(old_score / p),
      Rcpp::Named("original_objective") = std::exp(old_score),
      Rcpp::Named("backend") = "metal",
      Rcpp::Named("device") = evaluator.device_name(),
      Rcpp::Named("gpu_evaluations") = static_cast<double>(evaluator.evaluations()));
  }
}

Rcpp::List strict_opt_metal_rcpp(
    const arma::mat& qmc, arma::mat A, int p, int l, double delta,
    int p0, int p1, double nugget,
    const std::string& kernel_path) {
  @autoreleasepool {
    const double beta = 1.0 / std::pow(2.0 * l, 2.0);
    IncrementalStrictEvaluator evaluator(
      qmc, beta, nugget, static_cast<arma::uword>(l * (p + 1)),
      static_cast<arma::uword>(l), kernel_path);
    StrictResult current = strict_construct(A, p, delta, evaluator);

    if (p0 > 1) {
      for (int i = 0; i < l - 1; ++i) {
        for (int j = i + 1; j < l; ++j) {
          for (int k = 1; k < p0; ++k) {
            arma::mat candidate_A = A;
            const int col = p1 + k;
            std::swap(candidate_A(i, col), candidate_A(j, col));
            StrictResult candidate = strict_construct(candidate_A, p, delta, evaluator);
            double candidate_score = candidate.score;
            double current_score = current.score;
            if (candidate_score > current_score &&
                !gpu_close(candidate_score, current_score)) {
              A = std::move(candidate_A);
              current = std::move(candidate);
              current.score = candidate_score;
            }
          }
        }
      }
    }
    return Rcpp::List::create(
      Rcpp::Named("design") = current.design,
      Rcpp::Named("log_objective") = current.score,
      Rcpp::Named("normalized_objective") = std::exp(current.score / p),
      Rcpp::Named("original_objective") = std::exp(current.score),
      Rcpp::Named("backend") = "metal",
      Rcpp::Named("device") = evaluator.device_name(),
      Rcpp::Named("gpu_evaluations") = static_cast<double>(evaluator.evaluations()));
  }
}

Rcpp::List mofat_metal_info() {
  Rcpp::List result = metal_info_rcpp();
  if (!Rcpp::as<bool>(result["available"])) {
    result["device"] = R_NilValue;
    result["reason"] = "No Metal device is available in this process.";
  }
  return result;
}

Rcpp::List mofat_metal_q(const arma::mat& qmc, const arma::mat& design,
                         double alpha, double beta, double nugget,
                         const std::string& kernel_path) {
  Rcpp::List raw = mmse_metal_log_rcpp(qmc, design, alpha, beta, nugget,
                                       kernel_path);
  return Rcpp::List::create(
    Rcpp::Named("log_Q") = raw["log_objective"],
    Rcpp::Named("Q") = raw["original_objective"],
    Rcpp::Named("normalized_Q") = raw["normalized_objective"]);
}

Rcpp::List mofat_metal_opt(const arma::mat& qmc, const arma::mat& A,
                           int p, int l, double delta, int remainder,
                           int quotient, const std::string& structure,
                           double nugget, const std::string& kernel_path) {
  const int p0 = remainder < 0 ? p : remainder;
  const int p1 = quotient < 0 ? 0 :
    static_cast<int>(std::tgamma(static_cast<double>(l) + 1.0)) * quotient;
  Rcpp::List raw = structure == "strict" ?
    strict_opt_metal_rcpp(qmc, A, p, l, delta, p0, p1, nugget,
                          kernel_path) :
    stand_opt_metal_rcpp(qmc, A, p, l, delta, p0, p1, nugget,
                         kernel_path);
  return Rcpp::List::create(
    Rcpp::Named("design") = raw["design"],
    Rcpp::Named("log_Q") = raw["log_objective"],
    Rcpp::Named("device") = raw["device"]);
}

Rcpp::List mofat_metal_augment(const arma::mat& design,
                               const arma::mat& candidates, int n_aug,
                               double beta, double nugget,
                               const std::string& kernel_path) {
  Rcpp::List raw = soft_aug_metal_rcpp(design, candidates, n_aug, beta,
                                       nugget, kernel_path);
  return Rcpp::List::create(
    Rcpp::Named("design") = raw["design"],
    Rcpp::Named("device") = raw["device"]);
}
