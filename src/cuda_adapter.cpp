// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cuda_build_config.h"
#include "soft_gpu_api.h"

#include <cmath>
#include <string>

#if defined(__linux__)
#include "cuda_helper_api.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <dlfcn.h>

namespace {

using AbiVersionFunction = int (*)(void);
using InfoFunction = int (*)(int*, int*, int*, int*, int*, char*, std::size_t,
                             char*, std::size_t);
using SelfTestFunction = int (*)(char*, std::size_t);
using QFunction = int (*)(const double*, int, int, const double*, int,
                          double, double, double, double*, char*, std::size_t);
using OptimizeFunction = int (*)(const double*, int, const double*, int, int,
                                 double, int, int, int, double, double*,
                                 double*, char*, std::size_t);
using AugmentFunction = int (*)(const double*, int, const double*, int, int,
                                int, double, double, double*, char*,
                                std::size_t);

bool file_exists(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  return input.good();
}

std::string join_path(const std::string& left, const std::string& right) {
  if (left.empty()) return right;
  return left + (left.back() == '/' ? "" : "/") + right;
}

std::vector<double> row_major(const arma::mat& matrix) {
  std::vector<double> result(matrix.n_elem);
  for (arma::uword row = 0; row < matrix.n_rows; ++row)
    for (arma::uword col = 0; col < matrix.n_cols; ++col)
      result[row * matrix.n_cols + col] = matrix(row, col);
  return result;
}

arma::mat column_major_matrix(const std::vector<double>& values,
                              int rows, int cols) {
  arma::mat result(rows, cols);
  for (int row = 0; row < rows; ++row)
    for (int col = 0; col < cols; ++col)
      result(row, col) = values[static_cast<std::size_t>(row) * cols + col];
  return result;
}

std::string helper_filename() {
  return "mofat_cuda_helper.so";
}

std::string discover_helper_path() {
  const char* override_path = std::getenv("MOFAT_CUDA_HELPER");
  if (override_path && override_path[0] != '\0') return override_path;

  Rcpp::Function system_file = Rcpp::Environment::base_env()["system.file"];
  Rcpp::CharacterVector roots = system_file(
    "libs", Rcpp::Named("package") = "MOFAT");
  if (roots.size() < 1 || roots[0] == NA_STRING) return helper_filename();
  const std::string root = Rcpp::as<std::string>(roots[0]);
  return join_path(root, helper_filename());
}

class CudaLoader {
 public:
  CudaLoader()
      : attempted_(false), compiled_(false), available_(false), handle_(nullptr),
        abi_version_(NA_INTEGER), toolkit_version_(NA_INTEGER),
        runtime_version_(NA_INTEGER),
        driver_version_(NA_INTEGER), compute_major_(NA_INTEGER),
        compute_minor_(NA_INTEGER), abi_(nullptr), info_(nullptr),
        self_test_(nullptr), q_(nullptr), optimize_(nullptr), augment_(nullptr) {}

  ~CudaLoader() { close(); }

  void refresh() {
    close();
    attempted_ = false;
    compiled_ = false;
    available_ = false;
    abi_version_ = toolkit_version_ = runtime_version_ = driver_version_ =
      NA_INTEGER;
    compute_major_ = compute_minor_ = NA_INTEGER;
    path_.clear();
    device_.clear();
    reason_.clear();
  }

  void initialize() {
    if (attempted_) return;
    attempted_ = true;
    path_ = discover_helper_path();
    compiled_ = file_exists(path_);
    if (!compiled_) {
      const std::string configured_reason = MOFAT_CUDA_BUILD_REASON;
      reason_ = configured_reason.empty() ?
        "CUDA helper was not built. Install MOFAT from source with a supported "
        "CUDA Toolkit and host compiler." : configured_reason;
      return;
    }
    try {
      open();
      resolve_all();
      abi_version_ = abi_();
      if (abi_version_ != MOFAT_CUDA_HELPER_ABI_VERSION)
        throw std::runtime_error(
          "CUDA helper ABI mismatch: package expects " +
          std::to_string(MOFAT_CUDA_HELPER_ABI_VERSION) + " but helper provides " +
          std::to_string(abi_version_) + ".");
      char device[MOFAT_CUDA_DEVICE_CAPACITY] = {0};
      char error[MOFAT_CUDA_ERROR_CAPACITY] = {0};
      int status = info_(&toolkit_version_, &runtime_version_, &driver_version_,
                         &compute_major_, &compute_minor_, device,
                         sizeof(device), error, sizeof(error));
      check_status(status, error, "CUDA device query failed");
      device_ = device;
      status = self_test_(error, sizeof(error));
      check_status(status, error, "CUDA self-test failed");
      available_ = true;
    } catch (const std::exception& error) {
      available_ = false;
      reason_ = error.what();
      close();
    }
  }

  bool compiled() { initialize(); return compiled_; }
  bool available() { initialize(); return available_; }
  const std::string& path() { initialize(); return path_; }
  const std::string& device() { initialize(); return device_; }
  const std::string& reason() { initialize(); return reason_; }
  int abi_version() { initialize(); return abi_version_; }
  int toolkit_version() { initialize(); return toolkit_version_; }
  int runtime_version() { initialize(); return runtime_version_; }
  int driver_version() { initialize(); return driver_version_; }
  int compute_major() { initialize(); return compute_major_; }
  int compute_minor() { initialize(); return compute_minor_; }

  double q(const arma::mat& qmc, const arma::mat& design,
           double alpha, double beta, double nugget) {
    ensure_available();
    const std::vector<double> qmc_values = row_major(qmc);
    const std::vector<double> design_values = row_major(design);
    char error[MOFAT_CUDA_ERROR_CAPACITY] = {0};
    double log_q = 0.0;
    const int status = q_(
      qmc_values.data(), static_cast<int>(qmc.n_rows),
      static_cast<int>(qmc.n_cols), design_values.data(),
      static_cast<int>(design.n_rows), alpha, beta, nugget, &log_q,
      error, sizeof(error));
    check_status(status, error, "CUDA Q evaluation failed");
    return log_q;
  }

  arma::mat optimize(const arma::mat& qmc, const arma::mat& initial_a,
                     int p, int l, double delta, int remainder, int quotient,
                     const std::string& structure, double nugget,
                     double& log_q) {
    ensure_available();
    const std::vector<double> qmc_values = row_major(qmc);
    const std::vector<double> a_values = row_major(initial_a);
    std::vector<double> output(static_cast<std::size_t>(l) * (p + 1) * p);
    char error[MOFAT_CUDA_ERROR_CAPACITY] = {0};
    const int status = optimize_(
      qmc_values.data(), static_cast<int>(qmc.n_rows), a_values.data(),
      p, l, delta, remainder, quotient, structure == "strict" ? 1 : 0,
      nugget, output.data(), &log_q, error, sizeof(error));
    check_status(status, error, "CUDA design optimization failed");
    return column_major_matrix(output, l * (p + 1), p);
  }

  arma::mat augment(const arma::mat& design, const arma::mat& candidates,
                    int n_aug, double beta, double nugget) {
    ensure_available();
    const std::vector<double> design_values = row_major(design);
    const std::vector<double> candidate_values = row_major(candidates);
    std::vector<double> output(
      static_cast<std::size_t>(design.n_rows + n_aug) * design.n_cols);
    char error[MOFAT_CUDA_ERROR_CAPACITY] = {0};
    const int status = augment_(
      design_values.data(), static_cast<int>(design.n_rows),
      candidate_values.data(), static_cast<int>(candidates.n_rows),
      static_cast<int>(design.n_cols), n_aug, beta, nugget, output.data(),
      error, sizeof(error));
    check_status(status, error, "CUDA design augmentation failed");
    return column_major_matrix(output,
                               static_cast<int>(design.n_rows) + n_aug,
                               static_cast<int>(design.n_cols));
  }

 private:
  bool attempted_, compiled_, available_;
  void* handle_;
  int abi_version_, toolkit_version_, runtime_version_, driver_version_;
  int compute_major_, compute_minor_;
  std::string path_, device_, reason_;
  AbiVersionFunction abi_;
  InfoFunction info_;
  SelfTestFunction self_test_;
  QFunction q_;
  OptimizeFunction optimize_;
  AugmentFunction augment_;

  void ensure_available() {
    initialize();
    if (!available_)
      throw std::runtime_error(reason_.empty() ?
        "CUDA backend is unavailable." : reason_);
  }

  static void check_status(int status, const char* error,
                           const char* fallback) {
    if (status == MOFAT_CUDA_SUCCESS) return;
    if (error && error[0] != '\0') throw std::runtime_error(error);
    throw std::runtime_error(std::string(fallback) +
                             " with status " + std::to_string(status) + ".");
  }

  void open() {
    handle_ = dlopen(path_.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (!handle_) {
      const char* detail = dlerror();
      throw std::runtime_error("Cannot load CUDA helper '" + path_ + "': " +
                               (detail ? detail : "unknown loader error"));
    }
  }

  void close() {
    if (!handle_) return;
    dlclose(handle_);
    handle_ = nullptr;
    abi_ = nullptr;
    info_ = nullptr;
    self_test_ = nullptr;
    q_ = nullptr;
    optimize_ = nullptr;
    augment_ = nullptr;
  }

  void* symbol(const char* name) {
    dlerror();
    void* result = dlsym(handle_, name);
    const char* detail = dlerror();
    if (detail || !result)
      throw std::runtime_error("CUDA helper is missing symbol '" +
                               std::string(name) + "'.");
    return result;
  }

  void resolve_all() {
    abi_ = reinterpret_cast<AbiVersionFunction>(
      symbol("mofat_cuda_helper_abi_version"));
    info_ = reinterpret_cast<InfoFunction>(
      symbol("mofat_cuda_helper_info"));
    self_test_ = reinterpret_cast<SelfTestFunction>(
      symbol("mofat_cuda_helper_self_test"));
    q_ = reinterpret_cast<QFunction>(symbol("mofat_cuda_helper_q"));
    optimize_ = reinterpret_cast<OptimizeFunction>(
      symbol("mofat_cuda_helper_optimize"));
    augment_ = reinterpret_cast<AugmentFunction>(
      symbol("mofat_cuda_helper_augment"));
  }
};

CudaLoader& cuda_loader() {
  static CudaLoader loader;
  return loader;
}

Rcpp::RObject nullable_integer(int value) {
  return value == NA_INTEGER ? Rcpp::RObject(R_NilValue) :
    Rcpp::RObject(Rcpp::wrap(value));
}

}  // namespace

Rcpp::List mofat_cuda_info(bool refresh) {
  CudaLoader& loader = cuda_loader();
  if (refresh) loader.refresh();
  const bool compiled = loader.compiled();
  const bool available = loader.available();
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("compiled") = compiled,
    Rcpp::Named("available") = available,
    Rcpp::Named("device") = available ? Rcpp::wrap(loader.device()) : R_NilValue,
    Rcpp::Named("reason") = available ? R_NilValue : Rcpp::wrap(loader.reason()),
    Rcpp::Named("helper_abi") = nullable_integer(loader.abi_version()),
    Rcpp::Named("toolkit_version") = nullable_integer(loader.toolkit_version()),
    Rcpp::Named("runtime_version") = nullable_integer(loader.runtime_version()),
    Rcpp::Named("driver_version") = nullable_integer(loader.driver_version()),
    Rcpp::Named("compute_capability") = available ?
      Rcpp::wrap(std::to_string(loader.compute_major()) + "." +
                 std::to_string(loader.compute_minor())) : R_NilValue,
    Rcpp::Named("helper_path") = loader.path());
  return result;
}

Rcpp::List mofat_cuda_q(const arma::mat& qmc, const arma::mat& design,
                        double alpha, double beta, double nugget) {
  const double log_q = cuda_loader().q(qmc, design, alpha, beta, nugget);
  return Rcpp::List::create(
    Rcpp::Named("log_Q") = log_q,
    Rcpp::Named("Q") = std::exp(log_q),
    Rcpp::Named("normalized_Q") = std::exp(log_q / design.n_cols));
}

Rcpp::List mofat_cuda_opt(const arma::mat& qmc, const arma::mat& initial_a,
                          int p, int l, double delta, int remainder,
                          int quotient, const std::string& structure,
                          double nugget) {
  double log_q = 0.0;
  arma::mat design = cuda_loader().optimize(
    qmc, initial_a, p, l, delta, remainder, quotient, structure, nugget, log_q);
  return Rcpp::List::create(
    Rcpp::Named("design") = design,
    Rcpp::Named("log_Q") = log_q,
    Rcpp::Named("device") = cuda_loader().device());
}

Rcpp::List mofat_cuda_augment(const arma::mat& design,
                              const arma::mat& candidates, int n_aug,
                              double beta, double nugget) {
  arma::mat result = cuda_loader().augment(
    design, candidates, n_aug, beta, nugget);
  return Rcpp::List::create(
    Rcpp::Named("design") = result,
    Rcpp::Named("device") = cuda_loader().device());
}

#else

namespace {

std::string cuda_unavailable_reason() {
  const std::string configured_reason = MOFAT_CUDA_BUILD_REASON;
  return configured_reason.empty() ?
    "CUDA acceleration is supported only on Linux." : configured_reason;
}

[[noreturn]] void stop_cuda_unavailable() {
  Rcpp::stop(cuda_unavailable_reason());
}

}  // namespace

Rcpp::List mofat_cuda_info(bool refresh) {
  (void) refresh;
  return Rcpp::List::create(
    Rcpp::Named("compiled") = false,
    Rcpp::Named("available") = false,
    Rcpp::Named("device") = R_NilValue,
    Rcpp::Named("reason") = cuda_unavailable_reason(),
    Rcpp::Named("helper_abi") = R_NilValue,
    Rcpp::Named("toolkit_version") = R_NilValue,
    Rcpp::Named("runtime_version") = R_NilValue,
    Rcpp::Named("driver_version") = R_NilValue,
    Rcpp::Named("compute_capability") = R_NilValue,
    Rcpp::Named("helper_path") = R_NilValue);
}

Rcpp::List mofat_cuda_q(const arma::mat& qmc, const arma::mat& design,
                        double alpha, double beta, double nugget) {
  (void) qmc;
  (void) design;
  (void) alpha;
  (void) beta;
  (void) nugget;
  stop_cuda_unavailable();
}

Rcpp::List mofat_cuda_opt(const arma::mat& qmc, const arma::mat& initial_a,
                          int p, int l, double delta, int remainder,
                          int quotient, const std::string& structure,
                          double nugget) {
  (void) qmc;
  (void) initial_a;
  (void) p;
  (void) l;
  (void) delta;
  (void) remainder;
  (void) quotient;
  (void) structure;
  (void) nugget;
  stop_cuda_unavailable();
}

Rcpp::List mofat_cuda_augment(const arma::mat& design,
                              const arma::mat& candidates, int n_aug,
                              double beta, double nugget) {
  (void) design;
  (void) candidates;
  (void) n_aug;
  (void) beta;
  (void) nugget;
  stop_cuda_unavailable();
}

#endif
