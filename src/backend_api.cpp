// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "soft_gpu_api.h"

// These wrappers are always compiled. Optional platform implementations can
// replace the unavailable results without changing the R-level API.

// [[Rcpp::export]]
Rcpp::List soft_backend_capabilities_rcpp(bool refresh = false) {
  Rcpp::List metal;
#ifdef MOFAT_HAS_METAL
  try {
    metal = mofat_metal_info();
    metal["compiled"] = true;
    if (!metal.containsElementNamed("reason")) {
      metal["reason"] = R_NilValue;
    }
  } catch (std::exception& error) {
    metal = Rcpp::List::create(
      Rcpp::Named("compiled") = true,
      Rcpp::Named("available") = false,
      Rcpp::Named("device") = R_NilValue,
      Rcpp::Named("reason") = error.what());
  }
#else
  metal = Rcpp::List::create(
    Rcpp::Named("compiled") = false,
    Rcpp::Named("available") = false,
    Rcpp::Named("device") = R_NilValue,
    Rcpp::Named("reason") = "Metal backend was not compiled.");
#endif
  Rcpp::List cuda;
  try {
    cuda = mofat_cuda_info(refresh);
  } catch (std::exception& error) {
    cuda = Rcpp::List::create(
      Rcpp::Named("compiled") = false,
      Rcpp::Named("available") = false,
      Rcpp::Named("device") = R_NilValue,
      Rcpp::Named("reason") = error.what());
  }
  return Rcpp::List::create(
    Rcpp::Named("cpu") = Rcpp::List::create(
      Rcpp::Named("compiled") = true,
      Rcpp::Named("available") = true,
      Rcpp::Named("device") = "CPU",
      Rcpp::Named("reason") = R_NilValue),
    Rcpp::Named("metal") = metal,
    Rcpp::Named("cuda") = cuda);
}

// [[Rcpp::export]]
Rcpp::List soft_q_metal_rcpp(const arma::mat& qmc, const arma::mat& design,
                             double alpha, double beta, double nugget,
                             std::string kernel_path) {
  (void) qmc; (void) design; (void) alpha; (void) beta; (void) nugget;
  (void) kernel_path;
#ifdef MOFAT_HAS_METAL
  return mofat_metal_q(qmc, design, alpha, beta, nugget, kernel_path);
#else
  Rcpp::stop("Metal backend is unavailable.");
#endif
}

// [[Rcpp::export]]
Rcpp::List soft_q_cuda_rcpp(const arma::mat& qmc, const arma::mat& design,
                            double alpha, double beta, double nugget) {
  (void) qmc; (void) design; (void) alpha; (void) beta; (void) nugget;
  return mofat_cuda_q(qmc, design, alpha, beta, nugget);
}

// [[Rcpp::export]]
Rcpp::List soft_opt_metal_rcpp(const arma::mat& qmc, const arma::mat& A,
                               int p, int l, double delta, int remainder,
                               int quotient, std::string structure,
                               double nugget, std::string kernel_path) {
  (void) qmc; (void) A; (void) p; (void) l; (void) delta;
  (void) remainder; (void) quotient; (void) structure; (void) nugget;
  (void) kernel_path;
#ifdef MOFAT_HAS_METAL
  return mofat_metal_opt(qmc, A, p, l, delta, remainder, quotient,
                         structure, nugget, kernel_path);
#else
  Rcpp::stop("Metal backend is unavailable.");
#endif
}

// [[Rcpp::export]]
Rcpp::List soft_opt_cuda_rcpp(const arma::mat& qmc, const arma::mat& A,
                              int p, int l, double delta, int remainder,
                              int quotient, std::string structure,
                              double nugget) {
  (void) qmc; (void) A; (void) p; (void) l; (void) delta;
  (void) remainder; (void) quotient; (void) structure; (void) nugget;
  return mofat_cuda_opt(qmc, A, p, l, delta, remainder, quotient,
                        structure, nugget);
}

// [[Rcpp::export]]
Rcpp::List soft_augment_metal_backend_rcpp(const arma::mat& design,
                                            const arma::mat& candidates,
                                            int n_aug, double beta,
                                            double nugget,
                                            std::string kernel_path) {
  (void) design; (void) candidates; (void) n_aug; (void) beta;
  (void) nugget; (void) kernel_path;
#ifdef MOFAT_HAS_METAL
  return mofat_metal_augment(design, candidates, n_aug, beta, nugget,
                             kernel_path);
#else
  Rcpp::stop("Metal backend is unavailable.");
#endif
}

// [[Rcpp::export]]
Rcpp::List soft_augment_cuda_backend_rcpp(const arma::mat& design,
                                           const arma::mat& candidates,
                                           int n_aug, double beta,
                                           double nugget) {
  (void) design; (void) candidates; (void) n_aug; (void) beta;
  (void) nugget;
  return mofat_cuda_augment(design, candidates, n_aug, beta, nugget);
}
