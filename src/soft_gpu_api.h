#ifndef MOFAT_SOFT_GPU_API_H
#define MOFAT_SOFT_GPU_API_H

#include <RcppArmadillo.h>
#include <string>

#ifdef MOFAT_HAS_METAL
Rcpp::List mofat_metal_info();
Rcpp::List mofat_metal_q(const arma::mat& qmc, const arma::mat& design,
                         double alpha, double beta, double nugget,
                         const std::string& kernel_path);
Rcpp::List mofat_metal_opt(const arma::mat& qmc, const arma::mat& A,
                           int p, int l, double delta, int remainder,
                           int quotient, const std::string& structure,
                           double nugget, const std::string& kernel_path);
Rcpp::List mofat_metal_augment(const arma::mat& design,
                               const arma::mat& candidates, int n_aug,
                               double beta, double nugget,
                               const std::string& kernel_path);
#endif

Rcpp::List mofat_cuda_info(bool refresh = false);
Rcpp::List mofat_cuda_q(const arma::mat& qmc, const arma::mat& design,
                        double alpha, double beta, double nugget);
Rcpp::List mofat_cuda_opt(const arma::mat& qmc, const arma::mat& A,
                          int p, int l, double delta, int remainder,
                          int quotient, const std::string& structure,
                          double nugget);
Rcpp::List mofat_cuda_augment(const arma::mat& design,
                              const arma::mat& candidates, int n_aug,
                              double beta, double nugget);

#endif
