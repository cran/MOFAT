// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// Calculate MIM kernel between two points A[ia,] and B[ib,]: R(A[ia,], B[ib,]).
// beta = theta^2
inline double corr_idx_fast(const arma::mat& A, arma::uword ia,
                            const arma::mat& B, arma::uword ib,
                            const double alpha, const double beta) {
  const arma::uword p = A.n_cols;
  const double invbeta = 1.0 / beta;
  
  if (alpha == 1.0) {
    double c = 1.0;
    for (arma::uword k = 0; k < p; ++k) {
      const double d = A(ia, k) - B(ib, k);
      c *= 1.0 / (1.0 + (d * d) * invbeta);
    }
    return c;
  }
  
  const double negalpha = -alpha;
  double logc = 0.0;
  for (arma::uword k = 0; k < p; ++k) {
    const double d = A(ia, k) - B(ib, k);
    logc += negalpha * std::log1p((d * d) * invbeta);
  }
  return std::exp(logc);
}

// Return log(K) while retaining the alpha = 1 direct-product fast path.
// If the product underflows, recompute in log space so high-dimensional
// designs keep the stabilization used by the public Q calculation.
inline double log_corr_idx_fast(const arma::mat& A, arma::uword ia,
                                const arma::mat& B, arma::uword ib,
                                const double alpha, const double beta) {
  const arma::uword p = A.n_cols;
  const double invbeta = 1.0 / beta;

  if (alpha == 1.0) {
    double product = 1.0;
    for (arma::uword k = 0; k < p; ++k) {
      const double d = A(ia, k) - B(ib, k);
      product *= 1.0 / (1.0 + (d * d) * invbeta);
    }
    if (product > 0.0) return std::log(product);
  }

  const double negalpha = -alpha;
  double logc = 0.0;
  for (arma::uword k = 0; k < p; ++k) {
    const double d = A(ia, k) - B(ib, k);
    logc += negalpha * std::log1p((d * d) * invbeta);
  }
  return logc;
}

// Build a column-scaled correlation RHS. For alpha = 1, the common path
// works with direct products and needs only one logarithm per evaluation
// point. If any product underflows, the whole column is recomputed in log
// space so the stabilized behavior is preserved.
inline void scaled_rhs_fast(const arma::mat& design,
                            const arma::mat& evaluation,
                            const double alpha, const double beta,
                            arma::mat& rhs, arma::rowvec& scales) {
  const arma::uword n = design.n_rows;
  const arma::uword m = evaluation.n_rows;
  rhs.set_size(n, m);
  scales.set_size(m);

  for (arma::uword q = 0; q < m; ++q) {
    bool use_log_space = alpha != 1.0;
    double maximum = 0.0;
    if (!use_log_space) {
      for (arma::uword row = 0; row < n; ++row) {
        const double value = corr_idx_fast(
          design, row, evaluation, q, alpha, beta);
        rhs(row, q) = value;
        maximum = std::max(maximum, value);
        if (value == 0.0) use_log_space = true;
      }
    }

    if (use_log_space) {
      maximum = -std::numeric_limits<double>::infinity();
      for (arma::uword row = 0; row < n; ++row) {
        const double value = log_corr_idx_fast(
          design, row, evaluation, q, alpha, beta);
        rhs(row, q) = value;
        maximum = std::max(maximum, value);
      }
      scales[q] = maximum;
      for (arma::uword row = 0; row < n; ++row) {
        rhs(row, q) = std::exp(rhs(row, q) - maximum);
      }
    } else {
      scales[q] = std::log(maximum);
      for (arma::uword row = 0; row < n; ++row) rhs(row, q) /= maximum;
    }
  }
}

// Stable log-space Q values for a complete design. This is the CPU reference
// used by the public Q function, SOFT optimization, and augmentation.
// [[Rcpp::export]]
Rcpp::List Q_score_cpu_rcpp(const arma::mat& qmc, const arma::mat& X,
                            const double alpha, const double beta,
                            const double nugget = 1e-6)
{
  if (X.n_rows == 0 || X.n_cols == 0)
    Rcpp::stop("design must be a nonempty matrix.");
  if (qmc.n_rows == 0 || qmc.n_cols != X.n_cols)
    Rcpp::stop("qmc must be nonempty and have the same column count as design.");
  if (!X.is_finite() || !qmc.is_finite())
    Rcpp::stop("design and qmc must contain only finite values.");
  if (!(alpha > 0.0) || !std::isfinite(alpha))
    Rcpp::stop("alpha must be positive and finite.");
  if (!(beta > 0.0) || !std::isfinite(beta))
    Rcpp::stop("beta must be positive and finite.");
  if (!(nugget >= 0.0) || !std::isfinite(nugget))
    Rcpp::stop("nugget must be nonnegative and finite.");

  const arma::uword n = X.n_rows;
  const arma::uword m = qmc.n_rows;
  const arma::uword p = X.n_cols;

  arma::mat R(n, n, arma::fill::eye);
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      const double corr = corr_idx_fast(X, i, X, j, alpha, beta);
      R(i, j) = corr;
      R(j, i) = corr;
    }
  }
  R.diag() += nugget;

  arma::mat L;
  if (!arma::chol(L, R, "lower"))
    Rcpp::stop("Cholesky failed for the design correlation matrix.");

  arma::mat rhs;
  arma::rowvec scales;
  scaled_rhs_fast(X, qmc, alpha, beta, rhs, scales);
  arma::mat solved = arma::solve(
    arma::trimatl(L), rhs, arma::solve_opts::fast);
  arma::rowvec norms = arma::sum(arma::square(solved), 0);
  const double smallest = std::numeric_limits<double>::min();
  arma::vec log_values(m);
  for (arma::uword j = 0; j < m; ++j) {
    log_values[j] = 2.0 * scales[j] +
      std::log(std::max(norms[j], smallest));
  }

  const double log_q = log_values.min();
  return Rcpp::List::create(
    Rcpp::Named("log_values") = log_values,
    Rcpp::Named("log_Q") = log_q,
    Rcpp::Named("Q") = std::exp(log_q),
    Rcpp::Named("normalized_Q") = std::exp(log_q / p)
  );
}

namespace {

arma::mat shifted_block_cpu(const arma::mat& block, int factor,
                            double delta) {
  arma::mat result = block;
  for (arma::uword row = 0; row < result.n_rows; ++row) {
    result(row, factor) += result(row, factor) < 0.5 ? delta : -delta;
  }
  return result;
}

arma::mat standard_design_cpu(const arma::mat& A, int p, double delta) {
  arma::mat result(A.n_rows * (p + 1), A.n_cols);
  result.rows(0, A.n_rows - 1) = A;
  for (int factor = 0; factor < p; ++factor) {
    result.rows((factor + 1) * A.n_rows,
                (factor + 2) * A.n_rows - 1) =
      shifted_block_cpu(A, factor, delta);
  }
  return result;
}

class CpuFullEvaluator {
 public:
  CpuFullEvaluator(const arma::mat& evaluation, double alpha, double beta,
                   double nugget)
      : evaluation_(evaluation), alpha_(alpha), beta_(beta), nugget_(nugget) {}

  double score(const arma::mat& design) {
    const arma::uword n = design.n_rows;
    const arma::uword m = evaluation_.n_rows;
    R_.eye(n, n);
    for (arma::uword i = 0; i < n; ++i) {
      for (arma::uword j = i + 1; j < n; ++j) {
        const double value = corr_idx_fast(
          design, i, design, j, alpha_, beta_);
        R_(i, j) = value;
        R_(j, i) = value;
      }
    }
    R_.diag() += nugget_;
    if (!arma::chol(L_, R_, "lower")) {
      Rcpp::stop("CPU Cholesky failed during standard optimization.");
    }

    scaled_rhs_fast(
      design, evaluation_, alpha_, beta_, rhs_, scales_);
    solved_ = arma::solve(
      arma::trimatl(L_), rhs_, arma::solve_opts::fast);
    norms_ = arma::sum(arma::square(solved_), 0);
    const double smallest = std::numeric_limits<double>::min();
    double minimum = std::numeric_limits<double>::infinity();
    for (arma::uword q = 0; q < m; ++q) {
      const double value = 2.0 * scales_[q] +
        std::log(std::max(norms_[q], smallest));
      minimum = std::min(minimum, value);
    }
    return minimum;
  }

 private:
  const arma::mat& evaluation_;
  double alpha_, beta_, nugget_;
  arma::mat R_, L_, rhs_, solved_;
  arma::rowvec scales_, norms_;
};

class CpuIncrementalEvaluator {
 public:
  struct Slot {
    arma::mat block;
    arma::mat W;
    arma::mat G;
    arma::mat Z;
    arma::mat correlations;
    arma::rowvec scales;
    arma::rowvec scale_values;
    arma::rowvec norms;
    arma::rowvec log_values;
    double minimum = -std::numeric_limits<double>::infinity();
    bool valid = false;
  };

  CpuIncrementalEvaluator(const arma::mat& evaluation, double alpha,
                          double beta, double nugget,
                          arma::uword maximum_n, arma::uword block_size)
      : evaluation_(evaluation), alpha_(alpha), beta_(beta), nugget_(nugget),
        max_n_(maximum_n), block_size_(block_size), n_(0) {
    const arma::uword m = evaluation_.n_rows;
    X_.set_size(max_n_, evaluation_.n_cols);
    L_.zeros(max_n_, max_n_);
    Y_.set_size(max_n_, m);
    scales_.set_size(m);
    norms_.set_size(m);
    log_values_.set_size(m);
    old_factors_.set_size(m);
    for (Slot& slot : slots_) {
      slot.block.set_size(block_size_, evaluation_.n_cols);
      slot.W.set_size(max_n_, block_size_);
      slot.G.set_size(block_size_, block_size_);
      slot.Z.set_size(block_size_, m);
      slot.correlations.set_size(m, block_size_);
      slot.scales.set_size(m);
      slot.scale_values.set_size(m);
      slot.norms.set_size(m);
      slot.log_values.set_size(m);
    }
  }

  void initialize(const arma::mat& design) {
    if (design.n_rows == 0 || design.n_rows > max_n_ ||
        design.n_cols != evaluation_.n_cols) {
      Rcpp::stop("Invalid CPU incremental initialization dimensions.");
    }
    n_ = design.n_rows;
    X_.rows(0, n_ - 1) = design;

    active_R_.eye(n_, n_);
    for (arma::uword i = 0; i < n_; ++i) {
      for (arma::uword j = i + 1; j < n_; ++j) {
        const double value = corr_idx_fast(
          X_, i, X_, j, alpha_, beta_);
        active_R_(i, j) = value;
        active_R_(j, i) = value;
      }
    }
    active_R_.diag() += nugget_;
    if (!arma::chol(active_L_, active_R_, "lower")) {
      Rcpp::stop("CPU incremental initialization Cholesky failed.");
    }
    L_.zeros();
    L_.submat(0, 0, n_ - 1, n_ - 1) = active_L_;

    scaled_rhs_fast(
      design, evaluation_, alpha_, beta_, active_rhs_, scales_);
    scale_values_ = arma::exp(scales_);
    active_solved_ = arma::solve(
      arma::trimatl(active_L_), active_rhs_, arma::solve_opts::fast);
    Y_.rows(0, n_ - 1) = active_solved_;
    norms_ = arma::sum(arma::square(active_solved_), 0);
    update_log_values(scales_, norms_, log_values_);
    for (Slot& slot : slots_) slot.valid = false;
  }

  bool evaluate(const arma::mat& block, int slot_index) {
    Slot& slot = slots_[slot_index];
    if (block.n_rows != block_size_ ||
        block.n_cols != evaluation_.n_cols ||
        n_ + block_size_ > max_n_) {
      Rcpp::stop("Invalid CPU incremental candidate dimensions.");
    }
    slot.block = block;

    active_C_.set_size(n_, block_size_);
    active_D_.eye(block_size_, block_size_);
    for (arma::uword i = 0; i < n_; ++i) {
      for (arma::uword j = 0; j < block_size_; ++j) {
        active_C_(i, j) = corr_idx_fast(
          X_, i, slot.block, j, alpha_, beta_);
      }
    }
    for (arma::uword i = 0; i < block_size_; ++i) {
      for (arma::uword j = i + 1; j < block_size_; ++j) {
        const double value = corr_idx_fast(
          slot.block, i, slot.block, j, alpha_, beta_);
        active_D_(i, j) = value;
        active_D_(j, i) = value;
      }
    }
    active_D_.diag() += nugget_;
    active_L_ = L_.submat(0, 0, n_ - 1, n_ - 1);
    active_W_ = arma::solve(
      arma::trimatl(active_L_), active_C_, arma::solve_opts::fast);
    active_S_ = active_D_ - active_W_.t() * active_W_;
    if (!arma::chol(active_G_, active_S_, "lower")) {
      slot.valid = false;
      return false;
    }
    slot.W.rows(0, n_ - 1) = active_W_;
    slot.G = active_G_;

    const arma::uword m = evaluation_.n_rows;
    const bool use_direct_products = alpha_ == 1.0;
    if (use_direct_products) {
      const double invbeta = 1.0 / beta_;
      slot.correlations.ones();
      for (arma::uword row = 0; row < block_size_; ++row) {
        double* products = slot.correlations.colptr(row);
        for (arma::uword dim = 0; dim < evaluation_.n_cols; ++dim) {
          const double block_value = slot.block(row, dim);
          const double* evaluation_values = evaluation_.colptr(dim);
          for (arma::uword q = 0; q < m; ++q) {
            const double difference = block_value - evaluation_values[q];
            products[q] *= 1.0 /
              (1.0 + difference * difference * invbeta);
          }
        }
      }
    }

    slot.minimum = std::numeric_limits<double>::infinity();
    for (arma::uword q = 0; q < m; ++q) {
      bool use_log_space = !use_direct_products || scale_values_[q] == 0.0;
      double block_maximum = 0.0;
      if (!use_log_space) {
        for (arma::uword row = 0; row < block_size_; ++row) {
          const double value = slot.correlations(q, row);
          slot.Z(row, q) = value;
          block_maximum = std::max(block_maximum, value);
          if (value == 0.0) use_log_space = true;
        }
      }

      double scale = scales_[q];
      double old_factor;
      if (use_log_space) {
        for (arma::uword row = 0; row < block_size_; ++row) {
          const double value = log_corr_idx_fast(
            slot.block, row, evaluation_, q, alpha_, beta_);
          slot.Z(row, q) = value;
          scale = std::max(scale, value);
        }
        old_factor = std::exp(scales_[q] - scale);
        slot.scale_values[q] = std::exp(scale);
      } else {
        slot.scale_values[q] = std::max(scale_values_[q], block_maximum);
        old_factor = scale_values_[q] / slot.scale_values[q];
        scale = std::log(slot.scale_values[q]);
      }

      old_factors_[q] = old_factor;
      slot.norms[q] = norms_[q] * old_factor * old_factor;
      for (arma::uword row = 0; row < block_size_; ++row) {
        slot.Z(row, q) = use_log_space ?
          std::exp(slot.Z(row, q) - scale) :
          slot.Z(row, q) / slot.scale_values[q];
      }
      slot.scales[q] = scale;
    }

    if (block_size_ == 1) {
      const double diagonal = active_G_(0, 0);
      for (arma::uword q = 0; q < m; ++q) {
        double adjustment = 0.0;
        for (arma::uword column = 0; column < n_; ++column) {
          adjustment += active_W_(column, 0) * Y_(column, q);
        }
        const double value =
          (slot.Z(0, q) - adjustment * old_factors_[q]) / diagonal;
        slot.Z(0, q) = value;
        slot.norms[q] += value * value;
      }
    } else {
      active_scaled_solutions_ = Y_.rows(0, n_ - 1);
      active_scaled_solutions_.each_row() %= old_factors_;
      slot.Z -= active_W_.t() * active_scaled_solutions_;
      slot.Z = arma::solve(
        arma::trimatl(active_G_), slot.Z, arma::solve_opts::fast);
      slot.norms += arma::sum(arma::square(slot.Z), 0);
    }
    for (arma::uword q = 0; q < m; ++q) {
      slot.log_values[q] = 2.0 * slot.scales[q] + std::log(std::max(
        slot.norms[q], std::numeric_limits<double>::min()));
      slot.minimum = std::min(slot.minimum, slot.log_values[q]);
    }
    slot.valid = true;
    return true;
  }

  void evaluate_pair(const arma::mat& left, const arma::mat& right) {
    if (!evaluate(left, 0) || !evaluate(right, 1)) {
      Rcpp::stop("CPU strict incremental factorization failed.");
    }
  }

  void accept(int slot_index) {
    Slot& slot = slots_[slot_index];
    if (!slot.valid) Rcpp::stop("Cannot accept an invalid CPU candidate.");
    const arma::uword old_n = n_;
    const arma::uword m = evaluation_.n_rows;
    for (arma::uword q = 0; q < m; ++q) {
      const double factor = std::exp(scales_[q] - slot.scales[q]);
      for (arma::uword row = 0; row < old_n; ++row) {
        Y_(row, q) *= factor;
      }
      for (arma::uword row = 0; row < block_size_; ++row) {
        Y_(old_n + row, q) = slot.Z(row, q);
      }
    }
    X_.rows(old_n, old_n + block_size_ - 1) = slot.block;
    L_.submat(old_n, 0, old_n + block_size_ - 1, old_n - 1) =
      slot.W.rows(0, old_n - 1).t();
    L_.submat(old_n, old_n,
              old_n + block_size_ - 1,
              old_n + block_size_ - 1) = slot.G;
    scales_ = slot.scales;
    scale_values_ = slot.scale_values;
    norms_ = slot.norms;
    log_values_ = slot.log_values;
    n_ += block_size_;
  }

  double candidate_min(int slot_index) const {
    return slots_[slot_index].minimum;
  }
  double current_min() const { return log_values_.min(); }
  const arma::rowvec& log_values() const { return log_values_; }
  arma::mat design() const { return X_.rows(0, n_ - 1); }

 private:
  static void update_log_values(const arma::rowvec& scales,
                                const arma::rowvec& norms,
                                arma::rowvec& output) {
    const double smallest = std::numeric_limits<double>::min();
    output.set_size(scales.n_elem);
    for (arma::uword q = 0; q < scales.n_elem; ++q) {
      output[q] = 2.0 * scales[q] +
        std::log(std::max(norms[q], smallest));
    }
  }

  const arma::mat& evaluation_;
  double alpha_, beta_, nugget_;
  arma::uword max_n_, block_size_, n_;
  Slot slots_[2];
  arma::mat X_, L_, Y_;
  arma::rowvec scales_, scale_values_, norms_, log_values_, old_factors_;
  arma::mat active_R_, active_L_, active_rhs_, active_solved_;
  arma::mat active_C_, active_D_, active_W_, active_S_, active_G_;
  arma::mat active_scaled_solutions_;
};

struct CpuDesignResult {
  arma::mat design;
  double log_Q;
};

CpuDesignResult strict_design_cpu(
    const arma::mat& A, int p, double delta,
    CpuIncrementalEvaluator& evaluator) {
  std::vector<arma::mat> blocks;
  blocks.reserve(p + 1);
  blocks.push_back(A);
  blocks.push_back(shifted_block_cpu(A, 0, delta));
  evaluator.initialize(arma::join_cols(blocks[0], blocks[1]));
  for (int factor = 1; factor < p; ++factor) {
    const arma::mat left = shifted_block_cpu(blocks.front(), factor, delta);
    const arma::mat right = shifted_block_cpu(blocks.back(), factor, delta);
    evaluator.evaluate_pair(left, right);
    if (evaluator.candidate_min(0) >= evaluator.candidate_min(1)) {
      blocks.insert(blocks.begin(), left);
      evaluator.accept(0);
    } else {
      blocks.push_back(right);
      evaluator.accept(1);
    }
  }
  arma::mat design(A.n_rows * (p + 1), A.n_cols);
  for (std::size_t index = 0; index < blocks.size(); ++index) {
    design.rows(index * A.n_rows, (index + 1) * A.n_rows - 1) =
      blocks[index];
  }
  return CpuDesignResult{design, evaluator.current_min()};
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List soft_opt_cpu_rcpp(const arma::mat& qmc, const arma::mat& initial_A,
                             int p, int l, double delta, int remainder,
                             int quotient, const std::string& structure,
                             double nugget) {
  const double alpha = 1.0;
  const double beta = 1.0 / std::pow(2.0 * l, 2.0);
  const int p0 = remainder < 0 ? p : remainder;
  const int p1 = quotient < 0 ? 0 :
    static_cast<int>(std::tgamma(static_cast<double>(l) + 1.0)) * quotient;
  arma::mat A = initial_A;
  CpuDesignResult current;

  if (structure == "strict") {
    CpuIncrementalEvaluator evaluator(
      qmc, alpha, beta, nugget, l * (p + 1), l);
    current = strict_design_cpu(A, p, delta, evaluator);
    if (p0 > 1) {
      for (int i = 0; i < l - 1; ++i) {
        for (int j = i + 1; j < l; ++j) {
          for (int k = 1; k < p0; ++k) {
            arma::mat candidate_A = A;
            std::swap(candidate_A(i, p1 + k), candidate_A(j, p1 + k));
            CpuDesignResult candidate = strict_design_cpu(
              candidate_A, p, delta, evaluator);
            if (candidate.log_Q > current.log_Q) {
              A = std::move(candidate_A);
              current = std::move(candidate);
            }
          }
        }
      }
    }
  } else {
    CpuFullEvaluator evaluator(qmc, alpha, beta, nugget);
    current.design = standard_design_cpu(A, p, delta);
    current.log_Q = evaluator.score(current.design);
    if (p0 > 1) {
      for (int i = 0; i < l - 1; ++i) {
        for (int j = i + 1; j < l; ++j) {
          for (int k = 1; k < p0; ++k) {
            arma::mat candidate_A = A;
            std::swap(candidate_A(i, p1 + k), candidate_A(j, p1 + k));
            arma::mat candidate_design = standard_design_cpu(
              candidate_A, p, delta);
            const double candidate_score = evaluator.score(candidate_design);
            if (candidate_score > current.log_Q) {
              A = std::move(candidate_A);
              current.design = std::move(candidate_design);
              current.log_Q = candidate_score;
            }
          }
        }
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("design") = current.design,
    Rcpp::Named("log_Q") = current.log_Q,
    Rcpp::Named("device") = "CPU");
}

// [[Rcpp::export]]
Rcpp::List soft_augment_cpu_rcpp(const arma::mat& design,
                                 const arma::mat& candidates,
                                 int n_aug, double beta, double nugget) {
  CpuIncrementalEvaluator evaluator(
    candidates, 1.0, beta, nugget, design.n_rows + n_aug, 1);
  evaluator.initialize(design);
  std::vector<bool> eligible(candidates.n_rows, true);
  for (int step = 0; step < n_aug; ++step) {
    while (true) {
      arma::uword selected = candidates.n_rows;
      double best = std::numeric_limits<double>::infinity();
      const arma::rowvec& values = evaluator.log_values();
      for (arma::uword index = 0; index < candidates.n_rows; ++index) {
        if (eligible[index] && std::isfinite(values[index]) &&
            values[index] < best) {
          selected = index;
          best = values[index];
        }
      }
      if (selected == candidates.n_rows) {
        Rcpp::stop("No numerically eligible CPU candidate remains.");
      }
      const arma::mat point = candidates.row(selected);
      if (evaluator.evaluate(point, 0)) {
        evaluator.accept(0);
        eligible[selected] = false;
        break;
      }
      eligible[selected] = false;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("design") = evaluator.design(),
    Rcpp::Named("device") = "CPU");
}
