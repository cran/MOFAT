#include <metal_stdlib>
using namespace metal;

struct CorrParams {
  uint n;
  uint m;
  uint p;
  uint stride;
  float beta;
  float nugget;
  float alpha;
};

kernel void correlation_xx(
    device const float* X [[buffer(0)]],
    device float* R [[buffer(1)]],
    constant CorrParams& par [[buffer(2)]],
    uint2 gid [[thread_position_in_grid]]) {
  const uint col = gid.x;
  const uint row = gid.y;
  if (row >= par.n || col >= par.n) return;

  float log_value = 0.0f;
  for (uint k = 0; k < par.p; ++k) {
    const float d = X[row * par.p + k] - X[col * par.p + k];
    log_value -= par.alpha * log(1.0f + d * d / par.beta);
  }
  float value = exp(log_value);
  if (row == col) value += par.nugget;
  R[row * par.stride + col] = value;
}

kernel void log_correlation_xqmc(
    device const float* X [[buffer(0)]],
    device const float* qmc [[buffer(1)]],
    device float* log_rhs [[buffer(2)]],
    constant CorrParams& par [[buffer(3)]],
    uint2 gid [[thread_position_in_grid]]) {
  const uint q = gid.x;
  const uint row = gid.y;
  if (row >= par.n || q >= par.m) return;

  float value = 0.0f;
  for (uint k = 0; k < par.p; ++k) {
    const float d = X[row * par.p + k] - qmc[q * par.p + k];
    value -= par.alpha * log(1.0f + d * d / par.beta);
  }
  log_rhs[row * par.stride + q] = value;
}

kernel void scale_log_columns(
    device float* log_rhs [[buffer(0)]],
    device float* column_scale [[buffer(1)]],
    constant CorrParams& par [[buffer(2)]],
    uint q [[thread_position_in_grid]]) {
  if (q >= par.m) return;
  float maximum = -INFINITY;
  for (uint row = 0; row < par.n; ++row)
    maximum = max(maximum, log_rhs[row * par.stride + q]);
  column_scale[q] = maximum;
  for (uint row = 0; row < par.n; ++row) {
    const uint idx = row * par.stride + q;
    log_rhs[idx] = exp(log_rhs[idx] - maximum);
  }
}

kernel void log_mmse_values(
    device const float* scaled_rhs [[buffer(0)]],
    device const float* solution [[buffer(1)]],
    device const float* column_scale [[buffer(2)]],
    device float* log_values [[buffer(3)]],
    constant CorrParams& par [[buffer(4)]],
    uint gid [[thread_position_in_grid]]) {
  if (gid >= par.m) return;
  float sum = 0.0f;
  for (uint row = 0; row < par.n; ++row) {
    const uint idx = row * par.stride + gid;
    sum = fma(scaled_rhs[idx], solution[idx], sum);
  }
  log_values[gid] = 2.0f * column_scale[gid] + log(max(sum, FLT_MIN));
}

kernel void reduce_min_log_values(
    device const float* log_values [[buffer(0)]],
    device float* minimum [[buffer(1)]],
    constant CorrParams& par [[buffer(2)]],
    uint gid [[thread_position_in_grid]]) {
  if (gid != 0) return;
  float answer = INFINITY;
  for (uint q = 0; q < par.m; ++q) answer = min(answer, log_values[q]);
  minimum[0] = answer;
}

struct IncrementalParams {
  uint n;
  uint k;
  uint m;
  uint p;
  uint y_stride;
  float beta;
  float alpha;
};

kernel void strict_initialize(
    device const float* X [[buffer(0)]],
    device const float* qmc [[buffer(1)]],
    device const float* L [[buffer(2)]],
    device float* Y [[buffer(3)]],
    device float* scale_out [[buffer(4)]],
    device float* norm_out [[buffer(5)]],
    device float* log_values [[buffer(6)]],
    constant IncrementalParams& par [[buffer(7)]],
    uint q [[thread_position_in_grid]]) {
  if (q >= par.m) return;
  float scale = -INFINITY;
  for (uint row = 0; row < par.n; ++row) {
    float value = 0.0f;
    for (uint dim = 0; dim < par.p; ++dim) {
      const float d = X[row * par.p + dim] - qmc[q * par.p + dim];
      value -= par.alpha * log(1.0f + d * d / par.beta);
    }
    Y[row * par.y_stride + q] = value;
    scale = max(scale, value);
  }
  float norm = 0.0f;
  for (uint row = 0; row < par.n; ++row) {
    float value = exp(Y[row * par.y_stride + q] - scale);
    for (uint col = 0; col < row; ++col)
      value -= L[row * par.n + col] * Y[col * par.y_stride + q];
    value /= L[row * par.n + row];
    Y[row * par.y_stride + q] = value;
    norm = fma(value, value, norm);
  }
  scale_out[q] = scale;
  norm_out[q] = norm;
  log_values[q] = 2.0f * scale + log(max(norm, FLT_MIN));
}

kernel void strict_incremental_candidate(
    device const float* block [[buffer(0)]],
    device const float* qmc [[buffer(1)]],
    device const float* W [[buffer(2)]],
    device const float* G [[buffer(3)]],
    device const float* old_Y [[buffer(4)]],
    device const float* old_scale [[buffer(5)]],
    device const float* old_norm [[buffer(6)]],
    device float* new_scale [[buffer(7)]],
    device float* new_norm [[buffer(8)]],
    device float* new_Z [[buffer(9)]],
    device float* log_values [[buffer(10)]],
    constant IncrementalParams& par [[buffer(11)]],
    uint q [[thread_position_in_grid]]) {
  if (q >= par.m) return;

  float scale = old_scale[q];
  for (uint row = 0; row < par.k; ++row) {
    float value = 0.0f;
    for (uint dim = 0; dim < par.p; ++dim) {
      const float d = block[row * par.p + dim] - qmc[q * par.p + dim];
      value -= par.alpha * log(1.0f + d * d / par.beta);
    }
    new_Z[row * par.y_stride + q] = value;
    scale = max(scale, value);
  }

  const float old_factor = exp(old_scale[q] - scale);
  for (uint row = 0; row < par.k; ++row) {
    float value = exp(new_Z[row * par.y_stride + q] - scale);
    for (uint col = 0; col < par.n; ++col)
      value -= W[col * par.k + row] * old_Y[col * par.y_stride + q] * old_factor;
    for (uint col = 0; col < row; ++col)
      value -= G[row * par.k + col] * new_Z[col * par.y_stride + q];
    value /= G[row * par.k + row];
    new_Z[row * par.y_stride + q] = value;
  }

  float norm = old_norm[q] * old_factor * old_factor;
  for (uint row = 0; row < par.k; ++row) {
    const float z = new_Z[row * par.y_stride + q];
    norm = fma(z, z, norm);
  }
  new_scale[q] = scale;
  new_norm[q] = norm;
  log_values[q] = 2.0f * scale + log(max(norm, FLT_MIN));
}

kernel void strict_accept_candidate(
    device float* Y [[buffer(0)]],
    device float* scale [[buffer(1)]],
    device float* norm [[buffer(2)]],
    device const float* candidate_scale [[buffer(3)]],
    device const float* candidate_norm [[buffer(4)]],
    device const float* candidate_Z [[buffer(5)]],
    constant IncrementalParams& par [[buffer(6)]],
    uint q [[thread_position_in_grid]]) {
  if (q >= par.m) return;
  const float factor = exp(scale[q] - candidate_scale[q]);
  for (uint row = 0; row < par.n; ++row)
    Y[row * par.y_stride + q] *= factor;
  for (uint row = 0; row < par.k; ++row)
    Y[(par.n + row) * par.y_stride + q] = candidate_Z[row * par.y_stride + q];
  scale[q] = candidate_scale[q];
  norm[q] = candidate_norm[q];
}

struct AugmentParams {
  uint n;
  uint m;
  uint p;
  uint y_stride;
  float beta;
  float g;
  float alpha;
};

kernel void augment_initialize(
    device const float* X [[buffer(0)]],
    device const float* candidates [[buffer(1)]],
    device const float* L [[buffer(2)]],
    device float* Y [[buffer(3)]],
    device float* scale_out [[buffer(4)]],
    device float* norm_out [[buffer(5)]],
    device float* log_values [[buffer(6)]],
    constant AugmentParams& par [[buffer(7)]],
    uint q [[thread_position_in_grid]]) {
  if (q >= par.m) return;

  float scale = -INFINITY;
  for (uint row = 0; row < par.n; ++row) {
    float value = 0.0f;
    for (uint dim = 0; dim < par.p; ++dim) {
      const float d = X[row * par.p + dim] - candidates[q * par.p + dim];
      value -= par.alpha * log(1.0f + d * d / par.beta);
    }
    Y[row * par.y_stride + q] = value;
    scale = max(scale, value);
  }

  float norm = 0.0f;
  for (uint row = 0; row < par.n; ++row) {
    float value = exp(Y[row * par.y_stride + q] - scale);
    for (uint col = 0; col < row; ++col)
      value -= L[row * par.n + col] * Y[col * par.y_stride + q];
    value /= L[row * par.n + row];
    Y[row * par.y_stride + q] = value;
    norm = fma(value, value, norm);
  }
  scale_out[q] = scale;
  norm_out[q] = norm;
  log_values[q] = 2.0f * scale + log(max(norm, FLT_MIN));
}

kernel void augment_accept(
    device const float* point [[buffer(0)]],
    device const float* candidates [[buffer(1)]],
    device const float* w [[buffer(2)]],
    device float* Y [[buffer(3)]],
    device float* scale [[buffer(4)]],
    device float* norm [[buffer(5)]],
    device float* log_values [[buffer(6)]],
    constant AugmentParams& par [[buffer(7)]],
    uint q [[thread_position_in_grid]]) {
  if (q >= par.m) return;

  float log_r = 0.0f;
  for (uint dim = 0; dim < par.p; ++dim) {
    const float d = point[dim] - candidates[q * par.p + dim];
    log_r -= par.alpha * log(1.0f + d * d / par.beta);
  }
  const float new_scale = max(scale[q], log_r);
  const float old_factor = exp(scale[q] - new_scale);
  float rhs = exp(log_r - new_scale);
  for (uint row = 0; row < par.n; ++row)
    rhs -= w[row] * Y[row * par.y_stride + q] * old_factor;
  const float z = rhs / par.g;

  for (uint row = 0; row < par.n; ++row)
    Y[row * par.y_stride + q] *= old_factor;
  Y[par.n * par.y_stride + q] = z;
  const float new_norm = fma(z, z, norm[q] * old_factor * old_factor);
  scale[q] = new_scale;
  norm[q] = new_norm;
  log_values[q] = 2.0f * new_scale + log(max(new_norm, FLT_MIN));
}
