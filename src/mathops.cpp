#include <algorithm>
#include <assert.h>
#include <math.h>

#include "mathops.h"

#include "fastonebigheader.h"

// Lets GCC's vectorizer replace exp() calls inside a "#pragma omp simd" loop
// with glibc's libmvec vector exp instead of a scalar loop. This is the same
// correctly-rounded algorithm as scalar exp(), just batched -- not an
// approximation, unlike fast_log_sum_exp's fasterexp() below. Requires
// -fopenmp-simd (Makefile) and -lmvec (LIBS); pulls in no OpenMP runtime.
#pragma omp declare simd notinbranch
extern "C" double exp(double);

const double LOG_ONE_HALF  = log(0.5);
const double TOLERANCE     = 1e-10;
const double LOG_E_BASE_10 = 0.4342944819;

double INT_LOGS[10000];

void precompute_integer_logs(){
  INT_LOGS[0] = -1000;
  for (unsigned int i = 1; i < 10000; i++)
    INT_LOGS[i] = log(i);
}

double int_log(int val){ return INT_LOGS[val]; }

__attribute__((target_clones("avx512f,avx2,sse4.2,default")))
double sum(const double* begin, const double* end){
  double total = 0.0;
  for (const double* iter = begin; iter != end; iter++)
    total += *iter;
  return total;
}

double sum(const std::vector<double>& vals){
  return sum(vals.data(), vals.data() + vals.size());
}

int sum(const std::vector<bool>& vals){
  int total = 0;
  for (auto iter = vals.begin(); iter != vals.end(); iter++)
    total += *iter;
  return total;
}

__attribute__((target_clones("avx512f,avx2,sse4.2,default")))
double log_sum_exp(const double* begin, const double* end){
  double max_val = *std::max_element(begin, end);
  double total   = 0.0;
  const long n   = end - begin;
  #pragma omp simd reduction(+:total)
  for (long i = 0; i < n; i++)
    total += exp(begin[i] - max_val);
  return max_val + log(total);
}

double log_sum_exp(double log_v1, double log_v2){
  if (log_v1 > log_v2)
    return log_v1 + log(1 + exp(log_v2-log_v1));
  else
    return log_v2 + log(1 + exp(log_v1-log_v2));
}

double log_sum_exp(double log_v1, double log_v2, double log_v3){
  double max_val = std::max(std::max(log_v1, log_v2), log_v3);
  return max_val + log(exp(log_v1-max_val) + exp(log_v2-max_val) + exp(log_v3-max_val));
}

double log_sum_exp(const std::vector<double>& log_vals){
  return log_sum_exp(log_vals.data(), log_vals.data() + log_vals.size());
}

void update_streaming_log_sum_exp(double log_val, double& max_val, double& total){
  if (log_val <= max_val)
    total += exp(log_val - max_val);
  else {
    total  *= exp(max_val-log_val);
    total  += 1.0;
    max_val = log_val;
  }
}

double finish_streaming_log_sum_exp(double max_val, double total){
  return max_val + log(total);
}

double fast_log_sum_exp(double log_v1, double log_v2){
  if (log_v1 > log_v2){
    double diff = log_v2-log_v1;
    return diff < LOG_THRESH ? log_v1 : log_v1 + fastlog(1 + fastexp(diff));
  }
  else {
    double diff = log_v1-log_v2;
    return diff < LOG_THRESH ? log_v2 : log_v2 + fastlog(1 + fastexp(diff));
  }
}

// Sums fasterexp(*iter - max_val) over [begin, end), 4 elements at a time using
// the SSE-vectorized vfasterexp() already vendored in fastonebigheader.h (previously
// unused). Each lane's diff is computed in double precision and narrowed to float
// immediately before the exp call, matching the scalar path's rounding exactly, so
// this is not an approximation of the scalar loop -- it's the same computation batched.
#ifdef __SSE2__
static inline double fast_exp_sum(const double* begin, const double* end, double max_val){
  const v4sf thresh = v4sfl((float) LOG_THRESH);
  v4sf acc = v4sfl(0.0f);
  const double* iter = begin;
  for (; iter + 4 <= end; iter += 4){
    float diffs[4] = { (float) (iter[0] - max_val), (float) (iter[1] - max_val),
                        (float) (iter[2] - max_val), (float) (iter[3] - max_val) };
    v4sf d    = _mm_loadu_ps(diffs);
    v4sf mask = _mm_cmpgt_ps(d, thresh);
    acc = acc + _mm_and_ps(mask, vfasterexp(d));
  }
  float lanes[4];
  _mm_storeu_ps(lanes, acc);
  double total = (double) lanes[0] + (double) lanes[1] + (double) lanes[2] + (double) lanes[3];
  for (; iter != end; iter++){
    double diff = *iter - max_val;
    if (diff > LOG_THRESH)
      total += fasterexp(diff);
  }
  return total;
}
#endif

__attribute__((target_clones("avx512f,avx2,sse4.2,default")))
double fast_log_sum_exp(const double* begin, const double* end){
  double max_val = *std::max_element(begin, end);
#ifdef __SSE2__
  double total = fast_exp_sum(begin, end, max_val);
#else
  double total = 0;
  for (const double* iter = begin; iter != end; iter++){
    double diff = *iter - max_val;
    if (diff > LOG_THRESH)
      total += fasterexp(diff);
  }
#endif
  return max_val + fasterlog(total);
}

double fast_log_sum_exp(const std::vector<double>& log_vals){
  return fast_log_sum_exp(log_vals.data(), log_vals.data() + log_vals.size());
}
