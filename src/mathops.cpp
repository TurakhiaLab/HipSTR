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
//
// This declaration is what makes GCC emit calls to libmvec's vector exp
// (_ZGVbN2v_exp and friends), so it has to be gated on libmvec actually being
// available: glibc only gained it in 2.22, and declaring the simd variant
// without the library present fails the link on undefined references rather
// than falling back to scalar. The Makefile probes for it and defines
// HAVE_LIBMVEC (alongside -fopenmp-simd) only when both halves work. Without
// it the loop below stays scalar, which is slower but correct.
#ifdef HAVE_LIBMVEC
#pragma omp declare simd notinbranch
extern "C" double exp(double);
#endif

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

// -ffp-contract=off keeps every target_clones ISA variant numerically
// identical (no FMA fusion pulled in by -mavx512f/-mavx2 that isn't also
// present in the default clone) -- see fast_log_sum_exp below for the
// concrete case this was needed for.
__attribute__((target_clones("avx512f,avx2,sse4.2,default")))
__attribute__((optimize("-ffp-contract=off")))
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
__attribute__((optimize("-ffp-contract=off")))
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

// Sums fasterexp(*iter - max_val) over [begin, end), evaluating the exp four
// elements at a time with the SSE-vectorized vfasterexp() already vendored in
// fastonebigheader.h (previously unused). Each lane's diff is computed in
// double precision and narrowed to float immediately before the exp call, and
// vfasterexp() is bit-identical to scalar fasterexp() lane for lane, so the
// exp evaluation itself is exactly the scalar computation, batched.
//
// The accumulator, however, must stay double and must add in scalar order.
// An earlier version accumulated into a v4sf, i.e. in single precision and
// lane-wise, then widened the four lanes at the end. That was NOT equivalent:
// it differed from the scalar loop in the large majority of realistic inputs,
// by up to ~5e-7 relative -- single-precision epsilon. Small as that is, it is
// enough to flip near-degenerate alignment paths, and it changed the alleles
// discovered at a homopolymer locus (chr2:33759762 on NA12892) relative to
// upstream HipSTR. Adding each lane into a double in the original order keeps
// this bit-identical to the scalar path while retaining the vectorized exp;
// verified over 40k randomized arrays across a range of value spreads.
//
// The >LOG_THRESH decision is made in double precision, from the same ddiffN
// values used for the scalar tail below -- comparing the float-narrowed diff
// instead (as an earlier version of this did) can flip right at the boundary,
// since a diff just above LOG_THRESH in double can round to float and land at
// or below the float-narrowed threshold, silently dropping that term.
#ifdef __SSE2__
__attribute__((optimize("-ffp-contract=off")))
static inline double fast_exp_sum(const double* begin, const double* end, double max_val){
  double total = 0.0;
  const double* iter = begin;
  for (; iter + 4 <= end; iter += 4){
    double ddiff0 = iter[0] - max_val, ddiff1 = iter[1] - max_val;
    double ddiff2 = iter[2] - max_val, ddiff3 = iter[3] - max_val;
    float diffs[4] = { (float) ddiff0, (float) ddiff1, (float) ddiff2, (float) ddiff3 };
    float exps[4];
    _mm_storeu_ps(exps, vfasterexp(_mm_loadu_ps(diffs)));
    if (ddiff0 > LOG_THRESH) total += exps[0];
    if (ddiff1 > LOG_THRESH) total += exps[1];
    if (ddiff2 > LOG_THRESH) total += exps[2];
    if (ddiff3 > LOG_THRESH) total += exps[3];
  }
  for (; iter != end; iter++){
    double diff = *iter - max_val;
    if (diff > LOG_THRESH)
      total += fasterexp(diff);
  }
  return total;
}
#endif

// -ffp-contract=off: without it, the avx512f/avx2 clones let the compiler
// fuse multiply-adds inside this function (including in inlined callees
// like fasterexp()/vfasterexp()) that the default clone doesn't, so the
// SAME source can round differently depending on which ISA clone the CPU
// dispatches to at runtime -- confirmed by bisection to flip a stutter-block
// candidate's log-probability at a couple of homopolymer STR loci out of
// 1.5M genome-wide, changing which alleles get discovered as candidates.
__attribute__((target_clones("avx512f,avx2,sse4.2,default")))
__attribute__((optimize("-ffp-contract=off")))
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
