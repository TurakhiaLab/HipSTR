#ifndef HAP_ALIGNER_H_
#define HAP_ALIGNER_H_

#include <assert.h>
#include <string>
#include <vector>

#include "AlignmentData.h"
#include "AlignmentTraceback.h"
#include "../base_quality.h"
#include "Haplotype.h"
#include "StutterAlignerClass.h"

// Strided view into one channel (match/insert/deletion) of an interleaved DP
// matrix buffer laid out as [match0, insert0, deletion0, match1, insert1, ...].
// The alignment recursion's neighbor lookups repeatedly need two of the three
// channels for the same cell (e.g. match+deletion for the "north" and
// "northwest" neighbors), so interleaving keeps those values in the same
// cache line instead of three independently-strided arrays. Mimics enough of
// double*'s interface (indexing, dereference, +/-) to drop into existing
// pointer-style code unchanged.
class MatrixChannel {
 public:
  MatrixChannel() : base_(nullptr) {}
  explicit MatrixChannel(double* base) : base_(base) {}

  double& operator[](long idx) const { return base_[3*idx]; }
  double& operator*() const { return *base_; }
  MatrixChannel operator+(long n) const { return MatrixChannel(base_ + 3*n); }
  MatrixChannel operator-(long n) const { return MatrixChannel(base_ - 3*n); }
  MatrixChannel& operator+=(long n) { base_ += 3*n; return *this; }
  MatrixChannel& operator-=(long n) { base_ -= 3*n; return *this; }

 private:
  double* base_;
};

class HapAligner {
 private:
  Haplotype* fw_haplotype_;
  Haplotype* rev_haplotype_;
  std::vector<bool> realign_to_hap_;
  std::vector<HapBlock*> rev_blocks_;
  std::vector<int32_t> repeat_starts_;
  std::vector<int32_t> repeat_ends_;

  // Per-aligner scratch buffers reused by process_read. HapAligner instances are
  // not shared between worker threads, so these remove hot-path allocations
  // without adding synchronization.
  std::vector<double> base_log_wrong_buf_;
  std::vector<double> base_log_correct_buf_;
  // Interleaved [match, insert, deletion] triples, one per DP matrix cell.
  std::vector<double> l_matrix_buf_;
  std::vector<double> r_matrix_buf_;
  std::vector<int> l_best_artifact_size_buf_;
  std::vector<int> l_best_artifact_pos_buf_;
  std::vector<int> r_best_artifact_size_buf_;
  std::vector<int> r_best_artifact_pos_buf_;

  std::vector<double> block_probs_buf_;
  std::vector<double> artifact_log_prior_buf_;
  std::vector<double> aln_log_probs_buf_;
  StutterWorkspace ws;

  /**
   * Align the sequence contained in SEQ_0 -> SEQ_N using the recursion
   * 0 -> 1 -> 2 ... N
   **/
  void align_seq_to_hap(Haplotype* haplotype, bool reuse_alns,
			const char* seq_0, int seq_len,
			const double* base_log_wrong, const double* base_log_correct,
			MatrixChannel match_matrix, MatrixChannel insert_matrix, MatrixChannel deletion_matrix,
			int* best_artifact_size, int* best_artifact_pos, double& left_prob);

  /**
   * Compute the log-probability of the alignment given the alignment matrices for the left and right segments.
   * Stores the index of the haplotype position with which the seed base is aligned in the maximum likelihood alignment
   **/
  double compute_aln_logprob(int base_seq_len, int seed_base,
			     char seed_char, double log_seed_wrong, double log_seed_correct,
			     MatrixChannel l_match_matrix, MatrixChannel l_insert_matrix, MatrixChannel l_deletion_matrix, double l_prob,
			     MatrixChannel r_match_matrix, MatrixChannel r_insert_matrix, MatrixChannel r_deletion_matrix, double r_prob,
			     int& max_index);

  std::string retrace(Haplotype* haplotype, const char* read_seq, const double* base_log_correct,
		      int seq_len, int block_index, int base_index, int matrix_index, MatrixChannel l_match_matrix,
		      MatrixChannel l_insert_matrix, MatrixChannel l_deletion_matrix, int* best_artifact_size, int* best_artifact_pos,
		      AlignmentTrace& trace);

  void calc_best_seed_position(int32_t region_start, int32_t region_end,
			       int32_t& best_dist, int32_t& best_pos);


  // Private unimplemented copy constructor and assignment operator to prevent operations
  HapAligner(const HapAligner& other);
  HapAligner& operator=(const HapAligner& other);

 

 public:
  HapAligner(Haplotype* haplotype, std::vector<bool>& realign_to_haplotype){
    assert(realign_to_haplotype.size() == haplotype->num_combs());
    fw_haplotype_   = haplotype;
    rev_haplotype_  = haplotype->reverse(rev_blocks_);
    realign_to_hap_ = realign_to_haplotype;


    for (int i = 0; i < fw_haplotype_->num_blocks(); i++){
      HapBlock* block = fw_haplotype_->get_block(i);
      if (block->get_repeat_info() != NULL){
	repeat_starts_.push_back(block->start());
	repeat_ends_.push_back(block->end());
      }
    }
  }

  ~HapAligner(){
    for (unsigned int i = 0; i < rev_blocks_.size(); i++)
      delete rev_blocks_[i];
    rev_blocks_.clear();
    delete rev_haplotype_;
  }

  /**
   * Returns the 0-based index into the sequence string that should be used as the seed for alignment or -1 if no valid seed exists
   **/
  int calc_seed_base(const Alignment& alignment);

  /**
   * Align a half-open subrange of reads into the shared output arrays. Callers
   * assign disjoint ranges, so each worker writes only its own probability and
   * seed-position slots.
   */
  void process_reads_range(const std::vector<Alignment>& alignments,
			   int begin,
			   int end,
			   int init_read_index,
			   const BaseQuality* base_quality,
			   const std::vector<bool>& realign_read,
			   double* aln_probs,
			   int* seed_positions);

  void process_read(const Alignment& aln, int seed_base, const BaseQuality* base_quality, bool retrace_aln,
		    double* prob_ptr, AlignmentTrace& traced_aln);

  void process_reads(const std::vector<Alignment>& alignments, int init_read_index, const BaseQuality* base_quality, const std::vector<bool>& realign_read,
		     double* aln_probs, int* seed_positions);

  /*
    Retraces the Alignment's optimal alignment to the provided haplotype.
    Returns the result as a new Alignment relative to the reference haplotype
   */
  AlignmentTrace* trace_optimal_aln(const Alignment& orig_aln, int seed_base, int best_haplotype, const BaseQuality* base_quality);
};

#endif
