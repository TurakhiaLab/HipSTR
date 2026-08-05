#include <iomanip>
#include <chrono>
#include <iostream>
#include <memory>
#include <sstream>
#include <time.h>
#include <unordered_map>

//#include "sys/sysinfo.h"
//#include "sys/types.h"

#include "extract_indels.h"
#include "genotyper_bam_processor.h"

static double elapsed_seconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(
    std::chrono::steady_clock::now() - start
  ).count();
}

int parseLine(char* line){
  int i = strlen(line);
  while (*line < '0' || *line > '9') line++;
  line[i-3] = '\0';
  i = atoi(line);
  return i;
}

int getUsedPhysicalMemoryKB(){
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmRSS:", 6) == 0){
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  return result;
}

/*
  Left align BamAlignments in the provided vector and store those that successfully realign in the provided vector.
  Also extracts other information for successfully realigned reads into provided vectors.
 */
double GenotyperBamProcessor::left_align_reads(const RegionGroup& region_group, const std::string& chrom_seq, std::vector<BamAlnList>& alignments,
						     const std::vector< std::vector<double> >& log_p1, const std::vector< std::vector<double> >& log_p2,
						     std::vector< std::vector<double> >& filt_log_p1,  std::vector< std::vector<double> >& filt_log_p2,
						     std::vector<Alignment>& left_alns,
						     std::ostream& logger){
  auto left_aln_start = std::chrono::steady_clock::now();
  logger << "Left aligning reads" << std::endl;
  std::unordered_map<std::string, int> seq_to_alns;
  int32_t align_fail_count = 0, total_reads = 0;
  left_alns.clear(); filt_log_p1.clear(); filt_log_p2.clear();

  std::vector<bool> passes_region_filters; passes_region_filters.reserve(region_group.num_regions());
  for (unsigned int i = 0; i < alignments.size(); ++i){
    filt_log_p1.push_back(std::vector<double>());
    filt_log_p2.push_back(std::vector<double>());

    for (unsigned int j = 0; j < alignments[i].size(); ++j, ++total_reads){
      // Trim alignment if it extends very far upstream or downstream of the STR. For tractability, we limit it to 40bp
      alignments[i][j].TrimAlignment((region_group.start() > 40 ? region_group.start()-40 : 1), region_group.stop()+40);
      if (alignments[i][j].Length() == 0)
        continue;

      auto iter      = seq_to_alns.find(alignments[i][j].QueryBases());
      bool have_prev = (iter != seq_to_alns.end());
      if (have_prev)
        have_prev &= left_alns[iter->second].get_sequence().size() == alignments[i][j].QueryBases().size();

      if (!have_prev){
        left_alns.push_back(Alignment(alignments[i][j].Name()));
        if (alignments[i][j].MatchesReference())
          convertAlignment(alignments[i][j], chrom_seq, left_alns.back());
        else if (!realign(alignments[i][j], chrom_seq, left_alns.back())){
	  // Failed to realign read
          align_fail_count++;
          left_alns.pop_back();
          continue;
	}
	seq_to_alns[alignments[i][j].QueryBases()] = left_alns.size()-1;
      }
      else {
        // Reuse alignments if the sequence has already been observed and didn't lead to a soft-clipped alignment
        // Soft-clipping is problematic because it complicates base quality extration (but not really that much)
        Alignment& prev_aln = left_alns[iter->second];
        assert(prev_aln.get_sequence().size() == alignments[i][j].QueryBases().size());
	std::string bases = uppercase(alignments[i][j].QueryBases());
        Alignment new_aln(prev_aln.get_start(), prev_aln.get_stop(), alignments[i][j].IsReverseStrand(), alignments[i][j].Name(), alignments[i][j].Qualities(), bases, prev_aln.get_alignment());
        new_aln.set_cigar_list(prev_aln.get_cigar_list());
        left_alns.push_back(new_aln);
      }

      left_alns.back().check_CIGAR_string(); // Ensure alignment is properly formatted
      filt_log_p1[i].push_back(log_p1[i][j]);
      filt_log_p2[i].push_back(log_p2[i][j]);

      passes_region_filters.clear();
      BamProcessor::passes_filters(alignments[i][j], passes_region_filters);
      left_alns.back().set_hap_gen_info(passes_region_filters);
    }
  }

  double left_aln_time = elapsed_seconds(left_aln_start);
  if (align_fail_count != 0)
    logger << "Failed to left align " << align_fail_count << " out of " << total_reads << " reads" << std::endl;
  return left_aln_time;
}

StutterModel* GenotyperBamProcessor::learn_stutter_model(std::vector<BamAlnList>& alignments,
								 const std::vector< std::vector<double> >& log_p1s,
								 const std::vector< std::vector<double> >& log_p2s,
								 bool haploid, const std::vector<std::string>& rg_names, const Region& region,
								 std::ostream& logger, RegionResult* result){
  std::vector< std::vector<int> > str_bp_lengths(alignments.size());
  std::vector< std::vector<double> > str_log_p1s(alignments.size()), str_log_p2s(alignments.size());
  int inf_reads = 0;
  const int MAX_INF_READS = 10000;

  // Extract bp differences and phasing probabilities for each read if we need to train a stutter model
  for (unsigned int i = 0; i < alignments.size(); ++i){
    for (unsigned int j = 0; j < alignments[i].size(); ++j){
      int bp_diff;
      bool got_size = ExtractCigar(alignments[i][j].CigarData(), alignments[i][j].Position(), region.start()-region.period(), region.stop()+region.period(), bp_diff);
      if (got_size){
	if (bp_diff < -(int)(region.stop()-region.start()+1))
	  continue;
	inf_reads++;
	str_bp_lengths[i].push_back(bp_diff);
	if (log_p1s.size() == 0){
	  str_log_p1s[i].push_back(0); str_log_p2s[i].push_back(0); // Assign equal phasing LLs as no SNP info is available
	}
	else {
	  str_log_p1s[i].push_back(log_p1s[i][j]); str_log_p2s[i].push_back(log_p2s[i][j]);
	}
      }
    }
    if (inf_reads > MAX_INF_READS)
      break;
  }

  if (inf_reads < MIN_TOTAL_READS){
    logger << "Skipping locus with too few informative reads for stutter training: TOTAL=" << inf_reads << ", MIN=" << MIN_TOTAL_READS << std::endl;
    if (result != NULL) result->too_few_reads++;
    else too_few_reads_++;
    return NULL;
  }

  logger << "Building EM stutter model" << std::endl;
  EMStutterGenotyper length_genotyper(haploid, region.period(), str_bp_lengths, str_log_p1s, str_log_p2s, rg_names, 0);
  logger << "Training EM stutter model" << std::endl;
  bool trained = length_genotyper.train(MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE, false, logger);
  if (trained){
    if (output_stutter_models_){
      if (result != NULL){
	std::stringstream stutter_ss;
	length_genotyper.get_stutter_model()->write_model(region.chrom(), region.start(), region.stop(), stutter_ss);
	result->stutter_text += stutter_ss.str();
	result->has_stutter = true;
      }
      else
	length_genotyper.get_stutter_model()->write_model(region.chrom(), region.start(), region.stop(), stutter_model_out_);
    }
    if (result != NULL) result->num_em_converge++;
    else num_em_converge_++;
    StutterModel* stutter_model = length_genotyper.get_stutter_model()->copy();
    logger << "Learned stutter model " << *stutter_model;
    return stutter_model;
  }
  else {
    if (result != NULL) result->num_em_fail++;
    else num_em_fail_++;
    logger << "Stutter model training failed for locus " << region.chrom() << ":" << region.start() << "-" << region.stop()
		  << " with " << inf_reads << " informative reads" << std::endl;
    return NULL;
  }
}
/**
 * Begin of pipeline code
 * analyze_reads_and_phasing changed
 */

void GenotyperBamProcessor::analyze_reads_and_phasing(std::vector<BamAlnList>& alignments,
							      std::vector< std::vector<double> >& log_p1s,
							      std::vector< std::vector<double> >& log_p2s,
							      const std::vector<std::string>& rg_names, const RegionGroup& region_group, const std::string& chrom_seq){
  analyze_reads_and_phasing(alignments, log_p1s, log_p2s, rg_names, region_group, chrom_seq, false, NULL);
}

void GenotyperBamProcessor::analyze_reads_and_phasing(std::vector<BamAlnList>& alignments,
								      std::vector< std::vector<double> >& log_p1s,
								      std::vector< std::vector<double> >& log_p2s,
								      const std::vector<std::string>& rg_names,
								      const RegionGroup& region_group,
								      const std::string& chrom_seq,
								      bool too_many_reads,
								      RegionResult* result){
  std::ostringstream result_log;
  std::ostream& logger = (result != NULL ? result_log : selective_logger());

  int32_t total_reads = 0;
  for (unsigned int i = 0; i < alignments.size(); i++)
    total_reads += alignments[i].size();
  if (total_reads < MIN_TOTAL_READS){
    logger << "Skipping locus with too few reads: TOTAL=" << total_reads << ", MIN=" << MIN_TOTAL_READS << std::endl;
    if (result != NULL) {
      result->too_few_reads++;
      result->log_text += result_log.str();
    }
    else too_few_reads_++;
    return;
  }
  // Can't simply check the total number of reads because the bam processor may have stopped reading at the threshold and then removed PCR duplicates
  // Instead, we check this flag which it sets when too many reads are encountered during filtering
  if (too_many_reads){
    logger << "Skipping locus with too many reads: TOTAL=" << total_reads << ", MAX=" << MAX_TOTAL_READS << std::endl;
    if (result != NULL) {
      result->too_many_reads++;
      result->log_text += result_log.str();
    }
    else too_many_reads_++;
    return;
  }

  assert(alignments.size() == log_p1s.size() && alignments.size() == log_p2s.size() && alignments.size() == rg_names.size());
  bool haploid = (haploid_chroms_.find(region_group.chrom()) != haploid_chroms_.end());
  const std::vector<Region>& regions = region_group.regions();

  // Clip the reads using de Bruijn graph assembly principles
  //assembly_based_read_clipping(alignments, region_group, chrom_seq);

  // Learn the stutter model for each region
  std::vector<StutterModel*> stutter_models;
  auto locus_stutter_start = std::chrono::steady_clock::now();
  bool stutter_success = true;
  for (auto region_iter = regions.begin(); region_iter != regions.end(); region_iter++){
    StutterModel* stutter_model = NULL;
    if (def_stutter_model_ != NULL){
      logger << "Using default stutter model" << std::endl;
      stutter_model = def_stutter_model_->copy();
      stutter_model->set_period(region_iter->period());
    }
    else if (read_stutter_models_){
      // Attempt to extact model from dictionary
      auto model_iter = stutter_models_.find(*region_iter);
      if (model_iter != stutter_models_.end())
	stutter_model = model_iter->second->copy();
      else {
	logger << "WARNING: No stutter model found for " << region_iter->chrom() << ":" << region_iter->start() << "-" << region_iter->stop() << std::endl;
	if (result != NULL) result->num_missing_models++;
	else num_missing_models_++;
      }
    }
    else {
      // Learn stutter model using length-based EM algorithm
      stutter_model = learn_stutter_model(alignments, log_p1s, log_p2s, haploid, rg_names, *region_iter, logger, result);
    }
    stutter_models.push_back(stutter_model);
    stutter_success &= (stutter_model != NULL);
  }
  double locus_stutter_time = elapsed_seconds(locus_stutter_start);
  if (result != NULL) result->stutter_time += locus_stutter_time;
  else total_stutter_time_ += locus_stutter_time;

  // Genotype the regions, if requested
  auto locus_genotype_start = std::chrono::steady_clock::now();
  double locus_left_aln_time = 0;
  std::unique_ptr<VCF::VCFReader> local_ref_vcf;
  SeqStutterGenotyper* seq_genotyper = NULL;
  if (vcf_writer_.is_open() && stutter_success) {
    std::vector<Alignment> left_alignments;
    std::vector< std::vector<double> > filt_log_p1s, filt_log_p2s;
    locus_left_aln_time = left_align_reads(region_group, chrom_seq, alignments, log_p1s, log_p2s, filt_log_p1s,
		     filt_log_p2s, left_alignments, logger);
    if (result == NULL) total_left_aln_time_ += locus_left_aln_time;

    bool run_assembly = (REQUIRE_SPANNING == 0);
    VCF::VCFReader* ref_vcf = ref_vcf_;
    if (result != NULL && !ref_vcf_file_.empty()){
      local_ref_vcf.reset(new VCF::VCFReader(ref_vcf_file_));
      ref_vcf = local_ref_vcf.get();
    }
    seq_genotyper = new SeqStutterGenotyper(region_group, haploid, run_assembly, left_alignments, filt_log_p1s, filt_log_p2s, rg_names, chrom_seq,
						    stutter_models, ref_vcf, logger, READ_THREADS);

    if (seq_genotyper->genotype(MAX_TOTAL_HAPLOTYPES, MAX_FLANK_HAPLOTYPES, MIN_FLANK_FREQ, logger)) {
      bool pass = true;

      // If appropriate, recalculate the stutter model using the haplotype ML alignments,
      // realign the reads and regenotype the samples
      if (recalc_stutter_model_)
	pass = seq_genotyper->recompute_stutter_models(logger, MAX_TOTAL_HAPLOTYPES, MAX_FLANK_HAPLOTYPES, MIN_FLANK_FREQ, MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE);

      if (pass){
	if (result != NULL) result->num_genotype_success++;
	else num_genotype_success_++;
	if (result != NULL){
	  result->chrom = region_group.chrom();
	  result->pos = region_group.start();

	  std::stringstream viz_ss;
	  std::vector<SeqStutterGenotyper::BuiltVCFRecord> records;
	  seq_genotyper->build_vcf_records(samples_to_genotype_, chrom_seq, records, logger,
					   output_viz_, (VIZ_LEFT_ALNS == 1), &viz_ss);
	  for (size_t i = 0; i < records.size(); i++){
	    if (records[i].valid){
	      VCFRecord record;
	      record.chrom = records[i].chrom;
	      record.pos = records[i].pos;
	      record.text = records[i].text;
	      record.valid = records[i].valid;
	      result->vcf_records.push_back(record);
	    }
	  }
	  if (output_viz_){
	    result->viz_text = viz_ss.str();
	    result->has_viz = !result->viz_text.empty();
	  }
	}
	else {
	  seq_genotyper->write_vcf_record(samples_to_genotype_, chrom_seq, output_viz_, (VIZ_LEFT_ALNS == 1), viz_out_, &vcf_writer_, logger);
	}
      }
      else {
	if (result != NULL) result->num_genotype_fail++;
	else num_genotype_fail_++;
      }
    }
    else {
      if (result != NULL) result->num_genotype_fail++;
      else num_genotype_fail_++;
    }
  }
  double locus_genotype_time = elapsed_seconds(locus_genotype_start);
  if (result != NULL) {
    result->left_aln_time += locus_left_aln_time;
    result->genotype_time += locus_genotype_time;
  }
  else {
    locus_left_aln_time_ = locus_left_aln_time;
    locus_genotype_time_ = locus_genotype_time;
    locus_stutter_time_ = locus_stutter_time;
    total_genotype_time_ += locus_genotype_time_;
  }

  double bam_seek_time = (result != NULL ? result->bam_seek_time : locus_bam_seek_time());
  double read_filter_time = (result != NULL ? result->read_filter_time : locus_read_filter_time());
  double snp_phase_info_time = (result != NULL ? result->snp_phase_info_time : locus_snp_phase_info_time());
  logger << "Locus timing:"                                          << "\n"
		     << " BAM seek time       = " << bam_seek_time       << " seconds\n"
		     << " Read filtering      = " << read_filter_time    << " seconds\n"
		     << " SNP info extraction = " << snp_phase_info_time << " seconds\n"
		     << " Stutter estimation  = " << locus_stutter_time  << " seconds\n";
  if (stutter_success && vcf_writer_.is_open()){
    logger << " Genotyping          = " << locus_genotype_time       << " seconds\n";
    if (vcf_writer_.is_open()){
      assert(seq_genotyper != NULL);
      logger << "\t" << " Left alignment        = "  << locus_left_aln_time             << " seconds\n"
				 << "\t" << " Haplotype generation  = "  << seq_genotyper->hap_build_time()  << " seconds\n"
				 << "\t" << " Haplotype alignment   = "  << seq_genotyper->hap_aln_time()    << " seconds\n"
				 << "\t" << " Flank assembly        = "  << seq_genotyper->assembly_time()   << " seconds\n"
				 << "\t" << " Posterior computation = "  << seq_genotyper->posterior_time()  << " seconds\n"
				 << "\t" << " Alignment traceback   = "  << seq_genotyper->aln_trace_time()  << " seconds\n";

      if (result != NULL){
	result->hap_build_time += seq_genotyper->hap_build_time();
	result->hap_aln_time += seq_genotyper->hap_aln_time();
	result->assembly_time += seq_genotyper->assembly_time();
	result->posterior_time += seq_genotyper->posterior_time();
	result->aln_trace_time += seq_genotyper->aln_trace_time();
      }
      else {
	process_timer_.add_time("Left alignment",        locus_left_aln_time);
	process_timer_.add_time("Haplotype generation",  seq_genotyper->hap_build_time());
	process_timer_.add_time("Haplotype alignment",   seq_genotyper->hap_aln_time());
	process_timer_.add_time("Flank assembly",        seq_genotyper->assembly_time());
	process_timer_.add_time("Posterior computation", seq_genotyper->posterior_time());
	process_timer_.add_time("Alignment traceback",   seq_genotyper->aln_trace_time());
      }
    }
  }

  /*
  logger() << "Total memory in use = " << getUsedPhysicalMemoryKB() << " KB"
	   << std::endl;
  */

  logger << "\n";
  if (result != NULL)
    result->log_text += result_log.str();

  delete seq_genotyper;
  for (int i = 0; i < stutter_models.size(); i++)
    delete stutter_models[i];
}

void GenotyperBamProcessor::process_region_item(RegionWorkItem& item, RegionResult& result){
  // Move the fetch/filter output into the result first so ordered writing can
  // still emit BAM records even if genotyping records an early skip.
  result.region_idx = item.region_idx;
  result.log_text = item.log_text;
  result.bam_seek_time = item.bam_seek_time;
  result.read_filter_time = item.read_filter_time;
  result.snp_phase_info_time = item.snp_phase_info_time;
  result.passing_bam_records = std::move(item.passing_bam_records);
  result.filtered_bam_records = std::move(item.filtered_bam_records);
  analyze_reads_and_phasing(item.alignments, item.log_p1s, item.log_p2s,
			    item.rg_names, item.region_group, *item.chrom_seq,
			    item.too_many_reads, &result);
}

void GenotyperBamProcessor::write_region_result(const RegionResult& result) {
  // This runs on the serial pipeline stage, so it is safe to update aggregate
  // counters and write shared output streams here.
  total_bam_seek_time_    += result.bam_seek_time; 
  total_read_filter_time_ += result.read_filter_time;
  too_few_reads_ += result.too_few_reads;
  too_many_reads_ += result.too_many_reads;
  num_missing_models_ += result.num_missing_models;
  num_em_converge_ += result.num_em_converge;
  num_em_fail_ += result.num_em_fail;
  num_genotype_success_ += result.num_genotype_success;
  num_genotype_fail_ += result.num_genotype_fail;
  total_stutter_time_ += result.stutter_time;
  total_left_aln_time_ += result.left_aln_time;
  total_genotype_time_ += result.genotype_time;
  process_timer_.add_time("Left alignment",        result.left_aln_time);
  process_timer_.add_time("Haplotype generation",  result.hap_build_time);
  process_timer_.add_time("Haplotype alignment",   result.hap_aln_time);
  process_timer_.add_time("Flank assembly",        result.assembly_time);
  process_timer_.add_time("Posterior computation", result.posterior_time);
  process_timer_.add_time("Alignment traceback",   result.aln_trace_time);

  for (const auto& aln : result.passing_bam_records) {
    write_passing_alignment(const_cast<BamAlignment&>(aln), pass_writer_);
  }

  for (const auto& rec : result.filtered_bam_records) {
    write_filtered_alignment(const_cast<BamAlignment&>(rec.aln), rec.filter, filt_writer_);
  }

  if (!result.log_text.empty())
    selective_logger() << result.log_text;
  for (size_t i = 0; i < result.vcf_records.size(); ++i) {
    const auto& r = result.vcf_records[i];
    if (r.valid) {
      vcf_writer_.add_vcf_record(r.chrom, r.pos, r.text);
    }
  }
  if (result.has_viz)
    viz_out_ << result.viz_text;
  if (result.has_stutter)
    stutter_model_out_ << result.stutter_text;
}
