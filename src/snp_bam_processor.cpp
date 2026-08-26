#include <assert.h>
#include <chrono>
#include <time.h>

#include "snp_bam_processor.h"
#include "snp_phasing_quality.h"
#include "snp_tree.h"

static double elapsed_seconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(
    std::chrono::steady_clock::now() - start
  ).count();
}

void SNPBamProcessor::verify_vcf_chromosomes(const std::vector<std::string>& chroms){
  if (phased_snp_vcf_ == NULL)
    return;

  for (auto chrom_iter = chroms.begin(); chrom_iter != chroms.end(); chrom_iter++){
    std::string chrom = (*chrom_iter);
    if (!(phased_snp_vcf_->has_chromosome(chrom))){
      std::stringstream err_msg;
      err_msg << "No entries for chromosome " << chrom << " found in the SNP VCF file" << "\n"
	      << "\t" << "Please ensure that the chromosome names in your BED file match those in your SNP VCF file";
      full_logger() << "\n" << "ERROR: " << err_msg.str() << std::endl;

      std::vector<std::string> alt_names(1, "chr" + chrom);
      if (chrom.size() > 3 && chrom.substr(0, 3).compare("chr") == 0)
	alt_names.push_back(chrom.substr(3));

      // Prompt if simple changes to the chromosome name would solve the issue
      for (auto alt_iter = alt_names.begin(); alt_iter != alt_names.end(); alt_iter++)
	if (phased_snp_vcf_->has_chromosome(*alt_iter))
	  full_logger() << "\t" << "NOTE: Found chromosome " << (*alt_iter) << " in the VCF, but not chromosome " << chrom << std::endl;

      // Abort execution
      printErrorAndDie("Terminating HipSTR-MT as chromosomes in the region file are missing from the SNP VCF file. Please see the log for details");
    }
  }
}

void SNPBamProcessor::process_reads(std::vector<BamAlnList>& paired_strs_by_rg,
				    std::vector<BamAlnList>& mate_pairs_by_rg,
				    std::vector<BamAlnList>& unpaired_strs_by_rg,
				    const std::vector<std::string>& rg_names, const RegionGroup& region_group, const std::string& chrom_seq){
  std::vector<BamAlnList> alignments;
  std::vector< std::vector<double> > log_p1s, log_p2s;
  prepare_read_phasing(paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg,
		       rg_names, region_group, chrom_seq, alignments, log_p1s, log_p2s, selective_logger(), NULL);
  analyze_reads_and_phasing(alignments, log_p1s, log_p2s, rg_names, region_group, chrom_seq);
}

bool SNPBamProcessor::prepare_region_work_item(RegionWorkItem& item, std::ostream& logger){
  bool needs_serial_snp_state = phased_snp_vcf_ != NULL || haplotype_tracker_ != NULL;
  bool prepared = false;
  double phase_time = 0;
  if(needs_serial_snp_state){
    // VCFReader iteration and HaplotypeTracker updates are stateful, so keep
    // this preparation step serialized when those inputs are active.
    std::lock_guard<std::mutex> lock(snp_phase_mutex_);
    prepared = prepare_read_phasing(item.paired_strs_by_rg, item.mate_pairs_by_rg, item.unpaired_strs_by_rg,
                item.rg_names, item.region_group, *item.chrom_seq,
                item.alignments, item.log_p1s, item.log_p2s, logger, &phase_time);
  }
  else {
    prepared = prepare_read_phasing(item.paired_strs_by_rg, item.mate_pairs_by_rg, item.unpaired_strs_by_rg,
                item.rg_names, item.region_group, *item.chrom_seq,
                item.alignments, item.log_p1s, item.log_p2s, logger, &phase_time);
  }
  
  item.snp_phase_info_time = phase_time;
  return prepared;
}

bool SNPBamProcessor::prepare_read_phasing(std::vector<BamAlnList>& paired_strs_by_rg,
					   std::vector<BamAlnList>& mate_pairs_by_rg,
					   std::vector<BamAlnList>& unpaired_strs_by_rg,
					   const std::vector<std::string>& rg_names,
					   const RegionGroup& region_group,
					   const std::string& chrom_seq,
					   std::vector<BamAlnList>& alignments,
					   std::vector< std::vector<double> >& log_p1s,
					   std::vector< std::vector<double> >& log_p2s,
					   std::ostream& logger,
             double* phase_time_out){

  auto phase_start = std::chrono::steady_clock::now();
  int32_t local_match_count = 0;
  int32_t local_mismatch_count = 0;
  if (bams_from_10x_){
    assert(paired_strs_by_rg.size() == mate_pairs_by_rg.size() && paired_strs_by_rg.size() == unpaired_strs_by_rg.size());

    alignments.assign(paired_strs_by_rg.size(), BamAlnList());
    log_p1s.clear(); log_p2s.clear();
    int32_t phased_reads = 0, total_reads = 0;
    for (unsigned int i = 0; i < paired_strs_by_rg.size(); i++){
      alignments[i].insert(alignments[i].end(), paired_strs_by_rg[i].begin(),   paired_strs_by_rg[i].end());
      alignments[i].insert(alignments[i].end(), unpaired_strs_by_rg[i].begin(), unpaired_strs_by_rg[i].end());

      log_p1s.push_back(std::vector<double>());
      log_p2s.push_back(std::vector<double>());
      for (unsigned int j = 0; j < paired_strs_by_rg[i].size(); j++){
	total_reads++;
	int haplotype_1 = get_haplotype(paired_strs_by_rg[i][j]);
	int haplotype_2 = get_haplotype(mate_pairs_by_rg[i][j]);
	int haplotype = (haplotype_1 != haplotype_2 ? -1 : haplotype_1);
	if (haplotype != -1){
	  phased_reads++;
	  log_p1s[i].push_back(haplotype == 1 ? FROM_HAP_LL : OTHER_HAP_LL);
	  log_p2s[i].push_back(haplotype == 2 ? FROM_HAP_LL : OTHER_HAP_LL);
	}
	else {
	  log_p1s[i].push_back(0.0);
	  log_p2s[i].push_back(0.0);
	}
      }
      for (unsigned int j = 0; j < unpaired_strs_by_rg[i].size(); j++){
	total_reads++;
	int haplotype = get_haplotype(unpaired_strs_by_rg[i][j]);
	if (haplotype != -1){
	  phased_reads++;
	  log_p1s[i].push_back(haplotype == 1 ? FROM_HAP_LL : OTHER_HAP_LL);
	  log_p2s[i].push_back(haplotype == 2 ? FROM_HAP_LL : OTHER_HAP_LL);
	}
	else {
	  log_p1s[i].push_back(0.0);
	  log_p2s[i].push_back(0.0);
	}
      }
    }

    logger << "Phased SNPs add info for " << phased_reads << " out of " << total_reads << " reads" << std::endl;
    double elapsed = elapsed_seconds(phase_start);
    locus_snp_phase_info_time_ = elapsed;
    if (phase_time_out != NULL)
      *phase_time_out = elapsed;
    else
      total_snp_phase_info_time_ += elapsed;
    return true;
  }

  assert(paired_strs_by_rg.size() == mate_pairs_by_rg.size() && paired_strs_by_rg.size() == unpaired_strs_by_rg.size());
  
  alignments.assign(paired_strs_by_rg.size(), BamAlnList());
  log_p1s.clear(); log_p2s.clear();
  const std::vector<Region>& skip_regions = region_group.regions();
  bool got_snp_info = false;
  if (phased_snp_vcf_ != NULL){
    // If we are tracking SNP haplotypes for pedigree-based filtering, we need to update the haplotypes to the current position
    if (haplotype_tracker_ != NULL){
      std::set<std::string> sites_to_skip;
      haplotype_tracker_->advance(region_group.chrom(), region_group.start(), sites_to_skip);
    }

    std::vector<SNPTree*> snp_trees;
    std::map<std::string, unsigned int> sample_indices;      
    if (create_snp_trees(region_group.chrom(), (region_group.start() > MAX_MATE_DIST ? region_group.start()-MAX_MATE_DIST : 1), region_group.stop()+MAX_MATE_DIST,
			 skip_regions, SKIP_PADDING, phased_snp_vcf_, haplotype_tracker_, sample_indices, snp_trees, logger)){
      got_snp_info = true;
      std::set<std::string> bad_samples, good_samples;
      for (unsigned int i = 0; i < paired_strs_by_rg.size(); ++i){
	if (sample_indices.find(rg_names[i]) != sample_indices.end()){
	  good_samples.insert(rg_names[i]);
	  std::vector<double> log_p1, log_p2;
	  SNPTree* snp_tree = snp_trees[sample_indices[rg_names[i]]];
	  calc_het_snp_factors(paired_strs_by_rg[i], mate_pairs_by_rg[i], base_quality_, snp_tree, log_p1, log_p2, local_match_count, local_mismatch_count);
	  calc_het_snp_factors(unpaired_strs_by_rg[i], base_quality_, snp_tree, log_p1, log_p2, local_match_count, local_mismatch_count);
	  log_p1s.push_back(log_p1); log_p2s.push_back(log_p2);
	}
	else {
	  std::vector<double> log_p1, log_p2;
	  for (unsigned int j = 0; j < paired_strs_by_rg[i].size()+unpaired_strs_by_rg[i].size(); ++j){
	    log_p1.push_back(0); log_p2.push_back(0); // Assign equal phasing LLs as no SNP info is available
	  }
	  log_p1s.push_back(log_p1); log_p2s.push_back(log_p2);
	  bad_samples.insert(rg_names[i]);
	}
	
	// Copy alignments
	alignments[i].insert(alignments[i].end(), paired_strs_by_rg[i].begin(),   paired_strs_by_rg[i].end());
	alignments[i].insert(alignments[i].end(), unpaired_strs_by_rg[i].begin(), unpaired_strs_by_rg[i].end());
      }
      logger << "Found VCF info for " << good_samples.size() << " out of " << good_samples.size()+bad_samples.size() << " samples with STR reads" << std::endl;
    }
    else 
      logger << "Warning: Failed to construct SNP trees for " << region_group.chrom() << ":" << region_group.start() << "-" << region_group.stop() << std::endl;
    destroy_snp_trees(snp_trees);      
  }
  if (!got_snp_info){
    for (unsigned int i = 0; i < paired_strs_by_rg.size(); i++){
      // Copy alignments                                                                                                                                             
      alignments[i].insert(alignments[i].end(), paired_strs_by_rg[i].begin(),   paired_strs_by_rg[i].end());
      alignments[i].insert(alignments[i].end(), unpaired_strs_by_rg[i].begin(), unpaired_strs_by_rg[i].end());
      
      // Assign equal phasing LLs as no SNP info is available
      log_p1s.push_back(std::vector<double>(paired_strs_by_rg[i].size()+unpaired_strs_by_rg[i].size(), 0.0));
      log_p2s.push_back(std::vector<double>(paired_strs_by_rg[i].size()+unpaired_strs_by_rg[i].size(), 0.0));
    }
  }
  
  int phased_samples = 0, phased_reads = 0, total_reads = 0;
  for (unsigned int i = 0; i < alignments.size(); i++){
    bool sample_phased = false;
    for (unsigned int j = 0; j < alignments[i].size(); j++){
      sample_phased |= (log_p1s[i][j] != log_p2s[i][j]);
      phased_reads  += (log_p1s[i][j] != log_p2s[i][j]);
    }
    total_reads    += alignments[i].size();
    phased_samples += sample_phased;
  }

  logger << "Phased SNPs add info for " << phased_reads << " out of " << total_reads << " reads"
		     << " and " << phased_samples << " out of " << rg_names.size() <<  " samples" << std::endl;

  double elapsed = elapsed_seconds(phase_start);
  locus_snp_phase_info_time_ = elapsed;
  if(phase_time_out != NULL)
    *phase_time_out = elapsed;
  else
    total_snp_phase_info_time_ += elapsed;

  {
    std::lock_guard<std::mutex> lock(snp_stats_mutex_);
    match_count_ += local_match_count;
    mismatch_count_ += local_mismatch_count;
  }

  return true;
}

int SNPBamProcessor::get_haplotype(BamAlignment& aln) const {
  if (!aln.HasTag(HAPLOTYPE_TAG.c_str()))
    return -1;
  int64_t haplotype;
  if (!aln.GetIntTag(HAPLOTYPE_TAG.c_str(), haplotype))
    printErrorAndDie("Failed to extract haplotype tag");
  assert(haplotype == 1 || haplotype == 2);
  return (int)haplotype;
}

/*
** Exploratory function for analyzing data from 10X Genomics BAMs
** These BAMs contain haplotype tags, which can be used in place of the physical-phasing + VCF approach
** used in the standard process_reads function
 */
void SNPBamProcessor::process_10x_reads(std::vector<BamAlnList>& paired_strs_by_rg,
					std::vector<BamAlnList>& mate_pairs_by_rg,
					std::vector<BamAlnList>& unpaired_strs_by_rg,
					const std::vector<std::string>& rg_names, const RegionGroup& region_group,
					const std::string& chrom_seq){
  auto phase_start = std::chrono::steady_clock::now();
  assert(paired_strs_by_rg.size() == mate_pairs_by_rg.size() && paired_strs_by_rg.size() == unpaired_strs_by_rg.size());

  std::vector<BamAlnList> alignments(paired_strs_by_rg.size());
  std::vector< std::vector<double> > log_p1s, log_p2s;
  int32_t phased_reads = 0, total_reads = 0;
  for (unsigned int i = 0; i < paired_strs_by_rg.size(); i++){
    // Copy alignments
    alignments[i].insert(alignments[i].end(), paired_strs_by_rg[i].begin(),   paired_strs_by_rg[i].end());
    alignments[i].insert(alignments[i].end(), unpaired_strs_by_rg[i].begin(), unpaired_strs_by_rg[i].end());

    log_p1s.push_back(std::vector<double>());
    log_p2s.push_back(std::vector<double>());
    for (unsigned int j = 0; j < paired_strs_by_rg[i].size(); j++){
      total_reads++;
      int haplotype_1 = get_haplotype(paired_strs_by_rg[i][j]);
      int haplotype_2 = get_haplotype(mate_pairs_by_rg[i][j]);

      // If the two mate pairs don't have the same haplotype index, it's possible that
      // i)  One of them is unmapped (and therefore has a -1)
      // ii) They map to two different phase sets. This is essentially a phasing breakpoint and we
      //     probably want to avoid using phase information for these reads
      int haplotype;
      if (haplotype_1 != haplotype_2)
	haplotype = -1;
      else
	haplotype = haplotype_1;
      if (haplotype != -1){
	phased_reads++;
	log_p1s[i].push_back(haplotype == 1 ? FROM_HAP_LL : OTHER_HAP_LL);
	log_p2s[i].push_back(haplotype == 2 ? FROM_HAP_LL : OTHER_HAP_LL);
      }
      else {
	log_p1s[i].push_back(0.0);
	log_p2s[i].push_back(0.0);
      }
    }
    for (unsigned int j = 0; j < unpaired_strs_by_rg[i].size(); j++){
      total_reads++;
      int haplotype = get_haplotype(unpaired_strs_by_rg[i][j]);
      if (haplotype != -1){
	phased_reads++;
	log_p1s[i].push_back(haplotype == 1 ? FROM_HAP_LL : OTHER_HAP_LL);
	log_p2s[i].push_back(haplotype == 2 ? FROM_HAP_LL : OTHER_HAP_LL);
      }
      else {
	log_p1s[i].push_back(0.0);
	log_p2s[i].push_back(0.0);
      }
    }
  }

  selective_logger() << "Phased SNPs add info for " << phased_reads << " out of " << total_reads << " reads" << std::endl;
  locus_snp_phase_info_time_ = elapsed_seconds(phase_start);
  total_snp_phase_info_time_ += locus_snp_phase_info_time_;

  // Run any additional analyses using phasing probabilities
  analyze_reads_and_phasing(alignments, log_p1s, log_p2s, rg_names, region_group, chrom_seq);
}
