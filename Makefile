##
## Makefile for all executables
##

##USE THESE COMMANDS FOR CLEAN AND COMPILE:
## make clean-all
## make -j4
##
## For a build tuned with profile-guided optimization (recommended for
## production/deployment on the target machine):
## make pgo

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
CXXFLAGS= -O3 -g -flto=auto -D__STDC_LIMIT_MACROS -D_FILE_OFFSET_BITS=64 -std=c++20 -DMACOSX -pthread -Itaskflow  #-pedantic -Wunreachable-code -Weverything

## To create a static distribution file, run:
##   make static-dist
ifeq ($(STATIC),1)
LDFLAGS=-static 
else
LDFLAGS= 
endif

## Source code files, add new files to this list
SRC_COMMON  = src/base_quality.cpp src/error.cpp src/region.cpp src/stringops.cpp src/zalgorithm.cpp src/alignment_filters.cpp src/extract_indels.cpp src/mathops.cpp src/pcr_duplicates.cpp src/bam_io.cpp src/adapter_trimmer.cpp
SRC_HIPSTR  = src/hipstr_main.cpp src/bam_processor.cpp src/stutter_model.cpp src/snp_phasing_quality.cpp src/snp_tree.cpp src/em_stutter_genotyper.cpp src/seq_stutter_genotyper.cpp src/snp_bam_processor.cpp src/genotyper_bam_processor.cpp src/vcf_input.cpp src/read_pooler.cpp src/version.cpp src/haplotype_tracker.cpp src/pedigree.cpp src/vcf_reader.cpp src/genotyper.cpp src/directed_graph.cpp src/debruijn_graph.cpp src/fasta_reader.cpp src/vcf_writer.cpp
SRC_SEQALN  = src/SeqAlignment/HapAligner.cpp src/SeqAlignment/AlignmentModel.cpp src/SeqAlignment/AlignmentOps.cpp src/SeqAlignment/HapBlock.cpp src/SeqAlignment/NeedlemanWunsch.cpp src/SeqAlignment/Haplotype.cpp src/SeqAlignment/HaplotypeGenerator.cpp src/SeqAlignment/HTMLCreator.cpp src/SeqAlignment/AlignmentViz.cpp src/SeqAlignment/AlignmentTraceback.cpp src/SeqAlignment/StutterAlignerClass.cpp
SRC_DENOVO  = src/denovos/denovo_main.cpp src/error.cpp src/stringops.cpp src/version.cpp src/pedigree.cpp src/haplotype_tracker.cpp src/vcf_input.cpp src/denovos/denovo_scanner.cpp src/mathops.cpp src/vcf_reader.cpp src/denovos/denovo_allele_priors.cpp src/denovos/trio_denovo_scanner.cpp

SRC = $(SRC_COMMON) $(SRC_HIPSTR) $(SRC_SEQALN) $(SRC_DENOVO)

# For each CPP file, generate an object file
OBJ_COMMON  := $(SRC_COMMON:.cpp=.o)
OBJ_HIPSTR  := $(SRC_HIPSTR:.cpp=.o)
OBJ_SEQALN  := $(SRC_SEQALN:.cpp=.o)
OBJ_DENOVO  := $(SRC_DENOVO:.cpp=.o)
DEP         := $(SRC:.cpp=.d)

CEPHES_ROOT=lib/cephes
HTSLIB_ROOT=lib/htslib

# ====================================================================
# 1. DEFINE PATHS AND ARTIFACTS FOR THE LOCAL MIMALLOC BUILD
# ====================================================================
MIMALLOC_ROOT = lib/mimalloc
# Target the static library produced by mimalloc's build system
MIMALLOC_LIB  = $(MIMALLOC_ROOT)/build/libmimalloc.a

# Locally vendored libdeflate (CMake-only upstream build, no system package
# required -- see lib/libdeflate/build). htslib's config.h defines
# HAVE_LIBDEFLATE so bgzf.c etc. call into it instead of zlib for BGZF/gzip.
LIBDEFLATE_ROOT = lib/libdeflate
LIBDEFLATE_LIB  = $(LIBDEFLATE_ROOT)/build/libdeflate.a

LIBS = -L./ -lm -L$(HTSLIB_ROOT)/ -lz -lcurl -lcrypto -L$(CEPHES_ROOT)/ -llzma -lbz2 $(LIBDEFLATE_LIB) -Wl,--whole-archive $(MIMALLOC_LIB) -Wl,--no-whole-archive
INCLUDE   = -Ilib -Ilib/htslib -Itaskflow -I$(MIMALLOC_ROOT)/include -I$(LIBDEFLATE_ROOT)
CEPHES_LIB        = lib/cephes/libprob.a
HTSLIB_LIB        = $(HTSLIB_ROOT)/libhts.a

# ====================================================================
# 2. MAIN BUILD TARGETS
# ====================================================================
.PHONY: all
all: $(MIMALLOC_LIB) HipSTR DenovoFinder test/fast_ops_test test/haplotype_test test/read_vcf_alleles_test test/snp_tree_test test/vcf_snp_tree_test

# Create a tarball with static binaries
.PHONY: static-dist
static-dist:
	rm -f HipSTR
	$(MAKE) STATIC=1
	( VER="$$(git describe --abbrev=7 --dirty --always --tags)" ;\
	  DST="HipSTR-$${VER}-static-$$(uname -s)-$$(uname -m)" ; \
	  mkdir "$${DST}" && \
	        mkdir "$${DST}/scripts" && \
	        cp HipSTR VizAln VizAlnPdf README.md "$${DST}" && \
	        cp scripts/filter_haploid_vcf.py scripts/filter_vcf.py scripts/generate_aln_html.py scripts/html_alns_to_pdf.py "$${DST}/scripts" && \
	        tar -czvf "$${DST}.tar.gz" "$${DST}" && \
	        rm -r "$${DST}/" \
	    )

version:
	git describe --abbrev=7 --dirty --always --tags | awk '{print "#include \"version.h\""; print "const std::string VERSION = \""$$0"\";"}' > src/version.cpp

# ====================================================================
# PROFILE-GUIDED OPTIMIZATION (PGO)
# ====================================================================
# `make pgo` builds HipSTR twice: an instrumented pass trained on the
# bundled fixture in test/pgo/ (a 1Mb chr20 slice with a matching FASTA,
# 654 STR loci, and 13 subsetted BAMs -- the original 12-sample cluster
# plus a 20%-downsampled real single-sample BAM covering the full 1Mb
# span), then a final pass compiled against the resulting profile. This
# trains and rebuilds locally on whatever machine runs `make pgo`, so
# the result is tuned for that hardware rather than shipped as a
# prebuilt binary. The profile is regenerated from the fixture on every
# `make pgo`, so it can never go stale relative to the current source.
#
# The fixture was widened from an earlier 3-locus/~400kb version (whose
# BAMs, it turned out, also had reads rebased into a fake local
# coordinate frame that didn't match its own bundled SNP VCF -- the
# --snp-vcf training pass was silently exercising zero SNPs). 654 real
# chr20 loci across a real, coordinate-consistent coverage profile is a
# small, git-friendly fixture (~9MB total) that's actually representative.
# Override PGO_FASTA/PGO_BED/PGO_BAMS/PGO_SNP_VCF on the command line to
# train against a larger or different workload, e.g.:
#   make pgo PGO_FASTA=/path/ref.fa PGO_BED=/path/regions.bed \
#            PGO_BAMS="/path/a.bam /path/b.bam"
#
# CAUTION -- measured on this GCC 11.4 toolchain: `make pgo` currently
# ships a binary ~7% SLOWER than plain `make` on realistic multi-
# thousand-region workloads, regardless of fixture quality (verified
# with both the old 3-locus fixture and the new 654-locus one -- same
# regression either way) and regardless of -flto-partition (tried
# `=one` as well as the default). Isolating -flto from the comparison
# shows PGO alone is roughly neutral (within run-to-run noise); it's
# specifically the combination of -fprofile-use with -flto=auto in
# CXXFLAGS above that regresses, i.e. an LTO+PGO codegen interaction on
# this toolchain, not a training-data problem. Until that's resolved
# (a newer GCC, or a flag combination that avoids it), plain `make` is
# the faster build -- don't reach for `make pgo` by default here.
#
# -fprofile-update=atomic is required (not optional) because the
# Taskflow worker pool updates profile counters from multiple threads
# concurrently during the instrumented training run; without it,
# counters can be corrupted by races between workers.
PGO_DIR       = pgo-data
PGO_FASTA     = test/pgo/pgo_fixture.fa
PGO_BED       = test/pgo/pgo_fixture.bed
PGO_BAMS      = $(wildcard test/pgo/bams/*.pgo.bam)
PGO_SNP_VCF   = test/pgo/pgo_fixture_snps.vcf.gz
empty         :=
space         := $(empty) $(empty)
comma         := ,
PGO_BAM_LIST  := $(subst $(space),$(comma),$(PGO_BAMS))
PGO_TRAIN_VCF = $(PGO_DIR)/train.vcf.gz

.PHONY: pgo
pgo:
	rm -rf $(PGO_DIR)
	mkdir -p $(PGO_DIR)
	$(MAKE) clean
	$(MAKE) HipSTR CXXFLAGS="$(CXXFLAGS) -fprofile-generate=$(CURDIR)/$(PGO_DIR) -fprofile-update=atomic"
	mv HipSTR HipSTR.pgo-instrument
	./HipSTR.pgo-instrument --bams $(PGO_BAM_LIST) --fasta $(PGO_FASTA) --regions $(PGO_BED) --str-vcf $(PGO_TRAIN_VCF) --threads 1 --min-reads 10
	./HipSTR.pgo-instrument --bams $(PGO_BAM_LIST) --fasta $(PGO_FASTA) --regions $(PGO_BED) --str-vcf $(PGO_TRAIN_VCF) --threads 4 --min-reads 10
	./HipSTR.pgo-instrument --bams $(PGO_BAM_LIST) --fasta $(PGO_FASTA) --regions $(PGO_BED) --str-vcf $(PGO_TRAIN_VCF) --threads $(shell nproc) --min-reads 10
	./HipSTR.pgo-instrument --bams $(PGO_BAM_LIST) --fasta $(PGO_FASTA) --regions $(PGO_BED) --str-vcf $(PGO_TRAIN_VCF) --snp-vcf $(PGO_SNP_VCF) --threads 4 --min-reads 10
	rm -f HipSTR.pgo-instrument
	$(MAKE) clean
	$(MAKE) HipSTR CXXFLAGS="$(CXXFLAGS) -fprofile-use=$(CURDIR)/$(PGO_DIR) -fprofile-correction -Wno-coverage-mismatch -Wno-missing-profile"
	mv HipSTR HipSTR.pgo-tmp
	$(MAKE) clean
	mv HipSTR.pgo-tmp HipSTR
	$(MAKE) DenovoFinder
	@echo "[pgo] HipSTR rebuilt with profile-guided optimization ($(PGO_DIR)/)"

# Clean the generated files of the main project only
.PHONY: clean
clean:
	rm -f *~ src/*.o src/*.d src/*~ src/SeqAlignment/*~ src/SeqAlignment/*.o src/SeqAlignment/*.d src/denovos/*~ src/denovos/*.o src/denovos/*.d HipSTR DenovoFinder test/allele_expansion_test test/fast_ops_test test/haplotype_test test/read_vcf_alleles_test test/snp_tree_test test/vcf_snp_tree_test

# ====================================================================
# 3. ADD AUTOMATION TO CLEAN THE MIMALLOC ARTIFACTS
# ====================================================================
.PHONY: clean-all
clean-all: clean
	cd lib/htslib && $(MAKE) clean
	rm -f lib/cephes/*.o $(CEPHES_LIB)
	rm -rf $(MIMALLOC_ROOT)/build
	rm -rf $(LIBDEFLATE_ROOT)/build
	rm -rf pgo-data HipSTR.pgo-instrument

# Include auto-generated header dependencies when present.
-include $(DEP)

# ====================================================================
# 4. DEPENDENCY TRACKING: ENSURE HIPSTR REBUILDS IF MIMALLOC CHANGES
# ====================================================================
HipSTR: $(OBJ_COMMON) $(OBJ_HIPSTR) $(CEPHES_LIB) $(HTSLIB_LIB) $(MIMALLOC_LIB) $(LIBDEFLATE_LIB) $(OBJ_SEQALN)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $(filter-out $(MIMALLOC_LIB) $(LIBDEFLATE_LIB),$^) $(LIBS)

DenovoFinder: $(OBJ_DENOVO) $(HTSLIB_LIB) $(MIMALLOC_LIB) $(LIBDEFLATE_LIB)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $(filter-out $(MIMALLOC_LIB) $(LIBDEFLATE_LIB),$^) $(LIBS)
	
PhasingChecker: src/check_phasing.cpp src/region.cpp src/error.cpp src/haplotype_tracker.cpp src/version.cpp src/pedigree.cpp src/vcf_reader.cpp src/stringops.cpp $(HTSLIB_LIB)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/haplotype_test: test/haplotype_test.cpp src/SeqAlignment/Haplotype.cpp src/SeqAlignment/HapBlock.cpp src/SeqAlignment/NeedlemanWunsch.cpp src/error.cpp src/stringops.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/em_stutter_test: test/em_stutter_test.cpp src/em_stutter_genotyper.cpp src/genotyper_bam_processor.cpp src/error.cpp src/mathops.cpp src/stringops.cpp src/stutter_model.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/fast_ops_test: test/fast_ops_test.cpp src/mathops.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^

test/read_vcf_alleles_test: test/read_vcf_alleles_test.cpp src/error.cpp src/region.cpp src/vcf_input.cpp src/vcf_reader.cpp $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/snp_tree_test: src/snp_tree.cpp src/error.cpp test/snp_tree_test.cpp src/haplotype_tracker.cpp src/vcf_reader.cpp $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/vcf_snp_tree_test: test/vcf_snp_tree_test.cpp src/error.cpp src/snp_tree.cpp src/haplotype_tracker.cpp src/vcf_reader.cpp $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# Build each object file independently
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -MMD -MP -o $@ -c $<

# Rebuild CEPHES library if needed
$(CEPHES_LIB):
	cd lib/cephes && $(MAKE)

# Rebuild htslib library if needed. Needs libdeflate built first so its
# header/lib are present for HAVE_LIBDEFLATE (see lib/htslib/config.h).
$(HTSLIB_LIB): $(LIBDEFLATE_LIB)
	cd lib/htslib && $(MAKE) lib-static CPPFLAGS="-I$(CURDIR)/$(LIBDEFLATE_ROOT) -DHAVE_LIBDEFLATE"

# ====================================================================
# 5b. THE BUILD RECIPE FOR LIBDEFLATE (vendored; CMake-only upstream build)
# ====================================================================
$(LIBDEFLATE_LIB):
	mkdir -p $(LIBDEFLATE_ROOT)/build && cd $(LIBDEFLATE_ROOT)/build && \
	cmake -DCMAKE_BUILD_TYPE=Release -DLIBDEFLATE_BUILD_SHARED_LIB=OFF -DLIBDEFLATE_BUILD_GZIP=OFF .. && \
	$(MAKE) -j4

# ====================================================================
# 5. THE BUILD RECIPE FOR MIMALLOC (Handles Old System CMake Versions)
# ====================================================================
$(MIMALLOC_LIB):
	@echo "Checking system CMake version..."
	@CMAKE_BIN="cmake"; \
	CMAKE_VERSION=$$(cmake --version 2>/dev/null | head -n 1 | grep -oE '[0-9]+\.[0-9]+' | head -n 1); \
	MAJOR=$$(echo $$CMAKE_VERSION | cut -d. -f1); \
	MINOR=$$(echo $$CMAKE_VERSION | cut -d. -f2); \
	if [ -z "$$MAJOR" ] || [ $$MAJOR -lt 3 ] || { [ $$MAJOR -eq 3 ] && [ $$MINOR -lt 18 ]; }; then \
		echo "System CMake is too old or missing. Downloading local portable CMake..."; \
		mkdir -p $(MIMALLOC_ROOT)/cmake_local; \
		wget -qO- https://github.com/Kitware/CMake/releases/download/v3.26.4/cmake-3.26.4-linux-x86_64.tar.gz | tar -xzf - -C $(MIMALLOC_ROOT)/cmake_local --strip-components=1; \
		CMAKE_BIN="$$(pwd)/$(MIMALLOC_ROOT)/cmake_local/bin/cmake"; \
	fi; \
	echo "Building mimalloc using: $$CMAKE_BIN"; \
	mkdir -p $(MIMALLOC_ROOT)/build && cd $(MIMALLOC_ROOT)/build && \
	$$CMAKE_BIN -DMI_BUILD_SHARED=OFF -DMI_BUILD_OBJECT=OFF -DMI_BUILD_TESTS=OFF .. && \
	$(MAKE) -j4
