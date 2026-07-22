##
## Makefile for all executables
##

##USE THESE COMMANDS FOR CLEAN AND COMPILE:
## make clean-all
## make -j4

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
CXXFLAGS= -O3 -g -D__STDC_LIMIT_MACROS -D_FILE_OFFSET_BITS=64 -std=c++20 -DMACOSX -pthread -Itaskflow  #-pedantic -Wunreachable-code -Weverything

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

LIBS = -L./ -lm -L$(HTSLIB_ROOT)/ -lz -lcurl -lcrypto -L$(CEPHES_ROOT)/ -llzma -lbz2 -Wl,--whole-archive $(MIMALLOC_LIB) -Wl,--no-whole-archive
INCLUDE   = -Ilib -Ilib/htslib -Itaskflow -I$(MIMALLOC_ROOT)/include
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

# Include auto-generated header dependencies when present.
-include $(DEP)

# ====================================================================
# 4. DEPENDENCY TRACKING: ENSURE HIPSTR REBUILDS IF MIMALLOC CHANGES
# ====================================================================
HipSTR: $(OBJ_COMMON) $(OBJ_HIPSTR) $(CEPHES_LIB) $(HTSLIB_LIB) $(MIMALLOC_LIB) $(OBJ_SEQALN)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $(filter-out $(MIMALLOC_LIB),$^) $(LIBS)

DenovoFinder: $(OBJ_DENOVO) $(HTSLIB_LIB) $(MIMALLOC_LIB)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $(filter-out $(MIMALLOC_LIB),$^) $(LIBS)
	
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

# Rebuild htslib library if needed
$(HTSLIB_LIB):
	cd lib/htslib && $(MAKE)

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
