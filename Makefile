##
## Makefile for all executables
##

##USE THESE COMMANDS FOR CLEAN AND COMPILE:
## make clean-all
## make -j4
##
## `make pgo` is also available (see the PGO section below) but is currently
## measured slower than a plain build on this toolchain -- not recommended.

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
## -flto=auto enables link-time optimization across all translation units.
## -fopenmp-simd enables recognition of "#pragma omp simd"/"declare simd"
## (mathops.cpp uses it to vectorize log_sum_exp's reduction via libmvec's
## vector exp -- see -lmvec below). It does not pull in libgomp or any
## OpenMP runtime; Taskflow remains the only threading in this codebase.
## The leading $(CXXFLAGS) keeps whatever the environment already set, instead
## of discarding it. A plain `=` assignment takes precedence over the
## environment in GNU Make, which silently dropped the flags that packaging
## toolchains rely on -- conda-build passes `-isystem $PREFIX/include` and
## `-L$PREFIX/lib -Wl,-rpath,$PREFIX/lib` this way to point the build at its own
## zlib/bzip2/xz/libcurl/openssl. The environment's flags go first so that this
## project's own settings win where the two conflict (conda-build passes -O2;
## later flags beat earlier ones in GCC, so -O3 below has to come after it).
## `:=` is required here: a recursive `=` that refers to itself is an infinite
## recursion error. Overriding on the command line (`make CXXFLAGS=...`) still
## wins over both, as documented above.
CXXFLAGS := $(CXXFLAGS) -O3 -g -flto=auto -fopenmp-simd -D__STDC_LIMIT_MACROS -D_FILE_OFFSET_BITS=64 -std=c++20 -DMACOSX -pthread -Itaskflow  #-pedantic -Wunreachable-code -Weverything

## To create a static distribution file, run:
##   make static-dist
ifeq ($(STATIC),1)
LDFLAGS := -static $(LDFLAGS)
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
# required -- see lib/libdeflate/build). -DHAVE_LIBDEFLATE, passed to the
# htslib build below, makes bgzf.c etc. call into it instead of zlib for
# BGZF/gzip; htslib's own config.h fallback does not set this on its own.
LIBDEFLATE_ROOT = lib/libdeflate
LIBDEFLATE_LIB  = $(LIBDEFLATE_ROOT)/build/libdeflate.a

# glibc's libmvec (the vector math library backing log_sum_exp's -fopenmp-simd
# exp reduction) was only added in glibc 2.22, so probe for it rather than
# hardcoding -lmvec. conda-forge/Bioconda builds default to a glibc 2.17
# sysroot, where the library does not exist -- GCC correspondingly emits no
# vector-math calls there, making the link both impossible and unnecessary.
MVEC_LIB := $(shell echo 'int main(){return 0;}' | $(CXX) -x c++ - -lmvec -o /dev/null 2>/dev/null && echo -lmvec)

LIBS = -L./ -lm $(MVEC_LIB) -L$(HTSLIB_ROOT)/ -lz -lcurl -lcrypto -L$(CEPHES_ROOT)/ -llzma -lbz2 $(LIBDEFLATE_LIB) -Wl,--whole-archive $(MIMALLOC_LIB) -Wl,--no-whole-archive
INCLUDE   = -Ilib -Ilib/htslib -Itaskflow -I$(MIMALLOC_ROOT)/include -I$(LIBDEFLATE_ROOT)
CEPHES_LIB        = lib/cephes/libprob.a
HTSLIB_LIB        = $(HTSLIB_ROOT)/libhts.a

# ====================================================================
# 2. MAIN BUILD TARGETS
# ====================================================================
.PHONY: all
all: $(MIMALLOC_LIB) HipSTR-MT DenovoFinder test/fast_ops_test test/haplotype_test test/read_vcf_alleles_test test/snp_tree_test test/vcf_snp_tree_test

# Create a tarball with static binaries
.PHONY: static-dist
static-dist:
	rm -f HipSTR-MT
	$(MAKE) STATIC=1
	( VER="$$(git describe --abbrev=7 --dirty --always --tags)" ;\
	  DST="HipSTR-MT-$${VER}-static-$$(uname -s)-$$(uname -m)" ; \
	  mkdir "$${DST}" && \
	        mkdir "$${DST}/scripts" && \
	        cp HipSTR-MT VizAln VizAlnPdf README.md "$${DST}" && \
	        cp scripts/filter_haploid_vcf.py scripts/filter_vcf.py scripts/generate_aln_html.py scripts/html_alns_to_pdf.py "$${DST}/scripts" && \
	        tar -czvf "$${DST}.tar.gz" "$${DST}" && \
	        rm -r "$${DST}/" \
	    )

version:
	git describe --abbrev=7 --dirty --always --tags | awk '{print "#include \"version.h\""; print "const std::string VERSION = \""$$0"\";"}' > src/version.cpp

# ====================================================================
# PROFILE-GUIDED OPTIMIZATION (PGO)
# ====================================================================
# `make pgo` builds HipSTR-MT twice: an instrumented pass trained on the
# bundled fixture in test/pgo/ (a real 1Mb chr20 slice, 654 STR loci, 2
# subsetted sample BAMs), then a final pass compiled against the resulting
# profile. Training and both compiles run locally, so the result is
# tuned for whatever machine ran `make pgo` rather than shipped as a
# prebuilt binary, and the profile can't go stale relative to the source.
# Override PGO_FASTA/PGO_BED/PGO_BAMS/PGO_SNP_VCF to train against a
# different workload, e.g.:
#   make pgo PGO_FASTA=/path/ref.fa PGO_BED=/path/regions.bed \
#            PGO_BAMS="/path/a.bam /path/b.bam"
#
# CAUTION -- measured on GCC 11.4 and 13.1: `make pgo` currently ships a
# binary ~5-7% SLOWER than a plain build on multi-thousand-region
# workloads. Isolating -flto from the comparison shows PGO alone is
# roughly neutral; the regression is specific to combining -fprofile-use
# with -flto=auto in CXXFLAGS above. Until that toolchain interaction is
# resolved, plain `make` is the faster build.
#
# -fprofile-update=atomic is required: the Taskflow worker pool updates
# profile counters from multiple threads during the instrumented
# training run, and counters would otherwise race.
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
	$(MAKE) HipSTR-MT CXXFLAGS="$(CXXFLAGS) -fprofile-generate=$(CURDIR)/$(PGO_DIR) -fprofile-update=atomic"
	mv HipSTR-MT HipSTR-MT.pgo-instrument
	./HipSTR-MT.pgo-instrument --bams $(PGO_BAM_LIST) --fasta $(PGO_FASTA) --regions $(PGO_BED) --str-vcf $(PGO_TRAIN_VCF) --threads 1 --min-reads 10
	./HipSTR-MT.pgo-instrument --bams $(PGO_BAM_LIST) --fasta $(PGO_FASTA) --regions $(PGO_BED) --str-vcf $(PGO_TRAIN_VCF) --threads 4 --min-reads 10
	./HipSTR-MT.pgo-instrument --bams $(PGO_BAM_LIST) --fasta $(PGO_FASTA) --regions $(PGO_BED) --str-vcf $(PGO_TRAIN_VCF) --threads $(shell nproc) --min-reads 10
	./HipSTR-MT.pgo-instrument --bams $(PGO_BAM_LIST) --fasta $(PGO_FASTA) --regions $(PGO_BED) --str-vcf $(PGO_TRAIN_VCF) --snp-vcf $(PGO_SNP_VCF) --threads 4 --min-reads 10
	rm -f HipSTR-MT.pgo-instrument
	$(MAKE) clean
	$(MAKE) HipSTR-MT CXXFLAGS="$(CXXFLAGS) -fprofile-use=$(CURDIR)/$(PGO_DIR) -fprofile-correction -Wno-coverage-mismatch -Wno-missing-profile"
	mv HipSTR-MT HipSTR-MT.pgo-tmp
	$(MAKE) clean
	mv HipSTR-MT.pgo-tmp HipSTR-MT
	$(MAKE) DenovoFinder
	@echo "[pgo] HipSTR-MT rebuilt with profile-guided optimization ($(PGO_DIR)/)"

# Clean the generated files of the main project only
.PHONY: clean
clean:
	rm -f *~ src/*.o src/*.d src/*~ src/SeqAlignment/*~ src/SeqAlignment/*.o src/SeqAlignment/*.d src/denovos/*~ src/denovos/*.o src/denovos/*.d HipSTR-MT DenovoFinder test/allele_expansion_test test/fast_ops_test test/haplotype_test test/read_vcf_alleles_test test/snp_tree_test test/vcf_snp_tree_test

# ====================================================================
# 3. ADD AUTOMATION TO CLEAN THE MIMALLOC ARTIFACTS
# ====================================================================
.PHONY: clean-all
clean-all: clean
	cd lib/htslib && $(MAKE) clean
	rm -f lib/cephes/*.o $(CEPHES_LIB)
	rm -rf $(MIMALLOC_ROOT)/build
	rm -rf $(LIBDEFLATE_ROOT)/build
	rm -rf pgo-data HipSTR-MT.pgo-instrument

# Include auto-generated header dependencies when present.
-include $(DEP)

# ====================================================================
# 4. DEPENDENCY TRACKING: ENSURE HIPSTR REBUILDS IF MIMALLOC CHANGES
# ====================================================================
HipSTR-MT: $(OBJ_COMMON) $(OBJ_HIPSTR) $(CEPHES_LIB) $(HTSLIB_LIB) $(MIMALLOC_LIB) $(LIBDEFLATE_LIB) $(OBJ_SEQALN)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $(filter-out $(MIMALLOC_LIB) $(LIBDEFLATE_LIB),$^) $(LIBS)

DenovoFinder: $(OBJ_DENOVO) $(HTSLIB_LIB) $(MIMALLOC_LIB) $(LIBDEFLATE_LIB)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $(filter-out $(MIMALLOC_LIB) $(LIBDEFLATE_LIB),$^) $(LIBS)
	
PhasingChecker: src/check_phasing.cpp src/region.cpp src/error.cpp src/haplotype_tracker.cpp src/version.cpp src/pedigree.cpp src/vcf_reader.cpp src/stringops.cpp $(HTSLIB_LIB)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/haplotype_test: test/haplotype_test.cpp src/SeqAlignment/Haplotype.cpp src/SeqAlignment/HapBlock.cpp src/SeqAlignment/NeedlemanWunsch.cpp src/error.cpp src/stringops.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/em_stutter_test: test/em_stutter_test.cpp src/em_stutter_genotyper.cpp src/genotyper.cpp src/genotyper_bam_processor.cpp src/error.cpp src/mathops.cpp src/stringops.cpp src/stutter_model.cpp
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
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDE) -MMD -MP -o $@ -c $<

# Rebuild CEPHES library if needed
$(CEPHES_LIB):
	cd lib/cephes && $(MAKE)

# Rebuild htslib library if needed. Needs libdeflate built first so its
# header/lib are present for HAVE_LIBDEFLATE (see lib/htslib/config.h).
$(HTSLIB_LIB): $(LIBDEFLATE_LIB)
	cd lib/htslib && $(MAKE) lib-static CPPFLAGS="$(CPPFLAGS) -I$(CURDIR)/$(LIBDEFLATE_ROOT) -DHAVE_LIBDEFLATE"

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
