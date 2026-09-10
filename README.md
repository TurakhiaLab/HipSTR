# HipSTR-MT
**H**aplotype **i**nference and **p**hasing for **S**hort **T**andem **R**epeats  
![HipSTR icon!](https://raw.githubusercontent.com/tfwillems/HipSTR/master/img/HipSTR_icon_small.png)
[![DOI](https://zenodo.org/badge/1273590821.svg)](https://doi.org/10.5281/zenodo.22650334)

#### Author: Thomas Willems <hipstrtool@gmail.com>

#### Optimizer: Joachim Galil <joachimbgalil@gmail.com>, Turakhia Lab -- the Gymrek Lab asked for the performance work in this fork

#### License: GNU v2

*This document is written in ASD-STE100 Simplified Technical English.*

[Introduction](#introduction)  
[Requirements](#requirements)  
[Installation](#installation)  
[Testing](#testing)  
[Quick Start](#quick-start)       
[HipSTR-MT Changes](#hipstr-mt-changes)  
[Tutorial](#tutorial)  
[In-depth Usage](#in-depth-usage)  
[Data Requirements](#data-requirements)  
[Phasing](#phasing)     
[Speed](#speed)  
[Default Filtering](#default-filtering)  
[Call Filtering](#call-filtering)  
[Additional Usage Options](#additional-usage-options)		  
[File Formats](#file-formats)     
[FAQ](#faq)     
[Help](#help)       
[Citation](#citation)

## Introduction
Short tandem repeats [(STRs)](http://en.wikipedia.org/wiki/Microsatellite) are genomic sequences that contain many copies of a short motif. The genomes of most organisms contain STRs. STRs are important, because they change much more quickly than other genomic elements. Thus they give much information for genomic identification, ancestry inference, and genealogy.

But STRs are difficult to genotype. The repetitive sequence that causes the high rate of change also causes frequent alignment errors. These errors make the subsequent analyses difficult, and can cause a bias. PCR stutter errors also cause reads that contain more repeat copies, or fewer repeat copies, than the correct genotype.

**HipSTR-MT** decreases the effect of these errors and gives more reliable STR genotypes. HipSTR-MT does these four functions:

1. It learns a PCR stutter model for each locus with an [EM algorithm](http://en.wikipedia.org/wiki/Expectation-maximization_algorithm).
2. It finds candidate STR alleles in population-scale sequencing data.
3. It aligns the reads to the candidate alleles with a special hidden Markov model that includes the STR artifacts.
4. It uses phased SNP haplotypes to genotype and to phase the STRs.

We think that these four functions make **HipSTR-MT** the most reliable tool to genotype STRs from **Illumina** sequencing data.

## Requirements

- **Operating system: Linux (x86_64)**. The memory-usage report reads `/proc/self/status`, and the CPU-affinity detection uses `sched_getaffinity`. Only Linux has these two APIs, and this program has no alternative code for other systems. Nobody has built or run this program on macOS, BSD, or Windows. The unmodified upstream HipSTR, on which this fork is based, has the same limit.
- **make**
- **g++ with C++20 support** — Taskflow, the parallelization library that this fork adds, needs C++20.
- **CMake 3.18 or a subsequent version** — the builds of the local copies of mimalloc and libdeflate need CMake, because these two projects have only a CMake build system. If the system has no CMake, or has a version before 3.18, the Makefile gets a portable CMake automatically.
- **zlib** (dev headers)
- **libbz2** (dev headers)
- **liblzma** (dev headers)
- **libcurl and OpenSSL** (dev headers) — the default build configuration of the local htslib 1.24 always compiles the remote-file support that uses libcurl.

This repository contains its own copy of the HTSlib source in `lib/htslib`, and the build makes HTSlib from that source. The repository also contains its own copies of mimalloc, libdeflate, and Taskflow in `lib/`. Thus the system does not need a `libhts`/`htslib` package, and the build does not use one.

On Ubuntu and Debian, this command installs all the packages in the list above. It does not install make and g++, because these two programs are usually already installed.

    apt install make g++ cmake zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev

The CI workflow ([.github/workflows/ci.yml](.github/workflows/ci.yml)) runs this same command before each build. Thus the command stays correct for the compilation of this project.

## Installation
This repository contains the Taskflow headers and the necessary mimalloc source directly, in the same manner as `lib/htslib`. These libraries are not git submodules. Thus a usual clone is sufficient. You do not need `--recurse-submodules`, and you do not need to remember `git submodule update --init --recursive`.

    git clone https://github.com/TurakhiaLab/HipSTR-MT

To build the program, use Make:

    cd HipSTR-MT
    make

The command makes an executable file with the name **HipSTR-MT** in the current directory. As part of the same command, it also builds mimalloc. To see the detailed help, use this command:

    ./HipSTR-MT --help

The Makefile now writes compiler dependency files with `-MMD -MP`. Thus a change to a header in `src`, `src/SeqAlignment`, or `src/denovos` causes a new build of the related object files.

### Building with profile-guided optimization
    make pgo

First, `make pgo` compiles an instrumented `HipSTR-MT`. Then it trains that program with a fixture in this repository. Then it compiles the final `HipSTR-MT` again, with the profile from the training. The fixture (`test/pgo/`) contains a real cluster of chr20 STR loci of approximately 1 Mb, with 654 loci, a related FASTA file, and 2 subsetted sample BAM files. The training and the two compilations occur on your machine. Thus the binary is tuned for the machine that ran `make pgo`, and the binary does not contain a profile from different hardware. `make pgo` does not change `DenovoFinder`: the build makes `DenovoFinder` again with the usual flags. The correctness regression check uses this fixture and one more -- refer to [Testing](#testing).

**Caution: do not use `make pgo` at this time.** On GCC 11.4 and GCC 13.1, `make pgo` makes a binary that is 5-7% *slower* than the binary from the usual `make`. An unwanted interaction between `-fprofile-use` and `-flto=auto` causes this decrease in speed. Refer to the PGO section of the Makefile for more data. Use the usual `make` until this interaction has a solution.

### Conda/Bioconda
A draft recipe is in [recipe/](recipe/meta.yaml) for a subsequent submission to [bioconda-recipes](https://github.com/bioconda/bioconda-recipes). The recipe uses the source tarball of the applicable release. The recipe does not contain the sha256 checksum of that tarball, because a tarball cannot contain its own checksum. The GitHub Release notes of each tag give the checksum, and the recipe shows the command that calculates it again. You must add the checksum when you copy the recipe into a bioconda-recipes pull request. Nobody has sent the recipe to Bioconda yet. Until the recipe is available there, build the program from the source as above.

## Testing
    test/run_tests.sh

`test/run_tests.sh` is self-contained. It builds the components that it needs, and it uses only the fixture data in `test/`. It does not download data.

**Correctness regression** (`test/check_correctness.sh`). This is the primary check, and a JOSS-style review will most probably examine it first. The script runs `HipSTR-MT` against two fixtures in this repository. It uses more than one value of `--threads` (the default values are 1, 2, 4, and 8, but the script does not use more threads than the number of cores in the machine). For each fixture, it compares the output of each thread count against one golden VCF file with `test/compare_vcf_tolerant.py`. Each golden VCF file is the output of the unmodified upstream `gymrek-lab/HipSTR` on the same fixture. Thus this one comparison checks the two conditions on which the correctness claim of this fork is based:

1. The genotype calls do not change when the thread count changes.
2. The genotype calls agree with the unmodified upstream `gymrek-lab/HipSTR` -- refer to [Correctness](#hipstr-mt-changes) below. The most recent check found 0 discrepancies on the full tutorial trio of 599 loci, on a full-genome NA12891 run of 1,512,240 loci, and on a full-genome NA12892 run of 1,171,158 loci, at `--threads` 1 to 64.

The two fixtures test different risks:

- **`test/pgo/`** gives breadth. It is the same real chr20 fixture of approximately 1 Mb that `make pgo` uses for its training. It has 654 STR loci and 2 subsetted sample BAM files -- refer to [Building with profile-guided optimization](#building-with-profile-guided-optimization).
- **`test/homopolymer/`** gives depth. It has one locus only, but that locus is difficult. It is a homopolymer of 20 T bases (chr2:33759762 in hg19) from sample NA12892. At a locus of this type, two or more alignment paths have almost equal probability. Thus a very small floating-point difference can change which candidate alleles the program finds, and this changes the genotype call itself. Such an error is much worse than a small change to a statistic. The 654-locus fixture contains no locus of this type: it gave a correct result both before and after we corrected an error of exactly this type in `fast_exp_sum` (refer to the comments in `src/mathops.cpp`). A full-genome run of a second sample was necessary to find that error. This fixture is only 424 KB, and it makes an error of this type fail the test immediately. The fixture uses the coordinates of a 10 kb slice of chr2, not the coordinates of the chromosome, to keep the FASTA file small.

`compare_vcf_tolerant.py` fails if a genotype call changes, or if a derived statistic changes more than 0.1%. It permits a change of less than 0.1%, because the order of the floating-point summation causes this small noise. To use different thread counts, give them as arguments: `test/check_correctness.sh 1 16 64`.

**Unit tests** (`test/*_test.cpp`). You can build each test independently with `make test/<name>`.
- **`snp_tree_test`** has a true pass/fail assertion. It builds the same SNP set in two different ways: with a brute-force scan, and with the interval-tree structure that `snp_bam_processor.cpp` uses. Then it asserts that the results of the two queries agree.
- **`fast_ops_test`**, **`haplotype_test`**, and **`read_vcf_alleles_test`** are diagnostic tests. They do not assert. They print their output (tables of the approximation error, the generated haplotype sequences, and the parsed VCF alleles) for a manual examination. `read_vcf_alleles_test` needs the `test/input/1kg.chr1.imputed.vcf.gz` fixture from this repository (the run above includes this fixture).
- **`em_stutter_test`** and **`vcf_snp_tree_test`** build correctly, but the default run does not include them. The intended driver of `em_stutter_test` (`run_stutter_em_test.sh`) makes its input with an external STR-mutation simulator, and this repository does not contain that simulator. `vcf_snp_tree_test` needs a chr22:10-20Mb VCF file that this repository does not contain. Both tests accept a VCF path or a BED path as an argument, if you want to use your own data.

## Quick Start
The most generally applicable mode of HipSTR-MT processes **all the samples together**. Use this syntax:

```
./HipSTR-MT --bams          run1.bam,run2.bam,run3.bam,run4.bam
         --fasta         genome.fa
         --regions       str_regions.bed
         --str-vcf       str_calls.vcf.gz
         --threads       8
```

* **bams** :  a list of [BAM/CRAM](#bams) files, separated by commas. Make these files with [BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml), then sort them and index them with [samtools](http://www.htslib.org/).
* **regions** : a [BED](#str-bed) file that contains the coordinates of each applicable STR region. You can download BED files for different organisms, and also for humans, from [here](https://github.com/HipSTR-Tool/HipSTR-references/).
* **fasta** : a [FASTA file](https://en.wikipedia.org/wiki/FASTA_format) that contains the sequence of each chromosome in the BED file. The coordinates of this genome build must agree with the coordinates of the STR regions.
* **str-vcf** : the output path for the STR genotypes.
* **threads** : the number of worker threads in the Taskflow executor. If you do not give this option, HipSTR-MT selects a default value from the hardware. It uses the CPU allocation of the scheduler, the Linux CPU affinity, or `std::thread::hardware_concurrency()`. HipSTR-MT keeps four pipeline lines in operation for each worker, to hide the latency of the serial fetch stage and the serial write stage.

For each region in *str_regions.bed*, **HipSTR-MT** does these three steps:

1. It learns a stutter model for the locus.
2. It uses the stutter model and the haplotype-based alignment algorithm to genotype each person.
3. It writes the STR genotypes to *str_calls.vcf.gz*, a [bgzipped](http://www.htslib.org/doc/tabix.html) [VCF](#str-vcf) file. This VCF file contains the calls for each sample in the read groups of the BAM/CRAM files.

## HipSTR-MT Changes
HipSTR-MT is a performance fork of [gymrek-lab/HipSTR](https://github.com/gymrek-lab/HipSTR). All the data below is a comparison with that baseline.

**Correctness**: we compared the output of HipSTR-MT with the output of the HipSTR version of the Gymrek Lab. The comparator permits small differences: the genotype calls must agree exactly, but a float-type field can have a relative difference of less than 1e-3. The genotype calls and all the derived statistics now agree exactly (0 difference) on the tutorial trio of 599 loci, on a full-genome NA12891 run of 1,512,240 loci, and on a full-genome NA12892 run of 1,171,158 loci. These comparisons found two correctness bugs, and we corrected both. The first was in the target_clones dispatch -- refer to [Vectorization](#vectorization). The second was a single-precision accumulator in `fast_exp_sum`: it changed the alleles that the program found at one homopolymer locus of the 1,171,158 loci of NA12892. The `test/homopolymer` fixture now holds that locus -- refer to [Testing](#testing).

### Parallelization
- `bam_processor.*` replaces the single-region loop with a Taskflow pipeline of three stages: serial creation of the region tokens, parallel read filtration and genotyping, and serial output in the correct sequence. The worker threads process the regions in a different sequence, but the program writes them in the sequence of the BED file. Each pipeline line has its own `BamCramMultiReader` instance and its own `AdapterTrimmer` instance. Each line keeps its pass BAM records and its filter BAM records in a buffer, and does not write them immediately. All the lines share one cached FASTA chromosome sequence, and no line makes its own copy.
- `snp_bam_processor.*` moves the preparation of the SNP phasing into the work item of the pipeline. It adds two mutexes: `snp_stats_mutex_` for the aggregate counters, and `snp_phase_mutex_` for the shared reader state and phasing state. Thus a `--snp-vcf` run is safe across the worker threads.
- `genotyper_bam_processor.*` keeps the VCF, log, visualization, stutter, timing, and BAM output of each region in a `RegionResult`. It merges the counters, and it writes all the data in the BED sequence from the serial output stage.
- `seq_stutter_genotyper.*` adds `build_vcf_record` and `build_vcf_records`. These functions make the VCF line of a locus as a `BuiltVCFRecord` string, and do not write to a stream. Thus a worker thread can complete a region without the output lock.
- `hipstr_main.cpp` adds `--threads <num_threads>`. If you do not give this option, `default_thread_count()` selects a default value from the hardware, in this sequence: the environment variables for the CPU allocation of the scheduler (`SLURM_CPUS_PER_TASK`, `SLURM_CPUS_ON_NODE`, `PBS_NP`, `NSLOTS`, `OMP_NUM_THREADS`), then the Linux CPU affinity (`sched_getaffinity`), then `std::thread::hardware_concurrency()`. The pipeline keeps `4 * threads` region contexts in operation.

### Thread-safety fixes that this required
Two parts of the initial single-threaded code kept a mutable state. This state is safe with only one caller, but it is not safe when the worker threads share it.
- **`StutterAlignerClass`**: its scratch buffers (`ins_probs_`, `del_probs_`, `match_probs_`, `log_probs_`) were instance members, and each `load_read()` call made them again. The haplotype/block structure owns these objects and the threads share them. `HapAligner` is different, because each thread has its own `HapAligner`. Thus concurrent `load_read()` calls caused a race on the shared state. To correct this, we moved the buffers into a `StutterWorkspace` that the caller supplies (one workspace for each `HapAligner`). The methods of `StutterAlignerClass` are now `const`, and they do not touch a shared mutable state. This change also added memoization: if a subsequent `load_read()` call has the same arguments, the function does not calculate the result again. This condition is frequent when the program uses the same alignments for more than one candidate haplotype.
- **`bdtr` from Cephes** (the binomial CDF, which `compute_allele_bias` uses) keeps an internal state, and concurrent calls are not safe. `seq_stutter_genotyper.cpp` now protects that call with a `std::mutex`.

### Memory optimizations
- **`HapAligner`** uses the same scratch buffers for all its reads (the base-quality arrays, the DP matrices, and the buffers for the artifact size and the artifact position). It does not do `new[]`/`delete[]` for each read. The threads do not share a `HapAligner` instance. The two largest allocations for each read are the DP matrices for the match, the insertion, and the deletion, and each matrix is `O(read_len × haplotype_len)`. These matrices are now interleaved in one buffer (`MatrixChannel`, `[match0, insert0, deletion0, match1, ...]`) for better cache locality.
- **ASCII-only conversion of the case** replaces the locale-aware `toupper()`/`tolower()` in the loops that process each base. This change is in the `uppercase()` function of `stringops.cpp`, in the CIGAR-driven comparison of the bases in `AlignmentOps.cpp`, in `base_to_int()` in `NeedlemanWunsch.cpp`, and in the comparisons of the prefix, the suffix, and the end match in `zalgorithm.cpp` and `alignment_filters.cpp`.
- **`mathops.cpp`** adds a pointer-pair variant of `fast_log_sum_exp` (`const double* begin, const double* end`), together with the initial `vector<double>` variant. Thus some callers do not make a copy of a vector.
- **mimalloc** is linked by default (refer to Installation). It decreases the overhead of the allocator, because the program makes very many small allocations for each read and for each locus.
- **The chromosome cache** lets the threads share the chromosomes. The program uses the chromosomes in sequence. Thus, when all the threads move to the subsequent chromosome, the cache removes the chromosome that no thread uses. This decreases the memory usage.

### Vectorization
- **`-flto=auto`** enables the link-time optimization across all the translation units. In our benchmarks this is 7-8% faster and the output is identical.
- **The `sum`, `log_sum_exp`, and `fast_log_sum_exp` functions of `mathops.cpp`** are compiled with `__attribute__((target_clones("avx512f,avx2,sse4.2,default")))`. The compiler makes one copy of the function for each ISA in the list, and the program selects the best copy for the CPU at run time. This is portable across machines, and this build does not use `-march=native`. The user does not configure it.
- **The `exp()` reduction of `log_sum_exp`** is vectorized with libmvec from glibc, which is correctly rounded and is not an approximation. The code uses `#pragma omp simd` and `-fopenmp-simd`, which do not add an OpenMP runtime.
- **`fast_log_sum_exp`** processes 4 elements at the same time with `vfasterexp()`. This SSE-vectorized helper function was already in `fastonebigheader.h`, but the program did not use it.
- All three functions that use `target_clones` also have `__attribute__((optimize("-ffp-contract=off")))`. Without this attribute, the avx512f clone and the avx2 clone let the compiler fuse the multiply-add operations, but the default clone does not. Thus identical source code could round differently, only because of the ISA clone that the CPU selects at run time. A full-genome test showed this effect: at some homopolymer STR loci it changed the log-probability of a stutter-block candidate, and at one locus it changed the alleles that the program found as candidates. The attribute corrects this, and the program keeps the ISA dispatch for each machine.

### Dependency and I/O fixes
- **htslib is now version 1.24, and was version 1.9.** The `fai_retrieve()` function of the local 1.9 copy read the FASTA sequence one byte at a time (a `bgzf_getc()` call and a locale-aware `isgraph()` test for each byte). The chromosome load occurs in the necessary serial stage of the pipeline. Thus more worker threads did not decrease this cost: on the tutorial data it was approximately 66% of the minimum wall-clock time at high thread counts. Version 1.24 supplies the block-read implementation that upstream already corrected, and this repository does not have a local patch. The total FASTA load time across chr1-22 decreased from 7.17 s to 1.36 s. The wall time at `--threads 24` decreased from approximately 10.85 s to approximately 5.5-6.7 s. The newer source needed two small corrections for C++ compatibility: an explicit cast in `cram/cram_io.h` (an implicit conversion from `void*` is correct in C, but not in C++), and a missing `<unistd.h>` include in `bam_io.h` and `denovo_main.cpp` for `access()` and `F_OK`. The transitive includes of htslib 1.9 concealed these two problems.
- **libdeflate was linked, but the program did not use it.** The build rule of `HTSLIB_LIB` did not have `-DHAVE_LIBDEFLATE`. Thus the compiler did not include the libdeflate code paths of `bgzf.c`, and all the BGZF/BAM decompression used the system zlib. The build rule now has the define, which corrects this problem.

### New and restored CLI flags
- **`--threads <num_threads>`** — refer to Parallelization above.
- **`--lib-from-samp`** — get the library of each read from its sample name. Thus each read group does not need an `LB` tag.
- **`--output-hap-fields`** — write the additional `LFLANKS`, `RFLANKS`, `HQ`, `PHQ`, `LFGT`, and `RFGT` fields, which give data about the full assembled haplotypes.

## Tutorial
This [tutorial](https://hipstr-tool.github.io/HipSTR-tutorial/) shows how to apply HipSTR-MT to whole-genome sequencing data. In less than 10 minutes, it shows how to genotype approximately 600 STRs in a trio of persons with deep sequencing, and how to examine the results.

## In-depth Usage
**HipSTR-MT** has different usage options for sequencing data with a different number of samples and a different coverage. Most conditions are in one of these four categories:

1. 100 or more samples with low coverage (approximately 5x)
    * There are sufficient reads for the stutter estimation.
    * There are sufficient reads to find the candidate STR alleles.
    * [**Use the de novo stutter estimation and the STR calling with de novo allele generation**](#mode-1)
2. 20 or more samples with high coverage (approximately 30x)
    * There are sufficient reads for the stutter estimation.
    * There are sufficient reads to find the candidate STR alleles.
    * [**Use the de novo stutter estimation and the STR calling with de novo allele generation**](#mode-1)
3. Few samples with low coverage (approximately 5x)
    * There are insufficient reads for the stutter estimation.
    * There are insufficient reads to find the candidate STR alleles.
    * [**Use external or default stutter models and the STR calling with a reference panel**](#mode-3)
4. Few samples with high coverage (approximately 30x)
    * There are insufficient samples for the stutter estimation.
    * There are sufficient reads to find the candidate STR alleles.
    * [**Use external or default stutter models and the STR calling with de novo allele generation**](#mode-2)

<a id="mode-1"></a>

#### Mode 1: De novo stutter estimation + STR calling with de novo allele generation
This mode is the same as the mode in the **Quick Start** section, because it is applicable to most tasks. HipSTR-MT writes the STR genotypes in the bgzipped VCF format to *str_calls.vcf.gz*.

```
./HipSTR-MT --bams             run1.bam,run2.bam,run3.bam,run4.bam
         --fasta            genome.fa
         --regions          str_regions.bed
         --str-vcf          str_calls.vcf.gz
```

<a id="mode-2"></a>

#### Mode 2: External stutter models + STR calling with de novo allele generation
In this mode, the program does not learn the stutter models with the EM algorithm. It reads them from the **stutter-in** file. This is the only difference. For more data about the format of the stutter model file, refer to [below](#stutter-file).

```
./HipSTR-MT --bams             run1.bam,run2.bam,run3.bam,run4.bam
         --fasta            genome.fa
         --regions          str_regions.bed
         --stutter-in       ext_stutter_models.txt
         --str-vcf          str_calls.vcf.gz
```
If you do not have external stutter models for the **stutter-in** option, use **def-stutter-model**. This option uses one simple stutter model for all the loci (the help message of HipSTR-MT gives the parameters).

<a id="mode-3"></a>

#### Mode 3: External stutter models + STR calling with a reference panel
This mode is almost the same as mode 2, but you also supply a VCF file with known STR genotypes for each locus in the **str-vcf** option. With this option, **HipSTR-MT** does not find more candidate STR alleles in the BAM/CRAM files. Thus it is better to use a VCF file with STR genotypes from many different populations and persons.

```
./HipSTR-MT --bams             run1.bam,run2.bam,run3.bam,run4.bam
         --fasta            genome.fa
         --regions          str_regions.bed
         --stutter-in       ext_stutter_models.txt
         --ref-vcf          ref_strs.vcf.gz
         --str-vcf          str_calls.vcf.gz
```

If you do not have external stutter models for the **stutter-in** option, use **def-stutter-model**. This option uses one simple stutter model for all the loci (the help message of HipSTR-MT gives the parameters).

## Data Requirements
**HipSTR-MT** needs Illumina sequencing data to genotype STRs. The depth of the sequencing and the read length can be very different in different data sets. Thus this section gives the important factors to think about before you make data for a HipSTR-MT analysis.

STRs are repetitive. Thus a read that does not go fully across the repeat gives only a minimum value for the length of the repeat. This minimum value is useful, and **HipSTR-MT** uses it. But accurate and reliable STR genotypes need reads that go fully across the repetitive sequence (*spanning reads*). The number of reads that span an STR is a function of the read length, the sequencing depth, the length of the repeat, and other factors. The relation between these factors is complex, but [**Figure 2** in a recent review](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4254273/figure/F2/) by *Press et al.* shows it clearly. Longer reads and higher sequencing coverage increase the number of spanning reads. A longer repeat decreases the number of spanning reads, and thus a long STR is more difficult to genotype accurately. If the number of spanning reads is less than approximately 10, there is a risk that all the reads come from only 1 of the 2 copies of the chromosome. In that condition, it is not possible to call the two alleles of a heterozygous person correctly.

In our tests, 100 bp Illumina reads are sufficient for most of the STRs in the human genome. But an STR that is longer than 70 bp (for example a very long forensic STR) always needs longer reads. Most human STRs are shorter than 70 bp. This read length is also sufficient for most model organisms, if their repeats are not much longer than the human repeats.

The best minimum sequencing depth is a function of your analysis. If you examine how the STRs change, or if you [want to find de novo mutations](#de-novo-mutations), use a minimum coverage of 30x. This coverage lets HipSTR-MT give the necessary high specificity. If you only examine the allele frequencies of the STRs in a population, a coverage of 10x is usually sufficient. But at this coverage there are many genotype errors, in which the program calls a heterozygous genotype as a homozygous genotype.

**min-reads** is an important option for your type of sequencing data. HipSTR-MT does not genotype an STR if the number of reads for all the persons is less than *N*. The default value is 100, because this is a good minimum threshold to learn a stutter model before the genotyping. If you analyze very few samples (for example one trio of a mother, a father, and a child), decrease this threshold, because you seldom have 100 reads. In this condition, use options such as **--min-reads 15 --def-stutter-model**. The second option uses a default stutter model, because there are too few reads to learn an accurate model. But if you analyze many samples (for example more than ten 30x genomes, or more than thirty 10x genomes), do not change this parameter. In these conditions, a region with less than 100 reads can have a high GC content that is difficult for Illumina sequencing, can be difficult to map, or can be too long for your read length.

## Phasing
HipSTR-MT uses phased SNP haplotypes to phase the STR genotypes. To do this, it looks for pairs of reads in which the read with the STR, or its mate, overlaps a heterozygous SNP of the sample. In these conditions, the quality score of the overlapped base gives the probability that the read comes from each haplotype. If this data is not available, the program gives the read an equal probability for the two haplotypes. The genotyping model of HipSTR-MT includes these probabilities, and it writes phased genotypes. The *PQ* FORMAT field gives the quality of the phasing, which is the posterior probability of the phased genotype of the sample. For a homozygous genotype, this value is always equal to the *Q* FORMAT field, because the phasing is not applicable. For a heterozygous genotype, a *PQ* value that is almost equal to *Q* shows that one of the two phasings is much more probable. But if no read of a sample overlaps a heterozygous SNP, the two phasings are equally probable, and *PQ* is almost equal to *Q*/2. To use the physical phasing, give HipSTR-MT the **snp-vcf** option and a SNP VCF file with **phased** haplotypes. This schematic shows the concepts of the physical phasing model of HipSTR-MT:

![Phasing schematic!](https://raw.githubusercontent.com/tfwillems/HipSTR/master/img/phasing.png)

## Speed
We did a full hg19 genotyping run (the sample NA12891/ERR194160, 1,512,240 STR loci in the full genome, two Xeon Silver 4216 CPUs, 64 logical CPUs). HipSTR-MT at `--threads 64` completed the run in 24.6 minutes. The unmodified upstream HipSTR, with one thread, completed the same run in 13.3 hours. Thus HipSTR-MT is 32.5 times more quick, and its output is identical for each genotype (refer to [Correctness](#hipstr-mt-changes)). This value also includes the optimizations below that are not related to the parallelization (LTO, SIMD, htslib, and mimalloc): even at `--threads 1`, HipSTR-MT is measurably more quick than the initial program. The scaling is almost linear to 16 threads (94% efficiency), and it decreases at 64 threads (51% efficiency), because the SMT contention and the serial stages of the pipeline become more important. The benchmark harness of this fork makes thread-scaling plots with the full data for the run time, the speedup, the memory, and the CPU utilization.

HipSTR-MT has multithreading at the region level. Use `--threads N` to set the number of workers in the Taskflow executor. If you do not give `--threads`, the program selects a default value in this sequence: the environment variables for the CPU allocation of the scheduler (for example `SLURM_CPUS_PER_TASK`), then the Linux CPU affinity, then `std::thread::hardware_concurrency()`. The pipeline keeps `4 * N` region contexts in operation. Thus the worker threads can continue to genotype while the serial stages get the subsequent region or write the completed output.

These are the most important internal optimizations in this fork:

1. The program processes independent regions concurrently, but it writes the VCF, log, BAM, visualization, and stutter-model output in the sequence of the BED file.
2. The program caches each chromosome FASTA sequence one time and shares it with all the pipeline lines. Thus each line does not make a large copy of the contig.
3. Each `HapAligner` uses the same DP matrices for the haplotype alignment more than one time. Thus the program does not do millions of allocations in the hot path of the read alignment.
4. The build links mimalloc by default, which decreases the overhead of the allocator in the read processing and the haplotype processing.
5. htslib 1.24 replaces a FASTA reader that read one byte at a time in the serial stage of the pipeline. The new block-read implementation removes a bottleneck that limited the scaling at high thread counts (refer to [Dependency and I/O fixes](#dependency-and-io-fixes)).
6. The link-time optimization and the portable CPU dispatch at run time (refer to [Vectorization](#vectorization)) make the numeric hot path more quick. They do not need `-march=native` or a configuration for each machine.

For a large run, start with a `--threads` value that is near the number of physical cores available to the job. Then do a benchmark with a small set of representative regions.

Before the `--threads` option was available, the only method to parallelize the single-threaded HipSTR was to divide the work manually across more than one OS process. `--threads N` replaces that method on one machine. The two options below are still useful, but only to divide the work *across* different machines or jobs (for example an HPC array). They are not an alternative to `--threads` on one machine.

Option 1: Analyze each chromosome in parallel with the **--chrom** option. For example, **--chrom chr2** genotypes only the BED regions on chr2.

Option 2: Divide your BED file into *N* files and analyze the *N* files in parallel. This method is similar to option 1, but it can be more quick if *N* is much larger than the number of chromosomes.

## Default Filtering
HipSTR-MT sometimes filters the genotype of a sample automatically, and then writes a missing value in the VCF file. The program applies these filters when the data of a sample shows that it cannot make a reliable genotype. For each locus, the **log** file gives a summary of the number of filtered samples. If you give the **--output-filters** option, the VCF file has a FORMAT field with the name **FILTER** for each sample. In this field, *PASS* shows a satisfactory sample, and other values give the reason for the filtration.

**A sample with a PASS value still needs the additional variant filtration (refer to the next section). PASS shows only that no catastrophic problem occurred during the genotyping.** This table gives the possible reasons for the filtration:

| Filter | Explanation 
| :----- | :---------
| NO_READS                | No alignment was available for the sample at this STR. If reads overlap the STR in the BAM/CRAM file, the program possibly filtered them because of the read quality, the mapping uniqueness, or other reasons
| FLANK_ASSEMBLY_CYCLIC   | During the genotyping, HipSTR-MT tries to assemble the sequences before and after the STR (the *flanks*), to find possible SNPs that it must include. This assembly fails if the assembly graph contains a cycle, which causes this filter
| FLANK_ASSEMBLY_INDEL    | The assembly found an insertion or a deletion in the *flanks*. These indels are a problem for the model of HipSTR-MT, and thus the program does not genotype the sample
| FLANK_INDEL_FRAC        | When the genotyping is complete, HipSTR-MT finds the maximum-likelihood alignment of each read against the called alleles of its sample. If a large fraction of these alignments have an indel in the *flanks*, the alignments are probably incorrect, and thus the program ignores the genotype of the sample
| LOW_FREQUENCY_ALT_FLANK | The program collects the flanking sequences from the assembly of each sample, and makes all the candidate haplotypes from them. The number of haplotypes increases exponentially with the number of these sequences. Thus HipSTR-MT decreases the time and discards a flank that is present in only a few samples. If the data of a sample supports a low-frequency flank, the program does not genotype that sample. To change this frequency threshold, use the **--min-flank-freq** option 



## Call Filtering
**HipSTR-MT** decreases many of the most frequent causes of STR genotyping errors. But it is still very important to filter the VCF file and to discard the calls with a low quality. To help you with this task, the VCF output contains FORMAT fields and INFO fields that usually show a problem with a call. An INFO field gives the aggregate data for a locus, and some values show that you must possibly discard the full locus. A FORMAT field gives the data for one sample at one locus, and some values show that you must possibly discard the genotypes of some samples. The list below gives some of these fields and their meaning. You can also filter a VCF file with most of the fields below with the [dumpSTR](https://trtools.readthedocs.io/en/stable/source/dumpSTR.html) utility from the [TRTools package](https://trtools.readthedocs.io/en/stable/). The Gymrek lab maintains TRTools, which has utilities to filter, merge, and calculate statistics on the VCF files from HipSTR-MT and from other STR genotypers. Upstream recommends this method.

#### INFO fields:  
1. **DP**: the total depth, that is the number of informative reads for all the samples at the locus. To calculate the mean coverage for each sample, divide this value by the number of samples that have a genotype. A genotype with a low mean coverage is not reliable, because the reads possibly show only one of the two alleles of a heterozygous person.
2. **DSTUTTER**: the total number of reads at a locus that have a stutter artifact, as HipSTR-MT calculates it. If the fraction of reads with stutter (DSTUTTER/DP) is high, the genotypes of the locus are not reliable, because the reads frequently do not show the correct genotype. Three conditions can cause a high fraction of reads with stutter: too much PCR amplification, a duplicated locus that maps to one location in the genome, or a failure of HipSTR-MT to find sufficient candidate alleles.  
3. **DFLANKINDEL**: the total number of reads whose maximum-likelihood alignment contains an indel in the regions that flank the STR. A high fraction of reads with this artifact (DFLANKINDEL/DP) can be the result of a true indel in a region near the STR. But it can also occur if HipSTR-MT does not find sufficient candidate alleles. If the true alleles have a very different size from the candidate alleles, or are not multiples of the repeat unit, the program frequently aligns them as indels in the flanking sequences.

#### FORMAT fields:  
1. **Q**: the posterior probability of the genotype. In our experience this is the best indicator of the quality of the genotype of one sample, and we almost always use it to filter the calls.   
2. **DP**, **DSTUTTER**, and **DFLANKINDEL**: these fields are also available for each sample, with the same meaning as the INFO fields above. Use them in the same manner to find calls with a problem.  
3. **AB** and **FS**: the log10 of the p-value of the allele bias and of the Fisher strand bias. A large negative value shows that the bias is very improbable by chance. The two fields compare the read counts between the two haplotype *copies* of the sample (the two parental chromosomes, which the phased SNPs identify). They do not compare the read counts between different STR allele values. Thus they have a meaning only for a diploid sample with phasing data, and they can be non-zero for a homozygous STR call if the reads divide unequally between the two copies. For **AB**, an improbable division shows that the read counts for each haplotype copy do not agree with the predicted genotype. For **FS**, it shows a relation between the sequencing strand and the copy of each read, which shows that sequencing errors possibly increase one of the reported alleles.   

**What thresholds do we recommend for these fields?** The correct thresholds are a function of the quality of the sequencing data, the ploidy of the chromosome, and your subsequent analyses. As an initial value, dumpSTR options such as `--hipstr-min-call-Q 0.9 --hipstr-max-call-flank-indel 0.15 --hipstr-max-call-stutter 0.15` are satisfactory. Alternatively, this repository contains the initial filtration scripts in the **scripts** subdirectory. These scripts apply the same type of thresholds, and they do not need TRTools.

```
python scripts/filter_vcf.py  --vcf                   diploid_calls.vcf.gz
                              --min-call-qual         0.9
                              --max-call-flank-indel  0.15
                              --max-call-stutter      0.15
                  --min-call-allele-bias  -2
                  --min-call-strand-bias  -2
    
python scripts/filter_haploid_vcf.py  --vcf                   haploid_calls.vcf.gz
                                      --min-call-qual         0.9
                                      --max-call-flank-indel  0.15
                                      --max-call-stutter      0.15
```

The script writes the new VCF file to the standard output. For each sample, the script removes a call if one or more of these conditions are true: the posterior is less than 90%, more than 15% of the reads have a flank indel, or more than 15% of the reads have a stutter artifact. For a diploid VCF file, these filters also remove a genotype whose allele bias p-value or Fisher strand bias p-value is less than 0.01 (10^-2). The script replaces the call of each sample that fails these conditions with the missing symbol "." from the VCF specification. To see more filtration options, use one of these two commands:

```
python scripts/filter_vcf.py -h
python scripts/filter_haploid_vcf.py -h
```

## Additional Usage Options

| Option  | Description  
| :------- | :----------- 
| **viz-out**       aln_viz.gz     | Write a file with the alignments of each locus, for a visualization with VizAln or [VizAlnPdf](#aln-viz) <br> **Why? You want to see or examine the STR genotypes**
| **log**         log.txt               | Write the log data to the given file (Default = the standard error)  
| **threads** num_threads                | The number of worker threads in the Taskflow executor (Default = automatic) <br> **Why? You want to replace the hardware-aware default of HipSTR-MT.** This option replaces the method in the [Speed](#speed) section, in which you divide your BED file manually and run more than one initial-HipSTR process in parallel. `--threads` does the same operation internally, in one process, and it writes the output in the sequence of the BED file.
| **haploid-chrs**  list_of_chroms      | A list of chromosomes, separated by commas, that the program must process as haploid (Default = all chromosomes are diploid) <br> **Why? You analyze a haploid chromosome such as chrY**  
| **no-rmdup**                            | Do not remove the PCR duplicates. By default, the program removes them <br> **Why? Your sequencing data is for PCR-amplified regions**  
| **use-unpaired**                        | Use the unpaired reads for the genotyping (Default = False) <br> **Why? Your sequencing data contains only single-end reads**  
| **snp-vcf**    phased_snps.vcf.gz     | A bgzipped input VCF file with the phased SNP genotypes of the samples to genotype. The program uses these SNPs to phase the STRs physically <br> **Why? You have phased SNP genotypes**  
| **bam-samps**     list_of_read_groups | A list of samples, separated by commas, in the same sequence as the BAM files. <br> The program gives each read the sample of its file. By default, <br> each read must have an RG tag, and the SM field gives the sample <br> **Why? The RG tags of your BAM file do not have an SM field**  
| **bam-libs**      list_of_read_groups | A list of libraries, separated by commas, in the same sequence as the BAM files. <br> The program gives each read the library of its file. By default, <br> each read must have an RG tag, and the LB field gives the library <br> NOTE: This option is necessary when you give --bam-samps <br> **Why? The RG tags of your BAM file do not have an LB tag**  
| **lib-from-samp**                       | Give each read the library of its sample name. Thus each read group does not need an LB tag <br> **Why? The RG tags of your BAM file do not have an LB field, and one library for each sample is sufficient**  
| **def-stutter-model**                   | For each locus, use a stutter model with PGEOM=0.9 and UP=DOWN=0.05 for the in-frame artifacts, and PGEOM=0.9 and UP=DOWN=0.01 for the out-of-frame artifacts <br> **Why? You have too few samples for the stutter estimation, and you do not have stutter models**  
| **min-reads** num_reads                           | 	The minimum number of reads to genotype a locus (Default = 100) <br> **Why? Refer to the discussion [above](#data-requirements)**  
|**output-filters**                        | Write the reason for the filtration of each call to the VCF file (Default = False)
| **output-hap-fields**                    | Write the additional `LFLANKS`, `RFLANKS`, `HQ`, `PHQ`, `LFGT`, and `RFGT` FORMAT fields, which give the full assembled haplotypes of each sample and not only the STR part (Default = False) <br> **Why? You want to examine the alignment of the flanking sequence, and not only the STR genotype**


This list gives the most useful and most frequent options, but it is not complete. To see all the options, use this command:

    ./HipSTR-MT --help

<a id="aln-viz"></a>

## Alignment Visualization
To examine an STR call, it is very useful to see the reads that support it. The **viz-out** option of HipSTR-MT writes a compressed file with the alignments of each call. The **VizAln** command in the main directory of HipSTR-MT shows these alignments. Before you can see the alignments, you must index the file with tabix.
For example, if you ran HipSTR-MT with the option `--viz-out aln.viz.gz`, use this command:

    tabix -p bed aln.viz.gz

The command makes a [tabix](http://www.htslib.org/doc/tabix.html) index for the file. With this index, the program can extract the alignments of a locus quickly. Use this command only one time, after the file is complete.

You can then see the calls of the sample *NA12878* at the locus *chr1 3784267* with this command:

    ./VizAln aln.viz.gz chr1 3784267 NA12878

The command opens an image of the alignments in your browser automatically. The image can be similar to this example:
![Read more words!](https://raw.githubusercontent.com/HipSTR-Tool/HipSTR-tutorial/master/viz_NA12878.png)
The top bar is the reference sequence. The red text gives the name of the sample and its call at the locus. The other rows give the alignment of each read that the program used for the genotyping. In this example, 14 reads have an *8 bp deletion* and 14 reads have a *4 bp insertion*. Thus HipSTR-MT calls the genotype of this sample as *-8 | 4*.

To see all the calls at the same locus, use this command:

    ./VizAln aln.viz.gz chr1 3784267

There is also a similar script that makes a PDF file from these alignments, for a publication. This script accepts only one sample. To make the image above in a file with the name alignments.pdf, use this command:

    ./VizAlnPdf aln.viz.gz chr1 3784267 NA12878 alignments 1

NOTE: The **viz-out** file can become very large if you genotype thousands of loci or thousands of samples. Thus it can be better to run HipSTR-MT again with this option, and to use only the loci that you want to see.

## File Formats
<a id="bams"></a>

### BAM/CRAM files
HipSTR-MT needs [BAM/CRAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files from an indel-sensitive aligner. You must sort these files by position with the `samtools sort` command, and then index them with `samtools index`. To find the sample of a read, HipSTR-MT uses the read group data in the header lines of the BAM/CRAM file. Each *@*RG line must have an *ID* field, an *LB* field with the library, and an *SM* field with the sample. For example, if a BAM/CRAM file has this header line

    @RG     ID:RUN1 LB:ERR12345        SM:SAMPLE789

the program gives an alignment with this RG tag

    RG:Z:RUN1

the sample *SAMPLE789* and the library *ERR12345*. Thus HipSTR-MT can analyze BAM/CRAM files that contain more than one sample and more than one library. It can also analyze the reads of one sample that are in more than one file.

If your BAM/CRAM files do not have *RG* data, use the **bam-samps** and **bam-libs** flags to give the sample and the library of each file. But in this condition, a BAM/CRAM file can contain only one library and one read group. For example, this command

```
./HipSTR-MT --bams             run1.bam,run2.bam,run3.bam,run4.cram
         --fasta            genome.fa
         --regions          str_regions.bed
         --str-vcf          str_calls.vcf.gz
         --bam-samps        SAMPLE1,SAMPLE1,SAMPLE2,SAMPLE3
         --bam-libs         LIB1,LIB2,LIB3,LIB4
```

tells HipSTR-MT to give all the reads in the first two BAM files the sample *SAMPLE1*, all the reads in the third file the sample *SAMPLE2*, and all the reads in the last file the sample *SAMPLE3*.


HipSTR-MT can analyze BAM files and CRAM files in the same run. If your project has the two file types, HipSTR-MT does the CRAM decompression automatically. **Caution: when you analyze CRAM files, the file that you give to --fasta must be the same FASTA file that you used to make the CRAM files.** If it is different, the CRAM decompression usually fails, and the program can operate incorrectly.

<a id="str-bed"></a>

### STR region BED file
The BED file with the applicable STR regions is a tab-delimited file. It has 5 necessary columns and one optional column:

1. The name of the chromosome of the STR
2. The start position of the STR on its chromosome
3. The end position of the STR on its chromosome
4. The length of the motif, that is the number of bases in the repeat unit
5. The number of copies of the repeat unit in the reference allele

The 6th column is optional. It contains the name of the STR locus, which the program writes to the ID column of the VCF file.
This example file contains 5 STR loci.

**NOTE: The header of the table is only for the explanation. The BED file must not have a header.**

CHROM | START       | END         | MOTIF_LEN | NUM_COPIES | NAME
----  | ----        | ----        | ---       | ---        | ---
chr1  | 13784267    | 13784306    | 4         | 10         | GATA27E01
chr1  | 18789523    | 18789555    | 3         | 11         | ATA008
chr2  | 32079410    | 32079469    | 4         | 15         | AGAT117
chr17 | 38994441    | 38994492    | 4         | 12         | GATA25A04
chr17 | 55299940    | 55299992    | 4         | 13         | AAT245

We supply *BED* files with the STR loci of different organisms, and also of humans, [here](https://github.com/HipSTR-Tool/HipSTR-references/).

For other model organisms, use the same [method](https://github.com/HipSTR-Tool/HipSTR-references/blob/master/mouse/mouse_reference.md)
that we used to make the mouse BED file.

<a id="str-vcf"></a>

### VCF file
For more data about the VCF file format, refer to the [VCF spec](http://samtools.github.io/hts-specs/VCFv4.2.pdf).

#### INFO fields
The INFO fields contain aggregate statistics about each genotyped STR in the VCF file. The INFO fields of HipSTR-MT give the learned or supplied stutter model of the locus, the reference coordinates of the STR (START and END), the allele counts (AC), and the number of reads that the program used to genotype all the samples (DP).

FIELD | DESCRIPTION
----- | -----------
INFRAME_PGEOM  | Parameter of the geometric distribution of the step size for the in-frame changes
INFRAME_UP     | Probability that stutter causes an in-frame increase of the observed STR size
INFRAME_DOWN   | Probability that stutter causes an in-frame decrease of the observed STR size
OUTFRAME_PGEOM | Parameter of the geometric distribution of the step size for the out-of-frame changes
OUTFRAME_UP    | Probability that stutter causes an out-of-frame increase of the observed STR size
OUTFRAME_DOWN  | Probability that stutter causes an out-of-frame decrease of the observed STR size
BPDIFFS        | Base pair difference between each alternate allele and the reference allele
START          | Inclusive start coordinate of the repetitive part of the reference allele
END            | Inclusive end coordinate of the repetitive part of the reference allele
PERIOD         | Length of the STR motif
AN             | Total number of alleles in the called genotypes
REFAC          | Count of the reference allele
AC             | Counts of the alternate alleles
NSKIP          | Number of samples that the program did not genotype, because of different problems
NFILT          | Number of samples that the program genotyped and then filtered
DP             | Total number of reads that the program used to genotype all the samples
DSNP           | Total number of reads with SNP data
DSTUTTER       | Total number of reads with a stutter indel in the STR region
DFLANKINDEL    | Total number of reads with an indel in the regions that flank the STR

#### FORMAT fields
The FORMAT fields give data about the genotype of each sample at the locus. HipSTR-MT gives the most probable phased genotype (GT), the posterior probability of that genotype (PQ), and the unphased equivalent (Q). Other useful data is the number of reads that the program used for the genotype (DP), and the alignment artifacts of these reads (DSTUTTER and DFLANKINDEL).

FIELD     | DESCRIPTION
--------- | -----------
GT        | Genotype
GB        | Base pair differences between the genotype and the reference
Q         | Posterior probability of the unphased genotype
PQ        | Posterior probability of the phased genotype
DP        | Number of valid reads that the program used for the genotype of the sample
DSNP      | Number of reads with SNP phasing data
PDP       | Fraction of the reads that support each haploid genotype
GLDIFF    | Difference in likelihood between the reported genotype and the second most probable genotype
DSNP      | Total number of reads with SNP data
PSNP      | Number of reads with SNPs that support each haploid genotype
DSTUTTER  | Number of reads with a stutter indel in the STR region
DFLANKINDEL | Number of reads with an indel in the regions that flank the STR
AB        | log10 of the p-value of the allele bias. A value of 0 shows no bias, and a more negative value shows more bias. This field compares the read counts between the two haplotype *copies* of the sample (the chromosomes), and not between the STR allele values. Thus it can be negative for a homozygous STR genotype (the same allele on the two copies), if the phased SNPs near the STR identify the two copies and the reads divide unequally between them
FS        | log10 of the p-value of the strand bias from the exact test of Fisher. A value of 0 shows no bias, and a more negative value shows more bias. This field has the same condition for each haplotype copy as AB above: it can be negative for a homozygous STR genotype, if the reads of the two phased copies divide unequally between the two strands
DAB       | Number of reads that the program used for the calculation of the allele bias
ALLREADS  | Base pair difference in the Needleman-Wunsch alignment of each read
MALLREADS | Maximum-likelihood base pair difference of each read, from the haplotype alignments
GL        | log-10 likelihoods of the genotypes
PL        | Phred-scaled likelihoods of the genotypes

<a id="stutter-file"></a>

### Stutter model
To model the PCR stutter artifacts, we use three types of stutter artifact:

1. **In-frame changes**: these change the size of the STR in the read by a multiple of the repeat unit. For example, if the repeat motif is AGAT, an in-frame change causes a difference of -8, -4, 4, 8, and so on.
2. **Out-of-frame changes**: these change the size of the STR by a value that is not a multiple of the repeat unit. For example, if the repeat motif is AGAT, an out-of-frame change causes a difference of -3, -2, -1, 1, 2, 3, and so on.
3. **No stutter change**: the size of the STR in the read is the same as the size of the true STR.


A stutter model file contains the data for these three artifacts. It is a **tab-delimited BED-like** file with exactly 10 columns. All the columns are necessary: `StutterModel::read()` in `src/stutter_model.cpp` cannot parse the file if one column is not present, and this includes `PERIOD`. This is an example of such a file:

CHROM  | START       | END      | IGEOM | IDOWN | IUP   | OGEOM | ODOWN | OUP   | PERIOD
-----  | ----------- | -------- | ----  | ----  | ---   | ----  | ---   | ---   | ---
chr1   | 13784267    | 13784306 | 0.95  | 0.05  | 0.01  | 0.9   | 0.01  | 0.001 | 4
chr1   | 18789523    | 18789555 | 0.8   | 0.01  | 0.05  | 0.9   | 0.001 | 0.001 | 3
chr2   | 32079410    | 32079469 | 0.9   | 0.01  | 0.01  | 0.9   | 0.001 | 0.001 | 4
chr17  | 38994441    | 38994492 | 0.9   | 0.001 | 0.001 | 0.9   | 0.001 | 0.001 | 4
chr17  | 55299940    | 55299992 | 0.95  | 0.01  | 0.01  | 0.9   | 0.001 | 0.001 | 4

**NOTE: The header of the table is only for the explanation. The stutter file must not have a header.**


The stutter parameters have these definitions:

| VARIABLE | DESCRIPTION
| -------- | --------
| IDOWN    | Probability that an in-frame change decreases the size of the observed STR allele
| IUP      | Probability that an in-frame change increases the size of the observed STR allele
| ODOWN    | Probability that an out-of-frame change decreases the size of the observed STR allele
| OUP      | Probability that an out-of-frame change increases the size of the observed STR allele
| IGEOM    | Parameter of the geometric distribution of the step size for the in-frame changes
| OGEOM    | Parameter of the geometric distribution of the step size for the out-of-frame changes
| PERIOD   | Length of the STR motif

## FAQ
1. **Can I use HipSTR-MT if my data contains only single-end reads?**     
**Yes.** HipSTR-MT is designed for paired-end reads. It uses the data of the mate to filter the reads that are possibly aligned to an incorrect STR before the genotyping. Thus, by default, HipSTR-MT removes all the reads that have no mate. But if your data contains only single-end reads, give the **use-unpaired** option to prevent this filtration.  
2. **Can I use HipSTR-MT to analyze PCR-amplified reads?**     
**Yes.** HipSTR-MT is designed to analyze WGS data. Thus it finds the PCR duplicates and removes them before the genotyping. When you analyze PCR-amplified reads, HipSTR-MT identifies most of the reads as PCR duplicates, because they have exactly the same coordinates. To prevent this problem, give the **no-rmdup** option, which disables the removal of the duplicates for this type of data.
3. **Why are some of the STRs in my BED file not in the output VCF file?**    
HipSTR-MT genotypes a region only if the number of reads that overlap the STR is not less than **min-reads** and not more than **max-reads**. It then tries to learn the stutter model (if this is applicable), to make the haplotypes of the region, and to do the genotyping. If one of these steps fails, the program does not genotype the STR and continues with the subsequent region. The **log** file gives the reason for each failed region, and a summary of all the skipped regions at the end.      
4. **How can I use HipSTR-MT if I have too few samples to learn the stutter models, and I have no external models?**     
In this condition, run HipSTR-MT with **def-stutter-model**. This option disables the algorithm that learns the stutter models. HipSTR-MT then uses the same fixed stutter model for each locus. Use this option only when it is necessary, because the genotypes are more accurate with a specific model for each STR.
5. **Which sequencing platforms does HipSTR-MT support?**		
HipSTR-MT is designed to analyze **Illumina** sequencing data. Do not use it with PacBio data or Oxford Nanopore data, because their different error profiles cause problems.

## Help
If you have a problem with your analysis:      

    i.   Read the HipSTR-MT tutorial at https://hipstr-tool.github.io/HipSTR-tutorial
    ii.  Use ./HipSTR-MT --help for the data about each command line option
    iii. Send an email to hipstrtool@gmail.com

If you find a bug, or if you have a request for a function that is specific to this fork (the parallelization, the build, the performance, and so on):

     i.  Make an issue on GitHub (https://github.com/TurakhiaLab/HipSTR-MT)

If you have a question about the genotyping model or the algorithm, which this fork shares with the initial HipSTR:

     i.  Make an issue on GitHub (https://github.com/tfwillems/HipSTR)
    ii. Send an email to hipstrtool@gmail.com

## Citation
If HipSTR-MT was useful for your work, cite these publications:

- The initial HipSTR manuscript: **[Genome-wide profiling of heritable and de novo STR variations](https://www.nature.com/articles/nmeth.4267)**
- The software archive of this fork: Zenodo, [10.5281/zenodo.22650334](https://doi.org/10.5281/zenodo.22650334). This DOI always points to the most recent version. The Zenodo page of each release also gives a DOI for that one version.
- The Journal of Open Source Software (JOSS) paper of this fork, which gives the parallelization and the performance work. We will add the citation and the DOI here after the publication.

The [CITATION.cff](CITATION.cff) file gives the same data in a machine-readable format. The "Cite this repository" button on the GitHub page uses it.
