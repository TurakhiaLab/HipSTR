# HipSTR-MT
**H**aplotype **i**nference and **p**hasing for **S**hort **T**andem **R**epeats  
![HipSTR icon!](https://raw.githubusercontent.com/tfwillems/HipSTR/master/img/HipSTR_icon_small.png)

#### Author: Thomas Willems <hipstrtool@gmail.com>

#### Optimizer: Joachim Galil <joachimbgalil@gmail.com>, Turakhia Lab -- this fork's performance work was done at the Gymrek Lab's request

#### License: GNU v2

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
Short tandem repeats [(STRs)](http://en.wikipedia.org/wiki/Microsatellite) are highly repetitive genomic sequences comprised of repeated copies of an underlying motif. Prevalent in most organisms' genomes, STRs are of particular interest because they mutate much more rapidly than most other genomic elements. As a result, they're extremely informative for genomic identification, ancestry inference and genealogy.

Despite their utility, STRs are particularly difficult to genotype. The repetitive sequence responsible for their high mutability also results in frequent alignment errors that can complicate and bias downstream analyses. In addition, PCR stutter errors often result in reads that contain additional or fewer repeat copies than the true underlying genotype. 

**HipSTR-MT** was specifically developed to deal with these errors in the hopes of obtaining more robust STR genotypes. In particular, it accomplishes this by:

1. Learning locus-specific PCR stutter models using an [EM algorithm](http://en.wikipedia.org/wiki/Expectation-maximization_algorithm)
2. Mining candidate STR alleles from population-scale sequencing data
3. Employing a specialized hidden Markov model to align reads to candidate alleles while accounting for STR artifacts
4. Utilizing phased SNP haplotypes to genotype and phase STRs

In our opinion, all of these factors make **HipSTR-MT** the most reliable tool for genotyping STRs from **Illumina** sequencing data.

## Requirements
HipSTR-MT keeps the original HipSTR runtime requirements and adds a modern C++ compiler for Taskflow:

- make
- g++ with C++20 support
- zlib
- libhts
- libbz2
- liblzma
- libcurl and OpenSSL (dev headers) — vendored htslib 1.24's default build config always compiles in libcurl-backed remote-file support
- CMake 3.18+ for building the vendored mimalloc and libdeflate (both CMake-only upstream builds); a portable CMake is downloaded automatically if the system one is too old or missing

On Ubuntu 16+ systems, the system packages can be installed with:

    apt install make g++ zlib1g-dev libhts-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev cmake

## Installation
Taskflow's headers and mimalloc's build-relevant source are vendored directly in this repo (the same way `lib/htslib` already is) rather than pulled in as git submodules, so a plain clone is all you need — no `--recurse-submodules`, no `git submodule update --init --recursive` to remember:

    git clone https://github.com/TurakhiaLab/HipSTR-MT

To build, use Make:

    cd HipSTR-MT
    make

The command constructs an executable file called **HipSTR-MT** in the current directory and builds mimalloc automatically as part of that. View detailed help with:

    ./HipSTR-MT --help

The Makefile now emits compiler dependency files with `-MMD -MP`, so header changes in `src`, `src/SeqAlignment`, and `src/denovos` trigger the required object rebuilds.

### Building with profile-guided optimization
    make pgo

`make pgo` compiles an instrumented `HipSTR-MT`, trains it against a bundled fixture (`test/pgo/`: a real ~1Mb chr20 STR locus cluster, 654 loci, with matching FASTA and 2 subsetted sample BAMs), then recompiles the final `HipSTR-MT` using the resulting profile. Training and both compiles happen locally, so the binary is tuned for whatever machine ran `make pgo` rather than shipping a profile baked in on other hardware. `DenovoFinder` is unaffected — it's rebuilt afterward with normal flags. This same fixture also backs the correctness regression check -- see [Testing](#testing).

**Not currently recommended**: measured on GCC 11.4 and 13.1, `make pgo` produces a binary 5-7% *slower* than plain `make`, from an interaction between `-fprofile-use` and `-flto=auto`. See the Makefile's PGO section for details. Use plain `make` until that interaction is resolved.

## Testing
    test/run_tests.sh

Everything `test/run_tests.sh` runs is self-contained -- it builds what it needs and uses only fixture data already checked into `test/` (no external downloads).

**Correctness regression** (`test/check_correctness.sh`, the primary check -- this is the one test a JOSS-style review will most want to see): runs `HipSTR-MT` against the same real ~1Mb chr20 fixture `make pgo` trains against (654 STR loci, 2 subsetted sample BAMs; see [Building with profile-guided optimization](#building-with-profile-guided-optimization)) at multiple `--threads` values (default: 1, 2, 4, 8, capped to the machine's core count), diffing *every* thread count's output against the same committed golden VCF (`test/pgo/expected_output.vcf.gz`) using `test/compare_vcf_tolerant.py`. That single golden-file comparison checks both things this fork's correctness claim rests on at once: that genotyping calls don't depend on thread count, and that they still match unmodified upstream `gymrek-lab/HipSTR`, which is what the golden VCF was itself validated against (see [Correctness](#hipstr-mt-changes) above, most recently confirmed with 0 discrepancies across the full 599-locus tutorial trio and a full-genome NA12891 run of 1,512,240 loci at `--threads` 1 through 64). `compare_vcf_tolerant.py` hard-fails on any genotype call change or >0.1% drift in derived statistics, but allows the sub-0.1% float noise that's inherent to floating-point summation order. Run explicit thread counts with `test/check_correctness.sh 1 16 64`.

**Unit tests** (`test/*_test.cpp`, each independently buildable via `make test/<name>`):
- **`snp_tree_test`** has a real pass/fail assertion: builds the same SNP set two ways (brute-force scan and the interval-tree structure `snp_bam_processor.cpp` actually uses) and asserts their query results agree.
- **`fast_ops_test`**, **`haplotype_test`**, **`read_vcf_alleles_test`** are diagnostic -- they print output (approximation error tables, generated haplotype sequences, parsed VCF alleles) for manual inspection rather than asserting. `read_vcf_alleles_test` needs the bundled `test/input/1kg.chr1.imputed.vcf.gz` fixture (included in the run above).
- **`em_stutter_test`** and **`vcf_snp_tree_test`** build but aren't included in the default run: `em_stutter_test`'s intended driver (`run_stutter_em_test.sh`) generates its input via an external STR-mutation simulator not bundled in this repo, and `vcf_snp_tree_test` expects a chr22:10-20Mb VCF not included here. Both take a VCF/BED path as an argument if you want to point them at your own data.

## Quick Start
To run HipSTR-MT in its most broadly applicable mode, run it on **all samples concurrently** using the syntax:

```
./HipSTR-MT --bams          run1.bam,run2.bam,run3.bam,run4.bam
         --fasta         genome.fa
         --regions       str_regions.bed
         --str-vcf       str_calls.vcf.gz
         --threads       8
```

* **bams** :  a comma-separated list of [BAM/CRAM](#bams) files generated by [BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml) and sorted and indexed using [samtools](http://www.htslib.org/)
* **regions** : a [BED](#str-bed) file containing the coordinates for each STR region of interest. Download BED files for various organisms, including humans, from [here](https://github.com/HipSTR-Tool/HipSTR-references/)
* **fasta** : [FASTA file](https://en.wikipedia.org/wiki/FASTA_format) containing the sequence for each chromosome in the BED file. This build's coordinates must match those of the STR regions
* **str-vcf** : The output path for the STR genotypes
* **threads** : Number of Taskflow executor worker threads. If omitted, HipSTR-MT chooses a hardware-aware default from scheduler CPU allocations, Linux CPU affinity, or `std::thread::hardware_concurrency()`. Internally, HipSTR-MT keeps four pipeline lines in flight per worker to hide serial fetch/write latency.

For each region in *str_regions.bed*, **HipSTR-MT** will:

1. Learn a stutter model for each locus
2. Use the stutter model and haplotype-based alignment algorithm to genotype each individual
3. Output the resulting STR genotypes to *str_calls.vcf.gz*, a [bgzipped](http://www.htslib.org/doc/tabix.html) [VCF](#str-vcf) file. This VCF will contain calls for each sample in any of the BAM/CRAM files' read groups.

## HipSTR-MT Changes
HipSTR-MT is a performance fork of [gymrek-lab/HipSTR](https://github.com/gymrek-lab/HipSTR) and that baseline is the comparison for all claims below.

**Correctness**: verified against the outputs of the Gymrek Lab's version of HipSTR using a tolerant VCF comparator (exact match required on genotype calls; float-typed fields allowed under 1e-3 relative drift). Genotype calls and all derived statistics now match exactly (0 drift) on both the 599-locus tutorial trio and a full-genome NA12891 run (1,512,240 loci) -- see [Vectorization](#vectorization) for a correctness bug this surfaced and fixed in the target_clones dispatch.

### Parallelization
- `bam_processor.*` replaces the single-region loop with a three-stage Taskflow pipeline: serial region token creation, parallel read filtering/genotyping, and serial ordered output. Regions are processed out of order across worker threads but written in the original BED order. Each pipeline line gets its own `BamCramMultiReader` and `AdapterTrimmer` instance, buffers its pass/filter BAM records instead of writing them inline, and all lines share one cached FASTA chromosome sequence rather than each copying it.
- `snp_bam_processor.*` moves SNP phasing preparation into the pipeline work item and adds two mutexes (`snp_stats_mutex_` for aggregate counters, `snp_phase_mutex_` for the shared reader/phasing state) so `--snp-vcf` runs are safe across worker threads.
- `genotyper_bam_processor.*` buffers per-region VCF, log, visualization, stutter, timing, and BAM output into a `RegionResult`, merges counters, and flushes everything in BED order from the serial output stage.
- `seq_stutter_genotyper.*` adds `build_vcf_record`/`build_vcf_records`, which render a locus's VCF line into a `BuiltVCFRecord` string instead of writing directly to a stream, allowing a worker thread to finish a region without needing to hold the output lock.
- `hipstr_main.cpp` adds `--threads <num_threads>`. If omitted, `default_thread_count()` picks a hardware-aware default: scheduler CPU allocation env vars first (`SLURM_CPUS_PER_TASK`, `SLURM_CPUS_ON_NODE`, `PBS_NP`, `NSLOTS`, `OMP_NUM_THREADS`), then Linux CPU affinity (`sched_getaffinity`), then `std::thread::hardware_concurrency()`. The pipeline keeps `4 * threads` region contexts in flight.

### Thread-safety fixes this required
Two pieces of the original single-threaded code held mutable state that's safe when there's exactly one caller but isn't once worker threads share it:
- **`StutterAlignerClass`**: its scratch buffers (`ins_probs_`, `del_probs_`, `match_probs_`, `log_probs_`) were instance members, reallocated on every `load_read()` call. These objects are owned by the haplotype/block structure and shared across threads (unlike `HapAligner`, which is per-thread), so concurrent `load_read()` calls would have raced on that shared state. Fixed by moving the buffers into a caller-supplied `StutterWorkspace` (one per `HapAligner`) — `StutterAlignerClass`'s methods are now `const` and touch no shared mutable state. This also added memoization: a repeated `load_read()` call with identical arguments (common when reusing alignments across candidate haplotypes) now skips recomputation instead of redoing it.
- **Cephes' `bdtr`** (binomial CDF, used by `compute_allele_bias`) keeps internal state that isn't safe to call concurrently. `seq_stutter_genotyper.cpp` now wraps that call in a `std::mutex`.

### Memory optimizations
- **`HapAligner`** reuses per-aligner scratch buffers (base-quality arrays, DP matrices, artifact size/position buffers) across reads instead of `new[]`/`delete[]` on every one. `HapAligner` instances aren't shared between threads. The two largest per-read allocations (the match/insert/deletion DP matrices, `O(read_len × haplotype_len)` each) are interleaved into one buffer (`MatrixChannel`, `[match0, insert0, deletion0, match1, ...]`) for improved cache locality.
- **ASCII-only case conversion** replaces locale-aware `toupper()`/`tolower()` in hot per-base loops: `stringops.cpp`'s `uppercase()`, `AlignmentOps.cpp`'s CIGAR-driven base comparison, `NeedlemanWunsch.cpp`'s `base_to_int()`, and `zalgorithm.cpp`/`alignment_filters.cpp`'s prefix/suffix/end-match comparisons.
- **`mathops.cpp`** adds a pointer-pair overload of `fast_log_sum_exp` (`const double* begin, const double* end`) alongside the original `vector<double>` one, avoiding a vector copy at a couple of call sites.
- **mimalloc** is linked in by default (see Installation) to cut allocator overhead from the volume of small per-read/per-locus allocations.
- **chromosome cache** is used to share chromosomes across threads. Since the program uses the chromosomes in order, when a chromosome is no longer in use due to all threads migrating to the next one, it is removed from the shared cache, reducing memory footprint.

### Vectorization
- **`-flto=auto`** enables link-time optimization across all translation units (~7-8% faster in benchmarking, identical output).
- **`mathops.cpp`'s `sum`/`log_sum_exp`/`fast_log_sum_exp`** are compiled with `__attribute__((target_clones("avx512f,avx2,sse4.2,default")))`, which builds one copy per listed ISA and dispatches to the best one the CPU supports at runtime. This is portable across machines (unlike `-march=native`, which isn't used anywhere in this build) and requires no user configuration.
- **`log_sum_exp`**'s `exp()` reduction is vectorized via glibc's libmvec (correctly-rounded, not an approximation) using `#pragma omp simd` and `-fopenmp-simd`, which pulls in no OpenMP runtime.
- **`fast_log_sum_exp`** batches 4 elements at a time using `vfasterexp()`, an SSE-vectorized helper already vendored in `fastonebigheader.h` but previously unused.
- All three `target_clones`'d functions carry `__attribute__((optimize("-ffp-contract=off")))`. Without it, the avx512f/avx2 clones let the compiler fuse multiply-adds that the default clone doesn't, so identical source could round differently depending purely on which ISA clone the CPU dispatches to at runtime -- confirmed by full-genome testing to flip a stutter-block candidate's log-probability at a couple of homopolymer STR loci, changing which alleles got discovered as candidates at one of them. Fixed without giving up per-machine ISA dispatch.

### Dependency and I/O fixes
- **htslib upgraded 1.9 → 1.24.** The vendored 1.9 copy's `fai_retrieve()` read FASTA sequence one byte at a time (`bgzf_getc()` plus a locale-aware `isgraph()` check per byte). Since chromosome loading runs in the pipeline's mandatory serial stage, this cost didn't shrink with more worker threads — on the tutorial dataset it was ~66% of the wall-clock floor at high thread counts. 1.24 pulls in upstream's already-fixed block-read implementation instead of a local patch: total FASTA load time across chr1–22 dropped from 7.17s to 1.36s, and wall time at `--threads 24` dropped from ~10.85s to ~5.5–6.7s. Picking up the newer vendored source needed two small C++-compatibility fixes: an explicit cast in `cram/cram_io.h` (implicit `void*` conversion is valid C, not C++) and a missing `<unistd.h>` include in `bam_io.h`/`denovo_main.cpp` for `access()`/`F_OK`, both previously masked by htslib 1.9's transitive includes.
- **libdeflate was linked but never active.** `HTSLIB_LIB`'s build rule was missing `-DHAVE_LIBDEFLATE`, so `bgzf.c`'s libdeflate code paths never compiled in and all BGZF/BAM decompression silently used system zlib instead. Fixed by adding the define.

### New / restored CLI flags
- **`--threads <num_threads>`** — see Parallelization above.
- **`--lib-from-samp`** — assign each read's library from its sample name instead of requiring an `LB` tag on every read group.
- **`--output-hap-fields`** — writes extra `LFLANKS`/`RFLANKS`/`HQ`/`PHQ`/`LFGT`/`RFGT` fields about the full assembled haplotypes,

## Tutorial
To demonstrate how you can quickly apply HipSTR-MT to whole-genome sequencing datasets, we've built a simple [tutorial](https://hipstr-tool.github.io/HipSTR-tutorial/). In less than 10 minutes, this tutorial will teach you how to genotype ~600 STRs in a deeply sequenced trio of individuals and inspect the results.

## In-depth Usage
**HipSTR-MT** has a variety of usage options designed to accomodate scenarios in which the sequencing data varies in terms of the number of samples and the coverage. Most scenarios will fall into one of the following categories:

1. 100 or more low-coverage (~5x) samples
    * Sufficient reads for stutter estimation
    * Sufficient reads to detect candidate STR alleles
    * [**Use de novo stutter estimation + STR calling with de novo allele generation**](#mode-1)
2. 20 or more high-coverage (~30x) samples
    * Sufficient reads for stutter estimation
    * Sufficient reads to detect candidate STR alleles
    * [**Use de novo stutter estimation + STR calling with de novo allele generation**](#mode-1)
3. Handful of low-coverage  (~5x) samples
    * Insufficient reads for stutter estimation
    * Insufficient reads to detect candidate STR alleles
    * [**Use external/default stutter models + STR calling with a reference panel**](#mode-3)
4. Handful of high-coverage (~30x) samples
    * Insufficient samples for stutter estimation
    * Sufficient reads to detect candidate STR alleles
    * [**Use external/default stutter models + STR calling with de novo allele generation**](#mode-2)

<a id="mode-1"></a>

#### Mode 1: De novo stutter estimation + STR calling with de novo allele generation
This mode is identical to the one suggested in the **Quick Start** section as it suits most applications. HipSTR-MT will output the STR genotypes in bgzipped VCF format to *str_calls.vcf.gz* 

```
./HipSTR-MT --bams             run1.bam,run2.bam,run3.bam,run4.bam
         --fasta            genome.fa
         --regions          str_regions.bed
         --str-vcf          str_calls.vcf.gz
```

<a id="mode-2"></a>

#### Mode 2: External stutter models + STR calling with de novo allele generation
The sole difference in this mode is that we no longer learn stutter models using the EM algorithm but instead input them from the **stutter-in** file. For more details on the stutter model file format, see [below](#stutter-file).

```
./HipSTR-MT --bams             run1.bam,run2.bam,run3.bam,run4.bam
         --fasta            genome.fa
         --regions          str_regions.bed
         --stutter-in       ext_stutter_models.txt
         --str-vcf          str_calls.vcf.gz
```
If you don't have access to external stutter models for the **stutter-in** option, use **def-stutter-model**. This will use a simplistic stutter model for all loci (see the HipSTR-MT help message for specifics).

<a id="mode-3"></a>

#### Mode 3: External stutter models + STR calling with a reference panel
This mode is very similar to mode 2, except that we provide an additional VCF file containing known STR genotypes at each locus using the **str-vcf** option. **HipSTR-MT** will not identify any additional candidate STR alleles in the BAMs/CRAMs when this option is specified, so it's best to use a VCF that contains STR genotypes for a wide range of populations and individuals. 

```
./HipSTR-MT --bams             run1.bam,run2.bam,run3.bam,run4.bam
         --fasta            genome.fa
         --regions          str_regions.bed
         --stutter-in       ext_stutter_models.txt
         --ref-vcf          ref_strs.vcf.gz
         --str-vcf          str_calls.vcf.gz
```

If you don't have access to external stutter models for the **stutter-in** option, use **def-stutter-model**. This will use a simplistic stutter model for all loci (see the HipSTR-MT help message for specifics).

## Data Requirements
To genotype STRs, **HipSTR-MT** requires Illumina sequencing data. However, as the depth of sequencing and the read length in these datasets can vary dramatically, here we briefly describe key factors to consider before generating data for HipSTR-MT analyses.

Because of the repetitive nature of STRs, reads that do not fully extend across the repeat only provide a lower bound on its length. While this lower bound is informative and is leveraged by **HipSTR-MT**, obtaining accurate and robust STR genotypes requires reads that fully extend across the repetitive sequence (*i.e. spanning reads*). The number of reads that span an STR is a function of the read length, the sequencing depth, and the length of the repeat (as well as various other factors). The interplay between these factors is relatively complex, but [**Figure 2** in a recent review](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4254273/figure/F2/) by *Press et al.* nicely highlights these dependencies. As one would intuitively expect, using longer read lengths and higher sequencing coverage increases the number of spanning reads. Conversely, increasing the length of the repeat reduces the number of spanning reads, making it more difficult to accurately genotype long STRs. When the number of spanning reads approaches single digits, you statistically run the risk of observing reads from only 1 out of 2 chromosome copies, making it impossible to correctly call both alleles in a heterozygous individual. 

In our own analyses, we've found that 100 bp Illumina reads are sufficient to characterize the majority of STRs in the human genome. However, genotyping STRs that exceed 70bp (such as very long forensic STRs) invariably requires longer reads, but most human STRs are much shorter than this threshold. This read length will likely be sufficient for most model organisms unless their repeats are substantially longer than those in humans.

The optimal minimum sequencing depth for HipSTR-MT largely depends on your intended analyses. If you are interested in studying how STRs mutate or [want to identify de novo mutations](#de-novo-mutations), 30x coverage is an ideal minimum that allows HipSTR-MT to provide the required high degree of specificity. Conversely, if you are merely interested in studying the allele frequencies for various STRs in a population, 10x coverage will likely be sufficient. However, in this setting, there will likely be many genotyping errors in which heterozygous genotypes are miscalled as homozygotes.

Based on the nature of your sequencing data, one important HipSTR-MT option to consider is **min-reads**. HipSTR-MT uses the value of this parameter to skip any STRs where few than *N* reads are available for genotyping across all individuals. By default, this value is 100, as we've found that this is a good minimum threshold for learning stutter models prior to genotyping. If you're analyzing very few samples (e.g. a single mother-father-child trio), you may want to consider lowering this threshold as you will seldom have 100 reads. In this setting, it makes sense to use options like **--min-reads 15 --def-stutter-model**, where the latter option uses a default stutter model as too few reads are available for accurately inferring one. However, if you're analyzing many samples (e.g. more than ten 30x genomes or more than thirty 10x genomes), it likely doesn't make sense to change this parameter. In these settings, regions with fewer than 100 reads may have high GC content that is problematic for Illumina sequencing, may be difficult to map to, or may merely be too long for your chosen read length.  

## Phasing
HipSTR-MT utilizes phased SNP haplotypes to phase the resulting STR genotypes. To do so, it looks for pairs of reads in which the STR-containing read or its mate pair overlap a samples's heterozygous SNP. In these instances, the quality score for the overlapping base can be used to determine the likelihood that the read came from each haplotype. Alternatively, when this information is not available, we assign the read an equal likelihood of coming from either strand. These likelihoods are incorporated into the HipSTR-MT genotyping model which outputs phased genotypes. The quality of a phasing is reflected in the *PQ* FORMAT field, which provides the posterior probability of each sample's phased genotype. For homozygous genotypes, this value will always equal the *Q* FORMAT field as phasing is irrelevant. However, for heterozygous genotypes, if *PQ ~ Q*, it indicates that one of the two phasings is much more favorable. Alterneatively, if none of a sample's reads overlap heterozygous SNPs, both phasings will be equally probable and *PQ ~ Q/2*. To enable the use of physical phasing, supply HipSTR-MT with the **snp-vcf** option and a SNP VCF containing **phased** haplotypes. The schematic below outlines the concepts underlying HipSTR-MT's physical phasing model:

![Phasing schematic!](https://raw.githubusercontent.com/tfwillems/HipSTR/master/img/phasing.png)

## Speed
On a full hg19 genotyping run (sample NA12891/ERR194160, 1,512,240 STR loci genome-wide, dual Xeon Silver 4216, 64 logical CPUs), HipSTR-MT at `--threads 64` finishes in 24.6 minutes versus 13.3 hours for unmodified upstream HipSTR run single-threaded -- a 32.5x speedup, genotype-for-genotype identical output (see [Correctness](#hipstr-mt-changes)). That number already includes the non-parallel optimizations below (LTO, SIMD, htslib, mimalloc): even at `--threads 1`, HipSTR-MT alone is measurably faster than the original in serial. Scaling is close to linear through 16 threads (94% efficiency) and tapers off by 64 (51% efficiency) as SMT contention and serial pipeline stages start to dominate -- see the thread-scaling plots this fork's benchmark harness produces for the full runtime/speedup/memory/CPU-utilization breakdown.

HipSTR-MT has built-in region-level multithreading. Use `--threads N` to set the number of Taskflow executor workers. If `--threads` is omitted, the executable selects a default from scheduler CPU allocation variables such as `SLURM_CPUS_PER_TASK`, then Linux CPU affinity, then `std::thread::hardware_concurrency()`. The pipeline keeps `4 * N` region contexts in flight so worker threads can continue genotyping while serial stages fetch the next region or flush completed output.

The highest-value internal optimizations in this fork are:

1. Independent regions are processed concurrently while VCF, log, BAM, visualization, and stutter-model output are still written in the original BED order.
2. Chromosome FASTA sequences are cached once per chromosome and shared across all pipeline lines, avoiding large per-line contig copies.
3. Haplotype-alignment DP matrices are reused inside each `HapAligner`, avoiding millions of repeated allocations in the read-alignment hot path.
4. mimalloc is linked by default to reduce allocator overhead that remains in read and haplotype processing.
5. htslib 1.24 replaces a byte-at-a-time FASTA reader that ran in the pipeline's serial stage with a block-read implementation, removing a bottleneck that had capped scaling at higher thread counts (see [Dependency and I/O fixes](#dependency-and-io-fixes)).
6. Link-time optimization and portable runtime CPU dispatch (see [Vectorization](#vectorization)) speed up the numeric hot path without requiring `-march=native` or any per-machine tuning.

For larger runs, start with `--threads` near the number of physical cores available to the job and benchmark a small representative region set.

Before `--threads` existed, the only way to parallelize the original single-threaded HipSTR was to manually split work across multiple OS processes -- `--threads N` replaces that within a single machine/process. The two options below are still useful, but now specifically for distributing work *across* separate machines/jobs (e.g. an HPC array), not as a substitute for `--threads` on one machine:

Option 1: Analyze each chromosome in parallel using the **--chrom** option. For example, **--chrom chr2** will only genotype BED regions on chr2

Option 2: Split your BED file into *N* files and analyze each of the *N* files in parallel. This allows you to parallelize analyses in a manner similar to option 1 but can be used for increased speed if *N* is much greater than the number of chromosomes.

## Default Filtering
HipSTR-MT sometimes automatically filters genotypes on a per-sample basis and will report a missing value in the VCF file. These filters are applied when a sample's data suggests that HipSTR-MT will not be able to produce a reliable genotype. For each locus, a summary of the number of filtered samples is output in the **log** file. If you specify the **--output-filters** command line option, a FORMAT field called **FILTER** will be reported in the VCF for each sample, where *PASS* designates ok samples and other values indicate the reason for filtering. 

**Samples with a PASS value should still undergo additional variant filtering (see below), as this merely indicates that no catastrophic issues were encountered during the genotyping process**. The table below summarizes the potential filtering reasons:  

| Filter | Explanation 
| :----- | :---------
| NO_READS                | No alignments were available for the sample at the current STR. If reads overlap the STR in the BAM/CRAM, they may have been filtered due to read quality issues, mapping uniqueness or other reasons
| FLANK_ASSEMBLY_CYCLIC   | During the genotyping process, HipSTR-MT attempts to assemble the sequences upstream and downstream of the STR (*flank*) to identify any potential SNPs it should consider. This assembly process fails if the resulting assembly graph contains a cycle, resulting in this filter
| FLANK_ASSEMBLY_INDEL    |This filter is triggered if the assembly process identifies an insertion or deletion in the *flanks*. These indels are problematic for HipSTR-MT's model and thus it does not attempt to genotype the sample
| FLANK_INDEL_FRAC        | When genotyping is complete, HipSTR-MT determines the maximum-likelihood alignment of each read relative to its sample's called alleles. If a large fraction of the resulting alignments have indels in the *flanks*, it's a strong indicator that they're misaligned and the sample's genotype is therefore ignored
| LOW_FREQUENCY_ALT_FLANK | Flanking sequences identified by the assembly process in each sample are pooled together to generate all candidate haplotypes. As the number of haplotypes grows exponentially with the number of such sequences, HipSTR-MT conserves time by discarding flanks that are only present in a few samples. If a sample's data supports a low-frequency flank, it is not genotyped. To adjust this frequency cutoff, use the **--min-flank-freq** option 



## Call Filtering
Although **HipSTR-MT** mitigates many of the most common sources of STR genotyping errors, it's still extremely important to filter the resulting VCFs to discard low quality calls. To facilitate this process, the VCF output contains various FORMAT and INFO fields that are usually indicators of problematic calls. The INFO fields indicate the aggregate data for a locus and, if certain flags are raised, may suggest that the entire locus should be discarded. In contrast, FORMAT fields are available on a per-sample basis for each locus and, if certain flags are raised, suggest that some samples' genotypes should be discarded. The list below includes some of these fields and how they can be informative. The [dumpSTR](https://trtools.readthedocs.io/en/stable/source/dumpSTR.html) utility in the [TRTools package](https://trtools.readthedocs.io/en/stable/) — actively maintained by the Gymrek lab, with utilities for filtering, merging, and computing statistics on VCFs from HipSTR-MT and other STR genotypers — can also be used to filter VCFs using most of the fields below, and is the currently recommended approach upstream.

#### INFO fields:  
1. **DP**: Reports the total depth/number of informative reads for all samples at the locus. The mean coverage per-sample can obtained by dividing this value by the number of samples with non-missing genotypes. In general, genotypes with a low mean coverage are unreliable because the reads may only have captured one of the two alleles if an individual is heterozygous.
2. **DSTUTTER**: Reports the total number of reads at a locus with what HipSTR-MT thinks is a stutter artifact. If the total fraction of reads with stutter (DSTUTTER/DP) is high, genotypes for a locus will be unreliable because the reads frequently don't reflect the true underlying genotype. A high fraction of stutter-containing reads can be caused by too much PCR amplification, a duplicated locus that is mapping to a single location in the genome, or a failure of HipSTR-MT to identify sufficient candidate alleles.  
3. **DFLANKINDEL**: Reports the total number of reads for which the maximum likelihood alignment contains an indel in the regions flanking the STR. A high fraction of reads with this artifact (DFLANKINDEL/DP) can be caused by an actual indel in a region neighboring the STR. However, it can also arise if HipSTR-MT fails to identify sufficient candidate alleles. When these alleles are very different in size from the candidate alleles or are non-unit multiples, they're frequently aligned as indels in the flanking sequences.

#### FORMAT fields:  
1. **Q**: Reports the posterior probability of the genotype. We've found that this is the best indicator of quality of an individual sample's genotype and almost always use it to filter calls.   
2. **DP**, **DSTUTTER** and **DFLANKINDEL**: Identical to the INFO field case, these fields are also available for each sample and can be used in the same way to identify problematic individual calls.  
3. **AB** and **FS**: Quantify the log10 p-value of the allele bias and the Fisher strand bias, respectively. Large negative values indicate that the degree of bias observed is very unlikely to occur by random chance. Both compare read counts between the sample's two haplotype *copies* (the two parental chromosomes, distinguished via phased SNPs), not between distinct STR allele values -- so they're only meaningful for diploid samples with phasing information, and can be nonzero even for a homozygous STR call if reads split unevenly between the two copies. In the case of **AB**, an unlikely split suggests the observed read counts per haplotype copy are inconsistent with the predicted genotype. In the case of **FS**, it indicates a non-random association between sequencing strand and which copy each read is assigned to, suggesting sequencing errors may be inflating one of the reported alleles.   

**So what thresholds do we suggest for each of these fields?** The answer really depends on the quality of the sequencing data, the ploidy of the chromosome and the downstream applications. As a starting point, dumpSTR options like `--hipstr-min-call-Q 0.9 --hipstr-max-call-flank-indel 0.15 --hipstr-max-call-stutter 0.15` are a reasonable default. Alternatively, this repo still bundles the original filtering scripts in the **scripts** subdirectory, which apply the same kind of thresholds directly without a TRTools dependency:

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

The resulting VCF, which is printed to the standard output stream, will omit calls on a sample-by-sample basis in which any of the following conditions are met: i) the posterior < 90%, ii) more than 15% of reads have a flank indel or iii) more than 15% of reads have a stutter artifact. For the diploid VCF, these filters will also remove genotypes with an allele bias or Fisher strand bias p-value less than 0.01 (10^-2). Calls for samples failing these criteria will be replaced with a "." missing symbol as per the VCF specification. For more filtering options, type either

```
python scripts/filter_vcf.py -h
python scripts/filter_haploid_vcf.py -h
```

## Additional Usage Options

| Option  | Description  
| :------- | :----------- 
| **viz-out**       aln_viz.gz     | Output a file of each locus' alignments for visualization with VizAln or [VizAlnPdf](#aln-viz) <br> **Why? You want to visualize or inspect the STR genotypes**
| **log**         log.txt               | Output the log information to the provided file (Default = Standard error)  
| **threads** num_threads                | Number of Taskflow executor worker threads (Default = auto) <br> **Why? You want to override HipSTR-MT's hardware-aware default.** This replaces the workaround under [Speed](#speed) of manually splitting your BED file and running several original-HipSTR processes in parallel -- `--threads` does the same thing internally, in one process, with output still written in original BED order.
| **haploid-chrs**  list_of_chroms      | Comma separated list of chromosomes to treat as haploid (Default = all diploid) <br> **Why? You're analyzing a haploid chromosome like chrY**  
| **no-rmdup**                            | Don't remove PCR duplicates. By default, they'll be removed <br> **Why? Your sequencing data  is for PCR-amplified regions**  
| **use-unpaired**                        | Use unpaired reads when genotyping (Default = False) <br> **Why? Your sequencing data only contains single-ended reads**  
| **snp-vcf**    phased_snps.vcf.gz     | Bgzipped input VCF file containing phased SNP genotypes for the samples to be genotyped. These SNPs will be used to physically phase STRs<br> **Why? You have available phased SNP genotypes**  
| **bam-samps**     list_of_read_groups | Comma separated list of samples in same order as BAM files. <br> Assign each read the sample corresponding to its file. By default, <br> each read must have an RG tag and and the sample is determined from the SM field <br> **Why? Your BAM file RG tags don't have an SM field**  
| **bam-libs**      list_of_read_groups | Comma separated list of libraries in same order as BAM files. <br> Assign each read the library corresponding to its file. By default, <br> each read must have an RG tag and and the library is determined from the LB field <br> NOTE: This option is required when --bam-samps has been specified <br> **Why? Your BAM file RG tags don't have an LB tag**  
| **lib-from-samp**                       | Assign each read the library corresponding to its sample name instead of requiring an LB tag on every read group <br> **Why? Your BAM file RG tags don't have an LB field, and per-sample library granularity is fine**  
| **def-stutter-model**                   | For each locus, use a stutter model with PGEOM=0.9 and UP=DOWN=0.05 for in-frame artifacts and PGEOM=0.9 and UP=DOWN=0.01 for out-of-frame artifacts <br> **Why? You have too few samples for stutter estimation and don't have stutter models**  
| **min-reads** num_reads                           | 	Minimum total reads required to genotype a locus (Default = 100) <br> **Why? Refer to the discussion [above](#data-requirements)**  
|**output-filters**                        | Write why individual calls were filtered to the VCF (Default = False)
| **output-hap-fields**                    | Write extra `LFLANKS`/`RFLANKS`/`HQ`/`PHQ`/`LFGT`/`RFGT` FORMAT fields describing each sample's full assembled haplotypes, not just the STR portion (Default = False) <br> **Why? You want to inspect or debug the flanking-sequence alignment, not just the STR genotype**


This list is comprised of the most useful and frequently used additional options, but is not all encompassing. For a complete list of options, please type

    ./HipSTR-MT --help

<a id="aln-viz"></a>

## Alignment Visualization
When deciphering and inspecting STR calls, it's extremely useful to visualize the supporting reads. HipSTR-MT facilitates this through the **viz-out** option, which writes a compressed file containing alignments for each call that can be readily visualized using the **VizAln** command included in HipSTR-MT's main directory. If you're interested in visualizing alignments, you first need to index the file using tabix. 
For example, if you ran HipSTR-MT with the option `--viz-out aln.viz.gz`, you should use the command

    tabix -p bed aln.viz.gz

to generate a [tabix](http://www.htslib.org/doc/tabix.html) index for the file so that we can rapidly extract alignments for a locus of interest. This command only needs to be run once after the file has been generated. 

You could then visualize the calls for sample *NA12878* at locus *chr1 3784267* using the command

    ./VizAln aln.viz.gz chr1 3784267 NA12878

This command will automatically open a rendering of the alignments in your browser and might look something like:
![Read more words!](https://raw.githubusercontent.com/HipSTR-Tool/HipSTR-tutorial/master/viz_NA12878.png)
The top bar represents the reference sequence and the red text indicates the name of the sample and its associated call at the locus. The remaining rows indicate the alignment for each read used in genotyping. In this particular example, 14 reads have an *8bp deletion* and 14 reads have a *4bp insertion*. HipSTR-MT therefore genotypes this sample as *-8 | 4*

If we wanted to inspect all calls for the same locus, we could  use the command 

    ./VizAln aln.viz.gz chr1 3784267

To facilitate rendering these images for publications, we've also created a similar script that converts
these alignments into a PDF. This script can only be applied to one sample at a time, but the image above
can be generated in a file alignments.pdf as follows:

    ./VizAlnPdf aln.viz.gz chr1 3784267 NA12878 alignments 1

NOTE: Because the **viz-out** file can become fairly large if you're genotyping thousands of loci or thousands of samples, in some scenarios it may be best to rerun HipSTR-MT using this option on the subset of loci which you wish to visualize.

## File Formats
<a id="bams"></a>

### BAM/CRAM files
HipSTR-MT requires [BAM/CRAM](https://samtools.github.io/hts-specs/SAMv1.pdf) files produced by any indel-sensitive aligner. These files must have been sorted by position using the `samtools sort` command and then indexed using `samtools index`. To associate a read with its sample of interest, HipSTR-MT uses read group information in the BAM/CRAM header lines. These *@*RG lines must contain an *ID* field, an *LB* field indicating the library and an *SM* field indicating the sample. For example, if a BAM/CRAM contained the following header line

    @RG     ID:RUN1 LB:ERR12345        SM:SAMPLE789

an alignment with the RG tag 

    RG:Z:RUN1

will be associated with sample *SAMPLE789* and library *ERR12345*. In this manner, HipSTR-MT can analyze BAMs/CRAMs containing more than one sample and/or more than one library and can handle cases in which a single sample's reads are spread across multiple files.

Alternatively, if your BAM/CRAM files lack *RG* information, you can use the **bam-samps** and **bam-libs** flags to specify the sample and library associated with each file. In this setting, however, a BAM/CRAM can only contain a single library and a single read group. For example, the command

```
./HipSTR-MT --bams             run1.bam,run2.bam,run3.bam,run4.cram
         --fasta            genome.fa
         --regions          str_regions.bed
         --str-vcf          str_calls.vcf.gz
         --bam-samps        SAMPLE1,SAMPLE1,SAMPLE2,SAMPLE3
         --bam-libs         LIB1,LIB2,LIB3,LIB4
```

essentially tells HipSTR-MT to associate all the reads in the first two BAMS with *SAMPLE1*, all the reads in the third file with *SAMPLE2* and all the reads in the last BAM with *SAMPLE3*.


HipSTR-MT can analyze both BAM and CRAM files simultaneously, so if your project contains a mixture of these two file types, HipSTR-MT will automatically perform CRAM decompression as necessary. **When analyzing CRAM files, please ensure that the file provided to --fasta is the same FASTA file used during CRAM generation**. Otherwise, CRAM decompression will likely fail and bizarre behavior may occur.

<a id="str-bed"></a>

### STR region BED file
The BED file containing each STR region of interest is a tab-delimited file comprised of 5 required columns and one optional column: 

1. The name of the chromosome on which the STR is located
2. The start position of the STR on its chromosome
3. The end position of the STR on its chromosome
4. The motif length (i.e. the number of bases in the repeat unit)
5. The number of copies of the repeat unit in the reference allele

The 6th column is optional and contains the name of the STR locus, which will be written to the ID column in the VCF. 
Below is an example file which contains 5 STR loci 

**NOTE: The table header is for descriptive purposes. The BED file should not have a header**

CHROM | START       | END         | MOTIF_LEN | NUM_COPIES | NAME
----  | ----        | ----        | ---       | ---        | ---
chr1  | 13784267    | 13784306    | 4         | 10         | GATA27E01
chr1  | 18789523    | 18789555    | 3         | 11         | ATA008
chr2  | 32079410    | 32079469    | 4         | 15         | AGAT117
chr17 | 38994441    | 38994492    | 4         | 12         | GATA25A04
chr17 | 55299940    | 55299992    | 4         | 13         | AAT245

We've provided various *BED* files containing STR loci for different organisms, including humans, [here](https://github.com/HipSTR-Tool/HipSTR-references/)

For other model organisms, we recommend that you modify the [framework](https://github.com/HipSTR-Tool/HipSTR-references/blob/master/mouse/mouse_reference.md)
we used to build the mouse BED file.

<a id="str-vcf"></a>

### VCF file
For more information on the VCF file format, please see the [VCF spec](http://samtools.github.io/hts-specs/VCFv4.2.pdf). 

#### INFO fields
INFO fields contains aggregated statistics about each genotyped STR in the VCF. The INFO fields reported by HipSTR-MT primarily describe the learned/supplied stutter model for the locus, the STR's reference coordinates (START and END) and information about the allele counts (AC) and number of reads used to genotype all samples (DP).

FIELD | DESCRIPTION
----- | -----------
INFRAME_PGEOM  | Parameter for in-frame geometric step size distribution
INFRAME_UP     | Probability that stutter causes an in-frame increase in obs. STR size
INFRAME_DOWN   | Probability that stutter causes an in-frame decrease in obs. STR size
OUTFRAME_PGEOM | Parameter for out-of-frame geometric step size distribution
OUTFRAME_UP    | Probability that stutter causes an out-of-frame increase in obs. STR size
OUTFRAME_DOWN  | Probability that stutter causes an out-of-frame decrease in obs. STR size
BPDIFFS        | Base pair difference of each alternate allele from the reference allele
START          | Inclusive start coodinate for the repetitive portion of the reference allele
END            | Inclusive end coordinate for the repetitive portion of the reference allele
PERIOD         | Length of STR motif
AN             | Total number of alleles in called genotypes
REFAC          | Reference allele count
AC             | Alternate allele counts
NSKIP          | Number of samples not genotyped due to various issues
NFILT          | Number of samples that were originally genotyped but have since been filtered
DP             | Total number of reads used to genotype all samples
DSNP           | Total number of reads with SNP information
DSTUTTER       | Total number of reads with a stutter indel in the STR region
DFLANKINDEL    | Total number of reads with an indel in the regions flanking the STR

#### FORMAT fields
FORMAT fields contain information about the genotype for each sample at the locus. In addition to the most probable phased genotype (GT), HipSTR-MT reports information about the posterior likelihood of this genotype (PQ) and its unphased analog (Q). Other useful information reported are the number of reads that were used to determine the genotype (DP) and whether these had any alignment artifacts (DSTUTTER and DFLANKINDEL).

FIELD     | DESCRIPTION
--------- | -----------
GT        | Genotype
GB        | Base pair differences of genotype from reference
Q         | Posterior probability of unphased genotype
PQ        | Posterior probability of phased genotype
DP        | Number of valid reads used for sample's genotype
DSNP      | Number of reads with SNP phasing information
PDP       | Fractional reads supporting each haploid genotype
GLDIFF    | Difference in likelihood between the reported and next best genotypes
DSNP      | Total number of reads with SNP information
PSNP      | Number of reads with SNPs supporting each haploid genotype
DSTUTTER  | Number of reads with a stutter indel in the STR region
DFLANKINDEL | Number of reads with an indel in the regions flanking the STR
AB        | log10 of the allele bias pvalue, where 0 is no bias and more negative values are increasingly biased. This compares read counts between the sample's two haplotype *copies* (chromosomes), not between STR allele values -- so it can still be negative for a homozygous STR genotype (same allele on both copies) if nearby phased SNPs distinguish the two copies and reads split unevenly between them
FS        | log10 of the strand bias pvalue from Fisher's exact test, where 0 is no bias and more negative values are increasingly biased. Same per-haplotype-copy caveat as AB above: can be negative for a homozygous STR genotype if the two phased copies' reads have an uneven strand split
DAB       | Number of reads used in the allele bias calculation
ALLREADS  | Base pair difference observed in each read's Needleman-Wunsch alignment
MALLREADS | Maximum likelihood bp diff in each read based on haplotype alignments
GL        | log-10 genotype likelihoods
PL        | Phred-scaled genotype likelihoods

<a id="stutter-file"></a>

### Stutter model
To model PCR stutter artifacts, we assume that there are three types of stutter artifacts:

1. **In-frame changes**: Change the size of the STR in the read by multiples of the repeat unit. For instance, if the repeat motif is AGAT, in-frame changes could lead to differences of -8, -4, 4, 8, and so on. 
2. **Out-of-frame changes**: Change the size of the STR by non-multiples of the repeat unit. For instance, if the repeat motif is AGAT, out-of-frame changes could lead to differences of -3, -2, -1, 1, 2, 3 and so on. 
3. **No stutter change**: The size of the STR in the read is the same as the size of the underlying STR. 


Stutter model files contain the information necessary to model each of these artifacts in a **tab-delimited BED-like** format with exactly 10 columns (all required -- `StutterModel::read()` in `src/stutter_model.cpp` fails to parse the file if any are missing, including `PERIOD`). An example of such a file is as follows:

CHROM  | START       | END      | IGEOM | IDOWN | IUP   | OGEOM | ODOWN | OUP   | PERIOD
-----  | ----------- | -------- | ----  | ----  | ---   | ----  | ---   | ---   | ---
chr1   | 13784267    | 13784306 | 0.95  | 0.05  | 0.01  | 0.9   | 0.01  | 0.001 | 4
chr1   | 18789523    | 18789555 | 0.8   | 0.01  | 0.05  | 0.9   | 0.001 | 0.001 | 3
chr2   | 32079410    | 32079469 | 0.9   | 0.01  | 0.01  | 0.9   | 0.001 | 0.001 | 4
chr17  | 38994441    | 38994492 | 0.9   | 0.001 | 0.001 | 0.9   | 0.001 | 0.001 | 4
chr17  | 55299940    | 55299992 | 0.95  | 0.01  | 0.01  | 0.9   | 0.001 | 0.001 | 4

**NOTE: The table header is for descriptive purposes. The stutter file should not have a header**


Each of the stutter parameters is defined as follows:

| VARIABLE | DESCRIPTION
| -------- | --------
| IDOWN    | Probability that in-frame changes decrease the size of the observed STR allele
| IUP      |  Probability that in-frame changes increase the size of the observed STR allele
| ODOWN    | Probability that out-of-frame changes decrease the size of the observed STR allele
| OUP      |  Probability that out-of-frame changes increase the size of the observed STR allele
| IGEOM    | Parameter governing geometric step size distribution for in-frame changes
| OGEOM    | Paramter  governing geometric step size distribution for out-of-frame changes
| PERIOD   | Length of STR motif

## FAQ
1. **Can I run HipSTR-MT if my dataset only contains single-ended reads?**     
**Yes.** HipSTR-MT is designed for paired-end reads and uses mate pair information to filter reads that are potentially aligned to an incorrect STR prior to genotyping. By default, HipSTR-MT therefore removes all reads without mate pairs. However, if your dataset only contains single-ended reads, specify the **use-unpaired** option to avoid performing this filtering.  
2. **Can I use HipSTR-MT to analyze PCR-amplified reads?**     
**Yes.** As HipSTR-MT was designed to analyze WGS data, HipSTR-MT automatically identifies and filters out PCR duplicates prior to genotyping. When analyzing PCR-amplified reads, HipSTR-MT will label most reads as PCR duplicates as they share exactly the same coordinates. To overcome this issue, specify the **no-rmdup** option to disable duplicate removal when analyzing this type of data.
3. **Why are some of the STRs in my BED file not present in the output VCF?**    
HipSTR-MT only genotypes a region if at least **min-reads** and at most **max-reads** overlap the STR. It then attempts to learn the stutter model (if appropriate), build haplotypes for the region and perform genotyping. If any of these stages is unsuccessful, it skips the STR and continues on to the next region. The **log** file contains the failure reason for each failed region as well as an overall summary of why regions were skipped at the end of the log.      
4. **How can I run HipSTR-MT if I have too few samples to learn stutter models and don't have external ones?**     
In this scenario, you can run HipSTR-MT with **def-stutter-model**. Invoking this option will disable the algorithm it uses to learn stutter models. Instead, HipSTR-MT will use the same fixed stutter model to genotype every locus. We don't recommend using this option unless necessary, as genotypes are more accurate if you learn a specific model for each STR.
5. **What sequencing platforms does HipSTR-MT support?**		
HipSTR-MT was designed to analyze **Illumina** sequencing data. We do not recommend running it on PacBio or Oxford Nanopore data, as the difference in error profiles will be problematic 

## Help
If you're having trouble getting your analysis up and running:      

    i.   Check out the HipSTR-MT tutorial at https://hipstr-tool.github.io/HipSTR-tutorial
    ii.  Type ./HipSTR-MT --help for details about each command line option
    iii. Email us at hipstrtool@gmail.com

If you encounter a bug/issue or have a feature request specific to this fork (parallelization, build, performance, etc.):

     i.  File an issue on GitHub (https://github.com/TurakhiaLab/HipSTR-MT)

For questions about the underlying genotyping model/algorithm itself, shared with the original HipSTR:

     i.  File an issue on GitHub (https://github.com/tfwillems/HipSTR)
    ii. Email us at hipstrtool@gmail.com

## Citation
If you found HipSTR-MT useful, please cite both:

- The original HipSTR manuscript: **[Genome-wide profiling of heritable and de novo STR variations](https://www.nature.com/articles/nmeth.4267)**
- This fork's Journal of Open Source Software (JOSS) paper, describing the parallelization and performance work: citation and DOI to be added here once it's published. Until then, please cite this repository directly: https://github.com/TurakhiaLab/HipSTR-MT
