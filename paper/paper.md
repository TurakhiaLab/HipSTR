---
title: 'HipSTR-MT: A Multi-threaded Extension of HipSTR for Faster Short Tandem Repeat Genotyping'
tags:
  - C++
  - bioinformatics
  - short tandem repeats
  - STR genotyping
  - population genetics
  - forensic genomics
  - whole-genome sequencing
  - multithreading
  - task-parallel pipeline
  - Taskflow
  - high-performance computing
authors:
  - name: Joachim Galil
    orcid: 0009-0000-2060-2074
    affiliation: 1
  - name: Melissa Gymrek
    orcid: 0000-0002-6086-3903
    affiliation: 2
  - name: Yatish Turakhia
    orcid: 0000-0001-5600-2900
    corresponding: true
    email: yturakhia@ucsd.edu
    affiliation: 1
affiliations:
  - name: Department of Electrical and Computer Engineering, University of California San Diego, San Diego, CA 92093, USA
    index: 1
    ror: 0168r3w48
  - name: Department of Medicine, University of California San Diego, San Diego, CA 92093, USA
    index: 2
    ror: 0168r3w48
date: 7 September 2026
bibliography: paper.bib
---

# Summary

Short tandem repeats (STRs) are highly variable, repetitive regions of the genome that serve as important markers in population genetics, forensics, and disease research [@budowle2024; @uguen2024]. HipSTR [@willems2017] is a widely used tool for genotyping STRs from whole-genome sequencing data, using a specialized hidden Markov model to align reads to candidate alleles while accounting for PCR stutter artifacts and phased SNP haplotypes. However, its serial execution model limits throughput on modern multi-core hardware, particularly in large cohort studies genotyping hundreds of thousands of loci across many samples. HipSTR-MT restructures HipSTR's core loop into a three-stage concurrent pipeline built on the Taskflow task-graph framework [@huang2022], preserving HipSTR's full feature set while producing genotype-identical output. On the NA12891 sample dataset, HipSTR-MT achieves a 32.5× reduction in wall-clock runtime at 64 threads relative to the serial baseline, letting researchers use multi-core workstations and HPC nodes without changing existing HipSTR workflows.

# Statement of need

Genome-wide STR genotyping is computationally intensive: for each locus, HipSTR seeks and filters BAM/CRAM records, constructs candidate haplotypes via a de Bruijn graph, aligns every read against every candidate haplotype under a stutter emission model, and estimates genotype posteriors—all on a single thread—making genome-wide analysis of large cohorts impractical on one machine. Cluster-based parallelism (running separate samples as independent jobs) is a common workaround, but it does not reduce the per-sample latency of genotyping a single large sample or a targeted large-region analysis, and requires infrastructure not all users have. HipSTR-MT closes this gap with intra-run thread-level parallelism, giving users of standard multi-core workstations a drop-in performance improvement via a single added `--threads` flag. It is intended for population geneticists, forensic genomics researchers, and any group running HipSTR as part of a larger cohort-scale variant-calling pipeline.

# State of the Field

HipSTR is one of several tools for genome-wide STR genotyping from short-read data, alongside GangSTR and ExpansionHunter, which additionally use mate-pair distance to detect STR expansions longer than a single read [@oketch2024]. Independent benchmarking shows GangSTR and ExpansionHunter can outperform HipSTR on call rate and memory usage, but HipSTR remains competitive in genotyping accuracy for common STR loci, including the CODIS core STRs used in forensics [@oketch2024]. None of these tools natively parallelizes within a single sample; throughput on large cohorts is instead typically achieved through cluster-based job scheduling across many machines [@ahmad2022; @mirchandani2024], which increases aggregate throughput but not per-sample latency, and assumes multi-node HPC access. HipSTR-MT takes a complementary approach: intra-run thread-level parallelism so a single workstation can genotype a sample or region substantially faster, while leaving HipSTR's statistical model and output format unchanged.

![Overview of the HipSTR-MT genotyping pipeline. A serial dispatch stage (1) creates one token per genomic region and loads the needed chromosome. Regions are distributed across parallel worker lines (2), each performing BAM/CRAM seeking, read filtering, SNP-based phasing, haplotype construction, and HMM-based realignment. A serial output stage (3) collects results and writes VCF, log, and BAM output in the original region order regardless of worker completion order, guaranteeing output identical to the serial tool.\label{fig:pipeline}](figures/figure1_pipeline.png)

# Software Design

The unit of work is the genomic region: regions are independent under HipSTR's model and map cleanly onto the existing per-locus loop. Rather than a parallel loop over regions, HipSTR-MT uses the three-stage Taskflow pipeline in \autoref{fig:pipeline}: serial dispatch, parallel workers, serial ordered output. The serial stages bound achievable scaling, but they let the tool emit VCF, log, and BAM records in the original region order regardless of which worker finishes first, keeping output directly comparable with the serial tool.

Treating that comparability as a hard constraint drove two trade-offs. Multiply-add contraction is disabled (`-ffp-contract=off`), giving up vectorization headroom to keep results bit-identical across instruction sets, since silent numerical divergence would make the fork unusable as a drop-in replacement. Each worker also keeps four region contexts in flight, with independent reader and alignment state, so the work-stealing scheduler can hide I/O and memory latency; this raises peak memory from ~1.6 GB to ~8 GB at 64 threads, favoring wall-clock time on machines where cores are scarcer than RAM.

Thread safety required eliminating two pieces of shared mutable state: `StutterAlignerClass`'s scratch buffers, moved into a per-`HapAligner` workspace, and the non-reentrant Cephes `bdtr` function, now mutex-guarded. Further optimizations include a shared chromosome cache, mimalloc [@leijen2019], an htslib upgrade, and SIMD dispatch via compiler target clones — improving single-threaded performance as well. Three additive flags extend the CLI without breaking backward compatibility.

+-----------------+-----------------------------+------------------------------------------------------------+
| Category        | Change                      | Effect                                                     |
+=================+=============================+============================================================+
| Parallelization | Taskflow pipeline           | Reads, filters, and writes regions across three concurrent |
|                 |                             | stages.                                                    |
|                 +-----------------------------+------------------------------------------------------------+
|                 | SNP-phasing safety          | Mutexes protect shared phasing state across worker         |
|                 |                             | threads.                                                   |
|                 +-----------------------------+------------------------------------------------------------+
|                 | Output buffering            | Per-thread output buffers are flushed in original region   |
|                 |                             | order.                                                     |
|                 +-----------------------------+------------------------------------------------------------+
|                 | Lock-free VCF records       | Workers render VCF lines as text without holding the       |
|                 |                             | output lock.                                               |
|                 +-----------------------------+------------------------------------------------------------+
|                 | `--threads` flag            | Sets worker thread count; auto-detected from hardware if   |
|                 |                             | unset.                                                     |
+-----------------+-----------------------------+------------------------------------------------------------+
| Thread safety   | StutterAligner buffers      | Per-`HapAligner` scratch buffers eliminate a shared-state  |
|                 |                             | race.                                                      |
|                 +-----------------------------+------------------------------------------------------------+
|                 | Cephes `bdtr` mutex         | Guards a non-reentrant function used in allele-bias        |
|                 |                             | computation.                                               |
+-----------------+-----------------------------+------------------------------------------------------------+
| Memory          | HapAligner buffer reuse     | Reuses per-aligner scratch buffers across reads instead of |
|                 |                             | reallocating.                                              |
|                 +-----------------------------+------------------------------------------------------------+
|                 | ASCII case conversion       | Replaces locale-aware `toupper`/`tolower` in hot per-base  |
|                 |                             | loops.                                                     |
|                 +-----------------------------+------------------------------------------------------------+
|                 | `fast_log_sum_exp` overload | Pointer-pair variant avoids a vector copy at call sites.   |
|                 +-----------------------------+------------------------------------------------------------+
|                 | mimalloc                    | Reduces allocator overhead from small per-read/per-locus   |
|                 |                             | allocations.                                               |
|                 +-----------------------------+------------------------------------------------------------+
|                 | Chromosome cache + eviction | Chromosomes shared across workers instead of duplicating   |
|                 |                             | per thread. Evicted once there are no workers using it.    |
+-----------------+-----------------------------+------------------------------------------------------------+
| Vectorization   | `-flto=auto`                | Enables link-time optimization across translation units.   |
|                 +-----------------------------+------------------------------------------------------------+
|                 | `target_clones` dispatch    | Builds one code variant per ISA; dispatches to the best at |
|                 |                             | runtime.                                                   |
|                 +-----------------------------+------------------------------------------------------------+
|                 | libmvec `exp()`             | Vectorizes the exponential reduction with                  |
|                 |                             | correctly-rounded results.                                 |
|                 +-----------------------------+------------------------------------------------------------+
|                 | SSE batching                | Processes four elements per instruction via                |
|                 |                             | `vfasterexp()`.                                            |
|                 +-----------------------------+------------------------------------------------------------+
|                 | `-ffp-contract=off`         | Disables multiply-add fusion to keep output bit-identical  |
|                 |                             | across ISAs.                                               |
+-----------------+-----------------------------+------------------------------------------------------------+
| Dependency/IO   | htslib 1.9 -> 1.24          | Replaces byte-at-a-time FASTA reads with a block-read      |
|                 |                             | implementation.                                            |
|                 +-----------------------------+------------------------------------------------------------+
|                 | libdeflate enabled          | Activates a previously unused BGZF decompression path.     |
+-----------------+-----------------------------+------------------------------------------------------------+
| CLI flags       | `--lib-from-samp`           | Assigns library name from sample name when LB tags are     |
|                 |                             | absent.                                                    |
|                 +-----------------------------+------------------------------------------------------------+
|                 | `--output-hap-fields`       | Adds extra FORMAT fields describing full assembled         |
|                 |                             | haplotypes.                                                |
+-----------------+-----------------------------+------------------------------------------------------------+

: Summary of code-level changes introduced in HipSTR-MT, grouped by category (parallelization, thread safety, memory, vectorization, dependency/I/O, and CLI flags), with the effect of each change on behavior or performance. Full implementation detail for each entry is provided in the repository README.

# Performance and correctness

![Performance results from processing the full NA12891 genome for all 1.5M+ STRs. a) Runtime vs. thread count. b) Parallel speedup vs. thread count, with data points annotated with scaling efficiency. c) Peak memory vs. thread count. d) CPU utilization vs. thread count.\label{fig:performance}](figures/figure2_performance.png)

Benchmarks used the NA12891 sample (accession ERR194160) against the genome-wide hg19 STR panel (1,512,240 loci):

```bash
./HipSTR-MT --bams ERR194160.bam --fasta hg19.fa --regions hg19_all_strs.bed \
  --str-vcf output.vcf.gz --min-reads 25 --def-stutter-model --threads N
```

with $N \in \{1, 2, 4, 8, 16, 32, 64\}$, measured via `/usr/bin/time -v` on dual Intel Xeon Silver 4216 CPUs (64 logical CPUs), Ubuntu 22.04. As \autoref{fig:performance} shows, HipSTR-MT at 64 threads finishes in 24.6 minutes versus 13.3 hours for the unmodified serial baseline (a 32.5× reduction), already including the non-parallel optimizations above, so even single-threaded HipSTR-MT is measurably faster than upstream. Scaling is near-linear through 16 threads (94% efficiency), tapering to 51% by 64 threads as SMT contention and the pipeline's mandatory serial stages dominate; per-locus haplotype alignment and traceback work is small, memory-irregular, and imbalanced across loci, leaving some workers idle.

Output equivalence was verified at every thread count with a tolerant VCF comparator requiring exact genotype-call matches and allowing $\leq 10^{-3}$ relative drift in floating-point fields (GLDIFF, PDP) from SIMD/codegen-dependent summation order. All genotype records were identical across thread counts. Maximum RSS grows with in-flight regions: ~8 GB at 64 threads versus ~1.6 GB at 1 thread.

# Research Impact

HipSTR-MT is a recent release, so its significance rests on the established user base it serves together with the performance and equivalence results above. The original HipSTR has been cited in over 330 publications since 2017 [@willems2017] and underpins large-scale STR resources. It profiled autosomal STRs in 1,916 individuals from 479 Simons Simplex Collection family quads---a mean of 1.14 million STRs per sample---for the first genome-wide SNP + STR imputation reference panel [@saini2018], and it is one of four genotypers combined by EnsembleTR into a catalog of more than 1.7 million tandem repeat loci across 3,550 individuals from the 1000 Genomes Project and H3Africa cohorts [@ziaeijam2023]. It has also been repackaged by an independent group as a graphical front end for routine forensic casework [@frontanilla2026]. Because HipSTR-MT produces genotype-identical output and preserves the command-line interface apart from three additive flags, these workflows can adopt it without re-validating existing call sets.

# Acknowledgements

We thank Gymrek and Turakhia lab members at UCSD for feedback and for providing the motivation and test data for this project. HipSTR-MT is a performance fork of HipSTR, originally developed by Thomas Willems and maintained by the Gymrek lab (gymrek-lab/HipSTR). Our fork is available at <https://github.com/TurakhiaLab/HipSTR-MT> under GPL v2, and we track upstream developments and intend to maintain it as an ongoing open-source project.

# AI usage disclosure

Claude Code, Codex, and DeepSeek were used to help understand the original HipSTR codebase, write test benches and data-collection/debugging scripts, and assist with some code comments and small memory refactors. All resulting code and results were reviewed, edited, and validated by the human authors, who made all core design decisions and take full responsibility for the accuracy, originality, and ethical compliance of the submitted materials.

# References
