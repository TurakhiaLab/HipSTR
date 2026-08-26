#!/bin/bash
# Correctness regression check: builds HipSTR-MT (if needed), then runs it
# against the self-contained fixture in test/pgo/ (a real ~1Mb chr20 STR
# locus cluster, matching FASTA, and 2 subsetted sample BAMs -- the same
# fixture `make pgo` trains against) at several different --threads values,
# diffing every run's output against one committed golden VCF
# (test/pgo/expected_output.vcf.gz) using compare_vcf_tolerant.py.
#
# This is the primary correctness check for this fork, and specifically
# targets its core claim: that parallelizing HipSTR doesn't change what it
# genotypes. Diffing every thread count against the SAME golden VCF checks
# both things a reviewer would ask about in one pass -- (a) correctness vs.
# the unmodified upstream gymrek-lab/HipSTR tool the golden VCF was itself
# validated against (see README.md's Correctness section), and (b) that
# --threads 1/2/4/8 (whichever are requested) all agree with each other,
# since they all have to agree with the same golden file.
#
# Usage: test/check_correctness.sh [THREADS...]
#   (default thread counts: 1 2 4 8, capped to nproc if the machine has fewer)
# Exit status: 0 if every requested thread count's output matches the golden
# VCF within tolerance (no genotype changes, no >0.1% drift in derived
# statistics); nonzero otherwise. Suitable for CI.

set -uo pipefail
cd "$(dirname "$0")/.."

BIN=./HipSTR-MT
FIXTURE_DIR=test/pgo
GOLDEN="$FIXTURE_DIR/expected_output.vcf.gz"

if [ ! -x "$BIN" ]; then
    echo "[check_correctness] $BIN not found, building..."
    make -j"$(nproc)" HipSTR-MT || exit 1
fi

if [ "$#" -gt 0 ]; then
    THREAD_COUNTS=("$@")
else
    NPROC=$(nproc)
    THREAD_COUNTS=()
    for t in 1 2 4 8; do
        [ "$t" -le "$NPROC" ] && THREAD_COUNTS+=("$t")
    done
    # Always test at least --threads 1, even on a single-CPU runner.
    [ "${#THREAD_COUNTS[@]}" -eq 0 ] && THREAD_COUNTS=(1)
fi

OVERALL_STATUS=0
for T in "${THREAD_COUNTS[@]}"; do
    OUT="$FIXTURE_DIR/actual_output_t${T}.vcf.gz"
    LOG="$FIXTURE_DIR/actual_output_t${T}.log"

    echo "[check_correctness] --threads $T: running $BIN against the pgo fixture..."
    "$BIN" \
        --bams "$FIXTURE_DIR/bams/sampleA.pgo.bam,$FIXTURE_DIR/bams/sampleB.pgo.bam" \
        --fasta "$FIXTURE_DIR/pgo_fixture.fa" \
        --regions "$FIXTURE_DIR/pgo_fixture.bed" \
        --str-vcf "$OUT" \
        --threads "$T" --min-reads 10 \
        --log "$LOG"
    if [ $? -ne 0 ]; then
        echo "[check_correctness] --threads $T: FAIL -- HipSTR-MT itself exited non-zero (kept $OUT / $LOG for inspection)"
        OVERALL_STATUS=1
        continue
    fi

    echo "[check_correctness] --threads $T: comparing against golden output..."
    python3 test/compare_vcf_tolerant.py "$GOLDEN" "$OUT"
    STATUS=$?

    if [ "$STATUS" -eq 0 ]; then
        echo "[check_correctness] --threads $T: PASS"
        rm -f "$OUT" "$LOG"
    else
        echo "[check_correctness] --threads $T: FAIL -- see findings above (kept $OUT / $LOG for inspection)"
        OVERALL_STATUS=1
    fi
    echo
done

if [ "$OVERALL_STATUS" -eq 0 ]; then
    echo "[check_correctness] PASS across all thread counts tested: ${THREAD_COUNTS[*]}"
else
    echo "[check_correctness] FAIL"
fi
exit $OVERALL_STATUS
