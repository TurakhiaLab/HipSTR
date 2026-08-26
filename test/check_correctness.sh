#!/bin/bash
# Correctness regression check: builds HipSTR-MT (if needed), runs it against
# the self-contained fixture in test/pgo/ (a real ~1Mb chr20 STR locus
# cluster, matching FASTA, and 2 subsetted sample BAMs -- the same fixture
# `make pgo` trains against), and diffs the output against a committed
# golden VCF (test/pgo/expected_output.vcf.gz) using compare_vcf_tolerant.py.
#
# This is the primary correctness check for this fork: it validates that the
# whole read -> haplotype-alignment -> genotyping pipeline still produces the
# same calls, not just that individual functions behave in isolation. The
# golden VCF was itself validated against the unmodified upstream
# gymrek-lab/HipSTR tool -- see README.md's "Correctness" section -- so this
# script is a fast, no-external-data-required proxy for re-running that full
# comparison after a change.
#
# Usage: test/check_correctness.sh
# Exit status: 0 if the current build's output matches the golden VCF within
# tolerance (no genotype changes, no >0.1% drift in derived statistics);
# nonzero otherwise. Suitable for CI.

set -euo pipefail
cd "$(dirname "$0")/.."

BIN=./HipSTR-MT
FIXTURE_DIR=test/pgo
GOLDEN="$FIXTURE_DIR/expected_output.vcf.gz"
OUT="$FIXTURE_DIR/actual_output.vcf.gz"

if [ ! -x "$BIN" ]; then
    echo "[check_correctness] $BIN not found, building..."
    make -j"$(nproc)" HipSTR-MT
fi

echo "[check_correctness] Running $BIN against the pgo fixture..."
"$BIN" \
    --bams "$FIXTURE_DIR/bams/sampleA.pgo.bam,$FIXTURE_DIR/bams/sampleB.pgo.bam" \
    --fasta "$FIXTURE_DIR/pgo_fixture.fa" \
    --regions "$FIXTURE_DIR/pgo_fixture.bed" \
    --str-vcf "$OUT" \
    --threads 1 --min-reads 10 \
    --log "$FIXTURE_DIR/actual_output.log"

echo "[check_correctness] Comparing against golden output..."
set +e
python3 test/compare_vcf_tolerant.py "$GOLDEN" "$OUT"
STATUS=$?
set -e

rm -f "$OUT" "$FIXTURE_DIR/actual_output.log"

if [ "$STATUS" -eq 0 ]; then
    echo "[check_correctness] PASS"
else
    echo "[check_correctness] FAIL -- see findings above"
fi
exit $STATUS
