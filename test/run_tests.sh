#!/bin/bash
# Runs everything in test/ that's actually buildable and runnable without
# external data this repo doesn't ship. See README.md's "Testing" section
# for what each of these covers and which ones assert pass/fail vs. print
# output for manual inspection.
#
# Usage: test/run_tests.sh   (run from anywhere; cd's to the repo root)

set -uo pipefail
cd "$(dirname "$0")/.."

FAIL=0

echo "=== check_correctness.sh (genotyping regression vs. golden VCF) ==="
test/check_correctness.sh || FAIL=1

echo
echo "=== snp_tree_test (brute-force vs. interval-tree agreement, asserts) ==="
make -j"$(nproc)" test/snp_tree_test >/dev/null
( cd test && ./snp_tree_test ) || FAIL=1

echo
echo "=== read_vcf_alleles_test (diagnostic -- inspect output, no assertions) ==="
make -j"$(nproc)" test/read_vcf_alleles_test >/dev/null
( cd test && ./read_vcf_alleles_test input/chr1_regions.bed input/1kg.chr1.imputed.vcf.gz > /dev/null )
( cd test && ./read_vcf_alleles_test input/chr1_regions_v2.bed input/1kg.chr1.imputed.vcf.gz > /dev/null )

echo
echo "=== fast_ops_test / haplotype_test (diagnostic -- inspect output, no assertions) ==="
make -j"$(nproc)" test/fast_ops_test test/haplotype_test >/dev/null
( cd test && ./fast_ops_test > /dev/null )
( cd test && ./haplotype_test > /dev/null )

echo
if [ "$FAIL" -eq 0 ]; then
    echo "run_tests.sh: PASS"
else
    echo "run_tests.sh: FAIL"
fi
exit $FAIL
