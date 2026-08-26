#!/bin/bash
set -euo pipefail

# mimalloc/libdeflate's vendored CMake builds need to see conda-build's
# cross-compiling toolchain, not a system one -- plain `make` already picks
# up $CXX/$CC from the environment conda-build sets, so no extra flags
# needed here beyond that.
make -j"${CPU_COUNT}" HipSTR-MT

install -d "${PREFIX}/bin"
install -m 0755 HipSTR-MT "${PREFIX}/bin/HipSTR-MT"
