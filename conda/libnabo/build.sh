#!/bin/bash
set -euo pipefail

mkdir -p build
cd build

cmake ${CMAKE_ARGS} \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DLIBNABO_BUILD_TESTS=OFF \
  -DLIBNABO_BUILD_EXAMPLES=OFF \
  ..

ninja -j${CPU_COUNT}
ninja install
