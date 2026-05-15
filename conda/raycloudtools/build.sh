#!/bin/bash
set -euo pipefail

# --- Build raycloudtools ---
mkdir -p build-rct
pushd build-rct

cmake ${CMAKE_ARGS} \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DWITH_LAS=ON \
  -DWITH_QHULL=ON \
  -DWITH_TIFF=ON \
  -DWITH_TBB=ON \
  -DRAYCLOUD_BUILD_TESTS=OFF \
  -DRAYCLOUD_BUILD_DOXYGEN=OFF \
  ..

ninja -j${CPU_COUNT}
ninja install
popd

# --- Build treetools (separate CMake project; requires raycloudtools installed) ---
mkdir -p build-tree
pushd build-tree

cmake ${CMAKE_ARGS} \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -Draycloudtools_DIR=$PREFIX/lib/cmake/raycloudtools \
  -DWITH_TIFF=ON \
  -DTREE_BUILD_TESTS=OFF \
  -DTREE_BUILD_DOXYGEN=OFF \
  ../treetools

ninja -j${CPU_COUNT}
ninja install
popd
