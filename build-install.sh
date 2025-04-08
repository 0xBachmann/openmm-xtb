#!/bin/bash
set -e
mkdir -p build
cd build
source $(conda info --base)/etc/profile.d/conda.sh
conda activate QMEDS_OpenMM
cmake -DCMAKE_BUILD_TYPE=Release -DCONDA_ROOT=$CONDA_PREFIX ..
make -j8 all
make -j8 test
make -j8 install
make -j8 PythonInstall
