#!/bin/bash
set -e

# Initialize variables
skip_conda=false
conda_env="EDS_OpenMM"  # Default conda environment
shift_count=0

# Function to display help message
show_help() {
    echo "Usage: $0 [ENV_NAME | --no-conda OPENMM_ROOT XTB_ROOT INSTALL_ROOT PYTHON_ROOT]"
    echo ""
    echo "Build and install the OpenMM xtb Plugin."
    echo ""
    echo "Options:"
    echo "  --help, -h       Show this help message and exit"
    echo "  --no-conda       Skip conda environment activation"
    echo ""
    echo "Arguments:"
    echo "  ENV_NAME         Conda environment to activate (default: $conda_env"
    echo "                   Use --no-conda to skip conda activation"
    echo "  OPENMM_ROOT      Path to OpenMM installation (default: \$CONDA_PREFIX)"
    echo "  XTB_ROOT         Path to XTB installation (default: \$CONDA_PREFIX)"
    echo "  INSTALL_ROOT     Installation directory (default: \$CONDA_PREFIX)"
    echo "  PYTHON_ROOT      Path to Python installation (default: \$CONDA_PREFIX)"
    echo ""
    exit 0
}

# Check for help flag
for arg in "$@"; do
    if [ "$arg" = "--help" ] || [ "$arg" = "-h" ]; then
        show_help
    fi
done


# Parse first argument - either ENV_NAME or --no-conda
if [ $# -gt 0 ]; then
    if [ "$1" = "--no-conda" ]; then
        skip_conda=true
        shift_count=1
    elif [[ "$1" != -* ]]; then  # If it doesn't start with a dash, treat as ENV_NAME
        conda_env="$1"
        shift_count=1
    fi
fi

# Create and enter build directory
mkdir -p _build
cd _build

# Activate conda if not skipped
if [ "$skip_conda" = false ]; then
    if command -v conda &> /dev/null; then
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate "$conda_env"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to activate conda environment '$conda_env'"
            exit 1
        fi
    else
        echo "Warning: conda not found but --no-conda was not specified"
    fi
fi

# Shift arguments if needed to account for env name or --no-conda
args=("$@")
if [ $shift_count -gt 0 ]; then
    args=("${@:$((shift_count+1))}")
fi

# Set arguments with defaults
OPENMM_ROOT=${args[0]:-$CONDA_PREFIX}
XTB_ROOT=${args[1]:-$CONDA_PREFIX}
INSTALL_ROOT=${args[2]:-$CONDA_PREFIX}
PYTHON_ROOT=${args[3]:-$CONDA_PREFIX}

# Display configuration
echo "Build configuration:"
echo "OPENMM_ROOT: $OPENMM_ROOT"
echo "XTB_ROOT: $XTB_ROOT"
echo "INSTALL_ROOT: $INSTALL_ROOT"
echo "PYTHON_ROOT: $PYTHON_ROOT"
if [ "$skip_conda" = true ]; then
    echo "Conda activation: skipped"
else
    echo "Conda environment: $conda_env"
fi
echo ""

# Configure and build
cmake -DCMAKE_BUILD_TYPE=Release \
      -DOPENMM_ROOT=$OPENMM_ROOT \
      -DXTB_ROOT=$XTB_ROOT \
      -DINSTALL_ROOT=$INSTALL_ROOT \
      -DPYTHON_ROOT=$PYTHON_ROOT ..

make -j8 all
make -j8 test
make -j8 install
make -j8 PythonInstall