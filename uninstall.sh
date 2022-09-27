#!/bin/bash

# Bash strict mode
set -euo pipefail

marple_pst_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

source "$marple_pst_dir"/marple_pst_miniconda/bin/activate
for i in $(seq ${CONDA_SHLVL}); do conda deactivate ; done
conda env remove -n marple-pst

rm "$marple_pst_dir"/marple_pst_miniconda.sh
rm -rf "$marple_pst_dir"/marple_pst_miniconda
