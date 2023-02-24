#!/bin/bash

# Bash strict mode
set -euo pipefail

install_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
marple_pst_dir=$(dirname $install_dir)

if [ -f "$marple_pst_dir"/conda/bin/activate ] ; then
    source "$marple_pst_dir"/conda/bin/activate
    for i in $(seq ${CONDA_SHLVL}); do conda deactivate ; done
    conda env remove -n marple-pst
fi

# remove installations from v0.3.0+
[ -f $install_dir/Mambaforge.sh ] && rm $install_dir/Mambaforge.sh
[ -d "$marple_pst_dir"/conda ] && rm -rf "$marple_pst_dir"/conda
# remove installations from v0.2.0
[ -f $install_dir/marple_pst_miniconda.sh ] && rm $install_dir/marple_pst_miniconda.sh
[ -d $marple_pst_dir/marple_pst_miniconda ] && rm -rf $marple_pst_dir/marple_pst_miniconda

exit 0
