#!/bin/bash

# Bash strict mode
set -euo pipefail

install_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
marple_pst_dir=$(dirname $install_dir)
echo "Installing into $marple_pst_dir"

# is_conda_active=$(conda --version > /dev/null 2> /dev/null && echo "yes" || echo "no")
# if [ "$is_conda_active" = "yes" ]; then
#     echo "SUCCESS: nothing needed doing, conda is already installed and activated" 
#     exit 0
# fi

full_operating_system_name="$(uname -s)"
case "${full_operating_system_name}" in
    Linux*)     operating_system_name=Linux;;
    Darwin*)    operating_system_name=Mac;;
    *)          operating_system_name="unsupported"
esac

if [ "$operating_system_name" = "unsupported" ]; then
    echo "ERROR: only Linux and Mac are supported, try installing manually" >&2
    exit 1
fi

machine_name="$(uname -m)"
# if [ "$machine_name" != "x86_64" ]; then
#     echo "ERROR: only x86_64 is supported, try installing manually" >&2
#     exit 1
# fi

# for conda_directory in ~/minconda ~/.conda ~/.anaconda ; do
#     if [ -d $conda_directory ]; then
#         echo "ERROR: halted because conda may already be installed - $conda_directory already exists. You may need to activate it." >&2
#         exit 1
#     fi
# done

# for conda_file in ~/miniconda.sh ; do
#     if [ -f $conda_file ]; then
#         echo "ERROR: halted because conda may already be installed - $conda_file already exists. You may need to activate it." >&2
#         exit 1
#     fi
# done

echo "INFO: Downloading miniconda installer"
if [ "$operating_system_name" = "Mac" ] && [ "$machine_name" = "x86_64" ]; then
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o "$install_dir"/marple_pst_miniconda.sh
fi
if [ "$operating_system_name" = "Mac" ] && [ "$machine_name" = "arm64" ]; then
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o "$install_dir"/marple_pst_miniconda.sh
fi
if [ "$operating_system_name" = "Linux" ]; then
    curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o "$install_dir"/marple_pst_miniconda.sh
fi
bash "$install_dir"/marple_pst_miniconda.sh -b -p "$marple_pst_dir"/marple_pst_miniconda
source "$marple_pst_dir"/marple_pst_miniconda/bin/activate
conda --version
echo "SUCCESS: activate conda with this command:"
echo "source ~/marple_pst_miniconda/bin/activate"
# for i in $(seq ${CONDA_SHLVL}); do conda deactivate ; done

if [ "$operating_system_name" = "Mac" ] && [ "$machine_name" = "arm64" ]; then
    CONDA_SUBDIR=osx-64 conda env create --force -f env.yml
    conda activate marple-pst
    conda config --env --set subdir osx-64
else 
    conda env create --force -f env.yml
    conda activate marple-pst
fi

./test/run_tests.sh
